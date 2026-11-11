library(jsonlite)
library(tidyr)
library(ggplot2)
library(httr2)
library(dplyr)
library(purrr)
library(terra)
library(sf)
library(patchwork)

# ==============================================================================
# 1. Load District Boundaries (Moroto, Amudat, Napak)
# ==============================================================================
target_districts <- c("Moroto", "Amudat", "Napak")
cache_file       <- "uga_districts_adm2.gpkg"

if (file.exists(cache_file)) {
  cat("Loading cached Uganda district boundaries...\n")
  uga_adm2 <- st_read(cache_file, quiet = TRUE)
} else {
  endpoints <- c(
    "https://www.geoboundaries.org/api/current/gbHumanitarian/UGA/ADM2/",
    "https://www.geoboundaries.org/api/current/gbOpen/UGA/ADM3/"
  )
  
  uga_adm2 <- NULL
  for (ep in endpoints) {
    candidate <- tryCatch({
      gb_meta <- request(ep) |> req_perform() |> resp_body_json()
      st_read(gb_meta$gjDownloadURL, quiet = TRUE) |>
        st_make_valid() |>
        st_transform(4326)
    }, error = function(e) NULL)
    
    if (!is.null(candidate)) {
      clean_names <- trimws(sub("\\s+district$", "", candidate$shapeName, ignore.case = TRUE))
      if (all(toupper(target_districts) %in% toupper(clean_names))) {
        cat("Downloaded modern district boundaries from:\n ", ep, "\n")
        uga_adm2 <- candidate
        st_write(uga_adm2, cache_file, quiet = TRUE, delete_dsn = TRUE)
        break
      }
    }
  }
  
  if (is.null(uga_adm2)) {
    stop("Could not retrieve Moroto, Amudat, and Napak from geoBoundaries.")
  }
}

districts_sf <- uga_adm2 |>
  mutate(
    clean_name = trimws(sub("\\s+district$", "", shapeName, ignore.case = TRUE)),
    district   = tools::toTitleCase(tolower(clean_name))
  ) |>
  filter(toupper(clean_name) %in% toupper(target_districts)) |>
  group_by(district) |>
  summarise(.groups = "drop")

cat("Matched districts:", paste(districts_sf$district, collapse = ", "), "\n")

# ==============================================================================
# 2. Build ~9 km (0.08 deg) Query Grid Over the 3-District Bounding Box
# ==============================================================================
bbox <- st_bbox(districts_sf)
pad  <- 0.20

grid <- expand.grid(
  lat = round(seq(floor((bbox["ymin"] - pad) * 10) / 10,
                  ceiling((bbox["ymax"] + pad) * 10) / 10,
                  by = 0.08), 4),
  lon = round(seq(floor((bbox["xmin"] - pad) * 10) / 10,
                  ceiling((bbox["xmax"] + pad) * 10) / 10,
                  by = 0.08), 4)
) |>
  mutate(pt_id = row_number())

cat("Total ~9km grid points to query:", nrow(grid), "\n")

batch_size   <- 80
grid_batches <- split(grid, ceiling(seq_len(nrow(grid)) / batch_size))
cat("Created", length(grid_batches), "batches.\n")

# ==============================================================================
# 3. Query Open-Meteo Forecast API (ECMWF IFS HRES 9km Daily Totals)
# ==============================================================================
base_url <- "https://api.open-meteo.com/v1/forecast"

fetch_batch <- function(batch_df, batch_idx, total_batches) {
  cat(sprintf("Fetching batch %s of %s...\n", batch_idx, total_batches))
  
  req <- request(base_url) |>
    req_url_query(
      latitude      = paste(batch_df$lat, collapse = ","),
      longitude     = paste(batch_df$lon, collapse = ","),
      daily         = "precipitation_sum",
      models        = "ecmwf_ifs",
      timezone      = "Africa/Kampala",
      forecast_days = 9  # Guarantees full Tomorrow (Day +1) through Day +7
    ) |>
    req_retry(
      max_tries    = 5,
      backoff      = ~ .x * 3,
      is_transient = function(resp) resp_status(resp) %in% c(429, 500, 502, 503, 504)
    )
  
  resp <- req_perform(req)
  raw_data <- resp_body_json(resp, simplifyVector = FALSE)
  
  if (!is.null(names(raw_data)) && "daily" %in% names(raw_data)) {
    raw_data <- list(raw_data)
  }
  
  Sys.sleep(1)
  
  map2_dfr(raw_data, seq_len(nrow(batch_df)), function(loc, idx) {
    tibble(
      pt_id     = batch_df$pt_id[idx],
      ifs_lat   = loc$latitude,   # Native ~9km IFS cell centre
      ifs_lon   = loc$longitude,  # Native ~9km IFS cell centre
      date      = as.Date(unlist(loc$daily$time)),
      precip_mm = map_dbl(
        loc$daily$precipitation_sum,
        ~ if (is.null(.x)) NA_real_ else as.numeric(.x)
      )
    )
  })
}

daily_raw <- imap_dfr(grid_batches, ~ fetch_batch(.x, .y, length(grid_batches)))

# ==============================================================================
# 4. Compute Tomorrow (Day +1) & 7-Day Accumulations + District Trigger Metrics
# ==============================================================================
thresh_tom_mm    <- 20.0   # Immediate trigger: 20 mm tomorrow (00:00-24:00 EAT)
thresh_7d_mm     <- 50.0   # Outlook trigger:   50 mm in 7 days (Tomorrow - Day +7)
trigger_quantile <- 0.90   # Spatial P90 (top 10% of district grid cells)

# Anchor timestamps to East Africa Time (UTC+3)
now_eat         <- as.POSIXct(format(Sys.time(), tz = "Africa/Kampala", usetz = FALSE),
                              tz = "Africa/Kampala")
today_eat       <- as.Date(format(now_eat, "%Y-%m-%d"))
tomorrow_date   <- today_eat + 1
target_dates_7d <- seq(from = tomorrow_date, by = "1 day", length.out = 7)

missing_dates <- setdiff(as.character(target_dates_7d), as.character(unique(daily_raw$date)))
if (length(missing_dates) > 0) {
  stop("Missing forecast dates from API: ", paste(missing_dates, collapse = ", "))
}

label_tom_window <- sprintf("%s (00:00 – 24:00 EAT)", format(tomorrow_date, "%A, %d %b %Y"))
label_7d_window  <- sprintf("%s to %s (EAT)",
                            format(min(target_dates_7d), "%d %b %Y"),
                            format(max(target_dates_7d), "%d %b %Y"))

cat("Tomorrow Window:", label_tom_window, "\n")
cat("7-Day Window:   ", label_7d_window, "\n")

# Deduplicate daily data to unique native IFS cells for Tomorrow & 7-Day windows
daily_dedup <- daily_raw |>
  filter(date %in% target_dates_7d) |>
  group_by(ifs_lat, ifs_lon, date) |>
  summarise(precip_mm = mean(precip_mm, na.rm = TRUE), .groups = "drop")

grid_summary <- daily_dedup |>
  group_by(ifs_lat, ifs_lon) |>
  summarise(
    rain_tom_mm = sum(precip_mm[date == tomorrow_date], na.rm = TRUE),
    rain_7d_mm  = sum(precip_mm, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  mutate(cell_id = row_number())

# Helper to smoothly interpolate any gridded field (0.08 deg -> 0.01 deg)
smooth_field <- function(df, field_col, breaks, labels) {
  pts_vect <- vect(df, geom = c("ifs_lon", "ifs_lat"), crs = "EPSG:4326")
  r_coarse <- rast(
    ext(min(df$ifs_lon) - 0.04, max(df$ifs_lon) + 0.04,
        min(df$ifs_lat) - 0.04, max(df$ifs_lat) + 0.04),
    resolution = 0.08,
    crs = "EPSG:4326"
  )
  r_binned <- rasterize(pts_vect, r_coarse, field = field_col, fun = "mean")
  r_filled <- focal(r_binned, w = 3, fun = "mean", na.policy = "only")
  r_fine   <- rast(ext(r_filled), resolution = 0.01, crs = "EPSG:4326")
  r_smooth <- resample(r_filled, r_fine, method = "bilinear")
  r_smooth[r_smooth < 0] <- 0
  
  out_df <- as.data.frame(r_smooth, xy = TRUE, na.rm = TRUE)
  names(out_df) <- c("lon", "lat", "precip_mm")
  out_df |>
    mutate(
      precip_cat = factor(
        cut(precip_mm, breaks = breaks, labels = labels,
            include.lowest = TRUE, right = FALSE),
        levels = labels
      )
    )
}

# Colour bins: 20 mm sits cleanly on a break for Tomorrow; 50 mm sits on a break for 7d
breaks_tom <- c(0, 1, 2, 5, 10, 15, 20, 30, 50, Inf)
labels_tom <- c("0 – 1", "1 – 2", "2 – 5", "5 – 10", "10 – 15",
                "15 – 20", "20 – 30", "30 – 50", "> 50")

breaks_7d  <- c(0, 2, 5, 10, 25, 50, 75, 100, 150, Inf)
labels_7d  <- c("0 – 2", "2 – 5", "5 – 10", "10 – 25", "25 – 50",
                "50 – 75", "75 – 100", "100 – 150", "> 150")

df_map_tom <- smooth_field(grid_summary, "rain_tom_mm", breaks_tom, labels_tom)
df_map_7d  <- smooth_field(grid_summary, "rain_7d_mm",  breaks_7d,  labels_7d)

# District-level spatial P90, Mean & Max for BOTH windows
cells_sf <- st_as_sf(grid_summary, coords = c("ifs_lon", "ifs_lat"), crs = 4326)

district_stats <- map_dfr(c("Moroto", "Napak", "Amudat"), function(dist_name) {
  poly_sf    <- districts_sf |> filter(district == dist_name)
  inside_ids <- cells_sf$cell_id[st_intersects(cells_sf, poly_sf, sparse = FALSE)[, 1]]
  if (length(inside_ids) == 0) {
    inside_ids <- cells_sf$cell_id[st_nearest_feature(poly_sf, cells_sf)]
  }
  
  sub_cells <- grid_summary |> filter(cell_id %in% inside_ids)
  
  p90_tom  <- unname(quantile(sub_cells$rain_tom_mm, probs = trigger_quantile, na.rm = TRUE))
  mean_tom <- mean(sub_cells$rain_tom_mm, na.rm = TRUE)
  max_tom  <- max(sub_cells$rain_tom_mm, na.rm = TRUE)
  
  p90_7d   <- unname(quantile(sub_cells$rain_7d_mm, probs = trigger_quantile, na.rm = TRUE))
  mean_7d  <- mean(sub_cells$rain_7d_mm, na.rm = TRUE)
  max_7d   <- max(sub_cells$rain_7d_mm, na.rm = TRUE)
  
  tibble(
    district  = dist_name,
    n_cells   = length(inside_ids),
    p90_tom   = p90_tom,
    mean_tom  = mean_tom,
    max_tom   = max_tom,
    hit_tom   = p90_tom >= thresh_tom_mm,
    p90_7d    = p90_7d,
    mean_7d   = mean_7d,
    max_7d    = max_7d,
    hit_7d    = p90_7d >= thresh_7d_mm,
    triggered = (p90_tom >= thresh_tom_mm) | (p90_7d >= thresh_7d_mm)
  )
})

print(district_stats)

# ==============================================================================
# 5. Build Two-Section Dashboard (Row 1 = Tomorrow, Row 2 = 7-Day Outlook)
# ==============================================================================
map_xlim <- c(min(df_map_7d$lon), max(df_map_7d$lon))
map_ylim <- c(min(df_map_7d$lat), max(df_map_7d$lat))

label_pts <- districts_sf |>
  st_point_on_surface() |>
  st_coordinates() |>
  as_tibble() |>
  bind_cols(district = districts_sf$district)

build_precip_map <- function(df_raster, legend_labels, legend_title, map_title, map_sub) {
  ggplot() +
    geom_raster(data = df_raster, aes(x = lon, y = lat, fill = precip_cat)) +
    geom_sf(data = uga_adm2, fill = NA, colour = "white", linewidth = 0.25, alpha = 0.5) +
    geom_sf(data = districts_sf, fill = NA, colour = "white", linewidth = 1.35) +
    geom_sf(data = districts_sf, fill = NA, colour = "#0F172A", linewidth = 0.7) +
    geom_label(
      data          = label_pts,
      aes(x = X, y = Y, label = district),
      fontface      = "bold",
      size          = 3.2,
      fill          = alpha("white", 0.82),
      colour        = "#0F172A",
      label.size    = 0.25,
      label.padding = unit(0.16, "lines")
    ) +
    scale_fill_viridis_d(
      option    = "mako",
      direction = -1,
      name      = legend_title,
      limits    = legend_labels,
      drop      = FALSE
    ) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    scale_x_continuous(breaks = seq(33.8, 35.2, by = 0.4)) +
    labs(title = map_title, subtitle = map_sub, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 10.5) +
    theme(
      plot.title      = element_text(face = "bold", size = 12.5, colour = "#0F172A"),
      plot.subtitle   = element_text(size = 9.5, colour = "grey30"),
      panel.border    = element_rect(colour = "grey30", fill = NA, linewidth = 0.6),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 9.5),
      axis.text       = element_text(colour = "grey25", size = 8.5)
    ) +
    guides(fill = guide_legend(reverse = TRUE))
}

build_cards_panel <- function(cards_df, header_rule_text) {
  ggplot(cards_df) +
    annotate(
      "text", x = 0, y = 1.03, hjust = 0, vjust = 0,
      label = header_rule_text,
      size = 3.2, fontface = "bold", colour = "#334155"
    ) +
    geom_rect(
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax,
          fill = card_bg, colour = border_col),
      linewidth = 0.75
    ) +
    geom_text(
      aes(x = 0.035, y = ymax - 0.04, label = paste0(district, " District")),
      hjust = 0, vjust = 1, fontface = "bold", size = 4.0, colour = "#0F172A"
    ) +
    geom_label(
      aes(x = 0.965, y = ymax - 0.04, label = badge_lbl,
          fill = badge_bg, colour = badge_fg),
      hjust = 1, vjust = 1, fontface = "bold", size = 2.75,
      label.size = 0, label.padding = unit(0.22, "lines")
    ) +
    geom_text(
      aes(x = 0.035, y = ymax - 0.115, label = primary_line),
      hjust = 0, vjust = 1, family = "mono", fontface = "bold",
      size = 3.3, colour = "#0F172A"
    ) +
    geom_text(
      aes(x = 0.035, y = ymax - 0.19, label = secondary_line),
      hjust = 0, vjust = 1, family = "mono",
      size = 3.0, colour = "#475569"
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.05, 1.08), clip = "off") +
    theme_void()
}

# --- Row 1: Tomorrow (Day +1) Map + 3 District Cards ---
p_map_tom <- build_precip_map(
  df_raster     = df_map_tom,
  legend_labels = labels_tom,
  legend_title  = "Tomorrow\nRain (mm)",
  map_title     = "Immediate Forecast: Tomorrow (Day +1)",
  map_sub       = label_tom_window
)

cards_tom_df <- district_stats |>
  mutate(
    idx            = row_number(),
    ymax           = 1 - (idx - 1) * 0.31,
    ymin           = ymax - 0.27,
    badge_lbl      = if_else(hit_tom, "TOMORROW TRIGGERED", "BELOW 20 MM"),
    badge_bg       = if_else(hit_tom, "#DC2626", "#E2E8F0"),
    badge_fg       = if_else(hit_tom, "white", "#334155"),
    border_col     = if_else(hit_tom, "#DC2626", "#CBD5E1"),
    card_bg        = if_else(hit_tom, "#FEF2F2", "#F8FAFC"),
    primary_line   = sprintf("Tomorrow P90 (10%% area): %5.1f mm  [Trigger >= %.0f mm]",
                             p90_tom, thresh_tom_mm),
    secondary_line = sprintf("District Mean: %.1f mm   |   Peak 9km Cell: %.1f mm",
                             mean_tom, max_tom)
  )

p_cards_tom <- build_cards_panel(
  cards_tom_df,
  sprintf("Immediate Trigger Rule: Spatial P90 >= %.0f mm Tomorrow (%s)",
          thresh_tom_mm, format(tomorrow_date, "%d %b"))
)

# --- Row 2: 7-Day Cumulative Map + 3 District Cards ---
p_map_7d <- build_precip_map(
  df_raster     = df_map_7d,
  legend_labels = labels_7d,
  legend_title  = "7-Day Rain\n(mm)",
  map_title     = "7-Day Outlook: Tomorrow – Day +7",
  map_sub       = label_7d_window
)

cards_7d_df <- district_stats |>
  mutate(
    idx            = row_number(),
    ymax           = 1 - (idx - 1) * 0.31,
    ymin           = ymax - 0.27,
    badge_lbl      = if_else(hit_7d, "7-DAY TRIGGERED", "BELOW 50 MM"),
    badge_bg       = if_else(hit_7d, "#DC2626", "#E2E8F0"),
    badge_fg       = if_else(hit_7d, "white", "#334155"),
    border_col     = if_else(hit_7d, "#DC2626", "#CBD5E1"),
    card_bg        = if_else(hit_7d, "#FEF2F2", "#F8FAFC"),
    primary_line   = sprintf("7-Day P90 (10%% area):    %5.1f mm  [Trigger >= %.0f mm]",
                             p90_7d, thresh_7d_mm),
    secondary_line = sprintf("District Mean: %.1f mm   |   Peak 9km Cell: %.1f mm",
                             mean_7d, max_7d)
  )

p_cards_7d <- build_cards_panel(
  cards_7d_df,
  sprintf("7-Day Trigger Rule: Spatial P90 >= %.0f mm (Tomorrow through Day +7)", thresh_7d_mm)
)

# ==============================================================================
# 6. Assemble 2-Row Dashboard & Save PNG Artifact
# ==============================================================================
row_tom <- (p_map_tom + p_cards_tom) + plot_layout(widths = c(1.12, 1.38))
row_7d  <- (p_map_7d  + p_cards_7d)  + plot_layout(widths = c(1.12, 1.38))

final_plot <- (row_tom / row_7d) +
  plot_annotation(
    title    = "ECMWF IFS HRES (9 km): Karamoja Precipitation Trigger Monitor",
    subtitle = sprintf("Moroto, Napak & Amudat Districts, Uganda  |  Run: %s EAT",
                       format(now_eat, "%d %b %Y %H:%M")),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 16, colour = "#0F172A"),
      plot.subtitle = element_text(size = 11, colour = "grey30")
    )
  )

ggsave(
  filename = "karamoja_ifs_hres_trigger_dashboard.png",
  plot     = final_plot,
  width    = 13.5,
  height   = 13.2,
  dpi      = 300,
  bg       = "white"
)

# ==============================================================================
# 7. Write Daily RSS Trigger Feed & Landing Page for GitHub Pages + Teams
# ==============================================================================
any_triggered <- any(district_stats$triggered)
send_always   <- tolower(Sys.getenv("SEND_TEST_ITEM", "false")) == "true"

pages_url <- Sys.getenv("PAGES_BASE_URL")
run_url   <- Sys.getenv("GITHUB_RUN_URL")
link_url  <- if (nzchar(pages_url)) pages_url else if (nzchar(run_url)) run_url else "https://github.com"

escape_xml <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

if (any_triggered || send_always) {
  triggered_details <- district_stats |>
    filter(triggered) |>
    mutate(
      reason = case_when(
        hit_tom & hit_7d ~ paste0(district, " (Tomorrow & 7-Day)"),
        hit_tom          ~ paste0(district, " (Tomorrow >=20mm)"),
        hit_7d           ~ paste0(district, " (7-Day >=50mm)")
      )
    ) |>
    pull(reason) |>
    paste(collapse = ", ")
  
  headline <- if (any_triggered) {
    sprintf("WE HAVE TRIGGERED: %s", triggered_details)
  } else {
    sprintf("TEST CHECK - No Districts Triggered (Forecast for %s)",
            format(tomorrow_date, "%d %b"))
  }
  
  summary_lines <- district_stats |>
    mutate(
      status_lbl = case_when(
        hit_tom & hit_7d ~ "TRIGGERED [Tomorrow + 7d]",
        hit_tom          ~ "TRIGGERED [Tomorrow]",
        hit_7d           ~ "TRIGGERED [7-Day]",
        TRUE             ~ "Below triggers"
      ),
      line = sprintf(
        "%s: %s | Tomorrow (%s) P90 = %.1f mm (Mean %.1f) | 7d P90 = %.1f mm (Mean %.1f)",
        district, status_lbl, format(tomorrow_date, "%d %b"),
        p90_tom, mean_tom, p90_7d, mean_7d
      )
    ) |>
    pull(line) |>
    paste(collapse = " &#10; ")
  
  guid_val <- if (send_always && !any_triggered) {
    paste0("karamoja-test-", as.integer(Sys.time()))
  } else {
    paste0("karamoja-trigger-", as.character(today_eat))
  }
  
  rss_xml <- sprintf(
    paste0(
      '<?xml version="1.0" encoding="UTF-8"?>\n',
      '<rss version="2.0"><channel>\n',
      '  <title>Karamoja IFS HRES Trigger Monitor</title>\n',
      '  <link>%s</link>\n',
      '  <description>Daily ECMWF IFS HRES trigger alerts for Moroto, Napak and Amudat</description>\n',
      '  <item>\n',
      '    <guid isPermaLink="false">%s</guid>\n',
      '    <pubDate>%s</pubDate>\n',
      '    <title>%s</title>\n',
      '    <link>%s</link>\n',
      '    <description>%s</description>\n',
      '  </item>\n',
      '</channel></rss>\n'
    ),
    escape_xml(link_url),
    guid_val,
    format(Sys.time(), "%a, %d %b %Y %H:%M:%S GMT", tz = "GMT"),
    escape_xml(headline),
    escape_xml(link_url),
    summary_lines
  )
} else {
  rss_xml <- sprintf(
    paste0(
      '<?xml version="1.0" encoding="UTF-8"?>\n',
      '<rss version="2.0"><channel>\n',
      '  <title>Karamoja IFS HRES Trigger Monitor</title>\n',
      '  <link>%s</link>\n',
      '  <description>No active triggers today</description>\n',
      '</channel></rss>\n'
    ),
    escape_xml(link_url)
  )
}

dir.create("public", showWarnings = FALSE)
writeLines(rss_xml, "public/trigger_feed.xml")
file.copy("karamoja_ifs_hres_trigger_dashboard.png",
          "public/karamoja_ifs_hres_trigger_dashboard.png",
          overwrite = TRUE)

index_html <- sprintf(
  paste0(
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    "<title>Karamoja Trigger Monitor</title>",
    "<style>body{font-family:sans-serif;margin:24px;background:#F8FAFC;color:#0F172A;max-width:1200px;}",
    "img{width:100%%;height:auto;border:1px solid #CBD5E1;border-radius:6px;background:white;}</style></head>",
    "<body><h2>Karamoja Daily IFS HRES Trigger Dashboard</h2>",
    "<p>Last updated: %s EAT (%s UTC) | <a href='trigger_feed.xml'>RSS Trigger Feed</a></p>",
    "<img src='karamoja_ifs_hres_trigger_dashboard.png?v=%d' alt='Karamoja Forecast Dashboard'/>",
    "</body></html>"
  ),
  format(now_eat, "%Y-%m-%d %H:%M"),
  format(Sys.time(), "%H:%M", tz = "UTC"),
  as.integer(Sys.time())
)
writeLines(index_html, "public/index.html")

cat("Wrote public/trigger_feed.xml and public/index.html (triggered =", any_triggered, ")\n")
