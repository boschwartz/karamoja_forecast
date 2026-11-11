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
  summarise(geometry = st_union(geometry), .groups = "drop")

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
# 3. Query Open-Meteo Forecast API (ECMWF IFS HRES 9km)
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
      forecast_days = 9  # Today + 8 future days guarantees full Day +1 to Day +7
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
      req_lat   = batch_df$lat[idx],
      req_lon   = batch_df$lon[idx],
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

forecast_df <- imap_dfr(grid_batches, ~ fetch_batch(.x, .y, length(grid_batches)))

# ==============================================================================
# 4. Filter to Tomorrow through Day +7 (EAT) & Compute Trigger Metrics
# ==============================================================================
thresh_7d_mm     <- 50.0   # 50 mm in 7 days
thresh_1d_mm     <- 20.0   # 20 mm in any 1 day
trigger_quantile <- 0.90   # Spatial P90 (top 10% of district grid cells)

today_eat    <- as.Date(format(Sys.time(), tz = "Africa/Kampala", usetz = FALSE))
target_dates <- seq(from = today_eat + 1, by = "1 day", length.out = 7)

missing_dates <- setdiff(as.character(target_dates), as.character(unique(forecast_df$date)))
if (length(missing_dates) > 0) {
  stop("Missing forecast dates from API: ", paste(missing_dates, collapse = ", "))
}

cat(sprintf("7-Day Window (Day +1 to Day +7): %s to %s\n",
            min(target_dates), max(target_dates)))

# Deduplicate to unique native IFS cells so snapped points aren't double-counted
forecast_7d <- forecast_df |>
  filter(date %in% target_dates) |>
  group_by(ifs_lat, ifs_lon, date) |>
  summarise(precip_mm = mean(precip_mm, na.rm = TRUE), .groups = "drop") |>
  group_by(ifs_lat, ifs_lon) |>
  mutate(cell_id = cur_group_id()) |>
  ungroup()

grid_summary <- forecast_7d |>
  group_by(cell_id, ifs_lat, ifs_lon) |>
  summarise(
    total_7d_mm = sum(precip_mm, na.rm = TRUE),
    max_1d_mm   = max(precip_mm, na.rm = TRUE),
    .groups     = "drop"
  )

# Smoothly interpolate across the bounding box (0.08 deg -> 0.01 deg)
pts_vect <- vect(grid_summary, geom = c("ifs_lon", "ifs_lat"), crs = "EPSG:4326")

r_coarse <- rast(
  ext(min(grid_summary$ifs_lon) - 0.04, max(grid_summary$ifs_lon) + 0.04,
      min(grid_summary$ifs_lat) - 0.04, max(grid_summary$ifs_lat) + 0.04),
  resolution = 0.08,
  crs = "EPSG:4326"
)

r_binned <- rasterize(pts_vect, r_coarse, field = "total_7d_mm", fun = "mean")
r_filled <- focal(r_binned, w = 3, fun = "mean", na.policy = "only")
r_fine   <- rast(ext(r_filled), resolution = 0.01, crs = "EPSG:4326")
r_smooth <- resample(r_filled, r_fine, method = "bilinear")
r_smooth[r_smooth < 0] <- 0

precip_breaks <- c(0, 2, 5, 10, 25, 50, 75, 100, 150, Inf)
precip_labels <- c("0 – 2", "2 – 5", "5 – 10", "10 – 25", "25 – 50",
                   "50 – 75", "75 – 100", "100 – 150", "> 150")

df_map <- as.data.frame(r_smooth, xy = TRUE, na.rm = TRUE)
names(df_map) <- c("lon", "lat", "precip_mm")

df_map <- df_map |>
  mutate(
    precip_cat = factor(
      cut(precip_mm, breaks = precip_breaks, labels = precip_labels,
          include.lowest = TRUE, right = FALSE),
      levels = precip_labels
    )
  )

# District-level spatial P90 & Mean calculations
cells_sf <- st_as_sf(grid_summary, coords = c("ifs_lon", "ifs_lat"), crs = 4326)

district_stats <- map_dfr(c("Moroto", "Napak", "Amudat"), function(dist_name) {
  poly_sf    <- districts_sf |> filter(district == dist_name)
  inside_ids <- cells_sf$cell_id[st_intersects(cells_sf, poly_sf, sparse = FALSE)[, 1]]
  if (length(inside_ids) == 0) {
    inside_ids <- cells_sf$cell_id[st_nearest_feature(poly_sf, cells_sf)]
  }
  
  cells_7d <- grid_summary |> filter(cell_id %in% inside_ids)
  mean_7d  <- mean(cells_7d$total_7d_mm, na.rm = TRUE)
  p90_7d   <- unname(quantile(cells_7d$total_7d_mm, probs = trigger_quantile, na.rm = TRUE))
  max_7d   <- max(cells_7d$total_7d_mm, na.rm = TRUE)
  
  daily_df <- forecast_7d |>
    filter(cell_id %in% inside_ids) |>
    group_by(date) |>
    summarise(
      mean_1d = mean(precip_mm, na.rm = TRUE),
      p90_1d  = unname(quantile(precip_mm, probs = trigger_quantile, na.rm = TRUE)),
      max_1d  = max(precip_mm, na.rm = TRUE),
      .groups = "drop"
    )
  
  peak_day_idx <- which.max(daily_df$p90_1d)
  
  tibble(
    district     = dist_name,
    n_cells      = length(inside_ids),
    mean_7d      = mean_7d,
    trig_7d      = p90_7d,
    max_7d       = max_7d,
    mean_1d_peak = max(daily_df$mean_1d, na.rm = TRUE),
    trig_1d      = daily_df$p90_1d[peak_day_idx],
    trig_1d_date = format(daily_df$date[peak_day_idx], "%b %d"),
    max_1d       = max(daily_df$max_1d, na.rm = TRUE),
    hit_7d       = p90_7d >= thresh_7d_mm,
    hit_1d       = daily_df$p90_1d[peak_day_idx] >= thresh_1d_mm,
    triggered    = (p90_7d >= thresh_7d_mm) | (daily_df$p90_1d[peak_day_idx] >= thresh_1d_mm)
  )
})

print(district_stats)

# ==============================================================================
# 5. Build Single Bounding-Box Map + Right-Hand Trigger Summary Panel
# ==============================================================================
map_xlim <- c(min(df_map$lon), max(df_map$lon))
map_ylim <- c(min(df_map$lat), max(df_map$lat))

label_pts <- districts_sf |>
  st_point_on_surface() |>
  st_coordinates() |>
  as_tibble() |>
  bind_cols(district = districts_sf$district)

p_map <- ggplot() +
  geom_raster(data = df_map, aes(x = lon, y = lat, fill = precip_cat)) +
  geom_sf(data = uga_adm2, fill = NA, colour = "white", linewidth = 0.25, alpha = 0.5) +
  geom_sf(data = districts_sf, fill = NA, colour = "white", linewidth = 1.4) +
  geom_sf(data = districts_sf, fill = NA, colour = "#0F172A", linewidth = 0.75) +
  geom_label(
    data          = label_pts,
    aes(x = X, y = Y, label = district),
    fontface      = "bold",
    size          = 3.4,
    fill          = alpha("white", 0.82),
    colour        = "#0F172A",
    label.size    = 0.25,
    label.padding = unit(0.18, "lines")
  ) +
  scale_fill_viridis_d(
    option    = "mako",
    direction = -1,
    name      = "7-Day Rain\n(mm)",
    limits    = precip_labels,
    drop      = FALSE
  ) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  scale_x_continuous(breaks = seq(33.8, 35.2, by = 0.4)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(
    panel.border    = element_rect(colour = "grey30", fill = NA, linewidth = 0.6),
    panel.grid      = element_blank(),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 10),
    axis.text       = element_text(colour = "grey25", size = 9)
  ) +
  guides(fill = guide_legend(reverse = TRUE))

card_positions <- district_stats |>
  mutate(
    idx        = row_number(),
    ymax       = 1 - (idx - 1) * 0.31,
    ymin       = ymax - 0.28,
    badge_lbl  = if_else(triggered, "ACTION TRIGGERED", "BELOW TRIGGER"),
    badge_bg   = if_else(triggered, "#DC2626", "#E2E8F0"),
    badge_fg   = if_else(triggered, "white", "#334155"),
    border_col = if_else(triggered, "#DC2626", "#CBD5E1"),
    card_bg    = if_else(triggered, "#FEF2F2", "#F8FAFC"),
    body_text  = sprintf(
      paste0(
        "7-Day Total  (Trigger >= %.0f mm)\n",
        "  P90 (10%% area): %5.1f mm %s  |  Mean: %.1f mm (Max: %.1f)\n\n",
        "Max 1-Day    (Trigger >= %.0f mm)\n",
        "  P90 (10%% area): %5.1f mm (%s)%s |  Mean: %.1f mm"
      ),
      thresh_7d_mm,
      trig_7d, if_else(hit_7d, "[*]", "   "), mean_7d, max_7d,
      thresh_1d_mm,
      trig_1d, trig_1d_date, if_else(hit_1d, " [*]", "    "), mean_1d_peak
    )
  )

p_summary <- ggplot(card_positions) +
  annotate(
    "text", x = 0, y = 1.04, hjust = 0, vjust = 0,
    label = sprintf("Trigger Rule: P90 >= %.0f mm / 7d  OR  >= %.0f mm / 1d  ([*] = exceeded)",
                    thresh_7d_mm, thresh_1d_mm),
    size = 3.2, fontface = "bold", colour = "#334155"
  ) +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax,
        fill = card_bg, colour = border_col),
    linewidth = 0.7
  ) +
  geom_text(
    aes(x = 0.03, y = ymax - 0.035, label = paste0(district, " District")),
    hjust = 0, vjust = 1, fontface = "bold", size = 4.1, colour = "#0F172A"
  ) +
  geom_label(
    aes(x = 0.97, y = ymax - 0.035, label = badge_lbl,
        fill = badge_bg, colour = badge_fg),
    hjust = 1, vjust = 1, fontface = "bold", size = 2.8,
    label.size = 0, label.padding = unit(0.22, "lines")
  ) +
  geom_text(
    aes(x = 0.03, y = ymax - 0.095, label = body_text),
    hjust = 0, vjust = 1, family = "mono", size = 3.05,
    lineheight = 1.22, colour = "#1E293B"
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.05, 1.08), clip = "off") +
  theme_void()

# ==============================================================================
# 6. Combine Map + Summary Panel & Save PNG Artifact
# ==============================================================================
final_plot <- (p_map + p_summary) +
  plot_layout(widths = c(1.15, 1.35)) +
  plot_annotation(
    title    = "ECMWF IFS HRES (9 km): 7-Day Accumulated Precipitation & Trigger Monitor",
    subtitle = sprintf("Moroto, Napak & Amudat Districts, Uganda  |  Tomorrow – Day +7: %s to %s (EAT)",
                       format(min(target_dates), "%d %b %Y"),
                       format(max(target_dates), "%d %b %Y")),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15, colour = "#0F172A"),
      plot.subtitle = element_text(size = 11, colour = "grey30")
    )
  )

ggsave(
  filename = "karamoja_ifs_hres_trigger_dashboard.png",
  plot     = final_plot,
  width    = 13.5,
  height   = 7.8,
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
  triggered_names <- district_stats |>
    filter(triggered) |>
    pull(district) |>
    paste(collapse = ", ")
  
  headline <- if (any_triggered) {
    sprintf("WE HAVE TRIGGERED: %s (%s to %s)",
            triggered_names,
            format(min(target_dates), "%d %b"),
            format(max(target_dates), "%d %b"))
  } else {
    sprintf("TEST CHECK - No Districts Triggered (%s)", format(today_eat, "%d %b"))
  }
  
  summary_lines <- district_stats |>
    mutate(
      line = sprintf(
        "%s: %s | 7d P90 = %.1f mm (Mean %.1f) | Max 1d P90 = %.1f mm on %s",
        district,
        if_else(triggered, "TRIGGERED", "Below trigger"),
        trig_7d, mean_7d, trig_1d, trig_1d_date
      )
    ) |>
    pull(line) |>
    paste(collapse = " &#10; ")
  
  # Unique ID per day (or per manual test run) so Teams RSS picks it up cleanly
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

# Simple HTML page so clicking the link in Teams opens today's map right in your browser
index_html <- sprintf(
  paste0(
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    "<title>Karamoja Trigger Monitor</title>",
    "<style>body{font-family:sans-serif;margin:24px;background:#F8FAFC;color:#0F172A;}",
    "img{max-width:100%%;height:auto;border:1px solid #CBD5E1;border-radius:6px;}</style></head>",
    "<body><h2>Karamoja Daily IFS HRES Dashboard</h2>",
    "<p>Last updated: %s UTC | <a href='trigger_feed.xml'>RSS Trigger Feed</a></p>",
    "<img src='karamoja_ifs_hres_trigger_dashboard.png' alt='Karamoja Forecast Dashboard'/>",
    "</body></html>"
  ),
  format(Sys.time(), "%Y-%m-%d %H:%M", tz = "UTC")
)
writeLines(index_html, "public/index.html")

cat("Wrote public/trigger_feed.xml and public/index.html (triggered =", any_triggered, ")\n")