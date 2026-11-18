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
# 1. Load Boundaries: Uganda Districts (ADM2) + Kenya Sub-Counties (via GADM ADM3 Wards)
# ==============================================================================
fetch_geoboundaries <- function(endpoints, target_names, cache_path, strip_regex) {
  if (file.exists(cache_path)) {
    cat("Loading cached boundaries:", cache_path, "\n")
    return(st_read(cache_path, quiet = TRUE) |> st_make_valid() |> st_transform(4326))
  }
  
  for (ep in endpoints) {
    candidate <- tryCatch({
      gb_meta <- request(ep) |>
        req_options(connecttimeout = 30, ipresolve = 1L) |>
        req_retry(max_tries = 4, retry_on_failure = TRUE) |>
        req_perform() |>
        resp_body_json()
      st_read(gb_meta$gjDownloadURL, quiet = TRUE) |>
        st_make_valid() |>
        st_transform(4326)
    }, error = function(e) NULL)
    
    if (!is.null(candidate)) {
      clean <- trimws(sub(strip_regex, "", candidate$shapeName, ignore.case = TRUE))
      if (all(toupper(target_names) %in% toupper(clean))) {
        cat("Downloaded boundaries from:\n ", ep, "\n")
        st_write(candidate, cache_path, quiet = TRUE, delete_dsn = TRUE)
        return(candidate)
      }
    }
  }
  stop("Could not retrieve required boundaries for: ", paste(target_names, collapse = ", "))
}

# 1A. Uganda ADM2 Districts (Moroto, Napak, Amudat)
uga_targets <- c("Moroto", "Napak", "Amudat")
uga_adm2 <- fetch_geoboundaries(
  endpoints    = c("https://www.geoboundaries.org/api/current/gbHumanitarian/UGA/ADM2/",
                   "https://www.geoboundaries.org/api/current/gbOpen/UGA/ADM3/"),
  target_names = uga_targets,
  cache_path   = "uga_districts_adm2.gpkg",
  strip_regex  = "\\s+district$"
)

uga_targets_sf <- uga_adm2 |>
  mutate(
    clean_name = trimws(sub("\\s+district$", "", shapeName, ignore.case = TRUE)),
    area_name  = tools::toTitleCase(tolower(clean_name))
  ) |>
  filter(toupper(clean_name) %in% toupper(uga_targets)) |>
  group_by(area_name) |>
  summarise(.groups = "drop") |>
  mutate(country = "UG", unit_type = "District")

# 1B. Kenya Sub-Counties Dissolved from GADM 4.1 Level-3 Wards (Lokiriama, Loima, Pokot North)
ken_wards_cache <- "ken_wards_gadm41_adm3.gpkg"

if (file.exists(ken_wards_cache)) {
  cat("Loading cached Kenya ADM3 wards:", ken_wards_cache, "\n")
  ken_adm3 <- st_read(ken_wards_cache, quiet = TRUE) |> st_make_valid() |> st_transform(4326)
} else {
  cat("Downloading GADM 4.1 Kenya ADM3 wards (one-time ~480 KB zip)...\n")
  gadm_zip <- tempfile(fileext = ".zip")
  gadm_dir <- tempfile()
  request("https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_KEN_3.json.zip") |>
    req_options(connecttimeout = 30, ipresolve = 1L) |>
    req_timeout(seconds = 120) |>
    req_retry(max_tries = 5, retry_on_failure = TRUE) |>
    req_perform(path = gadm_zip)
  unzip(gadm_zip, exdir = gadm_dir)
  json_file <- list.files(gadm_dir, pattern = "\\.json$", full.names = TRUE)[1]
  ken_adm3  <- st_read(json_file, quiet = TRUE) |> st_make_valid() |> st_transform(4326)
  st_write(ken_adm3, ken_wards_cache, quiet = TRUE, delete_dsn = TRUE)
}

# Normalize ward names (removes spaces/slashes so "Lokiriama/Lorengippi" matches cleanly)
norm_str <- function(x) tolower(gsub("[^a-z0-9]", "", x))

# Exact ward -> sub-county mapping from GADM 4.1 KEN_3:
# - Lokiriama:   Kotaruk/Lobei + Lokiriama/Lorengippi (western border half of Loima)
# - Loima:       Loima + Turkwel (eastern half of Loima)
# - Pokot North: Alale + Kasei + Kiwawa (add "Kapchok", "Kodich", "Suam" if full Kacheliba needed)
lokiriama_wards   <- c("Kotaruk/Lobei", "Lokiriama/Lorengippi")
loima_wards       <- c("Loima", "Turkwel")
pokot_north_wards <- c("Alale", "Kasei", "Kiwawa")

ken_targets_sf <- ken_adm3 |>
  dplyr::mutate(
    area_name = dplyr::case_when(
      NAME_1 == "Turkana"   & NAME_2 == "Loima"     & NAME_3 %in% lokiriama_wards   ~ "Lokiriama",
      NAME_1 == "Turkana"   & NAME_2 == "Loima"     & NAME_3 %in% loima_wards       ~ "Loima",
      NAME_1 == "WestPokot" & NAME_2 == "Kacheliba" & NAME_3 %in% pokot_north_wards ~ "Pokot North",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(area_name)) |>
  dplyr::group_by(area_name) |>
  dplyr::summarise(.groups = "drop") |>
  dplyr::mutate(country = "KE", unit_type = "Sub-County")

if (nrow(ken_targets_sf) < 3) {
  stop("Expected 3 Kenyan sub-counties, matched: ", paste(ken_targets_sf$area_name, collapse = ", "))
}

# Standardize geometry column names before combining
st_geometry(uga_targets_sf) <- "geometry"
st_geometry(ken_targets_sf) <- "geometry"
st_geometry(uga_adm2)       <- "geometry"
st_geometry(ken_adm3)       <- "geometry"

target_order <- c("Moroto", "Napak", "Amudat", "Lokiriama", "Loima", "Pokot North")

areas_sf <- bind_rows(uga_targets_sf, ken_targets_sf) |>
  mutate(
    area_order = match(area_name, target_order),
    card_title = sprintf("%s %s (%s)", area_name, unit_type, country)
  ) |>
  arrange(area_order)

# Dissolve national borders cleanly for the UG/KE international boundary
country_borders_sf <- bind_rows(
  uga_adm2 |> summarise(.groups = "drop") |> mutate(adm0 = "UGA"),
  ken_adm3 |> summarise(.groups = "drop") |> mutate(adm0 = "KEN")
)

cat("Matched project areas:", paste(areas_sf$card_title, collapse = ", "), "\n")

# ==============================================================================
# 2. Build Native ~9 km (0.09 deg) Query Grid (Trimmed to 6 Areas + 15km Halo)
# ==============================================================================
bbox    <- st_bbox(areas_sf)
pad_lat <- 0.24
pad_lon <- 0.36

full_grid <- expand.grid(
  lat = round(seq(floor((bbox["ymin"] - pad_lat) * 10) / 10,
                  ceiling((bbox["ymax"] + pad_lat) * 10) / 10,
                  by = 0.09), 4),
  lon = round(seq(floor((bbox["xmin"] - pad_lon) * 10) / 10,
                  ceiling((bbox["xmax"] + pad_lon) * 10) / 10,
                  by = 0.09), 4)
)

sf_use_s2(FALSE)
mask_poly <- st_buffer(st_union(areas_sf), dist = 0.14)
sf_use_s2(TRUE)

grid_sf <- st_as_sf(full_grid, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
grid    <- full_grid[st_intersects(grid_sf, mask_poly, sparse = FALSE)[, 1], ] |>
  mutate(pt_id = row_number())

cat("Trimmed ~9km grid points to query:", nrow(grid), "\n")

batch_size   <- 60
grid_batches <- split(grid, ceiling(seq_len(nrow(grid)) / batch_size))
cat("Created", length(grid_batches), "batches.\n")

# ==============================================================================
# 3. Query Open-Meteo Forecast API (With SSL/IPv4 Timeout Fix & Local Dev Cache)
# ==============================================================================
base_url <- "https://api.open-meteo.com/v1/forecast"

fetch_batch <- function(batch_df, batch_idx, total_batches) {
  cat(sprintf("Fetching batch %s of %s...\n", batch_idx, total_batches))
  
  req <- request(base_url) |>
    req_headers(`User-Agent` = "Karamoja-CrossBorder-Trigger-Monitor/1.1 (github-actions)") |>
    req_url_query(
      latitude      = paste(batch_df$lat, collapse = ","),
      longitude     = paste(batch_df$lon, collapse = ","),
      daily         = "precipitation_sum",
      models        = "ecmwf_ifs",
      timezone      = "Africa/Kampala",  # UTC+3 (identical to Africa/Nairobi)
      forecast_days = 9
    ) |>
    req_options(
      connecttimeout = 30,  # Prevents 10s SSL handshake abort
      ipresolve      = 1L   # Forces IPv4 inside Docker container
    ) |>
    req_timeout(seconds = 60) |>
    req_retry(
      max_tries        = 6,
      retry_on_failure = TRUE,
      backoff          = ~ .x * 10 + runif(1, 0, 3),
      is_transient     = function(resp) resp_status(resp) %in% c(429, 500, 502, 503, 504)
    )
  
  resp <- req_perform(req)
  raw_data <- resp_body_json(resp, simplifyVector = FALSE)
  
  if (!is.null(names(raw_data)) && "daily" %in% names(raw_data)) {
    raw_data <- list(raw_data)
  }
  
  Sys.sleep(2.0)
  
  map2_dfr(raw_data, seq_len(nrow(batch_df)), function(loc, idx) {
    tibble(
      pt_id     = batch_df$pt_id[idx],
      ifs_lat   = loc$latitude,
      ifs_lon   = loc$longitude,
      date      = as.Date(unlist(loc$daily$time)),
      precip_mm = map_dbl(
        loc$daily$precipitation_sum,
        ~ if (is.null(.x)) NA_real_ else as.numeric(.x)
      )
    )
  })
}

raw_cache <- sprintf("daily_raw_subcounties_%s.rds", Sys.Date())
if (Sys.getenv("GITHUB_ACTIONS") != "true" && file.exists(raw_cache)) {
  cat("Loading locally cached API response from", raw_cache, "...\n")
  daily_raw <- readRDS(raw_cache)
} else {
  daily_raw <- imap_dfr(grid_batches, ~ fetch_batch(.x, .y, length(grid_batches)))
  if (Sys.getenv("GITHUB_ACTIONS") != "true") saveRDS(daily_raw, raw_cache)
}

# ==============================================================================
# 4. Compute Tomorrow (Day +1) & 7-Day Accumulations + Area Trigger Metrics
# ==============================================================================
thresh_tom_mm      <- 20.0   # Red trigger: >= 20.0 mm
orange_tom_mm      <- 15.0   # Orange card: > 15.0 mm
yellow_tom_mm      <- 10.0   # Yellow card: > 10.0 mm

thresh_7d_mm       <- 50.0   # Red trigger: >= 50.0 mm
orange_7d_mm       <- 37.0   # Orange card: > 37.0 mm
yellow_7d_mm       <- 25.0   # Yellow card: > 25.0 mm

trigger_quantile   <- 0.90   # Spatial P90 (top 10% of area grid cells)

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

mask_vect <- vect(st_as_sf(mask_poly))

smooth_field <- function(df, field_col, breaks, labels) {
  pts_vect <- vect(df, geom = c("ifs_lon", "ifs_lat"), crs = "EPSG:4326")
  r_coarse <- rast(
    ext(unname(bbox["xmin"]) - pad_lon, unname(bbox["xmax"]) + pad_lon,
        unname(bbox["ymin"]) - pad_lat, unname(bbox["ymax"]) + pad_lat),
    resolution = 0.09,
    crs = "EPSG:4326"
  )
  r_binned <- rasterize(pts_vect, r_coarse, field = field_col, fun = "mean")
  r_filled <- focal(r_binned, w = 3, fun = "mean", na.policy = "only")
  r_fine   <- rast(ext(r_filled), resolution = 0.015, crs = "EPSG:4326")
  r_smooth <- resample(r_filled, r_fine, method = "bilinear")
  r_smooth[r_smooth < 0] <- 0
  r_masked <- mask(r_smooth, mask_vect)
  
  out_df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
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

precip_breaks <- c(0, 2, 5, 10, 20, 35, 50, 75, 100, 150, Inf)
precip_labels <- c("0 – 2", "2 – 5", "5 – 10", "10 – 20", "20 – 35",
                   "35 – 50", "50 – 75", "75 – 100", "100 – 150", "> 150")

df_map_tom <- smooth_field(grid_summary, "rain_tom_mm", precip_breaks, precip_labels)
df_map_7d  <- smooth_field(grid_summary, "rain_7d_mm",  precip_breaks, precip_labels)

cells_sf <- st_as_sf(grid_summary, coords = c("ifs_lon", "ifs_lat"), crs = 4326)

area_stats <- map_dfr(seq_len(nrow(areas_sf)), function(i) {
  poly_sf    <- areas_sf[i, ]
  inside_ids <- cells_sf$cell_id[st_intersects(cells_sf, poly_sf, sparse = FALSE)[, 1]]
  if (length(inside_ids) == 0) {
    inside_ids <- cells_sf$cell_id[st_nearest_feature(poly_sf, cells_sf)]
  }
  
  sub_cells <- grid_summary |> filter(cell_id %in% inside_ids)
  
  p90_tom  <- round(unname(quantile(sub_cells$rain_tom_mm, probs = trigger_quantile, na.rm = TRUE)), 1)
  mean_tom <- mean(sub_cells$rain_tom_mm, na.rm = TRUE)
  max_tom  <- max(sub_cells$rain_tom_mm, na.rm = TRUE)
  
  p90_7d   <- round(unname(quantile(sub_cells$rain_7d_mm, probs = trigger_quantile, na.rm = TRUE)), 1)
  mean_7d  <- mean(sub_cells$rain_7d_mm, na.rm = TRUE)
  max_7d   <- max(sub_cells$rain_7d_mm, na.rm = TRUE)
  
  hit_tom <- p90_tom >= thresh_tom_mm
  hit_7d  <- p90_7d  >= thresh_7d_mm
  
  tibble(
    area_name  = poly_sf$area_name,
    country    = poly_sf$country,
    card_title = poly_sf$card_title,
    n_cells    = length(inside_ids),
    p90_tom    = p90_tom,
    mean_tom   = mean_tom,
    max_tom    = max_tom,
    hit_tom    = hit_tom,
    p90_7d     = p90_7d,
    mean_7d    = mean_7d,
    max_7d     = max_7d,
    hit_7d     = hit_7d,
    triggered  = hit_tom | hit_7d
  )
})

print(area_stats)

# ==============================================================================
# 5. Build Two-Section Cross-Border Dashboard (6 Traffic-Light Cards + Callouts)
# ==============================================================================
map_xlim <- c(unname(bbox["xmin"]) - pad_lon, unname(bbox["xmax"]) + pad_lon + 0.32)
map_ylim <- c(unname(bbox["ymin"]) - pad_lat, unname(bbox["ymax"]) + pad_lat)

# Uncrossed exterior label slots + short edge-touching pointer targets
callout_df <- tribble(
  ~area_name,    ~lbl_x, ~lbl_y, ~target_x, ~target_y,
  # Left side (Uganda Districts, North -> South)
  "Moroto",       33.98,   3.22,     34.48,     2.85,
  "Napak",        33.82,   2.30,     34.18,     2.38,
  "Amudat",       34.15,   1.48,     34.82,     1.75,
  # Right side (Kenya Sub-Counties, true North -> South order so lines never cross)
  "Loima",        35.92,   3.38,     35.42,     3.25,
  "Lokiriama",    35.95,   2.58,     35.24,     2.68,
  "Pokot North",  35.85,   1.68,     35.32,     1.98
)


build_precip_map <- function(df_raster, map_title, map_sub) {
  ggplot() +
    geom_raster(
      data        = df_raster,
      aes(x = lon, y = lat, fill = precip_cat),
      show.legend = TRUE
    ) +
    borders("world", xlim = c(32, 38), ylim = c(0, 6),
            colour = "grey35", linewidth = 0.45) +
    geom_sf(data = country_borders_sf, fill = NA, colour = "grey20",
            linewidth = 0.5, linetype = "dashed") +
    geom_sf(data = areas_sf, fill = NA, colour = "white", linewidth = 1.35) +
    geom_sf(data = areas_sf, fill = NA, colour = "#0F172A", linewidth = 0.7) +
    geom_segment(
      data = callout_df,
      aes(x = lbl_x, y = lbl_y, xend = target_x, yend = target_y),
      colour = "white", linewidth = 1.1
    ) +
    geom_segment(
      data = callout_df,
      aes(x = lbl_x, y = lbl_y, xend = target_x, yend = target_y),
      colour = "#0F172A", linewidth = 0.45
    ) +
    geom_point(
      data = callout_df,
      aes(x = target_x, y = target_y),
      size = 1.4, colour = "#0F172A"
    ) +
    geom_label(
      data          = callout_df,
      aes(x = lbl_x, y = lbl_y, label = area_name),
      fontface      = "bold",
      size          = 2.85,
      fill          = alpha("white", 0.92),
      colour        = "#0F172A",
      label.size    = 0.25,
      label.padding = unit(0.16, "lines")
    ) +
    scale_fill_viridis_d(
      option    = "mako",
      direction = -1,
      name      = "Rain (mm)",
      limits    = precip_labels,
      drop      = FALSE
    ) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    scale_x_continuous(breaks = seq(34.0, 36.0, by = 0.5)) +
    labs(title = map_title, subtitle = map_sub, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 10.5) +
    theme(
      plot.title       = element_text(face = "bold", size = 12.5, colour = "#0F172A"),
      plot.subtitle    = element_text(size = 9.5, colour = "grey30"),
      panel.background = element_rect(fill = "#F8FAFC", colour = NA),
      panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.6),
      panel.grid       = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 9.5),
      axis.text        = element_text(colour = "grey25", size = 8.5)
    ) +
    guides(fill = guide_legend(reverse = TRUE, override.aes = list(alpha = 1)))
}

# Compact 6-card vertical stack
build_cards_panel <- function(cards_df, header_rule_text) {
  ggplot(cards_df) +
    annotate(
      "text", x = 0, y = 1.022, hjust = 0, vjust = 0,
      label = header_rule_text,
      size = 3.15, fontface = "bold", colour = "#334155"
    ) +
    geom_rect(
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax,
          fill = card_bg, colour = border_col),
      linewidth = 0.72
    ) +
    geom_text(
      aes(x = 0.03, y = ymax - 0.016, label = card_title),
      hjust = 0, vjust = 1, fontface = "bold", size = 3.45, colour = "#0F172A"
    ) +
    geom_label(
      aes(x = 0.97, y = ymax - 0.014, label = badge_lbl,
          fill = badge_bg, colour = badge_fg),
      hjust = 1, vjust = 1, fontface = "bold", size = 2.45,
      label.size = 0, label.padding = unit(0.16, "lines")
    ) +
    geom_text(
      aes(x = 0.03, y = ymax - 0.060, label = primary_line),
      hjust = 0, vjust = 1, family = "mono", fontface = "bold",
      size = 2.9, colour = "#0F172A"
    ) +
    geom_text(
      aes(x = 0.03, y = ymax - 0.101, label = secondary_line),
      hjust = 0, vjust = 1, family = "mono",
      size = 2.75, colour = "#334155"
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.01, 1.065), clip = "off") +
    theme_void()
}

# --- Row 1: Tomorrow (Day +1) Map + 6 Area Cards ---
p_map_tom <- build_precip_map(
  df_raster = df_map_tom,
  map_title = "Immediate Forecast: Tomorrow (Day +1)",
  map_sub   = label_tom_window
)

cards_tom_df <- area_stats |>
  mutate(
    idx        = row_number(),
    ymax       = 1 - (idx - 1) * 0.165,
    ymin       = ymax - 0.146,
    badge_lbl  = if_else(hit_tom, "TOMORROW TRIGGERED", "BELOW 20 MM"),
    card_bg    = case_when(
      hit_tom                 ~ "#FEF2F2",
      p90_tom > orange_tom_mm ~ "#FFF7ED",
      p90_tom > yellow_tom_mm ~ "#FEFCE8",
      TRUE                    ~ "#F8FAFC"
    ),
    border_col = case_when(
      hit_tom                 ~ "#DC2626",
      p90_tom > orange_tom_mm ~ "#EA580C",
      p90_tom > yellow_tom_mm ~ "#CA8A04",
      TRUE                    ~ "#CBD5E1"
    ),
    badge_bg   = case_when(
      hit_tom                 ~ "#DC2626",
      p90_tom > orange_tom_mm ~ "#F97316",
      p90_tom > yellow_tom_mm ~ "#FDE047",
      TRUE                    ~ "#E2E8F0"
    ),
    badge_fg   = case_when(
      hit_tom                 ~ "white",
      p90_tom > orange_tom_mm ~ "white",
      p90_tom > yellow_tom_mm ~ "#1E293B",
      TRUE                    ~ "#334155"
    ),
    primary_line   = sprintf("Tomorrow P90 (10%% area): %5.1f mm  [Trigger >= %.0f mm]",
                             p90_tom, thresh_tom_mm),
    secondary_line = sprintf("Area Mean: %4.1f mm   |   Peak 9km Cell: %5.1f mm (%d cells)",
                             mean_tom, max_tom, n_cells)
  )

p_cards_tom <- build_cards_panel(
  cards_tom_df,
  sprintf("Immediate Trigger Rule: Spatial P90 >= %.1f mm Tomorrow (%s)",
          thresh_tom_mm, format(tomorrow_date, "%d %b"))
)

# --- Row 2: 7-Day Cumulative Map + 6 Area Cards ---
p_map_7d <- build_precip_map(
  df_raster = df_map_7d,
  map_title = "7-Day Outlook: Tomorrow – Day +7",
  map_sub   = label_7d_window
)

cards_7d_df <- area_stats |>
  mutate(
    idx        = row_number(),
    ymax       = 1 - (idx - 1) * 0.165,
    ymin       = ymax - 0.146,
    badge_lbl  = if_else(hit_7d, "7-DAY TRIGGERED", "BELOW 50 MM"),
    card_bg    = case_when(
      hit_7d                 ~ "#FEF2F2",
      p90_7d > orange_7d_mm  ~ "#FFF7ED",
      p90_7d > yellow_7d_mm  ~ "#FEFCE8",
      TRUE                   ~ "#F8FAFC"
    ),
    border_col = case_when(
      hit_7d                 ~ "#DC2626",
      p90_7d > orange_7d_mm  ~ "#EA580C",
      p90_7d > yellow_7d_mm  ~ "#CA8A04",
      TRUE                   ~ "#CBD5E1"
    ),
    badge_bg   = case_when(
      hit_7d                 ~ "#DC2626",
      p90_7d > orange_7d_mm  ~ "#F97316",
      p90_7d > yellow_7d_mm  ~ "#FDE047",
      TRUE                   ~ "#E2E8F0"
    ),
    badge_fg   = case_when(
      hit_7d                 ~ "white",
      p90_7d > orange_7d_mm  ~ "white",
      p90_7d > yellow_7d_mm  ~ "#1E293B",
      TRUE                   ~ "#334155"
    ),
    primary_line   = sprintf("7-Day P90 (10%% area):    %5.1f mm  [Trigger >= %.0f mm]",
                             p90_7d, thresh_7d_mm),
    secondary_line = sprintf("Area Mean: %4.1f mm   |   Peak 9km Cell: %5.1f mm (%d cells)",
                             mean_7d, max_7d, n_cells)
  )

p_cards_7d <- build_cards_panel(
  cards_7d_df,
  sprintf("7-Day Trigger Rule: Spatial P90 >= %.1f mm (Tomorrow through Day +7)",
          thresh_7d_mm)
)

# ==============================================================================
# 6. Assemble 2-Row Cross-Border Dashboard & Save PNG Artifact
# ==============================================================================
row_tom <- (p_map_tom + p_cards_tom) + plot_layout(widths = c(1.15, 1.35))
row_7d  <- (p_map_7d  + p_cards_7d)  + plot_layout(widths = c(1.15, 1.35))

final_plot <- (row_tom / row_7d) +
  plot_annotation(
    title    = "ECMWF IFS HRES (9 km): Karamoja–Turkana–Pokot Cross-Border Trigger Monitor",
    subtitle = sprintf(
      "Uganda (Moroto, Napak, Amudat) & Kenya Sub-Counties (Lokiriama, Loima, Pokot North)  |  Run: %s EAT",
      format(now_eat, "%d %b %Y %H:%M")
    ),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15.5, colour = "#0F172A"),
      plot.subtitle = element_text(size = 10.5, colour = "grey30")
    )
  )

ggsave(
  filename = "karamoja_ifs_hres_trigger_dashboard.png",
  plot     = final_plot,
  width    = 14.4,
  height   = 16.0,
  dpi      = 300,
  bg       = "white"
)

# ==============================================================================
# 7. Write Daily RSS Trigger Feed & Landing Page for GitHub Pages + Teams
# ==============================================================================
any_triggered <- any(area_stats$triggered)
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
  triggered_details <- area_stats |>
    filter(triggered) |>
    mutate(
      reason = case_when(
        hit_tom & hit_7d ~ sprintf("%s-%s (Tomorrow & 7d)", area_name, country),
        hit_tom          ~ sprintf("%s-%s (Tomorrow >=20mm)", area_name, country),
        hit_7d           ~ sprintf("%s-%s (7-Day >=50mm)", area_name, country)
      )
    ) |>
    pull(reason) |>
    paste(collapse = ", ")
  
  headline <- if (any_triggered) {
    sprintf("WE HAVE TRIGGERED: %s", triggered_details)
  } else {
    sprintf("TEST CHECK - No Areas Triggered (Forecast for %s)",
            format(tomorrow_date, "%d %b"))
  }
  
  summary_lines <- area_stats |>
    mutate(
      status_lbl = case_when(
        hit_tom & hit_7d ~ "TRIGGERED [Tomorrow + 7d]",
        hit_tom          ~ "TRIGGERED [Tomorrow]",
        hit_7d           ~ "TRIGGERED [7-Day]",
        TRUE             ~ "Below triggers"
      ),
      line = sprintf(
        "%s (%s): %s | Tomorrow (%s) P90 = %.1f mm (Mean %.1f) | 7d P90 = %.1f mm (Mean %.1f)",
        area_name, country, status_lbl, format(tomorrow_date, "%d %b"),
        p90_tom, mean_tom, p90_7d, mean_7d
      )
    ) |>
    pull(line) |>
    paste(collapse = " &#10; ")
  
  guid_val <- if (send_always && !any_triggered) {
    paste0("karamoja-turkana-test-", as.integer(Sys.time()))
  } else {
    paste0("karamoja-turkana-trigger-", as.character(today_eat))
  }
  
  rss_xml <- sprintf(
    paste0(
      '<?xml version="1.0" encoding="UTF-8"?>\n',
      '<rss version="2.0"><channel>\n',
      '  <title>Karamoja-Turkana-Pokot IFS HRES Trigger Monitor</title>\n',
      '  <link>%s</link>\n',
      '  <description>Daily ECMWF IFS HRES trigger alerts for Moroto, Napak, Amudat (UG) &amp; Lokiriama, Loima, Pokot North (KE)</description>\n',
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
      '  <title>Karamoja-Turkana-Pokot IFS HRES Trigger Monitor</title>\n',
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
    "<title>Karamoja–Turkana–Pokot Trigger Monitor</title>",
    "<style>body{font-family:sans-serif;margin:24px;background:#F8FAFC;color:#0F172A;max-width:1280px;}",
    "img{width:100%%;height:auto;border:1px solid #CBD5E1;border-radius:6px;background:white;}</style></head>",
    "<body><h2>Uganda–Kenya Cross-Border IFS HRES Trigger Dashboard</h2>",
    "<p>Areas: Moroto, Napak, Amudat (UG) &amp; Lokiriama, Loima, Pokot North Sub-Counties (KE) | ",
    "Last updated: %s EAT (%s UTC) | <a href='trigger_feed.xml'>RSS Trigger Feed</a></p>",
    "<img src='karamoja_ifs_hres_trigger_dashboard.png?v=%d' alt='Cross-Border Forecast Dashboard'/>",
    "</body></html>"
  ),
  format(now_eat, "%Y-%m-%d %H:%M"),
  format(Sys.time(), "%H:%M", tz = "UTC"),
  as.integer(Sys.time())
)
writeLines(index_html, "public/index.html")

cat("Wrote public/trigger_feed.xml and public/index.html (triggered =", any_triggered, ")\n")
