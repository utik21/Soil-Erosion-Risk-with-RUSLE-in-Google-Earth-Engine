# =============================================================================
# Soil Erosion Risk with RUSLE - analysis in R
# Upland sub-watershed east of Malang, East Java, likely part of the Amprong
# sub-watershed (upper Brantas basin), on the western slopes of the Tengger massif
# -----------------------------------------------------------------------------
# Self-initiated portfolio project using public data.
# Author: Utik Tri Wulan Cahya, PhD
#
# Input (exported from Google Earth Engine with RUSLE_soil_erosion_GEE.js):
#   RUSLE_layers.tif          A, R, K, LS, C, P, slope_pct, landcover, erosion_class (30 m)
#   RUSLE_area_by_class.csv   area per land cover x erosion class (cross-check)
#   RUSLE_watershed.geojson   watershed boundary (HydroBASINS level 9)
#
# Input files are read from data/ if that folder exists, otherwise from the
# project folder. Results are written to output/.
# Run:  Rscript RUSLE_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ggplot2)
  library(patchwork)
})

IN  <- if (dir.exists("data")) "data" else "."
inp <- function(f) file.path(IN, f)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
sink("output/analysis_output.txt", split = TRUE)
section <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"; nodata <- "#d9d8d3"
# Ordinal single-hue ramp for the five erosion classes (validated: monotone,
# visible steps, lightest step clears 2:1 on the light surface)
class_pal <- c("#eda070", "#e07a42", "#c85a24", "#9c3f16", "#652509")
class_lab <- c("Very light (< 15)", "Light (15-60)", "Moderate (60-180)",
               "Heavy (180-480)", "Very heavy (> 480)")
# Categorical slots 1-3 (validated for all pairs)
lc_pal <- c("Cropland" = "#2a78d6", "Grassland" = "#eb6834", "Tree cover" = "#1baf7a")

theme_pub <- function(base = 11) {
  theme_minimal(base_size = base, base_family = "DejaVu Sans") +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
          axis.title = element_text(colour = ink2),
          plot.title = element_text(face = "bold", size = base + 3),
          plot.subtitle = element_text(colour = ink2),
          plot.caption = element_text(colour = ink2, size = base - 2, hjust = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
          legend.position = "top", legend.justification = "left",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title.position = "plot", plot.caption.position = "plot")
}
theme_map <- function(base = 11) {
  theme_void(base_family = "DejaVu Sans") +
    theme(text = element_text(colour = ink),
          plot.title = element_text(face = "bold", size = base + 3, colour = ink),
          plot.subtitle = element_text(colour = ink2, size = base - 0.5),
          plot.caption = element_text(colour = ink2, size = base - 2.5, hjust = 0),
          legend.text = element_text(colour = ink, size = base - 1.5),
          legend.title = element_text(colour = ink2, size = base - 1),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(10, 12, 10, 12),
          plot.title.position = "plot", plot.caption.position = "plot")
}

# =============================================================================
section("1. DATA")
# =============================================================================
r   <- rast(inp("RUSLE_layers.tif"))
bnd <- st_transform(st_read(inp("RUSLE_watershed.geojson"), quiet = TRUE), 32749)
px_ha <- prod(res(r)) / 1e4

df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
df <- df[!is.na(df$landcover), ]
lc_names <- c("10" = "Tree cover", "20" = "Shrubland", "30" = "Grassland", "40" = "Cropland",
              "50" = "Built-up", "60" = "Bare / sparse", "80" = "Water", "90" = "Wetland")
df$lc <- factor(lc_names[as.character(df$landcover)],
                levels = c("Tree cover", "Grassland", "Cropland", "Bare / sparse", "Built-up", "Water"))
df$slope_cls <- cut(df$slope_pct, c(-Inf, 8, 15, 25, 40, Inf),
                    labels = c("0-8%", "8-15%", "15-25%", "25-40%", "> 40%"))
soil <- df[!is.na(df$A), ]

cat(sprintf("Watershed: HydroBASINS %s, %.1f km2 (catalogue area %.1f km2)\n",
            bnd$HYBAS_ID, nrow(df) * px_ha / 100, bnd$SUB_AREA))
cat(sprintf("Soil-bearing area analysed (built-up and water excluded): %.0f ha\n", nrow(soil) * px_ha))
lc_area <- aggregate(list(area_ha = rep(px_ha, nrow(df))), list(land_cover = df$lc), sum)
lc_area$share_pct <- round(100 * lc_area$area_ha / sum(lc_area$area_ha), 1)
lc_area$area_ha <- round(lc_area$area_ha)
cat("\nLand cover (ESA WorldCover 2021):\n"); print(lc_area, row.names = FALSE)
write.csv(lc_area, "output/tables/land_cover_area.csv", row.names = FALSE)

# Cross-check with the area table computed in Earth Engine
ee_tab <- read.csv(inp("RUSLE_area_by_class.csv"))
cat(sprintf("\nCross-check: Earth Engine area table %.0f ha vs raster %.0f ha (difference %.1f%%)\n",
            sum(ee_tab$area_ha), nrow(soil) * px_ha,
            100 * (sum(ee_tab$area_ha) / (nrow(soil) * px_ha) - 1)))

# =============================================================================
section("2. RUSLE FACTORS")
# =============================================================================
fsum <- function(x) c(min = min(x), p05 = unname(quantile(x, .05)), median = median(x),
                      mean = mean(x), p95 = unname(quantile(x, .95)), max = max(x))
factors <- do.call(rbind, lapply(c("R", "K", "LS", "C", "P", "slope_pct", "A"),
                                 function(b) data.frame(factor = b, t(signif(fsum(soil[[b]]), 3)))))
print(factors, row.names = FALSE)
write.csv(factors, "output/tables/rusle_factor_summary.csv", row.names = FALSE)
# Mean annual rainfall implied by R (inverse of Renard & Freimund 1994, P > 850 mm)
P_from_R <- function(R) (1.219 + sqrt(1.219^2 - 4 * 0.004105 * (587.8 - R))) / (2 * 0.004105)
cat(sprintf("\nImplied mean annual rainfall (CHIRPS 2014-2023): %.0f to %.0f mm\n",
            P_from_R(min(soil$R)), P_from_R(max(soil$R))))

# =============================================================================
section("3. SOIL LOSS AND EROSION CLASSES")
# =============================================================================
total_t <- sum(soil$A * px_ha)
cat(sprintf("Median soil loss: %.1f t/ha/yr | area-weighted mean: %.1f t/ha/yr\n",
            median(soil$A), total_t / (nrow(soil) * px_ha)))
cat(sprintf("Estimated gross soil loss for the watershed: %.2f million t/yr\n", total_t / 1e6))

cls <- aggregate(list(area_ha = rep(px_ha, nrow(soil))), list(class = soil$erosion_class), sum)
cls$label <- class_lab[cls$class]
cls$share_pct <- round(100 * cls$area_ha / sum(cls$area_ha), 1)
cls$area_ha <- round(cls$area_ha)
cat("\nArea by erosion hazard class (t/ha/yr):\n"); print(cls[, c("label", "area_ha", "share_pct")], row.names = FALSE)
write.csv(cls, "output/tables/area_by_erosion_class.csv", row.names = FALSE)

# Concentration: how much of the soil loss comes from the worst areas
o  <- sort(soil$A, decreasing = TRUE)
cs <- cumsum(o) / sum(o)
conc <- data.frame(area_share_pct = c(1, 5, 10, 20),
                   soil_loss_share_pct = round(100 * cs[ceiling(c(.01, .05, .10, .20) * length(o))], 1))
cat("\nConcentration of soil loss:\n"); print(conc, row.names = FALSE)
write.csv(conc, "output/tables/soil_loss_concentration.csv", row.names = FALSE)

# =============================================================================
section("4. SOIL LOSS BY LAND COVER AND SLOPE")
# =============================================================================
main_lc <- c("Tree cover", "Grassland", "Cropland")
s3 <- soil[soil$lc %in% main_lc, ]
s3$lc <- droplevels(s3$lc)

by_lc <- do.call(rbind, lapply(main_lc, function(l) {
  x <- s3[s3$lc == l, ]
  data.frame(land_cover = l, area_ha = round(nrow(x) * px_ha),
             median_C = signif(median(x$C), 2),
             median_A = round(median(x$A), 1),
             share_on_slopes_over_40pct = round(100 * mean(x$slope_pct > 40), 1),
             share_of_total_soil_loss_pct = round(100 * sum(x$A) / sum(soil$A), 1))
}))
cat("Summary by land cover:\n"); print(by_lc, row.names = FALSE)
write.csv(by_lc, "output/tables/soil_loss_by_land_cover.csv", row.names = FALSE)

grid <- aggregate(list(median_A = s3$A), list(land_cover = s3$lc, slope_class = s3$slope_cls), median)
cnt  <- aggregate(list(area_ha = rep(px_ha, nrow(s3))), list(land_cover = s3$lc, slope_class = s3$slope_cls), sum)
grid <- merge(grid, cnt)
grid$median_A <- round(grid$median_A, 1); grid$area_ha <- round(grid$area_ha)
grid <- grid[order(grid$land_cover, grid$slope_class), ]
cat("\nMedian soil loss (t/ha/yr) by land cover x slope class:\n"); print(grid, row.names = FALSE)
write.csv(grid, "output/tables/soil_loss_land_cover_x_slope.csv", row.names = FALSE)

cat("\nNo significance tests are reported: neighbouring 30 m pixels are spatially\n",
    "autocorrelated, so pixel-level p-values would be misleadingly small.\n", sep = "")

# Conservation priority: farmed or grazed land (cropland, grassland) on slopes
# > 25% with heavy or very heavy erosion (>= 180 t/ha/yr)
prio <- soil$lc %in% c("Cropland", "Grassland") & soil$slope_pct > 25 & soil$A >= 180
cat(sprintf("\nPriority area (cropland or grassland, slope > 25%%, soil loss >= 180 t/ha/yr): %.0f ha\n",
            sum(prio) * px_ha))
cat(sprintf("  = %.1f%% of the analysed area, producing %.1f%% of the estimated soil loss\n",
            100 * mean(prio), 100 * sum(soil$A[prio]) / sum(soil$A)))
print(round(tapply(rep(px_ha, sum(prio)), droplevels(soil$lc[prio]), sum)))

steep <- soil$slope_pct > 100
cat(sprintf("\nSlopes > 100%% (beyond the range RUSLE was developed for): %.1f%% of area, %.1f%% of soil loss\n",
            100 * mean(steep), 100 * sum(soil$A[steep]) / sum(soil$A)))

# =============================================================================
section("4b. SENSITIVITY: LIKELY 2023 FIRE SCAR")
# =============================================================================
# A compact patch near the caldera rim (north-east of the watershed) had very low
# NDVI in 2023 (median C about 0.25). It matches the area burned in the Bromo
# savanna fire of September 2023, so its C factor - and soil loss - reflect
# temporary post-fire conditions rather than normal land cover.
fire_box <- soil$x > 706500 & soil$x < 712000 & soil$y > 9115800 & soil$y < 9118500
fire <- fire_box & soil$C >= 0.1
fire_ha <- sum(fire) * px_ha
cat(sprintf("Patch with low 2023 cover (C >= 0.1) near the caldera rim: %.0f ha (%.1f%% of area)
",
            fire_ha, 100 * mean(fire)))
cat(sprintf("  median C %.2f | median soil loss %.0f t/ha/yr | %.1f%% of the estimated soil loss
",
            median(soil$C[fire]), median(soil$A[fire]), 100 * sum(soil$A[fire]) / sum(soil$A)))
cat("  land cover before the fire (WorldCover 2021, ha):\n")
print(round(tapply(rep(px_ha, sum(fire)), droplevels(soil$lc[fire]), sum)))
nf <- soil[!fire, ]
onf <- sort(nf$A, decreasing = TRUE)
prio_nf <- nf$lc %in% c("Cropland", "Grassland") & nf$slope_pct > 25 & nf$A >= 180
sens <- data.frame(
  measure = c("Median soil loss (t/ha/yr)", "Share of soil loss from the worst 10% of land (%)",
              "Grassland share of soil loss (%)", "Priority area (ha)", "Priority area share of soil loss (%)"),
  all_land = c(round(median(soil$A), 1), conc$soil_loss_share_pct[3],
               by_lc$share_of_total_soil_loss_pct[by_lc$land_cover == "Grassland"],
               round(sum(prio) * px_ha), round(100 * sum(soil$A[prio]) / sum(soil$A), 1)),
  without_fire_patch = c(round(median(nf$A), 1),
                         round(100 * sum(onf[1:ceiling(0.1 * length(onf))]) / sum(onf), 1),
                         round(100 * sum(nf$A[nf$lc == "Grassland"]) / sum(nf$A), 1),
                         round(sum(prio_nf) * px_ha), round(100 * sum(nf$A[prio_nf]) / sum(nf$A), 1)))
cat("\nResults with and without the patch:\n"); print(sens, row.names = FALSE)
write.csv(sens, "output/tables/sensitivity_fire_patch.csv", row.names = FALSE)

# =============================================================================
section("5. FIGURES")
# =============================================================================
ext <- ext(r)
sb_x0 <- ext$xmin + 14500; sb_y0 <- ext$ymin + 700       # scale bar in the empty area below the valley
scalebar <- list(
  annotate("rect", xmin = sb_x0, xmax = sb_x0 + 5000, ymin = sb_y0, ymax = sb_y0 + 260, fill = ink),
  annotate("rect", xmin = sb_x0 + 2500, xmax = sb_x0 + 5000, ymin = sb_y0, ymax = sb_y0 + 260, fill = "white", colour = ink, linewidth = 0.3),
  annotate("text", x = c(sb_x0, sb_x0 + 2500, sb_x0 + 5000), y = sb_y0 + 650, label = c("0", "2.5", "5 km"),
           size = 3, colour = ink),
  annotate("text", x = ext$xmax - 1200, y = ext$ymax - 700, label = "N", fontface = "bold", size = 4.5, colour = ink),
  annotate("segment", x = ext$xmax - 1200, xend = ext$xmax - 1200, y = ext$ymax - 3000, yend = ext$ymax - 1300,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"), colour = ink, linewidth = 0.6)
)

# --- Figure 1: erosion hazard class map --------------------------------------
df$class_f <- factor(ifelse(is.na(df$erosion_class), NA, class_lab[df$erosion_class]), levels = class_lab)
fig1 <- ggplot() +
  geom_raster(data = df[is.na(df$A), ], aes(x, y), fill = nodata) +
  geom_raster(data = df[!is.na(df$A), ], aes(x, y, fill = class_f)) +
  geom_sf(data = bnd, fill = NA, colour = ink, linewidth = 0.4) +
  scalebar +
  annotate("rect", xmin = 706500, xmax = 712000, ymin = 9115800, ymax = 9118500,
           fill = NA, colour = ink, linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 706300, y = 9117150, hjust = 1, size = 3.1, colour = ink, lineheight = 0.95,
           label = "Likely burned in the\nSeptember 2023 Bromo fire") +
  scale_fill_manual(values = setNames(class_pal, class_lab), name = "Soil loss\n(t/ha/yr)", drop = FALSE) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(title = "Erosion hazard concentrates on the steep upper slopes",
       subtitle = sprintf("RUSLE soil loss, %.0f km2 sub-watershed east of Malang (Amprong area). Grey: built-up and water (not assessed).",
                          nrow(df) * px_ha / 100),
       caption = paste0("Classes as used for erosion hazard mapping in Indonesia. West = lower watershed towards Malang; east = Tengger massif.\n",
                        "Data: CHIRPS, OpenLandMap, Copernicus DEM GLO-30, MERIT Hydro, Sentinel-2, ESA WorldCover. Self-initiated portfolio project using public data.")) +
  theme_map() + theme(legend.position = "right")
ggsave("output/figures/fig1_erosion_class_map.png", fig1, width = 11, height = 4.7, dpi = 300, bg = "white")

# --- Figure 2: the five factors ----------------------------------------------
blue_ramp <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
fmap <- function(col, title, lims = NULL, trans = "identity", labels = waiver(), breaks = waiver()) {
  d <- soil[, c("x", "y", col)]; names(d)[3] <- "v"
  if (!is.null(lims)) d$v <- pmin(pmax(d$v, lims[1]), lims[2])
  ggplot() +
    geom_raster(data = df[is.na(df$A), ], aes(x, y), fill = nodata) +
    geom_raster(data = d, aes(x, y, fill = v)) +
    geom_sf(data = bnd, fill = NA, colour = ink, linewidth = 0.3) +
    scale_fill_gradientn(colours = blue_ramp, limits = lims, trans = trans, name = NULL,
                         labels = labels, breaks = breaks) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title) + theme_map(10) +
    theme(plot.title = element_text(size = 10.5, face = "bold"),
          legend.key.height = unit(0.35, "cm"), legend.key.width = unit(1.4, "cm"),
          legend.position = "bottom", plot.margin = margin(4, 6, 4, 6))
}
f_R  <- fmap("R",  "R: rainfall erosivity (MJ mm/ha/h/yr)")
f_K  <- fmap("K",  "K: soil erodibility (t ha h/ha/MJ/mm)")
f_LS <- fmap("LS", "LS: slope length and steepness (capped at 40)", lims = c(0, 40))
f_C  <- fmap("C",  "C: cover management (from 2023 NDVI)", lims = c(0.0001, 1), trans = "log10",
             breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1"))
fig2 <- (f_R + f_K) / (f_LS + f_C) +
  plot_annotation(title = "Slope and vegetation cover, not rainfall or soil, drive the pattern",
                  subtitle = "RUSLE factor maps. K varies little at 250 m soil-map resolution; LS and C vary by orders of magnitude.",
                  caption = "P (support practice) = 1 except on cropland, where contour farming values by slope class were applied (0.5-1.0).",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, family = "DejaVu Sans", colour = ink),
                                plot.subtitle = element_text(colour = ink2, size = 10, family = "DejaVu Sans"),
                                plot.caption = element_text(colour = ink2, size = 8.5, hjust = 0, family = "DejaVu Sans")))
ggsave("output/figures/fig2_rusle_factors.png", fig2, width = 11, height = 5.6, dpi = 250, bg = "white")

# --- Figure 3: median soil loss by land cover x slope ------------------------
grid$lc_f <- factor(grid$land_cover, levels = rev(main_lc))
grid$cls  <- cut(grid$median_A, c(-Inf, 15, 60, 180, 480, Inf), labels = class_lab)
grid$lab  <- ifelse(grid$median_A >= 100, sprintf("%.0f", grid$median_A), sprintf("%.1f", grid$median_A))
grid$sub  <- sprintf("%s ha", format(grid$area_ha, big.mark = ","))
fig3 <- ggplot(grid, aes(slope_class, lc_f)) +
  geom_tile(aes(fill = cls), colour = "white", linewidth = 2) +
  geom_text(aes(label = lab, colour = as.integer(cls) >= 3), fontface = "bold", size = 4.6, nudge_y = 0.1) +
  geom_text(aes(label = sub, colour = as.integer(cls) >= 3), size = 3, nudge_y = -0.22) +
  scale_fill_manual(values = setNames(class_pal, class_lab), name = "Median soil loss (t/ha/yr)", drop = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink), guide = "none") +
  labs(title = "Cropland and grassland on steep slopes erode most",
       subtitle = "Median soil loss by land cover and slope class; small text = area of each combination",
       x = "Slope class", y = NULL,
       caption = "Tree cover keeps soil loss low even on slopes over 40%. Bare/sparse land (< 4 ha) is not shown.") +
  guides(fill = guide_legend(nrow = 1, title.position = "top")) +
  theme_pub() + theme(panel.grid = element_blank(), legend.position = "bottom",
                      legend.text = element_text(size = 9), legend.key.width = unit(0.5, "cm"),
                      axis.text = element_text(size = 10.5, colour = ink))
ggsave("output/figures/fig3_soil_loss_land_cover_slope.png", fig3, width = 9, height = 4.8, dpi = 300, bg = "white")

# --- Figure 4: concentration curve -------------------------------------------
cc <- data.frame(area = seq_along(o) / length(o) * 100, loss = cs * 100)[seq(1, length(o), 50), ]
fig4 <- ggplot(cc, aes(area, loss)) +
  geom_abline(slope = 1, intercept = 0, colour = ink2, linetype = "dashed", linewidth = 0.4) +
  geom_line(colour = "#c85a24", linewidth = 1.1) +
  annotate("point", x = 10, y = conc$soil_loss_share_pct[3], size = 3, colour = "#c85a24") +
  annotate("text", x = 12, y = conc$soil_loss_share_pct[3] - 4, hjust = 0, size = 3.8, colour = ink,
           label = sprintf("The worst 10%% of the land\nproduces %.0f%% of the soil loss", conc$soil_loss_share_pct[3])) +
  annotate("text", x = 62, y = 52, label = "if erosion were evenly spread", angle = 36, size = 3.1, colour = ink2) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, 20)) +
  coord_equal() +
  labs(title = "A small share of land produces most soil loss",
       subtitle = "Cumulative share of estimated soil loss, from the most eroding pixels down",
       x = "Share of analysed area (%)", y = "Share of soil loss (%)") +
  theme_pub() + theme(plot.margin = margin(10, 20, 10, 10))
ggsave("output/figures/fig4_soil_loss_concentration.png", fig4, width = 6.8, height = 6.4, dpi = 300, bg = "white")

# --- Cover image for Upwork (16:9) -------------------------------------------
cov_map <- fig1 + labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
cov_heat <- fig3 + labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL) +
  theme(legend.position = "none", axis.text = element_text(size = 9.5, colour = ink),
        plot.margin = margin(2, 0, 0, 0))
cover <- cov_map / cov_heat + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    title = "Soil Erosion Risk with RUSLE in Google Earth Engine",
    subtitle = "Amprong area, east of Malang, East Java | erosion hazard map and median soil loss (t/ha/yr) by land cover and slope",
    caption = "Self-initiated portfolio project using public data (CHIRPS, OpenLandMap, Copernicus DEM, MERIT Hydro, Sentinel-2, ESA WorldCover)",
    theme = theme(plot.title = element_text(face = "bold", size = 17, family = "DejaVu Sans", colour = ink),
                  plot.subtitle = element_text(size = 11, family = "DejaVu Sans", colour = ink2),
                  plot.caption = element_text(colour = ink2, size = 8.5, hjust = 0, family = "DejaVu Sans"),
                  plot.background = element_rect(fill = "white", colour = NA),
                  plot.margin = margin(14, 22, 10, 22)))
ggsave("output/figures/cover_upwork.png", cover, width = 10, height = 5.625, dpi = 200, bg = "white")

section("SESSION INFO")
print(sessionInfo())
sink()
message("Done.")
