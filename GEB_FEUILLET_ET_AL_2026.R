# ==============================================================================
# SPATIAL ANALYSIS OF TREELINE DYNAMICS IN THE PYRENEES
# Hypotheses testing: H1 (elevation), H2 (spatial non-stationarity), H3 (local factors)
# ==============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(sp)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(GWmodel)
  library(mgcv)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
  library(grid)
  library(rstatix)
  library(ggsignif)
  library(broom)
  library(stringr)
  library(colorspace)
})

# ==============================================================================
# DATA PREPARATION
# ==============================================================================

# Load data
df_analysis <- read.csv("data_for_analysis.csv")

# Configuration
CONFIG <- list(
  elev_threshold = 1942,  # Montane-Subalpine threshold (m)
  kernel = "gaussian",
  adaptive = TRUE,
  alpha_sig = 0.10
)

# Prepare datasets for analysis
df_base <- df_analysis %>%
  mutate(
    warming = scale(CLIM_tmax_diff_abs)[,1],
    elev_class = ifelse(DEM_elev_mn < CONFIG$elev_threshold,
                        "MONTANE", "SUBALPINE")
  ) %>%
  drop_na(warming, X, Y)

# Separate datasets by response variable and elevation zone
df_shift_m <- df_base %>% filter(elev_class == "MONTANE") %>% drop_na(DIAC_shift)
df_shift_s <- df_base %>% filter(elev_class == "SUBALPINE") %>% drop_na(DIAC_shift)
df_dens_m <- df_base %>% filter(elev_class == "MONTANE") %>% drop_na(DIAC_dens)
df_dens_s <- df_base %>% filter(elev_class == "SUBALPINE") %>% drop_na(DIAC_dens)

df_shift <- df_base %>% drop_na(DIAC_shift)
df_dens <- df_base %>% drop_na(DIAC_dens)

# ==============================================================================
# H1: WARMING × ELEVATION INTERACTION
# ==============================================================================

test_H1 <- function(df, response_var, family = gaussian()) {
  
  df_clean <- df %>%
    drop_na(DEM_elev_mn, CLIM_tmax_diff_abs, !!sym(response_var)) %>%
    mutate(
      y = .data[[response_var]],
      warming = scale(CLIM_tmax_diff_abs)[,1],
      elev = scale(DEM_elev_mn)[,1]
    )
  
  # M1: Main effects only
  m1_main <- mgcv::gam(y ~ warming + elev, 
                       data = df_clean, 
                       family = family, 
                       method = "REML")
  
  # M2: Linear interaction
  m2_linear <- mgcv::gam(y ~ warming + elev + warming:elev, 
                         data = df_clean, 
                         family = family, 
                         method = "REML")
  
  # M3: Non-linear interaction
  m3_nonlinear <- mgcv::gam(y ~ s(warming, k = 6) + s(elev, k = 6) + 
                              ti(warming, elev, k = 6),
                            data = df_clean, 
                            family = family, 
                            method = "REML")
  
  # Tests
  anova_linear <- anova(m1_main, m2_linear, test = "Chisq")
  p_linear <- anova_linear$`Pr(>Chi)`[2]
  
  anova_nonlinear <- anova(m1_main, m3_nonlinear, test = "Chisq")
  p_nonlinear <- anova_nonlinear$`Pr(>Chi)`[2]
  
  # AIC comparison
  aic_table <- tibble(
    Model = c("M1: Main effects only",
              "M2: + Linear interaction",
              "M3: + Non-linear interaction"),
    AIC = c(AIC(m1_main), AIC(m2_linear), AIC(m3_nonlinear)),
    dAIC = c(AIC(m1_main), AIC(m2_linear), AIC(m3_nonlinear)) - 
      min(c(AIC(m1_main), AIC(m2_linear), AIC(m3_nonlinear)))
  ) %>%
    arrange(AIC)
  
  # Determine interaction type
  if (p_linear < 0.05) {
    h1_result <- "linear_interaction"
  } else if (p_nonlinear < 0.05) {
    h1_result <- "nonlinear_interaction"
  } else {
    h1_result <- "no_interaction"
  }
  
  list(
    h1_supported = (p_linear < 0.05 | p_nonlinear < 0.05),
    interaction_type = h1_result,
    p_linear = p_linear,
    p_nonlinear = p_nonlinear,
    aic_comparison = aic_table,
    best_model = aic_table$Model[1],
    models = list(
      m1_main = m1_main,
      m2_linear = m2_linear,
      m3_nonlinear = m3_nonlinear
    )
  )
}

# Run H1 tests
h1_shift <- test_H1(df_shift, "DIAC_shift", family = gaussian())
h1_dens  <- test_H1(df_dens,  "DIAC_dens",  family = gaussian())

# Summary
h1_summary <- tibble(
  Variable = c("Elevational shift", "Forest infilling"),
  `H1 supported` = c(h1_shift$h1_supported, h1_dens$h1_supported),
  `Interaction type` = c(h1_shift$interaction_type, h1_dens$interaction_type),
  `p (linear)` = c(h1_shift$p_linear, h1_dens$p_linear),
  `p (non-linear)` = c(h1_shift$p_nonlinear, h1_dens$p_nonlinear),
  `Best model` = c(h1_shift$best_model, h1_dens$best_model)
)

write.csv(h1_summary, "results/H1_summary.csv", row.names = FALSE)

# ==============================================================================
# FIGURE H1: BOXPLOTS BY WARMING CLASS
# ==============================================================================

# Create warming classes
df_base_classes <- df_base %>%
  drop_na(CLIM_tmax_diff_abs, DEM_elev_mn) %>%
  mutate(
    warming_class_raw = cut(
      CLIM_tmax_diff_abs,
      breaks = 10,
      labels = FALSE,
      include.lowest = TRUE
    ),
    warming_class = ifelse(warming_class_raw == 10, 9, warming_class_raw)
  )

# Calculate mean elevation per warming class
df_colors <- df_base_classes %>%
  group_by(warming_class) %>%
  summarise(
    warming_mean = mean(CLIM_tmax_diff_abs, na.rm = TRUE),
    elev_mean = mean(DEM_elev_mn, na.rm = TRUE),
    n_total = n(),
    .groups = "drop"
  )

# elev_mid <- mean(range(df_colors$elev_mean, na.rm = TRUE))
elev_mid <- 1942

df_shift <- df_base_classes %>%
  drop_na(DIAC_shift) %>%
  left_join(df_colors, by = "warming_class")

df_dens <- df_base_classes %>%
  drop_na(DIAC_dens) %>%
  left_join(df_colors, by = "warming_class")

# Common extent for color scale
extent_global <- c(
  min(c(df_shift$elev_mean, df_dens$elev_mean)), 
  max(c(df_shift$elev_mean, df_dens$elev_mean))
)

pal_brokenblue_orange <- divergingx_hcl(
  n = 7,
  h1 = 210,  
  h2 = 60,   
  c1 = 65,   
  c2 = 80,   
  l1 = 55,  
  l2 = 70    
)

swatchplot(pal_brokenblue_orange)
pal_grad <- colorRampPalette(pal_brokenblue_orange)

P1 <- ggplot(df_shift, aes(x = factor(warming_class), y = DIAC_shift, fill = elev_mean)) +
  geom_boxplot(outlier.alpha = 0.4, outlier.size = 1.2, linewidth = 0.5) +
  annotate("text", x = 0.7, y = max(df_shift$DIAC_shift, na.rm = TRUE) * 0.95,
           label = "A", size = 6, fontface = "bold", hjust = 0, vjust = 1) +
  
  scale_fill_gradientn(
    colours = pal_grad(200),
    values = scales::rescale(c(extent_global[1], elev_mid, extent_global[2])),
    limits = extent_global,
    name = "Mean elevation (m)",
    labels = scales::comma
  ) +
  
  scale_x_discrete(labels = NULL) +
  labs(x = "", y = "Elevational shift (m)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.ticks.x = element_blank()
  )

P2 <- ggplot(df_dens, aes(x = factor(warming_class), y = DIAC_dens, fill = elev_mean)) +
  geom_boxplot(outlier.alpha = 0.4, outlier.size = 1.2, linewidth = 0.5) +
  annotate("text", x = 0.7, y = max(df_dens$DIAC_dens, na.rm = TRUE) * 0.95,
           label = "B", size = 6, fontface = "bold", hjust = 0, vjust = 1) +
  
  scale_fill_gradientn(
    colours = pal_grad(200),
    values = scales::rescale(c(extent_global[1], elev_mid, extent_global[2])),
    limits = extent_global,
    name = "Mean elevation (m)",
    labels = scales::comma,
    guide = guide_colorbar(
      barwidth = 1.5, barheight = 12,
      title.position = "top", title.hjust = 0.5
    )
  ) +
  
  scale_x_discrete(labels = paste0(round(df_colors$warming_mean, 1), "°C")) +
  labs(x = "Warming magnitude (°C)", y = "Forest infilling (%)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold")
  )

fig_h1 <- (P1 / P2) + plot_layout(guides = "collect")

ggsave("figures/Figure_H1.png", fig_h1, width = 10, height = 10, dpi = 300, bg = "white")

# ==============================================================================
# H2: SPATIAL NON-STATIONARITY (GWR)
# ==============================================================================

dir.create("results/GWR", showWarnings = FALSE)

run_gwr_robust <- function(df, response_var, label) {
  
  df$y <- df[[response_var]]
  df <- df %>% drop_na(y)
  
  spdf <- st_as_sf(df, coords = c("X", "Y"), crs = 2154) %>% as("Spatial")
  
  form <- y ~ warming
  
  # Optimal bandwidth
  bw <- bw.gwr(form, data = spdf, approach = "AIC",
               kernel = CONFIG$kernel, adaptive = CONFIG$adaptive)
  
  # Robust GWR
  gwrfit <- gwr.robust(form, data = spdf, bw = bw,
                       kernel = CONFIG$kernel, adaptive = CONFIG$adaptive)
  
  sdf <- gwrfit$SDF
  
  result <- data.frame(
    X = st_coordinates(st_as_sf(sdf))[, 1],
    Y = st_coordinates(st_as_sf(sdf))[, 2],
    beta_warming = sdf@data$warming,
    se_warming = sdf@data$warming_SE,
    t_warming = sdf@data$warming_TV,
    label = label,
    bw = bw,
    stringsAsFactors = FALSE
  )
  
  out_sf <- st_as_sf(result, coords = c("X", "Y"), crs = 2154)
  st_write(out_sf, paste0("results/GWR/", gsub(" ", "_", label), ".gpkg"), 
           delete_layer = TRUE, quiet = TRUE)
  
  return(result)
}

# Run GWR for all scenarios
gwr_shift_m <- run_gwr_robust(df_shift_m, "DIAC_shift", "SHIFT_MONTANE")
gwr_shift_s <- run_gwr_robust(df_shift_s, "DIAC_shift", "SHIFT_SUBALPINE")
gwr_dens_m <- run_gwr_robust(df_dens_m, "DIAC_dens", "DENS_MONTANE")
gwr_dens_s <- run_gwr_robust(df_dens_s, "DIAC_dens", "DENS_SUBALPINE")

# Calculate significance
crit_value <- qnorm(1 - CONFIG$alpha_sig/2)

gwr_shift_m$signif <- abs(gwr_shift_m$t_warming) > crit_value
gwr_shift_s$signif <- abs(gwr_shift_s$t_warming) > crit_value
gwr_dens_m$signif <- abs(gwr_dens_m$t_warming) > crit_value
gwr_dens_s$signif <- abs(gwr_dens_s$t_warming) > crit_value

# Summary table
summary_gwr <- data.frame(
  Scenario = c("SHIFT_MONTANE", "SHIFT_SUBALPINE", "DENS_MONTANE", "DENS_SUBALPINE"),
  N = c(nrow(gwr_shift_m), nrow(gwr_shift_s), nrow(gwr_dens_m), nrow(gwr_dens_s)),
  BW = c(unique(gwr_shift_m$bw), unique(gwr_shift_s$bw),
         unique(gwr_dens_m$bw), unique(gwr_dens_s$bw)),
  Beta_mean = round(c(
    mean(gwr_shift_m$beta_warming, na.rm=TRUE),
    mean(gwr_shift_s$beta_warming, na.rm=TRUE),
    mean(gwr_dens_m$beta_warming, na.rm=TRUE),
    mean(gwr_dens_s$beta_warming, na.rm=TRUE)
  ), 3),
  Beta_sd = round(c(
    sd(gwr_shift_m$beta_warming, na.rm=TRUE),
    sd(gwr_shift_s$beta_warming, na.rm=TRUE),
    sd(gwr_dens_m$beta_warming, na.rm=TRUE),
    sd(gwr_dens_s$beta_warming, na.rm=TRUE)
  ), 3),
  Perc_sig = round(c(
    100 * mean(gwr_shift_m$signif, na.rm=TRUE),
    100 * mean(gwr_shift_s$signif, na.rm=TRUE),
    100 * mean(gwr_dens_m$signif, na.rm=TRUE),
    100 * mean(gwr_dens_s$signif, na.rm=TRUE)
  ), 1)
)

write.csv(summary_gwr, "results/GWR/H2_summary.csv", row.names=FALSE)


# ==============================================================================
# FIGURE H2: SPATIAL PATTERNS OF GWR COEFFICIENTS
# ==============================================================================

all_results <- bind_rows(gwr_shift_m, gwr_shift_s, gwr_dens_m, gwr_dens_s)

extent <- list(
  xmin = min(all_results$X), xmax = max(all_results$X),
  ymin = min(all_results$Y), ymax = max(all_results$Y)
)

lim_shift <- 25
lim_dens <- 1.5

create_panel <- function(df, lims, show_x = TRUE, show_y = TRUE, 
                         show_legend = FALSE, panel_label = "") {
  
  pct_sig <- round(100 * mean(df$signif, na.rm=TRUE), 1)
  beta_mean <- round(mean(df$beta_warming, na.rm=TRUE), 3)
  
  df$X_km <- df$X / 1000
  df$Y_km <- df$Y / 1000
  
  neg <- sequential_hcl(100, h = 200, c = c(0,60), l = c(100,50))
  pos <- sequential_hcl(100, h = 50,  c = c(0,70), l = c(100,60))
  
  colors_lch <- c(
    rev(neg),     # teal → blanc
    pos           # blanc → orange-rouge
  )
  
  p <- ggplot(df, aes(X_km, Y_km)) +
    geom_point(data = filter(df, !signif), aes(color = beta_warming),
               size = 2.5, alpha = 0.6) +
    geom_point(data = filter(df, signif), aes(fill = beta_warming),
               size = 2.5, alpha = 0.9, shape = 21, color = "black", stroke = 0.8) +
    scale_color_gradientn(
      colours = colors_lch,
      values = scales::rescale(seq(-1, 1, length.out = length(colors_lch))),
      limits = c(-lims, lims),
      oob = scales::squish,
      name = expression(beta[warming]),
      guide = if(show_legend) "colorbar" else "none"
    ) +
    scale_fill_gradientn(
      colours = colors_lch,
      values = scales::rescale(seq(-1, 1, length.out = length(colors_lch))),
      limits = c(-lims, lims),
      oob = scales::squish,
      guide = "none"
    ) +
    coord_equal(xlim = c(extent$xmin, extent$xmax) / 1000,
                ylim = c(extent$ymin, extent$ymax) / 1000) +
    labs(subtitle = sprintf("β̄ = %s | %s%% with trend", beta_mean, pct_sig),
         x = if(show_x) "X (km)" else "", y = if(show_y) "Y (km)" else "") +
    theme_minimal(base_size = 10) +
    theme(
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = if(show_legend) "right" else "none",
      legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.3, "cm"),
      legend.title = element_text(size = 9), legend.text = element_text(size = 8),
      panel.grid = element_line(color = "gray90", linewidth = 0.3),
      axis.text.x = if(show_x) element_text(size = 8) else element_blank(),
      axis.text.y = if(show_y) element_text(size = 8) else element_blank(),
      axis.title.x = if(show_x) element_text(size = 9) else element_blank(),
      axis.title.y = if(show_y) element_text(size = 9) else element_blank()
    )
  
  if(panel_label != "") {
    p <- p + annotate("text", 
                      x = extent$xmax / 1000 - (extent$xmax - extent$xmin) / 1000 * 0.05,
                      y = extent$ymax / 1000 - (extent$ymax - extent$ymin) / 1000 * 0.05,
                      label = panel_label, size = 6, fontface = "bold",
                      hjust = 1, vjust = 1)
  }
  
  return(p)
}

p1 <- create_panel(gwr_shift_m, lim_shift, show_x = FALSE, show_y = TRUE, 
                   show_legend = FALSE, panel_label = "A")
p2 <- create_panel(gwr_shift_s, lim_shift, show_x = FALSE, show_y = FALSE, 
                   show_legend = TRUE, panel_label = "B")
p3 <- create_panel(gwr_dens_m, lim_dens, show_x = TRUE, show_y = TRUE, 
                   show_legend = FALSE, panel_label = "C")
p4 <- create_panel(gwr_dens_s, lim_dens, show_x = TRUE, show_y = FALSE, 
                   show_legend = TRUE, panel_label = "D")

top_montane <- textGrob("Montane", gp = gpar(fontsize = 14, fontface = "bold"))
top_alpine <- textGrob("Alpine", gp = gpar(fontsize = 14, fontface = "bold"))
left_shift <- textGrob("Elevational shift", rot = 90, 
                       gp = gpar(fontsize = 14, fontface = "bold"))
left_infill <- textGrob("Forest infilling", rot = 90, 
                        gp = gpar(fontsize = 14, fontface = "bold"))
vertical_line <- linesGrob(x = unit(c(0.5, 0.5), "npc"), y = unit(c(0, 1), "npc"),
                           gp = gpar(col = "gray70", lwd = 0.5))
bottom_caption <- textGrob(
  "Circled points indicate locations with detectable trends (|t-value| > 1.65, α = 0.10)",
  gp = gpar(fontsize = 10, col = "grey50")
)

fig_h2 <- grid.arrange(
  top_montane, top_alpine,
  left_shift, p1, vertical_line, p2,
  left_infill, p3, vertical_line, p4,
  bottom_caption,
  layout_matrix = rbind(
    c(NA,  1, 1,  NA,  2, 2),
    c( 3,  4, 4,   5,  6, 6),
    c( 7,  8, 8,   9, 10, 10),
    c(NA, 11, 11, 11, 11, 11)
  ),
  widths = c(0.25, 5, 0.25, 0.05, 5, 0.25),
  heights = c(0.25, 4.75, 4.75, 0.35)
)

ggsave("figures/Figure_H2.png", fig_h2, width = 12, height = 6, dpi = 300, bg = "white")


# ==============================================================================
# H3: LOCAL ENVIRONMENTAL FACTORS
# ==============================================================================

local_vars <- c(
  "DEM_slope_mn", "DEM_TPI_mn", "LITH_PSed", "LITH_PCry", "LITH_shannon",
  "DEM_TWI_mn", "CA_core_mn_15", "AGG_clumpy_15", "GP_grazed_wood",
  "FOREST_conif", "FOREST_decid", "FOREST_mixed"
)

df_base_h3 <- df_base %>%
  mutate(across(any_of(local_vars), ~scale(.x)[,1], .names = "{.col}_std"))

analyze_h3 <- function(gwr_results, df_base, label, response_type, elevation_zone) {
  
  df_analysis <- gwr_results %>%
    select(X, Y, t_warming, signif, beta_warming) %>%
    mutate(
      cluster = case_when(
        !signif ~ "NS",
        signif & beta_warming > 0 ~ "POSITIVE",
        signif & beta_warming < 0 ~ "NEGATIVE"
      ),
      cluster = factor(cluster, levels = c("NEGATIVE", "NS", "POSITIVE"))
    ) %>%
    inner_join(df_base %>% select(X, Y, ends_with("_std")), by = c("X", "Y")) %>%
    drop_na()
  
  if (nrow(df_analysis) < 50) return(NULL)
  
  local_vars_available <- names(df_analysis)[grepl("_std$", names(df_analysis))]
  
  # ANOVA
  anova_results <- map_df(local_vars_available, function(var) {
    formula_str <- paste(var, "~ cluster")
    anova_res <- aov(as.formula(formula_str), data = df_analysis)
    anova_summary <- summary(anova_res)[[1]]
    
    tibble(
      Variable = var,
      F_statistic = anova_summary$`F value`[1],
      p_value = anova_summary$`Pr(>F)`[1]
    )
  }) %>%
    mutate(signif = p_value < 0.05) %>%
    arrange(p_value)
  
  # Tukey HSD
  tukey_results_list <- NULL
  significant_vars <- anova_results %>% filter(signif) %>% pull(Variable)
  
  if (length(significant_vars) > 0) {
    tukey_results_list <- map(significant_vars, function(var) {
      formula_str <- paste(var, "~ cluster")
      anova_res <- aov(as.formula(formula_str), data = df_analysis)
      tukey_res <- TukeyHSD(anova_res)
      
      as.data.frame(tukey_res$cluster) %>%
        rownames_to_column("comparison") %>%
        as_tibble() %>%
        mutate(Variable = var, signif_tukey = `p adj` < 0.05) %>%
        select(Variable, comparison, diff, `p adj`, signif_tukey)
    })
    
    tukey_results <- bind_rows(tukey_results_list)
  } else {
    tukey_results <- NULL
  }
  
  write.csv(anova_results, paste0("results/H3_ANOVA_", gsub(" ", "_", label), ".csv"), 
            row.names = FALSE)
  if (!is.null(tukey_results)) {
    write.csv(tukey_results, paste0("results/H3_TUKEY_", gsub(" ", "_", label), ".csv"),
              row.names = FALSE)
  }
  
  list(
    data = df_analysis,
    anova = anova_results,
    tukey = tukey_results,
    label = label,
    response_type = response_type,
    elevation_zone = elevation_zone
  )
}

h3_shift_m <- analyze_h3(gwr_shift_m, df_base_h3 %>% filter(elev_class == "MONTANE"),
                         "SHIFT_MONTANE", "Elevational shift", "Montane treeline ecotones")
h3_shift_s <- analyze_h3(gwr_shift_s, df_base_h3 %>% filter(elev_class == "SUBALPINE"),
                         "SHIFT_SUBALPINE", "Elevational shift", "Alpine treeline ecotones")
h3_dens_m <- analyze_h3(gwr_dens_m, df_base_h3 %>% filter(elev_class == "MONTANE"),
                        "DENS_MONTANE", "Forest infilling", "Montane treeline ecotones")
h3_dens_s <- analyze_h3(gwr_dens_s, df_base_h3 %>% filter(elev_class == "SUBALPINE"),
                        "DENS_SUBALPINE", "Forest infilling", "Alpine treeline ecotones")

# ==============================================================================
# FIGURE H3: ENVIRONMENTAL PROFILES BY CLUSTER
# ==============================================================================

key_variables <- c("LITH_PCry", "DEM_slope_mn", "DEM_TWI_mn")

all_possible_vars <- c(
  "DEM_slope_mn", "DEM_TPI_mn", "LITH_PSed", "LITH_PCry", "LITH_shannon",
  "DEM_TWI_mn", "CA_core_mn_15", "AGG_clumpy_15", "GP_grazed_wood",
  "FOREST_conif", "FOREST_decid", "FOREST_mixed"
)

variable_labels <- c(
  "DEM_slope_mn"   = "Slope",
  "DEM_TPI_mn"     = "TPI",
  "LITH_PSed"      = "% Sedimentary",
  "LITH_PCry"      = "% Crystalline",
  "LITH_shannon"   = "Lithodiversity",
  "DEM_TWI_mn"     = "TWI",
  "CA_core_mn_15"  = "Mean F. patch",
  "AGG_clumpy_15"  = "Clumpiness F. patch",
  "GP_grazed_wood" = "Grazing",
  "FOREST_conif"   = "% Conifer",
  "FOREST_decid"   = "% Caducifolius",
  "FOREST_mixed"   = "% Mixed F."
)

set.seed(123)
base_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(all_possible_vars))
color_palette_global <- setNames(base_colors, all_possible_vars)
is_key_var <- setNames(all_possible_vars %in% key_variables, all_possible_vars)

create_boxplot_panel <- function(h3_result, max_vars = 12, add_title = TRUE) {
  
  if (is.null(h3_result)) return(NULL)
  
  top_vars <- h3_result$anova %>%
    filter(signif) %>%
    arrange(p_value) %>%
    slice_head(n = max_vars) %>%
    pull(Variable)
  
  if (length(top_vars) == 0) return(NULL)
  
  var_names <- str_remove(top_vars, "_std$")
  
  n_plots <- length(var_names)
  ncol <- min(4, ceiling(sqrt(n_plots)))
  
  data_plot <- h3_result$data %>%
    select(cluster, all_of(top_vars)) %>%
    pivot_longer(-cluster, names_to = "Variable", values_to = "Value") %>%
    mutate(
      Variable = str_remove(Variable, "_std$"),
      Variable = factor(Variable, levels = var_names),
      cluster = factor(cluster, levels = c("NEGATIVE", "NS", "POSITIVE"),
                       labels = c("NEG", "NS", "POS"))
    )
  
  y_min_global <- min(data_plot$Value, na.rm = TRUE)
  y_max_global <- max(data_plot$Value, na.rm = TRUE)
  y_range_global <- y_max_global - y_min_global
  y_limits <- c(y_min_global - 0.05 * y_range_global, 
                y_max_global + 0.3 * y_range_global)
  
  colors_vars <- color_palette_global[var_names]
  alpha_vars <- ifelse(is_key_var[var_names], 0.9, 0.4)
  
  tukey_annotations <- NULL
  if (!is.null(h3_result$tukey)) {
    tukey_annotations <- h3_result$tukey %>%
      filter(signif_tukey, Variable %in% top_vars) %>%
      mutate(
        Variable = str_remove(Variable, "_std$"),
        group1 = case_when(
          str_detect(comparison, "^NEGATIVE") ~ "NEG",
          str_detect(comparison, "^NS") ~ "NS",
          str_detect(comparison, "^POSITIVE") ~ "POS",
          TRUE ~ str_extract(comparison, "^[^-]+")
        ),
        group2 = case_when(
          str_detect(comparison, "NEGATIVE$") ~ "NEG",
          str_detect(comparison, "NS$") ~ "NS",
          str_detect(comparison, "POSITIVE$") ~ "POS",
          TRUE ~ str_extract(comparison, "[^-]+$")
        ),
        label = case_when(
          `p adj` < 0.001 ~ "***",
          `p adj` < 0.01 ~ "**",
          `p adj` < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
  }
  
  plot_list <- list()
  
  for (i in seq_along(var_names)) {
    var <- var_names[i]
    data_var <- data_plot %>% filter(Variable == var)
    
    is_key <- is_key_var[var]
    current_alpha <- alpha_vars[i]
    
    show_ylabel_this <- ((i - 1) %% ncol == 0)
    
    title_size <- if (is_key) 18 else 16
    title_face <- if (is_key) "bold" else "plain"
    
    p_var <- ggplot(data_var, aes(x = cluster, y = Value, fill = Variable)) +
      geom_boxplot(alpha = current_alpha, outlier.size = 1, linewidth = 0.7,
                   show.legend = FALSE) +
      scale_fill_manual(values = colors_vars) +
      coord_cartesian(ylim = y_limits) +
      labs(
        title = variable_labels[var],
        x = NULL,
        y = if (show_ylabel_this) "Standardized\nvalue" else NULL
      ) +
      theme_minimal(base_size = 18) +
      theme(
        plot.title = element_text(face = title_face, size = title_size, hjust = 0.5),
        
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        
        axis.title.y = if (show_ylabel_this) {
          element_text(size = 18, face = "bold", lineheight = 0.9)
        } else {
          element_blank()
        },
        
        panel.grid.major = element_line(linewidth = 0.5),
        panel.grid.minor = element_line(linewidth = 0.3)
      )
    
    if (!is.null(tukey_annotations)) {
      tukey_var <- tukey_annotations %>% filter(Variable == var)
      
      if (nrow(tukey_var) > 0) {
        comparisons_list <- tukey_var %>%
          rowwise() %>%
          mutate(comp = list(c(group1, group2))) %>%
          pull(comp)
        
        annotations <- tukey_var$label
        
        p_var <- p_var + 
          geom_signif(
            comparisons = comparisons_list,
            annotations = annotations,
            y_position = y_max_global + y_range_global * seq(
              0.1, 0.1 + 0.15 * (nrow(tukey_var) - 1),
              length.out = nrow(tukey_var)
            ),
            tip_length = 0.01,
            textsize = 6,
            vjust = 0.5
          )
      }
    }
    
    plot_list[[var]] <- p_var
  }
  
  plot_title <- paste(h3_result$response_type, "-", h3_result$elevation_zone)
  
  if (add_title) {
    p_combined <- wrap_plots(plot_list, ncol = ncol) +
      plot_annotation(
        title = plot_title,
        theme = theme(plot.title = element_text(face = "bold", size = 28, hjust = 0.5))
      )
  } else {
    p_combined <- wrap_plots(plot_list, ncol = ncol)
  }
  
  list(plot = p_combined, title = plot_title, n_vars = n_plots)
}





result_shift_m <- create_boxplot_panel(h3_shift_m, max_vars = 12, add_title = TRUE)
result_shift_s <- create_boxplot_panel(h3_shift_s, max_vars = 12, add_title = TRUE)
result_dens_m <- create_boxplot_panel(h3_dens_m, max_vars = 12, add_title = TRUE)
result_dens_s <- create_boxplot_panel(h3_dens_s, max_vars = 12, add_title = TRUE)

if (!is.null(result_shift_m)) {
  ggsave("figures/H3_SHIFT_MONTANE.png", result_shift_m$plot, 
         width = 16, height = 12, dpi = 300)
}
if (!is.null(result_shift_s)) {
  ggsave("figures/H3_SHIFT_SUBALPINE.png", result_shift_s$plot, 
         width = 16, height = 12, dpi = 300)
}
if (!is.null(result_dens_m)) {
  ggsave("figures/H3_DENS_MONTANE.png", result_dens_m$plot, 
         width = 16, height = 12, dpi = 300)
}
if (!is.null(result_dens_s)) {
  ggsave("figures/H3_DENS_SUBALPINE.png", result_dens_s$plot, 
         width = 16, height = 12, dpi = 300)
}

# Combined figure
if (!is.null(result_shift_m) && !is.null(result_shift_s) && 
    !is.null(result_dens_m) && !is.null(result_dens_s)) {
  
  n_vars_shift_m <- result_shift_m$n_vars
  n_vars_shift_s <- result_shift_s$n_vars
  n_vars_dens_m <- result_dens_m$n_vars
  n_vars_dens_s <- result_dens_s$n_vars
  
  result_shift_m_notitle <- create_boxplot_panel(h3_shift_m, max_vars = 12, add_title = FALSE)
  result_shift_s_notitle <- create_boxplot_panel(h3_shift_s, max_vars = 12, add_title = FALSE)
  result_dens_m_notitle <- create_boxplot_panel(h3_dens_m, max_vars = 12, add_title = FALSE)
  result_dens_s_notitle <- create_boxplot_panel(h3_dens_s, max_vars = 12, add_title = FALSE)
  
  P_shift_m_with_title <- wrap_elements(
    grid::textGrob(paste0("(A) ", result_shift_m$title), 
                   gp = grid::gpar(fontface = "bold", fontsize = 25))
  ) / result_shift_m_notitle$plot + plot_layout(heights = c(0.05, 0.95))
  
  P_shift_s_with_title <- wrap_elements(
    grid::textGrob(paste0("(B) ", result_shift_s$title), 
                   gp = grid::gpar(fontface = "bold", fontsize = 25))
  ) / result_shift_s_notitle$plot + plot_layout(heights = c(0.05, 0.95))
  
  P_dens_m_with_title <- wrap_elements(
    grid::textGrob(paste0("(C) ", result_dens_m$title), 
                   gp = grid::gpar(fontface = "bold", fontsize = 25))
  ) / result_dens_m_notitle$plot + plot_layout(heights = c(0.05, 0.95))
  
  P_dens_s_with_title <- wrap_elements(
    grid::textGrob(paste0("(D) ", result_dens_s$title), 
                   gp = grid::gpar(fontface = "bold", fontsize = 25))
  ) / result_dens_s_notitle$plot + plot_layout(heights = c(0.05, 0.95))
  
  max_vars_top <- max(n_vars_shift_m, n_vars_shift_s)
  max_vars_bottom <- max(n_vars_dens_m, n_vars_dens_s)
  height_ratio <- c(max_vars_top, max_vars_bottom)
  
  P_combined <- 
    (
      (P_shift_m_with_title | P_shift_s_with_title) /
        plot_spacer() /
        (P_dens_m_with_title | P_dens_s_with_title)
    ) +
    plot_layout(
      ncol = 1,
      heights = c(height_ratio[1], 0.15, height_ratio[2])
    )
  
  
  height_per_var <- 1.5
  total_height <- (max_vars_top + max_vars_bottom) * height_per_var + 2
  
  ggsave("figures/Figure_H3.png", P_combined, 
         width = 24, height = total_height, dpi = 300)
}
