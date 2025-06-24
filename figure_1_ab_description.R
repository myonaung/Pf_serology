#!/usr/bin/env Rscript

library(plotmaps)
library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

# --- Set working directory ---
setwd("~/Documents/2024/Jasmin/MAP_output/")

# --- Load Coordinates ---
coord_file <- "India_Bangladesh_Myanmar.coord"
coords <- read.table(coord_file, header = FALSE)
colnames(coords) <- c("long", "lat")
cat("Longitude range:", range(coords$long), "\n")
cat("Latitude range:", range(coords$lat), "\n")

# --- Run plot_maps() to produce visual output ---
pdf("India_Bangladesh_Myanmar/MAPS_output.pdf", width = 10, height = 10)
plot_maps(
  add.pts       = TRUE,
  add.graph     = TRUE,
  add.countries = TRUE,
  longlat       = TRUE,
  mcmcpath      = "India_Bangladesh_Myanmar",
  outpath       = "India_Bangladesh_Myanmar"
)
dev.off()
cat("✅ plot_maps() completed and MAPS_output.pdf written\n")

# --- Generate mrates-mean.txt and Nsizes-mean.txt if missing ---
mrates_file <- "India_Bangladesh_Myanmar/mrates-mean.txt"
Nsizes_file <- "India_Bangladesh_Myanmar/Nsizes-mean.txt"

# Migration rates
if (!file.exists(mrates_file)) {
  if (file.exists("India_Bangladesh_Myanmar/mcmcmrates.txt")) {
    mrates_trace <- read.table("India_Bangladesh_Myanmar/mcmcmrates.txt", header = FALSE)
    mrates_mean <- colMeans(mrates_trace)
    write.table(mrates_mean, file = mrates_file, row.names = FALSE, col.names = FALSE)
    cat("✅ mrates-mean.txt computed from mcmcmrates.txt\n")
  } else {
    cat("❌ mcmcmrates.txt not found — cannot compute mrates-mean.txt\n")
  }
}

# Population sizes
if (!file.exists(Nsizes_file)) {
  if (file.exists("India_Bangladesh_Myanmar/mcmcthetas.txt")) {
    Nsizes_trace <- read.table("India_Bangladesh_Myanmar/mcmcthetas.txt", header = FALSE)
    Nsizes_mean <- colMeans(Nsizes_trace)
    write.table(Nsizes_mean, file = Nsizes_file, row.names = FALSE, col.names = FALSE)
    cat("✅ Nsizes-mean.txt computed from mcmcthetas.txt\n")
  } else {
    cat("❌ mcmcthetas.txt not found — cannot compute Nsizes-mean.txt\n")# Minimal Seroprevalence Comparison Analysis
    # For publication purposes
    
    # Load required libraries
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    
    # Set seed for reproducibility
    set.seed(1234)
    
    # ============================================================================
    # DATA LOADING AND PREPARATION
    # ============================================================================
    
    # Load data - replace with your actual file paths
    cohort1_data <- read.csv("cohort1_antibody_data.csv")
    cohort2_data <- read.csv("cohort2_antibody_data.csv")
    seroprevalence_data <- read.csv("seroprevalence_thresholds.csv")
    
    # Add cohort identifiers
    cohort1_data <- cohort1_data %>% mutate(cohort = "Cohort1")
    cohort2_data <- cohort2_data %>% mutate(cohort = "Cohort2")
    
    # Combine datasets
    combined_data <- rbind(cohort1_data, cohort2_data)
    
    # ============================================================================
    # DATA TRANSFORMATION
    # ============================================================================
    
    # Reshape to long format for antibody analysis
    # Adjust column range to match your antibody columns
    antibody_long <- combined_data %>%
      pivot_longer(
        cols = MSP1:Pfs25,  # Adjust this range to match your antibody columns
        names_to = "antigen", 
        values_to = "antibody_value"
      )
    
    # Merge with seroprevalence thresholds
    antibody_long_sero <- antibody_long %>%
      left_join(seroprevalence_data, by = c("cohort", "antigen", "Timepoint"))
    
    # Set factor levels for consistent ordering
    antibody_long_sero$cohort <- factor(antibody_long_sero$cohort, 
                                        levels = c("Cohort1", "Cohort2"))
    
    # Order antigens by median antibody value in reference cohort
    antigen_order <- antibody_long_sero %>%
      filter(cohort == "Cohort1") %>%  # Use first cohort as reference
      group_by(antigen) %>%
      summarise(median_value = median(antibody_value, na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(median_value)) %>%
      pull(antigen)
    
    antibody_long_sero$antigen <- factor(antibody_long_sero$antigen, levels = antigen_order)
    
    # ============================================================================
    # FIGURE 1A: ANTIBODY LEVEL COMPARISON
    # ============================================================================
    
    create_antibody_plot <- function(data) {
      ggplot(data, aes(x = Timepoint, y = antibody_value, fill = cohort)) +
        geom_violin(position = position_dodge(width = 0.8), alpha = 0.8, trim = TRUE) +
        geom_boxplot(width = 0.4, position = position_dodge(width = 0.8), 
                     outlier.shape = NA, alpha = 0.6) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), 
                    alpha = 0.1, size = 0.15, color = "black") +
        facet_wrap(~ antigen, scales = "free_y") +
        scale_y_log10(labels = scales::number_format(accuracy = 0.01)) +
        scale_fill_manual(values = c("Cohort1" = "#1f78b4", "Cohort2" = "#33a02c")) +
        labs(
          x = "Timepoint", 
          y = "Log10 Relative Antibody Unit", 
          fill = "Cohort"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          strip.background = element_rect(color = "gray", size = 0.25),
          strip.text = element_text(color = "black", face = "bold")
        )
    }
    
    # Create Figure 1A
    fig1A <- create_antibody_plot(antibody_long_sero)
    
    # Save Figure 1A
    ggsave("Figure1A_antibody_levels.png", plot = fig1A, width = 12, height = 8, dpi = 300)
    ggsave("Figure1A_antibody_levels.pdf", plot = fig1A, width = 12, height = 8)
    
    # ============================================================================
    # STATISTICAL TESTING
    # ============================================================================
    
    # Mann-Whitney U test for baseline vs endpoint comparisons
    perform_statistical_tests <- function(data) {
      data %>%
        filter(Timepoint %in% c("Baseline", "Endpoint")) %>%
        group_by(cohort, antigen) %>%
        summarise(
          n_baseline = sum(Timepoint == "Baseline" & !is.na(antibody_value)),
          n_endpoint = sum(Timepoint == "Endpoint" & !is.na(antibody_value)),
          median_baseline = median(antibody_value[Timepoint == "Baseline"], na.rm = TRUE),
          median_endpoint = median(antibody_value[Timepoint == "Endpoint"], na.rm = TRUE),
          p_value = ifelse(
            n_baseline > 1 & n_endpoint > 1,
            wilcox.test(
              antibody_value[Timepoint == "Baseline"], 
              antibody_value[Timepoint == "Endpoint"], 
              alternative = "two.sided"
            )$p.value, 
            NA
          ),
          .groups = "drop"
        ) %>%
        mutate(
          significant = ifelse(is.na(p_value), FALSE, p_value < 0.05),
          change_direction = case_when(
            is.na(p_value) ~ "No data",
            median_endpoint > median_baseline ~ "Increase",
            median_endpoint < median_baseline ~ "Decrease",
            TRUE ~ "No change"
          )
        )
    }
    
    # Run statistical tests
    statistical_results <- perform_statistical_tests(antibody_long_sero)
    
    # Save statistical results
    write.csv(statistical_results, "statistical_test_results.csv", row.names = FALSE)
    
    # Print summary
    cat("Statistical Testing Summary:\n")
    print(table(statistical_results$cohort, statistical_results$significant))
    
    # ============================================================================
    # FIGURE 1B: SEROPREVALENCE CHANGE PLOT
    # ============================================================================
    
    # Calculate seroprevalence summary
    calculate_seroprevalence <- function(data) {
      data %>%
        group_by(cohort, Timepoint, antigen) %>%
        summarize(
          mean_seroprevalence = mean(seroprevalence, na.rm = TRUE),
          n_samples = n(),
          .groups = 'drop'
        )
    }
    
    seroprevalence_summary <- calculate_seroprevalence(antibody_long_sero)
    
    # Prepare data for change plot
    seroprevalence_change <- seroprevalence_summary %>%
      pivot_wider(names_from = Timepoint, values_from = mean_seroprevalence) %>%
      filter(!is.na(Baseline) & !is.na(Endpoint))
    
    # Get antigen order from reference cohort
    reference_order <- seroprevalence_change %>%
      filter(cohort == "Cohort1") %>%
      arrange(Baseline) %>%
      pull(antigen)
    
    # Create change plots for each cohort
    create_change_plot <- function(data, cohort_name, antigen_order) {
      cohort_data <- data %>%
        filter(cohort == cohort_name) %>%
        mutate(antigen = factor(antigen, levels = antigen_order))
      
      ggplot(cohort_data, aes(x = Baseline, xend = Endpoint, y = antigen)) +
        geom_segment(
          aes(xend = Endpoint),
          arrow = arrow(length = unit(0.18, "cm"), type = "closed", angle = 30), 
          linewidth = 0.5,
          color = "black"
        ) +
        geom_point(aes(x = Baseline), size = 2, color = "#1f78b4", alpha = 0.7) +
        geom_point(aes(x = Endpoint), size = 2, color = "#e31a1c", alpha = 0.7) +
        scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
        labs(
          title = paste("Cohort:", cohort_name), 
          x = "Seroprevalence", 
          y = "Antigen"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
          panel.grid.major.x = element_line(color = "gray90", size = 0.3),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
        )
    }
    
    # Create plots for each cohort
    cohort_plots <- list()
    for (cohort_name in unique(seroprevalence_change$cohort)) {
      cohort_plots[[cohort_name]] <- create_change_plot(
        seroprevalence_change, 
        cohort_name, 
        reference_order
      )
    }
    
    # Combine plots
    fig1B <- wrap_plots(cohort_plots, ncol = 2)
    
    # Save Figure 1B
    ggsave("Figure1B_seroprevalence_change.png", plot = fig1B, width = 12, height = 8, dpi = 300)
    ggsave("Figure1B_seroprevalence_change.pdf", plot = fig1B, width = 12, height = 8)
    
    # ============================================================================
    # SUMMARY STATISTICS
    # ============================================================================
    
    # Calculate summary statistics
    create_summary_stats <- function(sero_data, stat_data) {
      
      # Seroprevalence summary
      sero_summary <- sero_data %>%
        group_by(cohort, Timepoint) %>%
        summarise(
          mean_seroprevalence = round(mean(mean_seroprevalence, na.rm = TRUE), 3),
          median_seroprevalence = round(median(mean_seroprevalence, na.rm = TRUE), 3),
          n_antigens = n(),
          .groups = 'drop'
        )
      
      # Statistical summary
      stat_summary <- stat_data %>%
        group_by(cohort) %>%
        summarise(
          total_comparisons = n(),
          significant_increases = sum(significant & change_direction == "Increase", na.rm = TRUE),
          significant_decreases = sum(significant & change_direction == "Decrease", na.rm = TRUE),
          no_change = sum(!significant, na.rm = TRUE),
          .groups = 'drop'
        )
      
      return(list(
        seroprevalence = sero_summary,
        statistics = stat_summary
      ))
    }
    
    summary_stats <- create_summary_stats(seroprevalence_summary, statistical_results)
    
    # Print summaries
    cat("\n=== SEROPREVALENCE SUMMARY ===\n")
    print(summary_stats$seroprevalence)
    
    cat("\n=== STATISTICAL CHANGES SUMMARY ===\n")
    print(summary_stats$statistics)
    
    # Save summaries
    write.csv(summary_stats$seroprevalence, "seroprevalence_summary.csv", row.names = FALSE)
    write.csv(summary_stats$statistics, "statistical_summary.csv", row.names = FALSE)
    
    # ============================================================================
    # FINAL OUTPUTS
    # ============================================================================
    
    cat("\n=== ANALYSIS COMPLETE ===\n")
    cat("Generated figures:\n")
    cat("- Figure1A_antibody_levels.png/pdf\n")
    cat("- Figure1B_seroprevalence_change.png/pdf\n")
    cat("\nGenerated data files:\n")
    cat("- statistical_test_results.csv\n")
    cat("- seroprevalence_summary.csv\n")
    cat("- statistical_summary.csv\n")
    
    # Display plots
    print("Figure 1A: Antibody Level Comparison")
    print(fig1A)
    
    print("Figure 1B: Seroprevalence Change")
    print(fig1B)
  }
}

# --- Load country outlines ---
world <- ne_countries(scale = "medium", returnclass = "sf")

# --- Load Grid Positions ---
grid_file <- "India_Bangladesh_Myanmar.grid"
grid <- read.table(grid_file, header = FALSE)
colnames(grid) <- c("deme", "x", "y")

# --- Define bounding box with padding ---
x_buffer <- 2
y_buffer <- 2
xlim_custom <- range(grid$x) + c(-x_buffer, x_buffer)
ylim_custom <- range(grid$y) + c(-y_buffer, y_buffer)

# --- Plot Migration Rates ---
if (file.exists(mrates_file)) {
  mrates <- read.table(mrates_file, header = FALSE)
  names(mrates) <- "rate"
  mrates <- cbind(grid, rate = mrates$rate)
  
  p_mrates <- ggplot() +
    geom_sf(data = world, fill = "gray95", color = "gray70", size = 0.3) +
    geom_point(data = mrates, aes(x = x, y = y, color = rate), size = 3) +
    scale_color_viridis_c(option = "plasma", trans = "log10", name = "Migration rate") +
    coord_sf(xlim = xlim_custom, ylim = ylim_custom, expand = FALSE) +
    theme_minimal(base_size = 14) +
    labs(title = "MAPS Inferred Migration Rates", x = "Longitude", y = "Latitude")
  
  ggsave("India_Bangladesh_Myanmar/mrates-mean-ggplot.pdf", plot = p_mrates, width = 10, height = 6)
  cat("✅ Enhanced migration plot saved: mrates-mean-ggplot.pdf\n")
} else {
  cat("❌ Migration rate file not found: mrates-mean.txt\n")
}

# --- Plot Population Sizes ---
if (file.exists(Nsizes_file)) {
  Nsizes <- read.table(Nsizes_file, header = FALSE)
  names(Nsizes) <- "size"
  Nsizes <- cbind(grid, size = Nsizes$size)
  
  p_Nsizes <- ggplot() +
    geom_sf(data = world, fill = "gray95", color = "gray70", size = 0.3) +
    geom_point(data = Nsizes, aes(x = x, y = y, color = size), size = 3) +
    scale_color_viridis_c(option = "viridis", name = "Population size") +
    coord_sf(xlim = xlim_custom, ylim = ylim_custom, expand = FALSE) +
    theme_minimal(base_size = 14) +
    labs(title = "MAPS Inferred Population Sizes", x = "Longitude", y = "Latitude")
  
  ggsave("India_Bangladesh_Myanmar/Nsizes-mean-ggplot.pdf", plot = p_Nsizes, width = 10, height = 6)
  cat("✅ Enhanced population size plot saved: Nsizes-mean-ggplot.pdf\n")
} else {
  cat("❌ Population size file not found: Nsizes-mean.txt\n")
}
