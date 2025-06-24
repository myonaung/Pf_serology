# Minimal Composite Antigen Durability Score Analysis
# For publication purposes

# Load required libraries
library(brms)
library(ggplot2)
library(dplyr)
library(stringr)
library(posterior)
library(readr)
library(tibble)

# Set seed for reproducibility
set.seed(1234)

# ============================================================================
# DATA LOADING AND SETUP
# ============================================================================

# Load time scaling data
antibody_data <- read.csv("antibody_kinetics_data.csv")
time_sd <- sd(antibody_data$time_since_infection, na.rm = TRUE)

cat("Time SD for scaling:", round(time_sd, 3), "\n")

# Define model files - replace with your actual file paths
model_files <- list.files(pattern = "kinetics_model.*\\.rds$", full.names = TRUE)

cat("Found", length(model_files), "model files\n")

# ============================================================================
# ANTIGEN LABELING
# ============================================================================

# Create antigen label mapping (customize for your antigens)
create_antigen_labels <- function() {
  tibble(
    Antigen = c("MSP1", "MSP2", "MSP3", "AMA1", "EBA175", "RH5", 
                "LSA3", "TRAP", "MSP7", "MSP4", "Pfs25"),
    Antigen_label = c("MSP1", "MSP2", "MSP3", "AMA1", "EBA175", "RH5",
                      "LSA3", "TRAP", "MSP7", "MSP4", "Pfs25")
  )
}

antigen_labels <- create_antigen_labels()

# ============================================================================
# HALF-LIFE EXTRACTION FUNCTION
# ============================================================================

extract_half_life <- function(model_file, parameter, component_name, time_scaling) {
  
  tryCatch({
    # Load model
    model <- readRDS(model_file)
    
    # Extract posterior samples
    posterior <- as.data.frame(as_draws_df(model))
    
    # Check if parameter exists
    if (!(parameter %in% names(posterior))) {
      return(NULL)
    }
    
    # Calculate half-lives
    decay_rate <- posterior[[parameter]] / time_scaling
    half_life <- log(2) / decay_rate
    
    # Filter valid half-lives
    half_life_clean <- half_life[is.finite(half_life) & half_life > 0 & half_life < 10000]
    
    if (length(half_life_clean) < 10) {
      return(NULL)
    }
    
    # Calculate summary statistics
    median_hl <- median(half_life_clean)
    mean_hl <- mean(half_life_clean)
    ci <- quantile(half_life_clean, c(0.025, 0.975))
    
    # Extract antigen name from filename
    antigen_name <- model_file %>%
      basename() %>%
      str_remove("^kinetics_model_.*?_") %>%
      str_remove("_.*$") %>%
      str_remove("\\.rds$")
    
    return(tibble(
      Antigen = antigen_name,
      Component = component_name,
      Median_HL = round(median_hl, 1),
      Mean_HL = round(mean_hl, 1),
      CI_Lower = round(ci[1], 1),
      CI_Upper = round(ci[2], 1)
    ))
    
  }, error = function(e) {
    warning(paste("Error processing", basename(model_file), ":", e$message))
    return(NULL)
  })
}

# ============================================================================
# PROCESS ALL MODELS
# ============================================================================

# Define decay components to extract
decay_components <- list(
  "b_rl_Intercept" = "long_lived",
  "b_rs_Intercept" = "short_lived", 
  "b_ra_Intercept" = "intermediate"
)

# Process all model files
all_results <- list()

for (file in model_files) {
  cat("Processing:", basename(file), "\n")
  
  file_results <- list()
  
  for (param in names(decay_components)) {
    component_result <- extract_half_life(
      file, 
      param, 
      decay_components[[param]], 
      time_sd
    )
    
    if (!is.null(component_result)) {
      file_results[[decay_components[[param]]]] <- component_result
    }
  }
  
  if (length(file_results) > 0) {
    all_results[[basename(file)]] <- bind_rows(file_results)
  }
}

# Combine all results
half_life_summary <- bind_rows(all_results)

cat("Extracted half-lives for", length(unique(half_life_summary$Antigen)), "antigens\n")

# ============================================================================
# COMPOSITE SCORE CALCULATION
# ============================================================================

# Define component weights (customize based on your biological interpretation)
component_weights <- c(
  "short_lived" = 0.3,     # Weight for short-lived component
  "intermediate" = 0.1,    # Weight for intermediate component  
  "long_lived" = 0.6       # Weight for long-lived component
)

calculate_composite_scores <- function(hl_data, weights, labels) {
  
  # Normalize scores within each component
  normalized_scores <- hl_data %>%
    filter(Component %in% names(weights)) %>%
    group_by(Component) %>%
    mutate(
      normalized_score = (Median_HL - min(Median_HL, na.rm = TRUE)) / 
        (max(Median_HL, na.rm = TRUE) - min(Median_HL, na.rm = TRUE))
    ) %>%
    ungroup()
  
  # Calculate composite scores
  composite_scores <- normalized_scores %>%
    select(Antigen, Component, normalized_score) %>%
    left_join(labels, by = "Antigen") %>%
    mutate(
      Antigen_label = coalesce(Antigen_label, Antigen),
      weight = weights[Component]
    ) %>%
    group_by(Antigen, Antigen_label) %>%
    summarise(
      Composite_Score = sum(normalized_score * weight, na.rm = TRUE),
      n_components = n(),
      .groups = "drop"
    ) %>%
    filter(n_components >= 2) %>%  # Require at least 2 components
    arrange(Composite_Score) %>%
    mutate(Antigen_label = factor(Antigen_label, levels = unique(Antigen_label)))
  
  return(composite_scores)
}

# Calculate composite scores
composite_results <- calculate_composite_scores(
  half_life_summary, 
  component_weights, 
  antigen_labels
)

cat("Calculated composite scores for", nrow(composite_results), "antigens\n")

# ============================================================================
# VISUALIZATION
# ============================================================================

create_composite_plot <- function(score_data) {
  
  ggplot(score_data, aes(x = Composite_Score, y = Antigen_label)) +
    geom_col(
      fill = "#2166ac", 
      color = "white", 
      width = 0.7,
      alpha = 0.8
    ) +
    labs(
      title = "Antibody Durability Composite Scores",
      subtitle = "Higher scores indicate longer antibody persistence",
      x = "Normalized Durability Score",
      y = "Antigen"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "gray50")
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, max(score_data$Composite_Score) * 1.05)
    )
}

# Create plot
composite_plot <- create_composite_plot(composite_results)

# ============================================================================
# SAVE OUTPUTS
# ============================================================================

# Save plot
ggsave("composite_durability_scores.png", composite_plot, 
       width = 10, height = 8, dpi = 300)
ggsave("composite_durability_scores.pdf", composite_plot, 
       width = 10, height = 8)

# Save data
write.csv(composite_results, "composite_scores.csv", row.names = FALSE)
write.csv(half_life_summary, "half_life_summary.csv", row.names = FALSE)

# Save RDS for further analysis
saveRDS(composite_results, "composite_scores.rds")
saveRDS(half_life_summary, "half_life_summary.rds")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

# Print summary
cat("\n=== COMPOSITE SCORE SUMMARY ===\n")
cat("Total antigens analyzed:", nrow(composite_results), "\n")
cat("Score range:", round(min(composite_results$Composite_Score), 3), 
    "to", round(max(composite_results$Composite_Score), 3), "\n")

# Top and bottom performers
cat("\nTop 5 most durable antigens:\n")
top_antigens <- composite_results %>%
  arrange(desc(Composite_Score)) %>%
  head(5) %>%
  select(Antigen_label, Composite_Score)
print(top_antigens)

cat("\nTop 5 least durable antigens:\n")
bottom_antigens <- composite_results %>%
  arrange(Composite_Score) %>%
  head(5) %>%
  select(Antigen_label, Composite_Score)
print(bottom_antigens)

# Component summary
cat("\n=== COMPONENT WEIGHTS USED ===\n")
weight_summary <- data.frame(
  Component = names(component_weights),
  Weight = component_weights,
  Description = c("Short-lived decay", "Intermediate decay", "Long-lived decay")
)
print(weight_summary)

# Half-life ranges by component
cat("\n=== HALF-LIFE RANGES BY COMPONENT ===\n")
hl_ranges <- half_life_summary %>%
  group_by(Component) %>%
  summarise(
    n_antigens = n(),
    median_hl = round(median(Median_HL, na.rm = TRUE), 1),
    min_hl = round(min(Median_HL, na.rm = TRUE), 1),
    max_hl = round(max(Median_HL, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(hl_ranges)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Generated files:\n")
cat("- composite_durability_scores.png/pdf\n")
cat("- composite_scores.csv\n") 
cat("- half_life_summary.csv\n")
cat("- composite_scores.rds\n")
cat("- half_life_summary.rds\n")

# Display plot
print(composite_plot)