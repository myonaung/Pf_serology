# Minimal Variance Analysis for Antibody Responses
# For publication purposes

# Load required libraries
library(car)
library(relaimpo)
library(ggplot2)
library(dplyr)
library(reshape2)

# Set seed for reproducibility
set.seed(1234)

# Load data
# Replace with your actual file paths
infection_data <- read.csv("infection_data.csv")
antibody_data <- read.csv("antibody_data.csv")

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Merge infection and antibody data
analysis_data <- merge(infection_data, antibody_data, by = "sampleid")

# Create categorical variables
analysis_data <- analysis_data %>%
  mutate(
    # Infection timing categories
    infection_timing = case_when(
      time_since_last_infection <= 30 ~ "Very_Recent",
      time_since_last_infection <= 90 ~ "Recent", 
      time_since_last_infection <= 180 ~ "Moderate",
      TRUE ~ "Distant"
    ),
    
    # Infection frequency categories
    infection_frequency = case_when(
      number_of_infections == 1 ~ "Single",
      number_of_infections <= 3 ~ "Few",
      TRUE ~ "Many"
    ),
    
    # Parasitemia categories
    parasitemia_category = case_when(
      max_parasitemia == 0 ~ "None",
      max_parasitemia <= 500 ~ "Low",
      max_parasitemia <= 5000 ~ "Medium",
      TRUE ~ "High"
    )
  )

# Convert to ordered factors
analysis_data$infection_timing <- factor(analysis_data$infection_timing,
                                         levels = c("Very_Recent", "Recent", "Moderate", "Distant"), 
                                         ordered = TRUE)

analysis_data$infection_frequency <- factor(analysis_data$infection_frequency,
                                            levels = c("Single", "Few", "Many"), 
                                            ordered = TRUE)

analysis_data$parasitemia_category <- factor(analysis_data$parasitemia_category,
                                             levels = c("None", "Low", "Medium", "High"), 
                                             ordered = TRUE)

# ============================================================================
# VARIANCE ANALYSIS FUNCTION
# ============================================================================

run_variance_analysis <- function(data, infection_var, model_name) {
  
  cat(paste("Running analysis:", model_name, "\n"))
  
  # Identify antibody columns (adjust this to match your data structure)
  antibody_columns <- names(data)[grepl("^(MSP|AMA|EBA|RH|LSA|TRAP|PF)", names(data))]
  
  cat(paste("Found", length(antibody_columns), "antibody columns\n"))
  
  # Initialize results
  results <- data.frame(
    antibody_response = character(),
    age_contribution = numeric(),
    infection_contribution = numeric(),
    parasitemia_contribution = numeric(),
    coinfection_contribution = numeric(),
    total_r2 = numeric(),
    model_type = character(),
    stringsAsFactors = FALSE
  )
  
  # Analyze each antibody
  for (antibody in antibody_columns) {
    
    # Prepare model data
    model_data <- data %>%
      select(all_of(c(antibody, "age", infection_var, "parasitemia_category", "number_of_coinfections"))) %>%
      filter(complete.cases(.))
    
    if (nrow(model_data) < 20) {
      cat(paste("Skipping", antibody, "- insufficient data\n"))
      next
    }
    
    tryCatch({
      # Build model formula
      formula_str <- paste(antibody, "~ log10(age) +", infection_var, 
                           "+ parasitemia_category + number_of_coinfections")
      model_formula <- as.formula(formula_str)
      
      # Fit linear model
      model <- lm(model_formula, data = model_data)
      
      # Calculate relative importance
      rel_imp <- calc.relimp(model, type = "lmg")
      
      # Store results
      new_row <- data.frame(
        antibody_response = antibody,
        age_contribution = rel_imp$lmg["log10(age)"],
        infection_contribution = rel_imp$lmg[infection_var],
        parasitemia_contribution = rel_imp$lmg["parasitemia_category"],
        coinfection_contribution = rel_imp$lmg["number_of_coinfections"],
        total_r2 = summary(model)$r.squared,
        model_type = model_name,
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, new_row)
      
    }, error = function(e) {
      cat(paste("Error with", antibody, ":", e$message, "\n"))
    })
  }
  
  return(results)
}

# ============================================================================
# RUN ANALYSES
# ============================================================================

# Option 1: Infection timing + parasitemia
timing_results <- run_variance_analysis(analysis_data, "infection_timing", "Timing_Model")

# Option 2: Infection frequency + parasitemia  
frequency_results <- run_variance_analysis(analysis_data, "infection_frequency", "Frequency_Model")

# ============================================================================
# VISUALIZATION FUNCTION
# ============================================================================

create_variance_plot <- function(results_df, model_name) {
  
  if (nrow(results_df) == 0) {
    warning(paste("No data for", model_name))
    return(NULL)
  }
  
  # Prepare data for plotting
  plot_data <- results_df %>%
    select(antibody_response, age_contribution, infection_contribution, 
           parasitemia_contribution, coinfection_contribution, total_r2) %>%
    gather(key = "variable", value = "contribution", -antibody_response, -total_r2) %>%
    mutate(
      antibody_response = factor(antibody_response, 
                                 levels = results_df %>% 
                                   arrange(desc(total_r2)) %>% 
                                   pull(antibody_response))
    )
  
  # Create plot
  variance_plot <- ggplot(plot_data, aes(x = antibody_response, y = contribution, fill = variable)) +
    geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
    scale_fill_manual(
      values = c(
        "age_contribution" = "#1f77b4",
        "infection_contribution" = "#ff7f0e", 
        "parasitemia_contribution" = "#2ca02c",
        "coinfection_contribution" = "#d62728"
      ),
      labels = c(
        "age_contribution" = "Age",
        "infection_contribution" = "Infection Pattern",
        "parasitemia_contribution" = "Parasitemia Level", 
        "coinfection_contribution" = "Co-infections"
      )
    ) +
    labs(
      title = paste("Variance Explained:", model_name),
      x = "Antibody Response",
      y = "Proportion of Variance Explained",
      fill = "Variable"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  return(variance_plot)
}

# ============================================================================
# CREATE VISUALIZATIONS
# ============================================================================

# Generate plots
timing_plot <- create_variance_plot(timing_results, "Infection Timing Model")
frequency_plot <- create_variance_plot(frequency_results, "Infection Frequency Model")

# Save plots
if (!is.null(timing_plot)) {
  ggsave("timing_variance_plot.png", timing_plot, width = 12, height = 8, dpi = 300)
}

if (!is.null(frequency_plot)) {
  ggsave("frequency_variance_plot.png", frequency_plot, width = 12, height = 8, dpi = 300)
}

# ============================================================================
# MODEL COMPARISON
# ============================================================================

if (nrow(timing_results) > 0 && nrow(frequency_results) > 0) {
  
  # Combine results
  combined_results <- rbind(timing_results, frequency_results)
  
  # Create comparison plot
  comparison_plot <- ggplot(combined_results, aes(x = model_type, y = total_r2, fill = model_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Model Performance Comparison",
      x = "Model Type",
      y = "Total R² (Explained Variance)",
      caption = "Each point represents one antibody response"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  ggsave("model_comparison.png", comparison_plot, width = 10, height = 6, dpi = 300)
  
  # Print summary statistics
  cat("\n=== MODEL PERFORMANCE SUMMARY ===\n")
  performance_summary <- combined_results %>%
    group_by(model_type) %>%
    summarise(
      mean_r2 = round(mean(total_r2, na.rm = TRUE), 3),
      median_r2 = round(median(total_r2, na.rm = TRUE), 3),
      sd_r2 = round(sd(total_r2, na.rm = TRUE), 3),
      n_antibodies = n(),
      .groups = 'drop'
    )
  
  print(performance_summary)
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

# Save detailed results
write.csv(timing_results, "timing_model_results.csv", row.names = FALSE)
write.csv(frequency_results, "frequency_model_results.csv", row.names = FALSE)

if (exists("combined_results")) {
  write.csv(combined_results, "combined_variance_results.csv", row.names = FALSE)
}

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

# Function to summarize key findings
summarize_findings <- function(results, model_name) {
  if (nrow(results) == 0) return(NULL)
  
  cat(paste("\n=== SUMMARY FOR", model_name, "===\n"))
  
  # Overall model performance
  cat("Overall model performance:\n")
  cat(paste("  Mean R²:", round(mean(results$total_r2), 3), "\n"))
  cat(paste("  Range:", round(min(results$total_r2), 3), "-", round(max(results$total_r2), 3), "\n"))
  
  # Variable importance
  avg_contributions <- results %>%
    summarise(
      age = round(mean(age_contribution, na.rm = TRUE), 3),
      infection = round(mean(infection_contribution, na.rm = TRUE), 3),
      parasitemia = round(mean(parasitemia_contribution, na.rm = TRUE), 3),
      coinfection = round(mean(coinfection_contribution, na.rm = TRUE), 3)
    )
  
  cat("Average variable contributions:\n")
  cat(paste("  Age:", avg_contributions$age, "\n"))
  cat(paste("  Infection pattern:", avg_contributions$infection, "\n"))
  cat(paste("  Parasitemia:", avg_contributions$parasitemia, "\n"))
  cat(paste("  Co-infections:", avg_contributions$coinfection, "\n"))
  
  # Top antibodies by explained variance
  top_antibodies <- results %>%
    arrange(desc(total_r2)) %>%
    head(5) %>%
    select(antibody_response, total_r2)
  
  cat("Top 5 antibodies by explained variance:\n")
  print(top_antibodies)
}

# Generate summaries
summarize_findings(timing_results, "TIMING MODEL")
summarize_findings(frequency_results, "FREQUENCY MODEL")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files saved:\n")
cat("- timing_variance_plot.png\n")
cat("- frequency_variance_plot.png\n") 
cat("- model_comparison.png\n")
cat("- timing_model_results.csv\n")
cat("- frequency_model_results.csv\n")
cat("- combined_variance_results.csv\n")