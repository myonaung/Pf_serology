# Minimal Bayesian Antibody Kinetics Model
# For publication purposes

# Load required libraries
library(dplyr)
library(tidyr)
library(brms)

# Set multicore processing
options(mc.cores = parallel::detectCores())

# Set seed for reproducibility
set.seed(1234)

# Load data
# Replace with your actual file paths
antibody_data <- read.csv("antibody_data.csv")
strain_data <- read.csv("strain_data.csv")  # Contains P_strain (haplotype prevalance)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Function to get strain-specific proportion
get_p_strain <- function(antigen, strain_lookup) {
  strain_row <- strain_lookup[strain_lookup$Antigens == antigen, ]
  if (nrow(strain_row) == 0) {
    warning(paste("No strain data found for", antigen, ". Using default P_strain = 0.02"))
    return(0.02)
  }
  return(strain_row$P_strain_fraction)
}

# Conservative priors for model stability
get_priors <- function() {
  c(
    set_prior("normal(0, 0.5)", nlpar = "Abg"),        # Background antibody
    set_prior("normal(0, 0.5)", nlpar = "A0"),         # Initial amplitude
    set_prior("exponential(10)", nlpar = "rl", lb = 0), # Long-lived decay
    set_prior("exponential(5)", nlpar = "rs", lb = 0),  # Short-lived decay
    set_prior("exponential(1)", nlpar = "ra", lb = 0),  # Atypical decay
    set_prior("normal(0, 0.5)", nlpar = "beta"),       # Boost amplitude
    set_prior("beta(2, 2)", nlpar = "rho", lb = 0, ub = 1), # Proportion parameter
    set_prior("normal(0, 0.25)", nlpar = "tau0"),      # Temporal shift
    set_prior("exponential(1)", class = "nu")          # Student-t degrees of freedom
  )
}

# ============================================================================
# MAIN MODEL FITTING FUNCTION
# ============================================================================

fit_kinetics_model <- function(antigen, data, strain_lookup) {
  cat("Fitting model for antigen:", antigen, "\n")
  
  # Check if antigen exists
  if (!antigen %in% colnames(data)) {
    warning(paste("Antigen", antigen, "not found in dataset"))
    return(NULL)
  }
  
  # Get strain proportion
  p_strain_value <- get_p_strain(antigen, strain_lookup)
  cat("Using P_strain =", p_strain_value, "\n")
  
  # Prepare data
  model_data <- data %>%
    pivot_longer(
      cols = all_of(antigen),
      names_to = "antibody",
      values_to = "observed_antibody"
    ) %>%
    filter(!is.na(observed_antibody)) %>%
    # Remove extreme outliers
    filter(
      observed_antibody >= quantile(observed_antibody, 0.01, na.rm = TRUE),
      observed_antibody <= quantile(observed_antibody, 0.99, na.rm = TRUE)
    ) %>%
    mutate(
      # Scale variables for numerical stability
      observed_antibody_scaled = scale(observed_antibody)[, 1],
      timeNumeric = scale(time_since_infection)[, 1],
      # Create infection categories
      infection_category = factor(case_when(
        number_of_infections == 0 ~ "No Infections",
        number_of_infections <= 2 ~ "1-2 Infections",
        number_of_infections > 2 ~ ">2 Infections"
      ), levels = c("No Infections", "1-2 Infections", ">2 Infections")),
      P_strain = p_strain_value
    ) %>%
    filter(!is.na(observed_antibody_scaled), !is.na(timeNumeric))
  
  cat("Sample size:", nrow(model_data), "\n")
  
  if (nrow(model_data) < 50) {
    warning(paste("Insufficient data for", antigen))
    return(NULL)
  }
  
  # Fit Bayesian model
  model <- tryCatch({
    brm(
      # Multi-component decay model with strain-specific scaling
      bf(
        observed_antibody_scaled ~ 
          P_strain * (
            Abg +  # Background level
              A0 * exp(-rl * (timeNumeric - tau0)) +  # Long-lived component
              beta * (  # Boost components
                (1 - rho) * (exp(-rs * (timeNumeric - tau0)) - exp(-ra * (timeNumeric - tau0))) / (ra - rs) +
                  rho * (exp(-rl * (timeNumeric - tau0)) - exp(-ra * (timeNumeric - tau0))) / (ra - rl)
              )
          ),
        # Model components with infection category effects
        Abg ~ 1 + infection_category + (1|individual_id),
        A0 ~ 1 + (1|individual_id),
        beta ~ 1 + infection_category,
        tau0 ~ 1,
        rho ~ 1,
        rl ~ 1,
        rs ~ 1,
        ra ~ 1,
        nl = TRUE
      ),
      data = model_data,
      family = student(),  # Robust to outliers
      prior = get_priors(),
      
      # Convergence settings
      control = list(
        adapt_delta = 0.95,
        max_treedepth = 12
      ),
      
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      
      # Initialization for stability
      init = function() {
        list(
          rho = 0.5,
          Abg = 0,
          A0 = 0,
          rl = 0.01,
          rs = 0.1,
          ra = 0.05,
          beta = 0,
          tau0 = 0
        )
      },
      
      silent = TRUE
    )
  }, error = function(e) {
    warning(paste("Model fitting failed for", antigen, ":", e$message))
    return(NULL)
  })
  
  if (is.null(model)) return(NULL)
  
  # Check convergence
  max_rhat <- max(rhat(model), na.rm = TRUE)
  min_ess <- min(neff_ratio(model), na.rm = TRUE)
  
  cat("Max R-hat:", round(max_rhat, 4), "\n")
  cat("Min ESS ratio:", round(min_ess, 4), "\n")
  
  # Convergence criteria
  converged <- (max_rhat <= 1.05 && min_ess >= 0.1)
  
  if (converged) {
    cat("✓ Model converged successfully!\n")
    
    # Save model
    saveRDS(model, paste0("kinetics_model_", antigen, ".rds"))
    
    # Calculate half-lives
    half_lives <- calculate_half_lives(model, data, antigen, p_strain_value)
    
    return(list(model = model, half_lives = half_lives, converged = TRUE))
  } else {
    cat("⚠ Poor convergence detected\n")
    return(list(model = model, half_lives = NULL, converged = FALSE))
  }
}

# ============================================================================
# HALF-LIFE CALCULATION
# ============================================================================

calculate_half_lives <- function(model, original_data, antigen, p_strain) {
  # Get time scaling factor
  time_sd <- sd(original_data$time_since_infection, na.rm = TRUE)
  
  # Extract posterior samples
  posterior <- as.data.frame(as_draws_df(model))
  
  # Initialize results
  results <- data.frame(
    parameter = character(),
    mean_half_life_days = numeric(),
    median_half_life_days = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    decay_rate_per_day = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate for each decay component
  for (param in c("rl", "rs", "ra")) {
    param_col <- paste0("b_", param, "_Intercept")
    
    if (param_col %in% names(posterior)) {
      # Convert to original time units
      decay_rates <- posterior[[param_col]] / time_sd
      half_lives <- log(2) / decay_rates
      
      # Remove extreme values
      hl_clean <- half_lives[half_lives > 0 & half_lives < 5000 & is.finite(half_lives)]
      
      if (length(hl_clean) > 10) {
        ci <- quantile(hl_clean, probs = c(0.025, 0.975), na.rm = TRUE)
        
        results <- rbind(results, data.frame(
          parameter = paste0(param, "_component"),
          mean_half_life_days = mean(hl_clean, na.rm = TRUE),
          median_half_life_days = median(hl_clean, na.rm = TRUE),
          ci_lower = ci[1],
          ci_upper = ci[2],
          decay_rate_per_day = mean(decay_rates, na.rm = TRUE)
        ))
      }
    }
  }
  
  # Save results
  write.csv(results, paste0("half_life_results_", antigen, ".csv"), row.names = FALSE)
  
  return(results)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Example usage for multiple antigens
antigens_to_model <- c("MSP1", "AMA1", "MSP2")  # Replace with your antigens

results_summary <- data.frame(
  antigen = character(),
  converged = logical(),
  max_rhat = numeric(),
  min_ess = numeric(),
  stringsAsFactors = FALSE
)

# Fit models for each antigen
for (antigen in antigens_to_model) {
  cat("\n=== Processing", antigen, "===\n")
  
  result <- fit_kinetics_model(antigen, antibody_data, strain_data)
  
  if (!is.null(result)) {
    max_rhat <- max(rhat(result$model), na.rm = TRUE)
    min_ess <- min(neff_ratio(result$model), na.rm = TRUE)
    
    results_summary <- rbind(results_summary, data.frame(
      antigen = antigen,
      converged = result$converged,
      max_rhat = max_rhat,
      min_ess = min_ess
    ))
    
    if (result$converged && !is.null(result$half_lives)) {
      cat("Half-life results:\n")
      print(result$half_lives)
    }
  }
}

# Save summary
write.csv(results_summary, "model_convergence_summary.csv", row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Total models attempted:", length(antigens_to_model), "\n")
cat("Successfully converged:", sum(results_summary$converged, na.rm = TRUE), "\n")
print(results_summary)