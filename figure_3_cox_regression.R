# Minimal Cox Regression Analysis with Forest Plot
# For publication purposes

# Load required libraries
library(survival)
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(1234)

# Load and prepare data
# Replace with your actual file paths
metadata <- read.csv("metadata.csv")
ab_data <- read.csv("antibody_data.csv")

# Transform antibody data into tertiles (0, 1, 2)
ab_data[, 3:ncol(ab_data)] <- lapply(ab_data[, 3:ncol(ab_data)], function(col) {
  cut(col, breaks = quantile(col, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
      include.lowest = TRUE, labels = c(0, 1, 2))
})

# Merge data and prepare for analysis
baseline_data <- subset(ab_data, Timepoint == "Baseline")
analysis_data <- merge(metadata, baseline_data, by = "sampleid")
analysis_data$village <- as.factor(analysis_data$village)

# Define variables for analysis (adjust column indices as needed)
antibody_vars <- names(analysis_data)[17:41]

# Univariate Cox regression with adjustments
univariate_results <- lapply(antibody_vars, function(var) {
  formula <- as.formula(paste("Surv(event, status) ~", var, 
                              "+ pfbase + hb + log10(age) + pmbase + pobase"))
  coxph(formula, data = analysis_data)
})

# Extract and format results
extract_results <- function(model) {
  summary_model <- summary(model)
  c(
    HR = summary_model$coef[1, 2],
    Lower_CI = summary_model$conf.int[1, "lower .95"],
    Upper_CI = summary_model$conf.int[1, "upper .95"],
    p_value = summary_model$coef[1, 5]
  )
}

results_df <- do.call(rbind, lapply(univariate_results, extract_results)) %>%
  as.data.frame() %>%
  mutate(
    Variable = antibody_vars,
    HR_label = paste0(round(HR, 2), " (", round(Lower_CI, 2), " - ", round(Upper_CI, 2), ")"),
    PPE = ifelse(Lower_CI > 1 | Upper_CI < 1, (1 - HR) * 100, NA)
  ) %>%
  arrange(HR) %>%
  mutate(Variable = factor(Variable, levels = Variable))

# Create forest plot
forest_plot <- ggplot(results_df, aes(y = Variable, x = HR, xmin = Lower_CI, xmax = Upper_CI)) +
  geom_pointrange(color = "black", size = 0.4) +
  geom_point(aes(color = HR), size = 2) +
  scale_color_gradient(low = "skyblue", high = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_label(
    data = subset(results_df, !is.na(PPE)),
    aes(x = max(HR, na.rm = TRUE) + 1.5, y = Variable, 
        label = paste0("PPE: ", round(PPE, 1), "%"), fill = PPE),
    hjust = 1, vjust = 0.5, size = 3, color = "black",
    angle = 90, na.rm = TRUE
  ) +
  scale_fill_gradientn(colors = c("white", "lightyellow", "skyblue"), guide = "none") +
  coord_flip() +
  scale_x_continuous(limits = c(0.1, 3)) +
  labs(x = "aHR (95% CI)", y = "") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 10, hjust = 1),
    axis.text.x = element_text(size = 10, angle = 90),
    axis.title.x = element_text(margin = margin(t = 10)),
    plot.margin = margin(5, 5, 5, 30),
    legend.position = "none"
  )

# Display and save plot
print(forest_plot)
ggsave("forest_plot.svg", plot = forest_plot, width = 8, height = 4)