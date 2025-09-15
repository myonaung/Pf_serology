suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

# --- Paths ---
infile      <- "raw_data/high_transmission_relative_antibody_unit.csv"
outfile_svg <- "corr_heatmap_spearman_high_transmission.svg"
outfile_pdf <- "corr_heatmap_spearman_high_transmission.pdf"

# --- Load (no Endpoint filter) ---
dat_raw <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

# Drop known non-feature columns if present
drop_cols <- intersect(c("sampleid", "Timepoint"), names(dat_raw))
dat <- dat_raw[, setdiff(names(dat_raw), drop_cols), drop = FALSE]

# Keep only numeric features
is_num <- vapply(dat, is.numeric, logical(1))
dat <- dat[, is_num, drop = FALSE]
if (ncol(dat) < 2) stop("Need at least 2 numeric feature columns.")

# Drop zero-variance features (avoid NA/NaN in correlations)
zero_var <- vapply(dat, function(x) isTRUE(var(x, na.rm = TRUE) == 0), logical(1))
dat <- dat[, !zero_var, drop = FALSE]
if (ncol(dat) < 2) stop("All remaining features had zero variance.")

# --- Spearman correlation ---
cor_mat <- cor(dat, use = "pairwise.complete.obs", method = "spearman")

# --- Cluster order (same for rows & cols) ---
d   <- as.dist(1 - abs(cor_mat))
hc  <- hclust(d, method = "ward.D2")
ord <- hc$order
cor_ord <- cor_mat[ord, ord]

# --- Upper triangle only (keep diagonal) ---
cor_upper <- cor_ord
cor_upper[lower.tri(cor_upper, diag = FALSE)] <- NA

# --- Data-driven color limits & midpoint (from upper triangle values) ---
vals <- cor_upper[upper.tri(cor_upper, diag = TRUE)]
vals <- vals[!is.na(vals)]
# Round to nearest 0.05 to keep legend tidy
round5_down <- function(x) floor(x * 20) / 20
round5_up   <- function(x) ceiling(x * 20) / 20
lo  <- round5_down(min(vals))
hi  <- round5_up(max(vals))
mid <- median(vals)
brks <- pretty(c(lo, hi), n = 7)

# --- Long format ---
mm <- melt(cor_upper, varnames = c("Features_y", "Features_x"),
           value.name = "value", na.rm = TRUE)
lvl <- rownames(cor_ord)
mm$Features_x <- factor(mm$Features_x, levels = lvl)
mm$Features_y <- factor(mm$Features_y, levels = lvl)

# --- Plot ---
p <- ggplot(mm, aes(Features_x, Features_y, fill = value)) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_gradient2(
    low  = "#4575b4",
    mid  = "white",
    high = "#d73027",
    midpoint = mid,
    limits   = c(lo, hi),
    breaks   = brks,
    name     = expression("Spearman " * rho),
    oob      = scales::squish
  ) +
  coord_fixed() +
  labs(x = "Features", y = "Features") +
  theme_minimal(base_size = 16) +  # <- increase global base size (default is 11)
  theme(
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin  = margin(10, 10, 10, 10)
  )

# Show & save
print(p)
ggsave(outfile_svg, p, width = 10, height = 8, device = "svg", dpi = 600)
ggsave(outfile_pdf, p, width = 10, height = 8, device = "pdf", dpi = 600)

