#!/usr/bin/env Rscript

# Draw QQ and Manhattan plots from harmonized GWAS summary statistics.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse input summary stats and output image paths.
args <- parse_args()
require_args(args, c("stats", "qq", "manhattan"))


# Keep only valid P values for plotting.
stats <- read_tsv(args$stats)
if (!"p" %in% names(stats)) die("GWAS stats file is missing p column")
stats$p_num <- suppressWarnings(as.numeric(stats$p))
stats <- stats[is.finite(stats$p_num) & stats$p_num > 0 & stats$p_num <= 1, , drop = FALSE]


# Draw the QQ plot, preserving an empty image for empty results.
ensure_parent(args$qq)
png(args$qq, width = 1200, height = 1200, res = 150)
par(mar = c(5, 5, 2, 1))
if (!nrow(stats)) {
  plot.new()
  title("QQ plot")
} else {
  observed <- -log10(sort(stats$p_num))
  expected <- -log10(ppoints(length(observed)))
  limit <- max(c(expected, observed), na.rm = TRUE)
  plot(expected, observed, pch = 16, cex = 0.45, col = "#2f5d8c",
       xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
       xlim = c(0, limit), ylim = c(0, limit), main = "QQ plot")
  abline(0, 1, col = "#7a7a7a", lwd = 1.2)
}
invisible(dev.off())


# Start the Manhattan plot and handle missing coordinates cleanly.
ensure_parent(args$manhattan)
png(args$manhattan, width = 1800, height = 900, res = 150)
par(mar = c(5, 5, 2, 1))
if (!all(c("chrom", "pos") %in% names(stats)) || !nrow(stats)) {
  plot.new()
  title("Manhattan plot")
  invisible(dev.off())
  quit(status = 0)
}


# Normalize positions and chromosome ordering for the x-axis.
stats$pos_num <- suppressWarnings(as.numeric(stats$pos))
stats <- stats[is.finite(stats$pos_num), , drop = FALSE]
chrom_order <- unique(stats$chrom)
chrom_num <- suppressWarnings(as.numeric(chrom_order))
if (all(is.finite(chrom_num))) chrom_order <- chrom_order[order(chrom_num)]


# Build cumulative chromosome offsets for a standard Manhattan x-axis.
offset <- 0
x <- numeric(nrow(stats))
axis_at <- numeric(length(chrom_order))
for (i in seq_along(chrom_order)) {
  idx <- which(stats$chrom == chrom_order[[i]])
  x[idx] <- stats$pos_num[idx] + offset
  axis_at[[i]] <- mean(range(x[idx], na.rm = TRUE))
  offset <- max(x[idx], na.rm = TRUE)
}


# Draw points with alternating chromosome colors and a genome-wide line.
y <- -log10(stats$p_num)
colors <- rep(c("#2f5d8c", "#8a6f2a"), length.out = length(chrom_order))
plot(x, y, pch = 16, cex = 0.35, col = colors[match(stats$chrom, chrom_order)], xaxt = "n",
     xlab = "Chromosome", ylab = "-log10(P)", main = "Manhattan plot")
axis(1, at = axis_at, labels = chrom_order, cex.axis = 0.8)
abline(h = -log10(5e-8), col = "#9c2f2f", lty = 2)
invisible(dev.off())


# Report generated plot paths for logs.
cat("Wrote GWAS plots:", args$qq, args$manhattan, "\n")
