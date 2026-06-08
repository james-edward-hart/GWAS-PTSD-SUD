#!/usr/bin/env Rscript

# Draw the POP-MaD reference/study projected PC space for QC reports.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Cairo-backed devices may need a writable Fontconfig cache on restricted nodes.
font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) Sys.setenv(XDG_CACHE_HOME = font_cache)


args <- parse_args()
require_args(args, c("reference-pcs", "study-pcs", "assignments", "excluded", "out"))


reference <- read_tsv(args[["reference-pcs"]])
study <- read_tsv(args[["study-pcs"]])
assignments <- read_tsv(args$assignments)
excluded <- read_tsv(args$excluded)

pc_cols <- paste0("PC", 1:3)
available_pcs <- pc_cols[pc_cols %in% names(reference) & pc_cols %in% names(study)]
if (length(available_pcs) < 2) die("POP-MaD projection plot requires at least PC1 and PC2")
require_columns(reference, c("FID", "IID", "super_population", available_pcs), "reference PC file")
require_columns(study, c("FID", "IID", available_pcs), "study PC file")


as_numeric_pc <- function(rows, cols, label) {
  for (col in cols) {
    value <- suppressWarnings(as.numeric(rows[[col]]))
    invalid <- !is.finite(value)
    if (any(invalid)) {
      warning("dropping ", sum(invalid), " ", label, " row(s) with missing or non-finite ", col)
      rows <- rows[!invalid, , drop = FALSE]
      value <- suppressWarnings(as.numeric(rows[[col]]))
    }
    rows[[col]] <- value
  }
  rows
}

reference <- as_numeric_pc(reference, available_pcs, "reference")
study <- as_numeric_pc(study, available_pcs, "study")
if (!nrow(reference) || !nrow(study)) die("POP-MaD projection plot has no plottable reference or study rows")


key <- function(rows) paste(rows$FID, rows$IID, sep = "\t")
study_key <- key(study)
study$plot_status <- "unassigned"
study$plot_ancestry <- "Unassigned"

if (nrow(assignments)) {
  require_columns(assignments, c("FID", "IID", "ancestry"), "POP-MaD assignments")
  idx <- match(study_key, key(assignments))
  assigned <- !is.na(idx)
  study$plot_status[assigned] <- "assigned"
  study$plot_ancestry[assigned] <- assignments$ancestry[idx[assigned]]
}

if (nrow(excluded)) {
  require_columns(excluded, c("FID", "IID", "best_ancestry", "reason"), "POP-MaD excluded samples")
  idx <- match(study_key, key(excluded))
  dropped <- !is.na(idx)
  study$plot_status[dropped] <- "excluded"
  best <- excluded$best_ancestry[idx[dropped]]
  best[is.na(best) | !nzchar(best)] <- "Excluded"
  study$plot_ancestry[dropped] <- best
}

reference$plot_ancestry <- reference$super_population
reference$plot_ancestry[is.na(reference$plot_ancestry) | !nzchar(reference$plot_ancestry)] <- "Reference"
study$plot_ancestry[is.na(study$plot_ancestry) | !nzchar(study$plot_ancestry)] <- "Unassigned"
labels <- sort(unique(c(reference$plot_ancestry, study$plot_ancestry)))
labels <- labels[!is.na(labels) & nzchar(labels)]
base_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
  "#56B4E9", "#6A3D9A", "#8A8A8A", "#117733", "#882255")
palette <- rep(base_palette, length.out = length(labels))
names(palette) <- labels
for (label in c("Unassigned", "Excluded")) {
  if (label %in% names(palette)) palette[[label]] <- "#777777"
}
if ("Reference" %in% names(palette)) palette[["Reference"]] <- "#777777"


pad_range <- function(x, frac = 0.08) {
  rng <- range(x, na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span == 0) span <- max(abs(rng), 1)
  rng + c(-1, 1) * span * frac
}

open_png <- function(path) {
  ensure_parent(path)
  png_args <- list(filename = path, width = 2400, height = 1800, res = 300)
  if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo-png"
  do.call(png, png_args)
}

draw_legend <- function(labels) {
  legend("topright", legend = labels, pch = 16, col = palette[labels],
    title = "Ancestry label", bty = "n", cex = 0.72, title.adj = 0)
  legend("bottomright", legend = c("Reference", "Study assigned", "Study excluded"),
    pch = c(16, 24, 4), col = c("#555555", "#111111", "#444444"),
    pt.bg = c(NA, "#FFFFFF", NA), bty = "n", cex = 0.72)
}

draw_2d <- function() {
  all_x <- c(reference$PC1, study$PC1)
  all_y <- c(reference$PC2, study$PC2)
  par(mar = c(4.4, 4.8, 2.8, 0.8), family = "sans", las = 1)
  plot(all_x, all_y, type = "n", xlab = "PC1", ylab = "PC2", bty = "n",
    xlim = pad_range(all_x), ylim = pad_range(all_y))
  grid(col = "#E7E7E7", lwd = 0.7)
  points(reference$PC1, reference$PC2, pch = 16, cex = 0.36,
    col = adjustcolor(palette[reference$plot_ancestry], alpha.f = 0.28))
  assigned <- study$plot_status == "assigned"
  if (any(assigned)) {
    points(study$PC1[assigned], study$PC2[assigned], pch = 24, cex = 0.72,
      col = "#111111", bg = palette[study$plot_ancestry[assigned]])
  }
  if (any(!assigned)) {
    points(study$PC1[!assigned], study$PC2[!assigned], pch = 4, cex = 0.74,
      col = "#444444", lwd = 1.1)
  }
  title("POP-MaD projected PC space")
  draw_legend(labels)
}

draw_3d <- function() {
  all_x <- c(reference$PC1, study$PC1)
  all_y <- c(reference$PC2, study$PC2)
  all_z <- c(reference$PC3, study$PC3)
  xlim <- pad_range(all_x)
  ylim <- pad_range(all_y)
  zlim <- pad_range(all_z)

  par(mar = c(2.4, 2.2, 3.0, 0.8), family = "sans", las = 1, xpd = NA)
  pmat <- persp(
    x = xlim,
    y = ylim,
    z = matrix(zlim[[1]], nrow = 2, ncol = 2),
    xlim = xlim,
    ylim = ylim,
    zlim = zlim,
    theta = 35,
    phi = 24,
    expand = 0.75,
    ticktype = "detailed",
    xlab = "PC1",
    ylab = "PC2",
    zlab = "PC3",
    col = NA,
    border = "#D8D8D8",
    main = "POP-MaD projected PC space"
  )

  ref_order <- order(reference$PC3)
  ref_xy <- trans3d(reference$PC1[ref_order], reference$PC2[ref_order], reference$PC3[ref_order], pmat)
  points(ref_xy, pch = 16, cex = 0.32,
    col = adjustcolor(palette[reference$plot_ancestry[ref_order]], alpha.f = 0.25))

  assigned <- study$plot_status == "assigned"
  if (any(assigned)) {
    idx <- which(assigned)[order(study$PC3[assigned])]
    xy <- trans3d(study$PC1[idx], study$PC2[idx], study$PC3[idx], pmat)
    points(xy, pch = 24, cex = 0.78, col = "#111111", bg = palette[study$plot_ancestry[idx]])
  }
  if (any(!assigned)) {
    idx <- which(!assigned)[order(study$PC3[!assigned])]
    xy <- trans3d(study$PC1[idx], study$PC2[idx], study$PC3[idx], pmat)
    points(xy, pch = 4, cex = 0.78, col = "#444444", lwd = 1.1)
  }
  draw_legend(labels)
}


open_png(args$out)
on.exit(invisible(dev.off()), add = TRUE)
if (length(available_pcs) >= 3) {
  draw_3d()
} else {
  draw_2d()
}

cat("Wrote POP-MaD projection plot:", args$out, "\n")
