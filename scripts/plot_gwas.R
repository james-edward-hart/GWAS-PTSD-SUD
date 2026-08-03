#!/usr/bin/env Rscript

# Draw QQ and Manhattan plots from harmonized GWAS summary statistics.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Cairo-backed devices may need a writable Fontconfig cache on restricted nodes.
font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) Sys.setenv(XDG_CACHE_HOME = font_cache)


# Parse input summary stats and output image paths.
args <- parse_args(defaults = list(
  "config" = "",
  "manhattan-pdf" = "",
  "metrics-out" = "",
  "top-hits-out" = "",
  "large-plot-threshold" = "Inf"
))
require_args(args, c("stats", "qq", "manhattan"))

analysis_name <- ""
if (nzchar(args$config)) {
  config <- load_config(args$config)
  analysis_name <- config$project$analysis_name %||% ""
}


# Native regenie files can contain leading metadata and dozens of columns. Read
# only the fields required for summaries and plots, directly from disk.
first_column <- function(columns, aliases) {
  hit <- aliases[aliases %in% columns]
  if (length(hit)) hit[[1]] else ""
}

stats_header <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) die("GWAS stats file is empty: ", path)
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  skipped <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) die("GWAS stats file has no header: ", path)
    if (!nzchar(trimws(line)) || startsWith(line, "##")) {
      skipped <- skipped + 1L
      next
    }
    sep <- if (grepl("\t", line, fixed = TRUE)) "\t" else ""
    columns <- strsplit(trimws(line), if (nzchar(sep)) "\t" else "[[:space:]]+")[[1]]
    return(list(columns = columns, skipped = skipped, sep = sep))
  }
}

read_stats <- function(path, include_top_hits = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) die("R package 'data.table' is required")
  header <- stats_header(path)
  sources <- list(
    p = first_column(header$columns, c("p", "P", "Pval", "LOG10P", "log10p")),
    chrom = first_column(header$columns, c("chrom", "CHROM", "#CHROM", "Chr")),
    pos = first_column(header$columns, c("pos", "GENPOS", "POS", "Pos")),
    variant_id = first_column(header$columns, c("variant_id", "ID", "Name")),
    trait = first_column(header$columns, c("trait", "Trait")),
    ancestry = first_column(header$columns, c("ancestry", "Ancestry")),
    build = first_column(header$columns, c("build", "Build")),
    effect = if (include_top_hits) first_column(header$columns, c("beta_or_log_or", "BETA", "Effect")) else "",
    se = if (include_top_hits) first_column(header$columns, c("SE", "se")) else ""
  )
  if (!nzchar(sources$p)) die("GWAS stats file is missing p or LOG10P column")
  selected <- unique(unname(unlist(sources, use.names = FALSE)))
  selected <- selected[nzchar(selected)]
  fread_args <- list(
    input = path,
    skip = header$skipped,
    select = selected,
    header = TRUE,
    data.table = FALSE,
    showProgress = FALSE,
    na.strings = c("NA", "NaN", ".")
  )
  if (nzchar(header$sep)) fread_args$sep <- header$sep
  rows <- tryCatch(
    do.call(data.table::fread, fread_args),
    error = function(err) die("could not read GWAS stats file ", path, ": ", conditionMessage(err))
  )

  raw_p <- suppressWarnings(as.numeric(rows[[sources$p]]))
  out <- data.frame(
    p_num = if (sources$p %in% c("LOG10P", "log10p")) 10 ^ -raw_p else raw_p,
    stringsAsFactors = FALSE
  )
  for (name in setdiff(names(sources), "p")) {
    source <- sources[[name]]
    if (nzchar(source)) out[[name]] <- rows[[source]]
  }
  out
}

large_plot_threshold <- suppressWarnings(as.numeric(args[["large-plot-threshold"]]))
if (length(large_plot_threshold) != 1L || is.na(large_plot_threshold) || large_plot_threshold <= 0) {
  die("large-plot-threshold must be a positive number or Inf")
}

stats <- read_stats(args$stats, include_top_hits = nzchar(args[["top-hits-out"]]))
total_variants <- nrow(stats)

format_lambda <- function(value) {
  if (is.finite(value)) sprintf("%.3f", value) else "NA"
}

plot_title <- function(df, analysis_name = "") {
  pieces <- character()
  for (column in c("trait", "ancestry", "build")) {
    if (column %in% names(df)) {
      value <- as.character(df[[column]])
      value <- unique(value[nzchar(value)])
      value <- value[!is.na(value) & value != "NA"]
      if (length(value) == 1) pieces <- c(pieces, value)
    }
  }
  label <- if (length(pieces)) paste(pieces, collapse = " / ") else "GWAS"
  if (nzchar(analysis_name)) paste(analysis_name, label, sep = " / ") else label
}


top_hits <- function(df, max_rows = 10L) {
  columns <- c("chrom", "pos", "variant_id", "effect", "se")
  if (!nrow(df)) {
    return(data.frame(chrom = character(), pos = character(), variant_id = character(),
      effect = character(), se = character(), p = numeric(), stringsAsFactors = FALSE))
  }
  ord <- order(df$p_num, na.last = NA)
  ord <- ord[seq_len(min(length(ord), max_rows))]
  values <- lapply(columns, function(column) {
    if (column %in% names(df)) as.character(df[[column]][ord]) else rep("NA", length(ord))
  })
  names(values) <- columns
  values$p <- df$p_num[ord]
  as.data.frame(values, stringsAsFactors = FALSE)
}


gwas_title <- plot_title(stats, analysis_name)
valid <- is.finite(stats$p_num) & stats$p_num > 0 & stats$p_num <= 1
stats <- stats[valid, , drop = FALSE]
valid_p_values <- nrow(stats)
hit_rows <- top_hits(stats)
stats <- stats[intersect(c("p_num", "chrom", "pos", "variant_id"), names(stats))]


# Draw the QQ plot, preserving an empty image for empty results.
lambda_gc <- genomic_lambda(stats$p_num)
sorted_p <- sort(stats$p_num)
qq_indices <- seq_along(sorted_p)
qq_thinned <- length(sorted_p) > large_plot_threshold
if (qq_thinned) {
  qq_indices <- sort(unique(c(seq.int(1L, length(sorted_p), by = 2L), which(sorted_p <= 1e-5), length(sorted_p))))
}
ppoints_a <- if (length(sorted_p) <= 10L) 3 / 8 else 1 / 2
expected_p <- (qq_indices - ppoints_a) / (length(sorted_p) + 1 - 2 * ppoints_a)
observed <- -log10(sorted_p[qq_indices])
expected <- -log10(expected_p)
ensure_parent(args$qq)
png(args$qq, width = 1200, height = 1200, res = 150)
par(mar = c(5, 5, 3.5, 1))
if (!nrow(stats)) {
  plot.new()
  title(paste0(gwas_title, "\nQQ plot (lambda GC = ", format_lambda(lambda_gc), ")"))
} else {
  limit <- max(c(expected, observed), na.rm = TRUE)
  plot(expected, observed, pch = 16, cex = 0.45, col = "#2f5d8c",
       xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
       xlim = c(0, limit), ylim = c(0, limit),
       main = paste0(gwas_title, "\nQQ plot (lambda GC = ", format_lambda(lambda_gc), ")"))
  abline(0, 1, col = "#7a7a7a", lwd = 1.2)
}
invisible(dev.off())


# Draw a publication-quality Manhattan plot while keeping base R dependencies.
open_manhattan_device <- function(path, format) {
  ensure_parent(path)
  if (format == "png") {
    png_args <- list(filename = path, width = 3600, height = 1800, res = 300)
    if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo-png"
    do.call(png, png_args)
  } else if (format == "pdf") {
    pdf(path, width = 12, height = 6, useDingbats = FALSE, colormodel = "srgb")
  } else {
    die("unsupported Manhattan plot format: ", format)
  }
}

clean_label <- function(value) {
  value <- gsub("^chr", "", as.character(value), ignore.case = TRUE)
  value <- trimws(value)
  value[value %in% c("", "NA", ".", "0")] <- NA_character_
  value
}

chrom_key <- function(value) {
  upper <- toupper(as.character(value))
  mapped <- suppressWarnings(as.numeric(upper))
  mapped[upper == "X"] <- 23
  mapped[upper == "Y"] <- 24
  mapped[upper %in% c("XY", "PAR")] <- 25
  mapped[upper %in% c("M", "MT", "MITO")] <- 26
  missing <- !is.finite(mapped)
  mapped[missing] <- 1000 + match(upper[missing], sort(unique(upper[missing])))
  mapped
}

lead_hit_labels <- function(df, max_labels = 8, window_bp = 500000) {
  if (!"variant_id" %in% names(df)) df$variant_id <- ""
  candidates <- df[df$p_num <= 5e-8, , drop = FALSE]
  if (!nrow(candidates)) return(integer())
  candidates <- candidates[order(candidates$p_num, na.last = NA), , drop = FALSE]
  selected <- integer()
  for (row_idx in seq_len(nrow(candidates))) {
    row <- candidates[row_idx, , drop = FALSE]
    if (length(selected)) {
      same_locus <- df$chrom_label[selected] == row$chrom_label &
        abs(df$pos_num[selected] - row$pos_num) <= window_bp
      if (any(same_locus, na.rm = TRUE)) next
    }
    selected <- c(selected, as.integer(rownames(row)))
    if (length(selected) >= max_labels) break
  }
  selected
}

point_size <- function(n) {
  if (n >= 1e6) return(0.14)
  if (n >= 5e5) return(0.18)
  if (n >= 1e5) return(0.23)
  if (n >= 2e4) return(0.30)
  0.44
}

draw_empty_manhattan <- function(message) {
  par(mar = c(4.4, 5.2, 2.4, 0.8), xaxs = "i", yaxs = "i", las = 1, family = "sans")
  plot.new()
  title(paste0(gwas_title, "\nManhattan plot"))
  text(0.5, 0.5, message, col = "#4D4D4D", cex = 0.95)
}

prepare_manhattan <- function(df, threshold = Inf) {
  empty <- function() list(df = NULL, eligible_variants = 0L, thinned = FALSE)
  if (!all(c("chrom", "pos") %in% names(df))) return(empty())
  df$pos_num <- suppressWarnings(as.numeric(df$pos))
  df$chrom_label <- clean_label(df$chrom)
  df <- df[is.finite(df$pos_num) & !is.na(df$chrom_label), , drop = FALSE]
  eligible_variants <- nrow(df)
  if (!eligible_variants) return(empty())
  thinned <- eligible_variants > threshold
  if (thinned) {
    keep <- order(df$p_num, seq_len(eligible_variants), na.last = NA)
    keep <- keep[seq_len(ceiling(eligible_variants / 2))]
    df <- df[keep, , drop = FALSE]
  }
  df$neg_log10_p <- -log10(df$p_num)
  df <- df[order(chrom_key(df$chrom_label), df$pos_num, df$p_num), , drop = FALSE]
  rownames(df) <- seq_len(nrow(df))

  chrom_order <- unique(df$chrom_label)
  chrom_span <- tapply(df$pos_num, df$chrom_label, max, na.rm = TRUE)
  gap <- max(1000000, max(chrom_span, na.rm = TRUE) * 0.005)
  offset <- 0
  df$plot_x <- NA_real_
  axis_at <- numeric(length(chrom_order))
  names(axis_at) <- chrom_order
  boundaries <- numeric(max(length(chrom_order) - 1, 0))

  for (i in seq_along(chrom_order)) {
    chrom <- chrom_order[[i]]
    idx <- which(df$chrom_label == chrom)
    df$plot_x[idx] <- df$pos_num[idx] + offset
    axis_at[[i]] <- mean(range(df$plot_x[idx], na.rm = TRUE))
    offset <- max(df$plot_x[idx], na.rm = TRUE)
    if (i < length(chrom_order)) {
      boundaries[[i]] <- offset + gap / 2
      offset <- offset + gap
    }
  }

  list(
    df = df,
    chrom_order = chrom_order,
    axis_at = axis_at,
    boundaries = boundaries,
    eligible_variants = eligible_variants,
    thinned = thinned
  )
}

draw_manhattan <- function(plot_data) {
  if (is.null(plot_data$df)) {
    draw_empty_manhattan("No variants with valid chromosome, position, and P value")
    return(invisible(NULL))
  }

  df <- plot_data$df
  if (!"variant_id" %in% names(df)) df$variant_id <- ""
  title <- gwas_title
  genomewide <- -log10(5e-8)
  suggestive <- -log10(1e-5)
  label_idx <- lead_hit_labels(df)
  ymax_data <- max(df$neg_log10_p, genomewide, na.rm = TRUE)
  ymax <- ceiling(ymax_data + if (length(label_idx)) max(1.4, 0.12 * ymax_data) else 0.8)
  xmax <- max(df$plot_x, na.rm = TRUE)
  colors <- rep(c("#1F5A85", "#B66D2D"), length.out = length(plot_data$chrom_order))
  point_colors <- colors[match(df$chrom_label, plot_data$chrom_order)]
  significant <- df$p_num <= 5e-8

  par(mar = c(4.4, 5.2, 2.4, 0.8), xaxs = "i", yaxs = "i", las = 1, family = "sans")
  plot(
    df$plot_x, df$neg_log10_p,
    type = "n", xaxt = "n", yaxt = "n", bty = "n",
    xlab = "Chromosome", ylab = expression(-log[10](italic(P))),
    xlim = c(0, xmax), ylim = c(0, ymax)
  )

  y_ticks <- pretty(c(0, ymax), n = 6)
  y_ticks <- y_ticks[y_ticks >= 0 & y_ticks <= ymax]
  abline(h = y_ticks, col = "#E7E7E7", lwd = 0.7)
  if (length(plot_data$boundaries)) abline(v = plot_data$boundaries, col = "#F0F0F0", lwd = 0.7)

  points(df$plot_x, df$neg_log10_p, pch = 16, cex = point_size(nrow(df)), col = point_colors)
  if (any(significant)) {
    points(df$plot_x[significant], df$neg_log10_p[significant], pch = 16, cex = point_size(nrow(df)) * 1.8, col = "#CC3311")
  }

  abline(h = suggestive, col = "#777777", lwd = 0.8, lty = 3)
  abline(h = genomewide, col = "#CC3311", lwd = 1.0, lty = 2)
  axis(1, at = plot_data$axis_at, labels = plot_data$chrom_order, cex.axis = 0.78, lwd = 0, lwd.ticks = 0.8, col.ticks = "#4D4D4D")
  axis(2, at = y_ticks, labels = y_ticks, cex.axis = 0.78, lwd = 0, lwd.ticks = 0.8, col.ticks = "#4D4D4D")
  mtext(title, side = 3, adj = 0, line = 0.7, font = 2, cex = 0.95)
  text(xmax * 0.995, genomewide + 0.04, "Genome-wide", adj = c(1, 0), cex = 0.58, col = "#CC3311")
  text(xmax * 0.995, suggestive + 0.04, "Suggestive", adj = c(1, 0), cex = 0.58, col = "#555555")

  if (length(label_idx)) {
    label_y <- pmin(df$neg_log10_p[label_idx] + 0.45, ymax - 0.4)
    raw_labels <- as.character(df$variant_id[label_idx])
    raw_labels[!nzchar(raw_labels) | raw_labels %in% c("NA", ".")] <- paste0(df$chrom_label[label_idx], ":", df$pos_num[label_idx])
    labels <- ifelse(nchar(raw_labels) > 28, paste0(substr(raw_labels, 1, 25), "..."), raw_labels)
    segments(df$plot_x[label_idx], df$neg_log10_p[label_idx] + 0.08, df$plot_x[label_idx], label_y - 0.08, col = "#4D4D4D", lwd = 0.5)
    text(df$plot_x[label_idx], label_y, labels, pos = 3, offset = 0.12, cex = 0.52, col = "#222222", xpd = NA)
  }
}

render_manhattan <- function(path, format, plot_data) {
  open_manhattan_device(path, format)
  on.exit(invisible(dev.off()), add = TRUE)
  draw_manhattan(plot_data)
}

manhattan_data <- prepare_manhattan(stats, large_plot_threshold)
render_manhattan(args$manhattan, "png", manhattan_data)
if (nzchar(args[["manhattan-pdf"]])) {
  render_manhattan(args[["manhattan-pdf"]], "pdf", manhattan_data)
}


metrics <- data.frame(
  metric = c(
    "total_variants",
    "valid_p_value_variants",
    "lambda_gc",
    "genomewide_significant_variants",
    "suggestive_variants",
    "qq_eligible_variants",
    "qq_points_plotted",
    "manhattan_eligible_variants",
    "manhattan_points_plotted",
    "large_plot_threshold",
    "plot_thinning_applied"
  ),
  value = c(
    total_variants,
    valid_p_values,
    if (is.finite(lambda_gc)) sprintf("%.15g", lambda_gc) else "NA",
    sum(stats$p_num <= 5e-8),
    sum(stats$p_num <= 1e-5),
    length(sorted_p),
    length(qq_indices),
    manhattan_data$eligible_variants,
    if (is.null(manhattan_data$df)) 0L else nrow(manhattan_data$df),
    if (is.finite(large_plot_threshold)) format(large_plot_threshold, scientific = FALSE, trim = TRUE) else "Inf",
    if (qq_thinned || manhattan_data$thinned) "True" else "False"
  ),
  stringsAsFactors = FALSE
)
if (nzchar(args[["metrics-out"]])) write_tsv(metrics, args[["metrics-out"]])
if (nzchar(args[["top-hits-out"]])) write_tsv(hit_rows, args[["top-hits-out"]])

if (qq_thinned || manhattan_data$thinned) {
  cat(
    "Applied large-plot fallback:", length(qq_indices), "of", length(sorted_p), "QQ points and",
    if (is.null(manhattan_data$df)) 0L else nrow(manhattan_data$df), "of",
    manhattan_data$eligible_variants, "Manhattan points plotted; exact summaries use all valid variants.\n"
  )
}

# Report generated plot paths for logs.
cat("Wrote GWAS plots:", args$qq, args$manhattan, args[["manhattan-pdf"]], "\n")
