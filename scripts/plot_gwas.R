#!/usr/bin/env Rscript

# Draw association plots and exact summaries from harmonized statistics.

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
  "trait-summary" = "",
  "mac-qq" = "",
  "effect-frequency" = "",
  "large-plot-threshold" = "Inf"
))
require_args(args, c("stats", "qq", "manhattan"))

analysis_name <- ""
config <- list()
if (nzchar(args$config)) {
  config <- load_config(args$config)
  analysis_name <- config$project$analysis_name %||% ""
}
rare_mode <- any(nzchar(c(args[["trait-summary"]], args[["mac-qq"]], args[["effect-frequency"]])))
if (rare_mode) {
  require_args(args, c("config", "trait-summary", "mac-qq", "effect-frequency", "metrics-out", "top-hits-out"))
}


# Native regenie files can contain leading metadata and dozens of columns. Read
# only the fields required for summaries and plots, directly from disk.
first_column <- function(columns, aliases) {
  hit <- aliases[aliases %in% columns]
  if (length(hit)) hit[[1]] else ""
}

stats_header <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) die("GWAS stats file is empty: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
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

read_stats <- function(path, include_top_hits = FALSE, rare = FALSE) {
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
    effect = if (include_top_hits || rare) first_column(header$columns, c("beta_or_log_or", "BETA", "Effect")) else "",
    se = if (include_top_hits && !rare) first_column(header$columns, c("SE", "se")) else "",
    lci = if (rare) first_column(header$columns, c("LCI_Effect", "LCI_effect")) else "",
    uci = if (rare) first_column(header$columns, c("UCI_Effect", "UCI_effect")) else "",
    aaf = if (rare) first_column(header$columns, "AAF") else "",
    model = if (rare) first_column(header$columns, "Model") else "",
    num_cases = if (rare) first_column(header$columns, "Num_Cases") else "",
    cases_ref = if (rare) first_column(header$columns, "Cases_Ref") else "",
    cases_het = if (rare) first_column(header$columns, "Cases_Het") else "",
    cases_alt = if (rare) first_column(header$columns, "Cases_Alt") else "",
    num_controls = if (rare) first_column(header$columns, "Num_Controls") else "",
    controls_ref = if (rare) first_column(header$columns, "Controls_Ref") else "",
    controls_het = if (rare) first_column(header$columns, "Controls_Het") else "",
    controls_alt = if (rare) first_column(header$columns, "Controls_Alt") else ""
  )
  if (!nzchar(sources$p)) die("GWAS stats file is missing p or LOG10P column")
  if (rare) {
    required <- c("chrom", "pos", "variant_id", "trait", "effect", "lci", "uci", "aaf",
      "num_cases", "cases_ref", "cases_het", "cases_alt", "num_controls",
      "controls_ref", "controls_het", "controls_alt")
    missing <- required[!nzchar(unlist(sources[required]))]
    if (length(missing)) die("rare-variant HTP file is missing field(s): ", paste(missing, collapse = ", "))
  }
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

stats <- read_stats(args$stats, include_top_hits = nzchar(args[["top-hits-out"]]), rare = rare_mode)
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


numeric_column <- function(df, name) suppressWarnings(as.numeric(df[[name]]))


add_rare_fields <- function(df, valid_p, summary, genotype_mode) {
  if (nrow(summary) != 1L || !summary$trait_type[[1]] %in% c("bt", "qt")) {
    die("rare-variant plotting requires one active BT or QT trait summary")
  }
  if (identical(summary$skipped[[1]], "True")) die("cannot plot a skipped rare-variant trait")
  observed_traits <- unique(as.character(df$trait[!is.na(df$trait) & nzchar(df$trait)]))
  if (length(observed_traits) && !identical(observed_traits, as.character(summary$trait[[1]]))) {
    die("rare-variant HTP trait does not match its singleton Phase 2 summary")
  }

  df$aaf <- numeric_column(df, "aaf")
  df$maf <- pmin(df$aaf, 1 - df$aaf)
  num_cases <- numeric_column(df, "num_cases")
  num_controls <- numeric_column(df, "num_controls")
  num_controls[!is.finite(num_controls)] <- 0
  called_n <- num_cases + num_controls

  required_frequency <- valid_p & (
    !is.finite(df$aaf) | df$aaf < 0 | df$aaf > 1 | !is.finite(called_n) | called_n <= 0
  )
  if (any(required_frequency)) {
    die("valid rare-variant HTP rows contain invalid AAF or sample counts")
  }

  if (identical(genotype_mode, "hardcall")) {
    counts <- do.call(cbind, lapply(
      c("cases_ref", "cases_het", "cases_alt", "controls_ref", "controls_het", "controls_alt"),
      function(name) numeric_column(df, name)
    ))
    if (identical(summary$trait_type[[1]], "qt")) counts[, 4:6] <- 0
    invalid_counts <- valid_p & apply(
      counts, 1L, function(row) any(!is.finite(row) | row < 0 | row != floor(row))
    )
    if (any(invalid_counts)) die("valid hardcall HTP rows contain invalid genotype counts")
    called_from_counts <- rowSums(counts[, c(1, 2, 3, 4, 5, 6), drop = FALSE])
    expected_called <- if (identical(summary$trait_type[[1]], "bt")) called_n else num_cases
    if (any(valid_p & called_from_counts != expected_called)) {
      die("rare-variant HTP genotype counts do not equal their reported sample totals")
    }
    alt_count <- counts[, 2] + 2 * counts[, 3] + counts[, 5] + 2 * counts[, 6]
    df$mac <- pmin(alt_count, 2 * called_from_counts - alt_count)
    mac_definition <- "exact hardcall MAC"
  } else if (identical(genotype_mode, "dosage")) {
    # Dosage AAF can imply a non-integer expected allele count.
    df$mac <- 2 * called_n * df$maf
    mac_definition <- "dosage-based expected MAC"
  } else {
    die("rare-variant plotting requires remeta.genotype_mode hardcall or dosage")
  }
  if (any(valid_p & (!is.finite(df$mac) | df$mac <= 0))) {
    die("valid rare-variant HTP rows contain an invalid minor allele count")
  }

  df$mac_bin <- cut(
    df$mac,
    breaks = c(-Inf, 1, 5, 10, 20, Inf),
    labels = c("1", "2-5", "6-10", "11-20", ">20"),
    right = TRUE
  )
  list(df = df, trait_type = summary$trait_type[[1]], mac_definition = mac_definition)
}


top_hits <- function(df, rare = FALSE, max_rows = 10L) {
  columns <- if (rare) {
    c("chrom", "pos", "variant_id", "effect", "lci", "uci", "aaf", "maf", "mac")
  } else {
    c("chrom", "pos", "variant_id", "effect", "se")
  }
  if (!nrow(df)) {
    values <- setNames(replicate(length(columns), character(), simplify = FALSE), columns)
    values$p <- numeric()
    return(as.data.frame(values, stringsAsFactors = FALSE))
  }
  ord <- order(df$p_num, seq_len(nrow(df)), na.last = NA)
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
rare_info <- NULL
if (rare_mode) {
  trait_summary <- read_tsv(args[["trait-summary"]])
  genotype_mode <- tolower(as.character(config$remeta$genotype_mode %||% ""))
  rare_info <- add_rare_fields(stats, valid, trait_summary, genotype_mode)
  stats <- rare_info$df
  gwas_title <- paste0(gwas_title, " / target-region single variants")
}
stats <- stats[valid, , drop = FALSE]
valid_p_values <- nrow(stats)
hit_rows <- top_hits(stats, rare = rare_mode)
if (rare_mode) names(hit_rows)[names(hit_rows) == "variant_id"] <- "cpra"
plot_columns <- c("p_num", "chrom", "pos", "variant_id")
if (rare_mode) plot_columns <- c(plot_columns, "effect", "aaf", "maf", "mac", "mac_bin")
stats <- stats[intersect(plot_columns, names(stats))]


# Preserve full-rank expectations when a very large plot is thinned.
qq_plot_indices <- function(sorted_p, threshold, force = FALSE) {
  indices <- seq_along(sorted_p)
  thinned <- force || length(sorted_p) > threshold
  if (thinned) {
    indices <- sort(unique(c(
      seq.int(1L, length(sorted_p), by = 2L), which(sorted_p <= 1e-5)
    )))
  }
  list(indices = indices, thinned = thinned)
}


qq_coordinates <- function(sorted_p, indices = seq_along(sorted_p)) {
  a <- if (length(sorted_p) <= 10L) 3 / 8 else 1 / 2
  expected_p <- (indices - a) / (length(sorted_p) + 1 - 2 * a)
  list(expected = -log10(expected_p), observed = -log10(sorted_p[indices]))
}


# Draw the QQ plot, preserving an empty image for empty results.
lambda_gc <- genomic_lambda(stats$p_num)
sorted_p <- sort(stats$p_num)
qq_selection <- qq_plot_indices(sorted_p, large_plot_threshold)
qq_indices <- qq_selection$indices
qq_thinned <- qq_selection$thinned
qq <- qq_coordinates(sorted_p, qq_indices)
observed <- qq$observed
expected <- qq$expected
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


mac_levels <- c("1", "2-5", "6-10", "11-20", ">20")
mac_colors <- setNames(c("#5E3C99", "#3288BD", "#66C2A5", "#E6AB02", "#D53E4F"), mac_levels)
mac_counts <- setNames(integer(length(mac_levels)), mac_levels)
mac_qq_eligible <- mac_qq_plotted <- effect_eligible <- effect_plotted <- 0L
mac_qq_thinned <- effect_thinned <- FALSE


mac_display_labels <- function(mac_definition) {
  if (grepl("expected", mac_definition, fixed = TRUE)) {
    c("<=1", ">1-5", ">5-10", ">10-20", ">20")
  } else {
    mac_levels
  }
}


draw_mac_qq <- function(df, path, threshold, mac_definition) {
  ensure_parent(path)
  png(path, width = 1400, height = 1200, res = 150)
  on.exit(invisible(dev.off()), add = TRUE)
  par(mar = c(5, 5, 3.5, 1))

  # Trigger fallback from the full diagnostic, while retaining each stratum's
  # own expected ranks and every suggestive point.
  thin_all <- nrow(df) > threshold
  strata <- lapply(mac_levels, function(label) {
    p <- sort(df$p_num[as.character(df$mac_bin) == label])
    selection <- qq_plot_indices(p, threshold, force = thin_all)
    list(label = label, n = length(p), selection = selection, xy = qq_coordinates(p, selection$indices))
  })
  names(strata) <- mac_levels
  populated <- strata[vapply(strata, function(x) x$n > 0L, logical(1))]
  if (!length(populated)) {
    plot.new()
    title(paste0(gwas_title, "\nMAC-stratified QQ plot"))
    text(0.5, 0.5, "No variants with valid P value and MAC", col = "#4D4D4D")
  } else {
    limit <- max(unlist(lapply(populated, function(x) c(x$xy$expected, x$xy$observed))), na.rm = TRUE)
    plot(NA, xlim = c(0, limit), ylim = c(0, limit),
      xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
      main = paste0(gwas_title, "\nMAC-stratified QQ plot"))
    abline(0, 1, col = "#7A7A7A", lwd = 1.2)
    for (entry in populated) {
      points(entry$xy$expected, entry$xy$observed, pch = 16, cex = 0.38,
        col = adjustcolor(mac_colors[[entry$label]], alpha.f = 0.72))
    }
    display <- setNames(mac_display_labels(mac_definition), mac_levels)
    legend("topleft",
      legend = vapply(populated, function(x) paste0("MAC ", display[[x$label]], " (n=", x$n, ")"), character(1)),
      col = mac_colors[vapply(populated, `[[`, character(1), "label")], pch = 16,
      bty = "n", cex = 0.78)
  }
  strata
}


draw_effect_frequency <- function(df, path, trait_type, threshold, mac_definition) {
  effect <- numeric_column(df, "effect")
  y <- if (identical(trait_type, "bt")) suppressWarnings(log2(effect)) else effect
  eligible <- is.finite(df$maf) & df$maf > 0 & df$maf <= 0.5 & is.finite(y)
  if (identical(trait_type, "bt")) eligible <- eligible & effect > 0
  plot_df <- df[eligible, , drop = FALSE]
  plot_df$effect_plot <- y[eligible]
  eligible_n <- nrow(plot_df)
  thinned <- eligible_n > threshold
  if (thinned) {
    ordered <- order(plot_df$maf, plot_df$p_num, na.last = NA)
    keep <- sort(unique(c(ordered[seq.int(1L, length(ordered), by = 2L)], which(plot_df$p_num <= 1e-5))))
    plot_df <- plot_df[keep, , drop = FALSE]
  }

  ensure_parent(path)
  png(path, width = 1400, height = 1200, res = 150)
  on.exit(invisible(dev.off()), add = TRUE)
  par(mar = c(5, 5, 3.5, 1))
  plot_title_text <- paste0(gwas_title, "\nEffect versus cohort minor-allele frequency")
  if (!nrow(plot_df)) {
    plot.new()
    title(plot_title_text)
    text(0.5, 0.5, "No variants with valid effect and allele frequency", col = "#4D4D4D")
  } else {
    colors <- adjustcolor(mac_colors[as.character(plot_df$mac_bin)], alpha.f = 0.55)
    ylab <- if (identical(trait_type, "bt")) "log2(ALT-allele odds ratio)" else "ALT-allele beta"
    plot(plot_df$maf, plot_df$effect_plot, log = "x", pch = 16,
      cex = point_size(nrow(plot_df)), col = colors,
      xlab = "Cohort minor-allele frequency (log scale)", ylab = ylab, main = plot_title_text)
    abline(h = 0, col = "#777777", lwd = 1, lty = 2)
    present <- mac_levels[mac_levels %in% as.character(plot_df$mac_bin)]
    display <- setNames(mac_display_labels(mac_definition), mac_levels)
    legend("topright", legend = paste("MAC", display[present]), col = mac_colors[present], pch = 16,
      bty = "n", cex = 0.78)
  }
  list(eligible = eligible_n, plotted = nrow(plot_df), thinned = thinned)
}


if (rare_mode) {
  mac_counts[] <- as.integer(table(factor(as.character(stats$mac_bin), levels = mac_levels)))
  strata <- draw_mac_qq(stats, args[["mac-qq"]], large_plot_threshold, rare_info$mac_definition)
  mac_qq_eligible <- sum(vapply(strata, `[[`, integer(1), "n"))
  mac_qq_plotted <- sum(vapply(strata, function(x) length(x$selection$indices), integer(1)))
  mac_qq_thinned <- any(vapply(strata, function(x) x$selection$thinned, logical(1)))
  effect_result <- draw_effect_frequency(
    stats, args[["effect-frequency"]], rare_info$trait_type, large_plot_threshold, rare_info$mac_definition
  )
  effect_eligible <- effect_result$eligible
  effect_plotted <- effect_result$plotted
  effect_thinned <- effect_result$thinned
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
    if (qq_thinned || manhattan_data$thinned || mac_qq_thinned || effect_thinned) "True" else "False"
  ),
  stringsAsFactors = FALSE
)
if (rare_mode) {
  metrics <- rbind(metrics, data.frame(
    metric = c(
      "rare_variant_mode", "genotype_mode", "mac_definition",
      "mac_bin_1", "mac_bin_2_5", "mac_bin_6_10", "mac_bin_11_20", "mac_bin_gt20",
      "mac_qq_eligible_variants", "mac_qq_points_plotted",
      "effect_frequency_eligible_variants", "effect_frequency_points_plotted"
    ),
    value = c(
      "True", genotype_mode, rare_info$mac_definition, unname(mac_counts),
      mac_qq_eligible, mac_qq_plotted, effect_eligible, effect_plotted
    ),
    stringsAsFactors = FALSE
  ))
}
if (nzchar(args[["metrics-out"]])) write_tsv(metrics, args[["metrics-out"]])
if (nzchar(args[["top-hits-out"]])) write_tsv(hit_rows, args[["top-hits-out"]])

if (qq_thinned || manhattan_data$thinned || mac_qq_thinned || effect_thinned) {
  cat(
    "Applied large-plot fallback:", length(qq_indices), "of", length(sorted_p), "QQ points and",
    if (is.null(manhattan_data$df)) 0L else nrow(manhattan_data$df), "of",
    manhattan_data$eligible_variants, "Manhattan points plotted; exact summaries use all valid variants.\n"
  )
  if (rare_mode) {
    cat(
      "Rare diagnostics:", mac_qq_plotted, "of", mac_qq_eligible,
      "MAC-QQ points and", effect_plotted, "of", effect_eligible,
      "effect-frequency points plotted.\n"
    )
  }
}

# Report generated plot paths for logs.
plot_paths <- c(args$qq, args$manhattan, args[["manhattan-pdf"]])
if (rare_mode) plot_paths <- c(plot_paths, args[["mac-qq"]], args[["effect-frequency"]])
cat("Wrote association plots:", plot_paths, "\n")
