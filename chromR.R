#!/usr/bin/env Rscript

# Charge library
#------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
})

# Parse manual arguments
#------------------------------------------------------------------------

parse_args_manual <- function(args) {
  
  # Initialize default values
  opts <- list(
    gff_file         = NULL,
    feature_file     = NULL,     # formerly fill_file
    feature_mode     = NULL,     # "fill" | "keyword"
    keywords         = NULL,
    strict_filter    = FALSE,
    max_chromosomes  = Inf,
    palette          = NULL,
    colors           = NULL,
    line_plot        = FALSE,
    accumulated_plot = FALSE,
    density_mode     = "hist",
    interactive      = FALSE,
    window_mode      = FALSE,
    window_size      = NULL,
    gene_list        = NULL,
    min_genes        = 2,
    facet_plot   = FALSE,
    facet_mode   = "density",
    facet_scales = "fixed"
  )
  
  i <- 1
  while (i <= length(args)) {
    flag <- args[i]


    #if (flag %in% c("-g","--gff_file")) {
    if (flag == "-g" || flag == "--gff_file") {
      opts$gff_file <- args[i + 1]; i <- i + 1
    } else if (flag == "ff" || flag == "--fill_file") {
      opts$feature_file <- args[i + 1]
      opts$feature_mode <- "fill"
      i <- i + 1
    } else if (flag == "--keywords" || flag == "-k") {
      opts$feature_mode <- "keyword"
      opts$keywords <- c()
      j <- i + 1
      while (j <= length(args) && !grepl("^--", args[j])) {
        opts$keywords <- c(opts$keywords, args[j])
        j <- j + 1
      }
      i <- j - 1
    } else if (flag == "--strict" || flag == "-s") {
      opts$strict_filter <- TRUE
    } else if (flag == "--number" || flag == "-n") {
      opts$max_chromosomes <- as.integer(args[i+1])
      i <- i + 1
    } else if (flag == "--palette" || flag == "-p") {
      opts$palette <- args[i+1]; i <- i + 1
    } else if (flag == "--line_plot" || flag == "-lp") {
      opts$line_plot <- TRUE
    } else if (flag == "--accumulated_plot" || flag == "-ap") {
      opts$accumulated_plot <- TRUE
    } else if (flag == "--density_mode" || flag == "-dm") {
      opts$density_mode <- args[i+1]; i <- i + 1
    }  else if (flag == "--facet_plot") {
      opts$facet_plot <- TRUE
    } else if (flag == "--facet_mode") {
      opts$facet_mode <- args[i+1]; i <- i + 1
    } else if (flag == "--facet_scales") {
      opts$facet_scales <- args[i+1]; i <- i + 1
    } else if (flag == "--interactive" || flag == "-int") {
      opts$interactive <- TRUE
    } else if (flag == "--window_count" || flag == "-wc") {
      opts$window_mode <- TRUE
    } else if (flag == "--window_size") {
      opts$window_size <- as.integer(args[i+1]); i <- i + 1
    } else if (flag == "--genes" || flag == "-gene") {
      opts$gene_list <- c()
      j <- i + 1
      while (j <= length(args) && !grepl("^--", args[j])) {
        opts$gene_list <- c(opts$gene_list, args[j])
        j <- j + 1
      }
      i <- j - 1
    } else if (flag == "--min_genes" || flag == "-mg") {
      opts$min_genes <- as.integer(args[i+1]); i <- i + 1
    } else if (flag == "--colors" || flag == "-c") {  # added colors option
      opts$colors <- c()
      j <- i + 1
      while (j <= length(args) && !grepl("^--", args[j])) {
        opts$colors <- c(opts$colors, args[j])
        j <- j + 1
      }
      i <- j - 1
    }

    i <- i + 1 # this increments the main loop counter
  }
  
  return(opts)
}

# Arguments validation and processing
#------------------------------------------------------------------------
validate_args <- function(opts) {

  if (is.null(opts$gff_file) || !file.exists(opts$gff_file))
    stop("Valid --gff file is required")

  if (is.null(opts$feature_mode))
    stop("You must specify --fill_file or --keywords")

  if (opts$feature_mode == "keyword" &&
      length(opts$keywords) %% 2 != 0)
    stop("--keywords must be attribute/type pairs")

  if (opts$window_mode) {
    if (is.null(opts$window_size))
      stop("--window_count requires --window_size")
    if (is.null(opts$gene_list))
      stop("--window_count requires --genes")
  }
}

# Data loading
#------------------------------------------------------------------------

load_gff <- function(file) {read_tsv(file, comment = "#", col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"), show_col_types = FALSE)
}

compute_chrom_limits <- function(gff) {
  gff %>%
    group_by(seqid) %>%
    summarise(
      chrom_start = min(start),
      chrom_end   = max(end),
      chrom_length = chrom_end - chrom_start,
      .groups = "drop"
    ) %>%
    arrange(chrom_length) %>%
    mutate(seqid = factor(seqid, levels = seqid))
}

# Feature extraction (unified)
#------------------------------------------------------------------------

extract_fill_features <- function(file) {
  read_tsv(file, col_names = FALSE, show_col_types = FALSE) %>%
    setNames(c("seqid","start","end","category")) %>%
    mutate(
      mid_position = (start + end) / 2,
      source_type  = "fill"
    )
}

extract_keyword_features <- function(gff, keywords) {
  attrs <- keywords[seq(1, length(keywords), 2)] # odd indices: attributes
  types <- keywords[seq(2, length(keywords), 2)] # even indices: types

  bind_rows(lapply(seq_along(attrs), function(i) { # for each attribute/type pair
    gff %>%
      filter(
        grepl(attrs[i], attributes, ignore.case = TRUE),
        grepl(types[i], type, ignore.case = TRUE)
      ) %>%
      mutate(
        # Force category to be the family name (attrs[i])
        category = attrs[i], # assign category
        mid_position = (start + end) / 2, # midpoint ubication
        source_type = "keyword" # mark source type
      )
  })) # end bind_rows
}

build_features <- function(gff, chrom_limits, opts) {

  df <- if (opts$feature_mode == "fill") { # fill mode
    extract_fill_features(opts$feature_file) # load fill features
  } else { # keyword mode
    extract_keyword_features(gff, opts$keywords) # load keyword features
  }

  df <- df %>%
    filter(seqid %in% chrom_limits$seqid) # keep only relevant chromosomes

  if (opts$strict_filter) {
    keep <- unique(df$seqid) # chromosomes with features
    df <- df %>% filter(seqid %in% keep) # strict filtering
  }

  df # return features
}

# Plotting functions
#------------------------------------------------------------------------

plot_chromosomes <- function(chrom_limits) {
  ggplot(chrom_limits) +
    geom_segment(
      aes(x = chrom_start, xend = chrom_end,
          y = seqid, yend = seqid),
      color = "grey50") +
    theme_minimal() +
    theme(
      # Force white background for the panel
      panel.background = element_rect(fill = "white", color = NA),
      # Force white background for the plot area
      plot.background = element_rect(fill = "white", color = NA))
}

add_feature_points <- function(p, features, colors) {
  p +
    geom_point(
      data = features,
      aes(x = mid_position, y = seqid, color = category),
      size = 1.4
    ) +
    scale_color_manual(values = colors)
}

add_feature_segments <- function(p, features, colors) {
  p +
    geom_segment(
      data = features,
      aes(x = start, xend = end,
          y = seqid, yend = seqid,
          color = category),
      size = 1.2
    ) +
    scale_color_manual(values = colors)
}

  plot_accumulated <- function(
    features,
    chrom_limits,
    colors,
    mode = "hist",
    binwidth = 0.05,
    # bins = bins,
    alpha = 0.6
  ) {

  df <- features %>%
    left_join(
      chrom_limits %>% select(seqid, chrom_start, chrom_length),
      by = "seqid"
    ) %>%
    mutate(relative_pos = (mid_position - chrom_start) / chrom_length)

  if (mode %in% c("overlay", "density")) {

    ggplot(df, aes(x = relative_pos, color = category, fill = category)) +
      geom_density(alpha = alpha) +
      # geom_density(aes(y = after_stat(scaled)), alpha = alpha) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      coord_cartesian(xlim = c(0, 1))

  } else {  # hist por defecto

    ggplot(df, aes(x = relative_pos, fill = category)) +
      geom_histogram(
        binwidth = binwidth,
        # bins = bins,
        alpha = alpha,
        position = "identity"
      ) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      coord_cartesian(xlim = c(0, 1))
  }
}

plot_faceted_accumulated <- function(
  features,
  chrom_limits,
  colors,
  geom = c("density", "hist"),
  scales = c("fixed", "free"),
  binwidth = 0.05,
  alpha = 0.6
) {

  geom   <- match.arg(geom)
  scales <- match.arg(scales)

    df <- features %>%
      left_join(
        chrom_limits %>% select(seqid, chrom_start, chrom_length),
        by = "seqid"
      ) %>%
      mutate(relative_pos = (mid_position - chrom_start) / chrom_length)

    # ---- N per category (for annotation) ----
    n_df <- df %>%
    count(category) %>%
    mutate(label = paste0("n = ", n))

    p <- ggplot(df, aes(x = relative_pos, fill = category, color = category))

    if (geom == "density") {

  p <- p +
    geom_density(
      # aes(y = after_stat(count)),
      alpha = alpha,
      linewidth = 0.8
    ) +
      geom_text(
        data = n_df,
        aes(
          x = Inf,
          y = Inf,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 1.1,
        vjust = 1.3,
        size = 3.2,
        color = "black"
      )

  ylab <- "Number of features"
} else {

      p <- p +
        geom_histogram(
          binwidth = binwidth,
          alpha = alpha,
          position = "identity",
          linewidth = 0.8
        ) +
        geom_text(
          data = n_df,
          aes(
            x = Inf,
            y = Inf,
            label = label
          ),
          inherit.aes = FALSE,
          hjust = 1.1,
          vjust = 1.3,
          size = 3.2,
          color = "black"
        )

      ylab <- "Count"
    }

    p +
      facet_wrap(~ category, scales = scales) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      coord_cartesian(xlim = c(0, 1)) +
      labs(
        x = "Relative chromosomal position",
        y = ylab
      ) +
      theme_minimal() +
      theme(panel.grid = element_blank(), strip.background = element_blank()) # remove grid and facet background for cleaner look
  }

# Windows count mode functions
#------------------------------------------------------------------------

run_window_mode <- function(opts, features) {
  message("→ Running Clustering Mode (Distance based)")
    
  target_genes <- if(!is.null(opts$gene_list)) opts$gene_list else unique(features$category)
    
  clusters <- features %>%
    filter(category %in% target_genes) %>%
    arrange(seqid, start) %>%
    group_by(seqid) %>%
    mutate(
      dist = start - lag(end),
      is_new = is.na(dist) | dist > opts$window_size,
      cluster_id = cumsum(is_new)
    ) %>%
    group_by(seqid, cluster_id) %>%
    summarise(
      start = min(start), end = max(end),
      total_features = n(),
      genes = paste(unique(category), collapse=","), .groups = "drop"
    ) %>%
    filter(total_features >= opts$min_genes)

  write_tsv(clusters, "genomic_clusters_report.tsv")
  message("-> Found ", nrow(clusters), " clusters. Saved to genomic_clusters_report.tsv")

  return(clusters)
}

# Main execution
#------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opts <- parse_args_manual(args)
validate_args(opts)

# Load GFF data
gff_data <- load_gff(opts$gff_file)
# Compute chromosome limits
chrom_limits <- compute_chrom_limits(gff_data)
# Build features based on mode
features <- build_features(gff_data, chrom_limits, opts)

# Aplly filtering based on arguments
#------------------------------------------------------------------------
# First: only apply --strict filter
if (opts$strict_filter) {
  chrom_limits <- chrom_limits %>%
    filter(seqid %in% unique(features$seqid))
}

# Second: only apply number of chromosome
if (!is.infinite(opts$max_chromosomes)) {
  chrom_limits <- chrom_limits %>%
    arrange(desc(chrom_length)) %>%
    head(opts$max_chromosomes)
} 

# third: apply order file if provided
# the most bigger chromosome is the first one, and the smaller is the last one
# Re-order data frames the last one is most smaller chromosome, and the first one is the most bigger chromosome
chrom_limits <- chrom_limits %>%
  arrange(chrom_length)

features <- features %>%
  filter(seqid %in% chrom_limits$seqid)
# apply factor levels to ensure correct ordering in plots
chrom_limits$seqid <- factor(chrom_limits$seqid, levels = chrom_limits$seqid)
features$seqid     <- factor(features$seqid,     levels = levels(chrom_limits$seqid))

if (opts$window_mode) {
  clusters <- run_window_mode(opts, features)

  if (nrow(clusters) > 0) {
    message("→ Generating cluster plot...")
    # Create base plot for clusters (similar to chromosome plot)
    clusters$seqid <- factor(clusters$seqid, levels = levels(chrom_limits$seqid))

    # p_clusters <- plot_clusters(clusters)
    # ggsave("genomic_clusters_map.png", p_clusters, width = 10, height = 8, dpi = 300)
  }

  quit(save = "no")
}

features <- build_features(gff_data, chrom_limits, opts)
# unique categories for color mapping
categories <- unique(features$category)

if (!is.null(opts$colors)) {
    # Use custom colors provided by the user
    colors <- opts$colors
    # If fewer colors than categories are provided, R will recycle colors automatically
} else if (!is.null(opts$palette)) {
    # Use a palette from RColorBrewer
    colors <- brewer.pal(max(3, length(categories)), opts$palette) # Ensure at least 3 colors for better visualization
} else {
    # Default color set
    colors <- brewer.pal(max(3, length(categories)), "Set1")
}

colors <- head(colors, length(categories)) # Trim colors to match number of categories
names(colors) <- categories # Name colors by category

# Create base plot and add features
# ------------------------------------------------------------------------
base_plot <- plot_chromosomes(chrom_limits)
point_plot <- add_feature_points(base_plot, features, colors)

# Save the base plot with points
ggsave("gene_positions_plot.pdf", point_plot, width = 8, height = 10)
ggsave("gene_positions_plot.png", point_plot, width = 8, height = 10, dpi = 900)

if (opts$line_plot) {
  line_plot <- add_feature_segments(base_plot, features, colors)
  ggsave("gene_positions_line_plot.pdf", line_plot, width = 8, height = 10)
}

if (opts$accumulated_plot) {
  acc_plot <- plot_accumulated(features, chrom_limits, colors, opts$density_mode)
  
  ggsave("gene_distribution_accumulated_plot.pdf", acc_plot, width = 8, height = 5)
  ggsave("gene_distribution_accumulated_plot.png", acc_plot, width = 8, height = 5, dpi = 800)
}

if (opts$facet_plot) {

  facet_plot <- plot_faceted_accumulated(
    features      = features,
    chrom_limits  = chrom_limits,
    colors        = colors,
    geom          = opts$facet_mode,
    scales        = opts$facet_scales
  )

  fname <- paste0(
    "gene_distribution_facet_",
    opts$facet_mode, "_",
    opts$facet_scales
  )

  ggsave(paste0(fname, ".pdf"), facet_plot, width = 8, height = 6)
  ggsave(paste0(fname, ".png"), facet_plot, width = 8, height = 6, dpi = 800)
}