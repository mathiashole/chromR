#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# init variable values
gff_file <- NULL
keyword_pairs <- NULL
layout_id <- NULL
pseudo_data <- NULL
strict <- FALSE
number <- Inf  # Default to Inf, meaning no limit (all chromosomes)
colors_input <- NULL
palette_name <- NULL # select palette color
line_plot <- FALSE
table_format <- NULL
fill_file <- NULL
file_format <- NULL
order_file <- NULL
accumulated_plot <- FALSE
density_mode <- "overlay" # density mode for accumulated plot: "overlay" (default), "stack", "facet"
summary_args <- FALSE
interactive_plot <- FALSE # interactive plot
# try new arguments
window_count_mode <- FALSE
window_size <- NULL
gene_list <- NULL
min_genes <- 2   # default

# Parse arguments manually
for (i in seq_along(args)) {
  if (args[i] == "--gff_file" || args[i] == "-g") {
    gff_file <- args[i + 1]
  } else if (args[i] == "--keywords" || args[i] == "-k") {
    keyword_pairs <- args[(i + 1):length(args)]
  } else if (args[i] == "--number" || args[i] == "-n") {
    number <- as.integer(args[i + 1])
  } else if (args[i] == "--strict" || args[i] == "-s") {
    strict <- TRUE
  } else if (args[i] == "--layout" || args[i] == "-l") {
    layout_id <- args[(i + 1):length(args)]
  } else if (args[i] == "--line_plot" || args[i] == "-lp") {
    line_plot <- TRUE
  } else if (args[i] == "--table" || args[i] == "-tab") {
    table_format <- args[i + 1]
  } else if (args[i] == "--colors" || args[i] == "-c") {
    colors_input <- args[(i + 1):length(args)]
  } else if (args[i] == "--palette" || args[i] == "-p") {
    palette_name <- args[i + 1]
  } else if (args[i] == "--fill_file" || args[i] == "-ff") {
    fill_file <- args[i + 1]
  } else if (args[i] == "--format" || args[i] == "-f") {
    file_format <- args[i + 1]
  } else if (args[i] == "--order_file" || args[i] == "-of") {
    order_file <- args[i + 1]
  } else if (args[i] == "--accumulated_plot" || args[i] == "-ap") {
    accumulated_plot <- TRUE
  } else if (args[i] == "--density_mode" || args[i] == "-dm") {
    density_mode <- args[i + 1]
  } else if (args[i] == "--interactive" || args[i] == "-int") {
    interactive_plot <- TRUE
  } else if (args[i] == "--summary" || args[i] == "-sm") {
    summary_args <- TRUE
  } else   if (args[i] == "--window_count") {
    window_count_mode <- TRUE
  } else if (args[i] == "--window_size") {
    window_size <- as.integer(args[i + 1])
  }
}

# Validate arguments
if (!file.exists(gff_file)) {
  stop("The provided GFF file does not exist.")
}

# Validate and parse colors
if (!is.null(colors_input)) {
  # Ensure the colors are either valid hex codes or named colors
  valid_colors <- colors_input %in% colors() | grepl("^#[A-Fa-f0-9]{6}$", colors_input) ## NEED DEBUG HEX COLORS
  if (all(valid_colors)) {
    custom_colors <- colors_input
  } else {
    stop("Invalid color format detected. Use valid R color names (e.g., 'red') or hex colors like #1f77b4.")
  }
} else if (!is.null(palette_name)) {
  # Option: select Palette of RColorBrewer
  available_palettes <- rownames(RColorBrewer::brewer.pal.info)
  if (!(palette_name %in% available_palettes)) {
    stop(paste("Invalid palette name. Use one from RColorBrewer:", paste(available_palettes, collapse = ", ")))
  }
  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  custom_colors <- RColorBrewer::brewer.pal(max_colors, palette_name)
} else {
  # Default colors if no colors are passed
  custom_colors <- RColorBrewer::brewer.pal(8, "Set1")
}

# Validate output format
if (!is.null(table_format) && !table_format %in% c("csv", "tsv")) {
  stop("Invalid format specified for --table. Use 'csv' or 'tsv'.")
}

#validate window count mode
if (window_count_mode) {
  if (is.null(window_size)) stop("You must provide --window_size.")
  if (is.null(gene_list)) stop("You must provide --genes.")
}

# Function to read fill file
read_fill_file <- function(fill_file, format = NULL) {
  if (is.null(format)) {
    stop("You must specify the format with --format (csv or tsv)")
  }
  format <- tolower(format)
  if (!format %in% c("csv", "tsv")) {
    stop("Invalid format. Use 'csv' or 'tsv'")
  }
  
  if (format == "csv") {
    data <- read_csv(fill_file, col_names = FALSE, show_col_types = FALSE)
  } else {
    data <- read_tsv(fill_file, col_names = FALSE, show_col_types = FALSE)
  }

  # Validate minimum structure
  if (ncol(data) < 4) {
    stop("The file must have at least 4 columns: contig, start, end, category")
  }
  
  # Assign names to the relevant columns
  colnames(data)[1:4] <- c("seqid", "start", "end", "category")
  return(data)
}

read_order_file <- function(order_file, format = NULL) {
  if(!file.exists(order_file)) {
    stop("The provided ORDER file does not exist")
  }

  format <- tolower(format)
  if (!format %in% c("csv", "tsv")) {
    stop("Invalid format. Use 'csv' or 'tsv'")
  }
  
  if (format == "csv") {
    custom_order <- read_csv(order_file, col_names = FALSE, show_col_types = FALSE)
  } else {
    custom_order <- read_tsv(order_file, col_names = FALSE, show_col_types = FALSE)
  }

  return(custom_order)
}

# charge library
library(dplyr)
library(readr)
library(ggplot2)

# Load GFF file
gff_data <- read_tsv(gff_file, comment = "#", col_names = FALSE)
# Rename GFF columns
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Calculate mimits of chromosome (start and end) per group
# Calculate chromosome limits
chrom_limits <- gff_data %>%
  group_by(seqid) %>%
  summarize(
    chrom_start = min(start),
    chrom_end = max(end),
    chrom_length = chrom_end - chrom_start,  # Calcualte chrmosome length
    .groups = "drop"
  ) # %>%
  # arrange(chrom_length) %>% # Order longest to smaller
  # mutate(seqid = factor(seqid, levels = seqid))

if (!is.null(order_file)) {
  custom_order <- read_order_file(order_file)
  chrom_limits <- chrom_limits %>%
    mutate(seqid = factor(seqid, levels = custom_order)) %>%
    arrange(seqid)  # Order as defined
} else {
  chrom_limits <- chrom_limits %>%
    arrange(chrom_length) %>%
    mutate(seqid = factor(seqid, levels = seqid)) # default order
}

# ========================================================================
# Need add enviroment in fill file and filter keyword
# and need clean code
if (!is.null(fill_file)) {
  # Modo fill_file
  if (!file.exists(fill_file)) {
    stop("The specified fill file does not exist")
  }
  
  fill_data <- read_fill_file(fill_file, format = file_format) %>%
    filter(seqid %in% chrom_limits$seqid) %>%
    mutate(
      mid_position = (start + end) / 2,
      seqid = factor(seqid, levels = levels(chrom_limits$seqid))
    )

  accum_data <- fill_data %>%
    left_join(chrom_limits %>% select(seqid, chrom_start, chrom_length), by = "seqid") %>%
    mutate(relative_pos = (mid_position - chrom_start) / chrom_length)
  
  # Apply strict keyword filtering if enabled
  if (strict && !is.null(fill_data)) {
    keyword_contigs <- unique(fill_data$seqid)
    chrom_limits <- chrom_limits %>% filter(seqid %in% keyword_contigs)
    fill_data <- fill_data %>% filter(seqid %in% chrom_limits$seqid)
  }

  # Apply filtering based on 'number' argument
  chrom_limits <- chrom_limits %>%
    arrange(desc(chrom_length)) %>%
    { if (!is.infinite(number)) slice_head(., n = number) else . } %>%
    arrange(chrom_length) %>%
    mutate(seqid = factor(seqid, levels = seqid))

  fill_data <- fill_data %>% # Apply filtering based on 'number' and id
    filter(seqid %in% chrom_limits$seqid)

  # Get unique category attributes
  unique_category <- unique(fill_data$category)
  # Create a named COLOR vector
  color_mapping <- setNames(custom_colors[seq_along(unique_category)], unique_category)

# ---------------------------------------------------------------------------
# Just try fill information of hover mode
fill_data <- fill_data %>%
  left_join(chrom_limits %>% select(seqid, chrom_start, chrom_length), by = "seqid")


  fill_data <- fill_data %>%
    mutate(
      hover_text = paste0(
        "<b>Category:</b> ", category, "<br>",
        "<b>Chromosome:</b> ", seqid, "<br>",
        "<b>Start:</b> ", start, "<br>",
        "<b>End:</b> ", end, "<br>",
        "<b>Length:</b> ", end - start, "<br>",
        "<b>Relative Pos:</b> ", round((mid_position - chrom_start) / chrom_length, 3)
      )
    )

# ---------------------------------------------------------------------------
  
  # Create the base plot object
  plot <- ggplot() +
    geom_segment(
      data = chrom_limits,
      aes(x = chrom_start, xend = chrom_end, y = seqid, yend = seqid),
      color = "gray50", size = 0.8, alpha = 0.8
    ) #+
    # geom_point(
    #   data = fill_data,
    #   aes(x = mid_position, y = seqid, color = factor(category)),
    #   size = 1.5
    # ) +
    # labs(color = "Category")
      # Points for genes
    point_plot <- plot +
      geom_point(
      data = fill_data,
      # aes(x = mid_position, y = seqid, color = factor(category)),
      aes(x = mid_position, y = seqid, color = factor(category), text = hover_text),
      size = 1.5
    ) +
    scale_color_manual(values = color_mapping) +  # Apply custom colors
    # Finalize plot aesthetics
    labs(
      x = "Position on Chromosome",
      y = "Chromosome",
      color = "Categories",  # category label is necesary
      title = "Genomic Features by Category"
    ) +
    theme_minimal()
  
} else {
    # Keywords mode
    if (!is.null(keyword_pairs)) {
    # Split keyword pairs into `attributes` and `type`
    keywords_attr <- keyword_pairs[seq(1, length(keyword_pairs), by = 2)]  # Odd indices: attributes
    keywords_type <- keyword_pairs[seq(2, length(keyword_pairs), by = 2)]  # Even indices: types

    # Get unique keyword attributes
    unique_keywords <- unique(keywords_attr)
    # Create a named COLOR vector
    color_mapping <- setNames(custom_colors[seq_along(unique_keywords)], unique_keywords)

    # ========================================================================
    # Filter GFF data using the keyword pairs
    filtered_data <- lapply(seq_along(keywords_attr), function(i) {
      gff_data %>%
        filter(
          grepl(keywords_attr[i], attributes),  # Filter by attribute keyword
          grepl(keywords_type[i], type)        # Filter by type keyword
        ) %>%
        mutate(
          mid_position = (start + end) / 2,  # Calculate midpoint
          keyword_attr = keywords_attr[i],  # Add attribute keyword as a column
          keyword_type = keywords_type[i],  # Add type keyword as a column
          seqid = factor(seqid)             # Ensure seqid is a factor
        )
    }) %>%
      bind_rows()  # Combine all filtered data

    accum_data <- filtered_data %>%
      left_join(chrom_limits %>% select(seqid, chrom_start, chrom_length), by = "seqid") %>%
      mutate(relative_pos = (mid_position - chrom_start) / chrom_length)

    # Apply strict keyword filtering if enabled
    if (strict && !is.null(filtered_data)) {
      keyword_contigs <- unique(filtered_data$seqid)
      chrom_limits <- chrom_limits %>% filter(seqid %in% keyword_contigs)
      filtered_data <- filtered_data %>% filter(seqid %in% chrom_limits$seqid)
    }

    # Apply filtering based on 'number' argument
    chrom_limits <- chrom_limits %>%
      arrange(desc(chrom_length)) %>%
      { if (!is.infinite(number)) slice_head(., n = number) else . } %>%
      arrange(chrom_length) %>%
      mutate(seqid = factor(seqid, levels = seqid))

    filtered_data <- filtered_data %>% # Apply filtering based on 'number' and id
      filter(seqid %in% chrom_limits$seqid)

# ---------------------------------------------------------------------------
# Just try fill information of hover mode
filtered_data <- filtered_data %>%
  left_join(chrom_limits %>% select(seqid, chrom_start, chrom_length), by = "seqid")

  filtered_data <- filtered_data %>%
    mutate(
      hover_text = paste0(
        "<b>Keyword:</b> ", keyword_attr, "<br>",
        "<b>Type:</b> ", keyword_type, "<br>",
        "<b>Chromosome:</b> ", seqid, "<br>",
        "<b>Start:</b> ", start, "<br>",
        "<b>End:</b> ", end, "<br>",
        "<b>Length:</b> ", end - start, "<br>",
        "<b>Mid Pos:</b> ", mid_position
      )
    )

# ---------------------------------------------------------------------------

    #make plot
    plot <- ggplot() +
      # Chromosome lines
      geom_segment(
        data = chrom_limits,
        aes(x = chrom_start, xend = chrom_end, y = seqid, yend = seqid),
        color = "gray50", size = 0.8, alpha = 0.8
      )# +

      # Points for genes
      point_plot <- plot +
        geom_point(
        data = filtered_data,
        # aes(x = mid_position, y = seqid, color = factor(keyword_attr)),
        aes(x = mid_position, y = seqid, color = factor(keyword_attr), text = hover_text),
        size = 1.5
      ) +
      scale_color_manual(values = color_mapping) +  # Apply custom colors
      # Finalize plot aesthetics
      labs(
        x = "Position on Chromosome",
        y = "Chromosome",
        color = "Keywords",  # keyword label
        title = "Gene Positions by Keywords on Chromosomes"
      ) +
      theme_minimal() #+
    #   # theme_classic()

  } else {
    stop("Either --keywords or --fill_file must be specified")
  }
}

# Process pseudogenes if file is provided
if (!is.null(layout_id)) {
  pseudogenes <- read_tsv(layout_id, col_names = FALSE)
  pseudo_pattern <- paste(pseudogenes[[1]], collapse = "|")
  
  if (!is.null(fill_file)) {
    pseudo_data <- fill_data %>%
      filter(grepl(pseudo_pattern, category)) %>%
      mutate(seqid = factor(seqid, levels = chrom_limits$seqid))
  } else {
    pseudo_data <- filtered_data %>%
      filter(grepl(pseudo_pattern, attributes)) %>%
      mutate(seqid = factor(seqid, levels = chrom_limits$seqid))
  }
  
  point_plot <- plot +
    geom_point(
      data = pseudo_data,
      aes(x = mid_position, y = seqid),
      color = "azure4", size = 1.6
    )
}

# Add gene representation
if (line_plot) {
  if(!is.null(fill_file)){
    line_plot <- plot +
      geom_segment(
        data = fill_data,
        aes(x = start, xend = end, y = seqid, yend = seqid, color = factor(category)),
        size = 1.2
      ) +
    scale_color_manual(values = color_mapping) +  # Apply custom colors
    # Finalize plot aesthetics
    labs(
      x = "Position on Chromosome",
      y = "Chromosome",
      color = "Categories",  # category label is necesary
      title = "Genomic Features by Category"
    ) +
    theme_minimal()

  } else {
  # Add genes as lines
  line_plot <- plot +
    geom_segment(
      data = filtered_data,
      aes(x = start, xend = end, y = seqid, yend = seqid, color = factor(keyword_attr)),
      size = 1.2
    ) +
  scale_color_manual(values = color_mapping) +  # Apply custom colors
  # Finalize plot aesthetics
  labs(
    x = "Position on Chromosome",
    y = "Chromosome",
    color = "Keywords",  # keyword label
    title = "Gene Positions by Keywords on Chromosomes"
  ) +
  theme_minimal()
  }
  
  if (!is.null(pseudo_data)) {
  line_plot <- line_plot +
    geom_segment(
      data = pseudo_data,
      aes(x = start, xend = end, y = seqid, yend = seqid),
      color = "azure4", size = 1.6
    )
  }

  # Save plot picture
  line_plot_file <- "gene_positions_line_plot.pdf"
  line_plot_file_png <- "gene_positions_line_plot.png"  # PNG format
  ggsave(line_plot_file_png, plot = line_plot, width = 8, height = 10, dpi = 600)
  ggsave(line_plot_file, plot = line_plot, width = 8, height = 10)
  cat("Plot saved to:", line_plot_file, "\n")

}

if (accumulated_plot) {
  # determine which column holds the category/key
  category_var <- if (!is.null(fill_file)) "category" else "keyword_attr"

  # acc_plot <- ggplot(accum_data, aes(x = relative_pos, color = !!sym(category_var))) +
  #   stat_ecdf(size = 1) +
  #   scale_color_manual(values = color_mapping) +
  #   labs(
  #     x = "Relative position",
  #     y = "Cumulative proportion of genes",
  #     color = "Category"
  #   ) +
  #   theme_minimal()
  
  accum_data[[category_var]] <- as.factor(accum_data[[category_var]])

  # Validate density_mode
  density_mode <- tolower(density_mode)
  if (!density_mode %in% c("overlay", "stack", "facet")) {
    warning("Invalid density_mode. Falling back to 'overlay'.")
    density_mode <- "overlay"
  }

if (density_mode == "overlay") {
    # overlay = density curves (semi-transparent) + ECDF option
    acc_plot <- ggplot(accum_data, aes_string(x = "relative_pos", color = category_var, fill = category_var)) +
      geom_density(alpha = 0.35, adjust = 1) +
      labs(
        x = "Relative position",
        y = "Density (normalized across chromosomes)",
        color = ifelse(!is.null(fill_file), "Category", "Keyword"),
        fill = ifelse(!is.null(fill_file), "Category", "Keyword"),
        title = "Normalized distribution (overlay) across chromosomes"
      ) +
      scale_color_manual(values = color_mapping) +
      scale_fill_manual(values = color_mapping) +
      theme_minimal()
  } else if (density_mode == "stack") {
    # stacked = densities stacked (like ridge/stacked area)
    acc_plot <- ggplot(accum_data, aes_string(x = "relative_pos", fill = category_var)) +
      geom_density(position = "stack", alpha = 0.9) +
      labs(
        x = "Relative position",
        y = "Stacked density (proportion)",
        fill = ifelse(!is.null(fill_file), "Category", "Keyword"),
        title = "Normalized distribution (stacked) across chromosomes"
      ) +
      scale_fill_manual(values = color_mapping) +
      theme_minimal()
  } else { # facet
    # facet = one density per facet (one row per chromosome or per category)
    # If you want per-chromosome facets use seqid; here we facet by category to compare shapes separately.
    acc_plot <- ggplot(accum_data, aes_string(x = "relative_pos")) +
      geom_density(fill = "grey60", alpha = 0.6) +
      facet_wrap(as.formula(paste("~", category_var)), scales = "free_y", ncol = 1) +
      labs(
        x = "Relative position",
        y = "Density",
        title = "Normalized distribution (faceted by category)"
      ) +
      theme_minimal() +
      theme(strip.text = element_text(size = 8))
  }

  # ggplot(accum_data, aes(x = relative_pos, fill = !!sym(category_var))) +
  #   geom_density(alpha = 0.6) +
  #   scale_fill_manual(values = color_mapping) +
  #   labs(
  #     x = "Relative position on chromosome (0-1)",
  #     y = "Density",
  #     fill = ifelse(!is.null(fill_file), "Category", "Keyword"),
  #     title = "Normalized gene distribution across chromosomes"
  #   ) +
  #   theme_minimal()

  ggsave("gene_distribution_accumulated_plot.pdf", acc_plot, width = 8, height = 5)
  ggsave("gene_distribution_accumulated_plot.png", acc_plot, width = 8, height = 5, dpi = 600)
  cat("Accumulated plot saved as gene_distribution_accumulated_plot.pdf/png\n")
}

output_prefix <- "Feature_chromosome"

if (summary_args) {
  summary_data <- if (!is.null(fill_file)) {
    fill_data %>%
      group_by(seqid, category) %>%
      summarise(count = n(), .groups = "drop")
  } else if (!is.null(keyword_pairs)) {
    filtered_data %>%
      group_by(seqid, keyword_attr) %>%
      summarise(count = n(), .groups = "drop")
  }

  write_tsv(summary_data, paste0(output_prefix, "_summary.tsv"))
  cat("Summary table saved to", paste0(output_prefix, "_summary.tsv"), "\n")
}

# Need stadistic information of position count 
# Save tables in the specified format without column names
if (!is.null(table_format)) {
  if (table_format == "csv") {
    write_csv(chrom_limits, "chrom_limits.csv", col_names = FALSE)
    # write_csv(filtered_data, "filtered_data.csv", col_names = FALSE)
    if (!is.null(fill_file)) {
      write_csv(fill_data, "fill_data.csv", col_names = FALSE)
    } else if (!is.null(keyword_pairs)) {
      write_csv(filtered_data, "filtered_data.csv", col_names = FALSE)
    }
  } else if (table_format == "tsv") {
    write_tsv(chrom_limits, "chrom_limits.tsv", col_names = FALSE)
    # write_tsv(filtered_data, "filtered_data.tsv", col_names = FALSE)
    if (!is.null(fill_file)) {
      write_tsv(fill_data, "fill_data.tsv", col_names = FALSE)
    } else if (!is.null(keyword_pairs)) {
      write_tsv(filtered_data, "filtered_data.tsv", col_names = FALSE)
    }
  } else {
    stop("Unsupported table format specified.")
  }
  cat("Tables saved in", table_format, "format without column names.\n")
}

# Save plot picture
plot_file <- "gene_positions_plot.pdf"
plot_file_png <- "gene_positions_plot.png"  # PNG format
ggsave(plot_file_png, plot = point_plot, width = 8, height = 10, dpi = 600)
ggsave(plot_file, plot = point_plot, width = 8, height = 10)
cat("Plot saved to:", plot_file, "\n")

if (interactive_plot) {
  suppressPackageStartupMessages({
    library(plotly)
    library(htmlwidgets)
  })

  interactive_file <- "gene_positions_plot.html"
  p_interactive <- ggplotly(point_plot, tooltip = "text")
  saveWidget(p_interactive, interactive_file, selfcontained = TRUE)
  cat("Interactive plot saved to:", interactive_file, "\n")

}

