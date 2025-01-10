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
line_plot <- FALSE
table_format <- NULL
# export_csv <- FALSE
# export_tsv <- FALSE

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
    break
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
} else {
  # Default colors if no colors are passed
  custom_colors <- RColorBrewer::brewer.pal(8, "Set1")
}

# Validate output format
if (!is.null(table_format) && !table_format %in% c("csv", "tsv")) {
  stop("Invalid format specified for --table. Use 'csv' or 'tsv'.")
}

# Split keyword pairs into `attributes` and `type`
keywords_attr <- keyword_pairs[seq(1, length(keyword_pairs), by = 2)]  # Odd indices: attributes
keywords_type <- keyword_pairs[seq(2, length(keyword_pairs), by = 2)]  # Even indices: types

# Get unique keyword attributes
unique_keywords <- unique(keywords_attr)
# Create a named COLOR vector
color_mapping <- setNames(custom_colors[seq_along(unique_keywords)], unique_keywords)

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
  )  %>%
  arrange(chrom_length)  # Order longest to smaller

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

# Apply strict keyword filtering if enabled
if (strict && !is.null(filtered_data)) {
  keyword_contigs <- unique(filtered_data$seqid)
  chrom_limits <- chrom_limits %>%
    filter(seqid %in% keyword_contigs)
  filtered_data <- filtered_data %>%
    filter(seqid %in% chrom_limits$seqid)
}

# Apply filtering based on 'number' argument
chrom_limits <- chrom_limits %>%
  arrange(desc(chrom_length)) %>%
  {
    if (!is.infinite(number)) slice_head(., n = number) else . 
  } %>%
  arrange(chrom_length) %>%
  mutate(seqid = factor(seqid, levels = seqid))

filtered_data <- filtered_data %>% # Apply filtering based on 'number' and id
  filter(seqid %in% chrom_limits$seqid)

# Process pseudogenes if file is provided
if (!is.null(layout_id)) {
  pseudogenes <- read_tsv(layout_id, col_names = FALSE)
  pseudo_pattern <- paste(pseudogenes[[1]], collapse = "|")
  
  pseudo_data <- filtered_data %>%
    filter(grepl(pseudo_pattern, attributes)) %>%
    mutate(
      mid_position = (start + end) / 2,
      seqid = factor(seqid, levels = chrom_limits$seqid)
    )
 }

# make plot
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
    aes(x = mid_position, y = seqid, color = factor(keyword_attr)),
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
  # theme_classic()

if (!is.null(pseudo_data)) {
  point_plot <- point_plot +
    geom_point(
      data = pseudo_data,
      aes(x = mid_position, y = seqid),
      color = "azure4", size = 1.6
    )
}

# Add gene representation
if (line_plot) {
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
  # ggsave(line_plot_file, width = 10, height = 6)
  line_plot_file_png <- "gene_positions_line_plot.png"  # PNG format
  ggsave(line_plot_file_png, plot = line_plot, width = 8, height = 10, dpi = 600)
  ggsave(line_plot_file, plot = line_plot, width = 8, height = 10)
  cat("Plot saved to:", line_plot_file, "\n")

}  

# # Export tables if requested
# if (export_csv) {
#   write_csv(chrom_limits, "chrom_limits.csv", col_names = FALSE)
#   write_csv(filtered_data, "filtered_data.csv", col_names = FALSE)
#   cat("Tables exported as CSV files: chrom_limits.csv and filtered_data.csv\n")
# }

# if (export_tsv) {
#   write_tsv(chrom_limits, "chrom_limits.tsv", col_names = FALSE)
#   write_tsv(filtered_data, "filtered_data.tsv", col_names = FALSE)
#   cat("Tables exported as TSV files: chrom_limits.tsv and filtered_data.tsv\n")
# }

# Save tables in the specified format without column names
if (!is.null(table_format)) {
  if (table_format == "csv") {
    write_csv(chrom_limits, "chrom_limits.csv", col_names = FALSE)
    write_csv(filtered_data, "filtered_data.csv", col_names = FALSE)
  } else if (table_format == "tsv") {
    write_tsv(chrom_limits, "chrom_limits.tsv", col_names = FALSE)
    write_tsv(filtered_data, "filtered_data.tsv", col_names = FALSE)
  }
  cat("Tables saved in", table_format, "format without column names.\n")
}

# Save plot picture
plot_file <- "gene_positions_plot.pdf"
# ggsave(plot_file, width = 10, height = 6)
plot_file_png <- "gene_positions_plot.png"  # PNG format
ggsave(plot_file_png, plot = point_plot, width = 8, height = 10, dpi = 600)
ggsave(plot_file, plot = point_plot, width = 8, height = 10)
cat("Plot saved to:", plot_file, "\n")
