#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# init variable values
gff_file <- NULL
keyword_pairs <- NULL
list_id <- NULL
pseudo_data <- NULL
strict <- FALSE
number <- Inf  # Default to Inf, meaning no limit (all chromosomes)

# Parse arguments manually
for (i in seq_along(args)) {
  if (args[i] == "--gff_file" || args[i] == "-g") {
    gff_file <- args[i + 1]
  } else if (args[i] == "--keywords" || args[i] == "-k") {
    keyword_pairs <- args[(i + 1):length(args)]
    #break
  } else if (args[i] == "--number" || args[i] == "-n") {
    number <- as.integer(args[i + 1])
  } else if (args[i] == "--strict" || args[i] == "-s") {
    strict <- TRUE
  } else if (args[i] == "--list" || args[i] == "-l") {
    list_id <- args[(i + 1):length(args)]
    break
  }
}

# Validate arguments
if (!file.exists(gff_file)) {
  stop("The provided GFF file does not exist.")
}

# Split keyword pairs into `attributes` and `type`
keywords_attr <- keyword_pairs[seq(1, length(keyword_pairs), by = 2)]  # Odd indices: attributes
keywords_type <- keyword_pairs[seq(2, length(keyword_pairs), by = 2)]  # Even indices: types

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

# Order seqid
chrom_limits <- chrom_limits %>%
  arrange(chrom_length) %>%
  slice_head(n = number) %>% # Apply filtering based on 'number' argument  ####### NEED DEBUGGING ########
  mutate(seqid = factor(seqid, levels = seqid))  # Reordenar levels

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

filtered_data <- filtered_data %>% # Apply filtering based on 'number' and id
  filter(seqid %in% chrom_limits$seqid)

# Apply strict keyword filtering if enabled
if (strict && !is.null(filtered_data)) {
  keyword_contigs <- unique(filtered_data$seqid)
  chrom_limits <- chrom_limits %>%
    filter(seqid %in% keyword_contigs)
  filtered_data <- filtered_data %>%
    filter(seqid %in% chrom_limits$seqid)
}

# Process pseudogenes if file is provided
if (!is.null(list_id)) {
  pseudogenes <- read_tsv(list_id, col_names = FALSE)
  pseudo_pattern <- paste(pseudogenes[[1]], collapse = "|")
  
  pseudo_data <- filtered_data %>%
    filter(grepl(pseudo_pattern, attributes)) %>%
    mutate(
      mid_position = (start + end) / 2,
      seqid = factor(seqid, levels = chrom_limits$seqid)
    )
 } #else {
#   pseudo_data <- NULL
# }

# make plot
plot <- ggplot() +
  # Chromosome lines
  geom_segment(
    data = chrom_limits,
    aes(x = chrom_start, xend = chrom_end, y = seqid, yend = seqid),
    color = "gray50", size = 0.8
  ) +
  # Points for genes
  geom_point(
    data = filtered_data,
    aes(x = mid_position, y = seqid, color = factor(keyword_attr)),
    size = 1.5
  ) +
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
  plot <- plot +
    geom_point(
      data = pseudo_data,
      aes(x = mid_position, y = seqid),
      color = "azure4", size = 1.6
    )
}

# Save plot picture
plot_file <- "gene_positions_plot.pdf"
# ggsave(plot_file, width = 10, height = 6)
plot_file_png <- "gene_positions_plot.png"  # PNG format
ggsave(plot_file_png, plot = plot, width = 8, height = 10, dpi = 600)
ggsave(plot_file, plot = plot, width = 8, height = 10)
cat("Plot saved to:", plot_file, "\n")
