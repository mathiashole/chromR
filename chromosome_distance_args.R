#!/usr/bin/env Rscript

# charge library
library(dplyr)
library(readr)
library(ggplot2)

# read arguments
args <- commandArgs(trailingOnly = TRUE)

# arguments validation
if (length(args) < 3) {
  stop("Usage: Rscript plot_dgf1_genes.R <gff_file> <keyword1> <keyword2>")
}

# global variable
gff_file <- args[1]
keyword1 <- args[2]
keyword2 <- args[3]

# charge gff file
gff_data <- read_tsv(gff_file, comment = "#", col_names = FALSE)

# change column name of gff
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Calcular los límites del cromosoma (inicio y final) por grupo
chrom_limits <- gff_data %>%
  group_by(seqid) %>%
  summarize(
    chrom_start = min(start),
    chrom_end = max(end),
    .groups = "drop"
  )

# sequence filter both keywords
filtered_data <- gff_data %>%
  filter(grepl(keyword1, attributes)) %>%
  filter(grepl(keyword2, type)) %>%
  mutate(
    mid_position = (start + end) / 2  # calculate mean position
  )

# Combinar los datos de límites del cromosoma con los datos filtrados
result <- filtered_data %>%
  left_join(chrom_limits, by = "seqid") %>%
  mutate(
    distance_to_chrom_start = start - chrom_start,  # Distancia del inicio del gen al inicio del cromosoma
    distance_to_chrom_end = chrom_end - end        # Distancia del final del gen al final del cromosoma
  )

filtered_result <- result %>%
  select(seqid, start, distance_to_chrom_start, distance_to_chrom_end)

# output_file <- "gene_positions.csv"
# write_csv(filtered_result, output_file)
# cat("Results saved to:", output_file, "\n")

# make plot
ggplot() +
  # Add line of chromosome
  geom_segment(
    data = chrom_limits,
    aes(x = chrom_start, xend = chrom_end, y = seqid, yend = seqid),
    color = "gray50", size = 0.8
  ) +
  # Add gene point to chromosome line
  geom_point(
    data = filtered_data,
    aes(x = mid_position, y = seqid),
    color = "red", size = 1.5
  ) +
  # fix label
  labs(
    x = "Position on Chromosome",
    y = "Chromosome",
    title = paste("Positions of", keyword1, "Genes on Chromosomes")
  ) +
  theme_minimal()

# Save plot picture
plot_file <- "gene_positions_plot.pdf"
# ggsave(plot_file, width = 10, height = 6)
ggsave(plot_file,width = 8, height = 10)
cat("Plot saved to:", plot_file, "\n")
