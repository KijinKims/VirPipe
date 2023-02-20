#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

require(MetaComp)

# load kraken result
assignment <- load_edge_assignment(args[1], type = "kraken")

# generate family-level heatmap
plot_merged_assignment(assignment, "family", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "family.svg", sep = ""))
# generate genus-level heatmap
plot_merged_assignment(assignment, "genus", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "genus.svg", sep = ""))
# generate species-level heatmap
plot_merged_assignment(assignment, "species", row_limit = 10, min_row_abundance = 1.0e-10, plot_title = args[2], filename = paste(args[2], "_", "species.svg", sep = ""))
