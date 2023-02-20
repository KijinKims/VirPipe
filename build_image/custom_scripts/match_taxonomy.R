args = commandArgs(trailingOnly=TRUE)

require(data.table)
require(dplyr)
require(taxonomizr)
require(DT)
require(htmlwidgets)

data <- as.data.frame(fread(args[1]))

uniq_taxid <- unique(data$"TAX_ID")

# find out all taxons of appearing taxid
taxons <- as.data.frame(getTaxonomy(uniq_taxid, args[3]))
taxons <- tibble::rownames_to_column(taxons, "TAX_ID")
taxons$TAX_ID <- as.integer(taxons$TAX_ID)

# match taxon name to data
taxon_mapped <- inner_join(data, taxons, by = "TAX_ID")

# export table
# write.table(taxon_mapped, file = args[2], row.names=FALSE, col.names = TRUE, sep="\t", quote = FALSE)

header_order <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "TAX_ID", "REF_ID", "REF_TITLE", "QUERY_ID", "REF_ID", "EVALUE", "BITSCORE", "PER_IDENT", "QUERY_LEN", "ALN_LEN", "MISMATCH", "GAPOPEN", "QUERY_START", "QUERY_END", "REF_START", "REF_END")


ordered <- taxon_mapped[, header_order]
dt <- datatable(ordered, 
                filter = 'top', 
                extensions = 'Buttons', 
                options = list(autoWidth = TRUE,
                               paging = TRUE,
                               searching = TRUE,
                               fixedColumns = TRUE,
                               ordering = TRUE,
                               dom = 'iftprB',
                               buttons = c('copy', 'csv', 'excel')),
                               class = 'display')
saveWidget(dt, args[2], selfcontained=FALSE)
