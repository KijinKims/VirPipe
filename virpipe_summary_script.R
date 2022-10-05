# A script to generate the pipeline result report using R markdown.
# 
# Author: Kijin Kim (skkujin@gmail.com)
###############################################################

suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("knitr"))

# create parser object
p <- arg_parser("Virpipe summary report generating script")
  
# --- Script parameter parsing ---

p <- add_argument(p, "--prefix", help="Sample prefix")

p <- add_argument(p, "--markdown", help="Rmarkdown template path")

p <- add_argument(p, "--qc_html", help="QC html file path")

p <- add_argument(p, "--map_summary", help="Reference map summary file path")

p <- add_argument(p, "--tax_classify_directory", help="Directory where taxonomic classifcation result is saved")

p <- add_argument(p, "--tax_classify_kraken_html", help="Kraken result html file path")

p <- add_argument(p, "--assembly_directory", help="Directory where assembly result is saved")

p <- add_argument(p, "--blastn_table", help="Blastn result file path")

p <- add_argument(p, "--blastn_html", help="Blastn result in functional table")

p <- add_argument(p, "--megablast_table", help="Megablast result file path")

p <- add_argument(p, "--megablast_html", help="Megablast result in functional table")

p <- add_argument(p, "--zoonotic_rank_predictions", help="Zoonotic rank result predictions file path")

# get command line options,
# if options not found on command line then set defaults automatically generated with prefix, 
args <- parse_args(p)
prefix <- args$prefix

# qc data
if(is.na(args$qc_html)){
  qc_link <- paste(prefix, "/qc/", prefix, ".qc.html", sep="")
} else{
  qc_link <- args$qc_html
}

# reference map data

if(is.na(args$map_summary)){
  map_summary_file <- paste(prefix, "/map/", prefix, ".map.filtered.tsv", sep="")
} else{
  map_summary_file <- args$map_summary
}

if(is.na(args$tax_classify_directory)){
  tax_classify_directory <- paste(prefix, "/tax_classify/", sep = "")
} else{
  tax_classify_directory <- args$tax_classify_directory
}

if(is.na(args$tax_classify_kraken_html)){
  tax_classify_kraken_link <- paste(prefix, "/tax_classify/", prefix, ".kraken.html", sep = "")
} else{
  tax_classify_kraken_link <- args$tax_classify_kraken_html
}

if(is.na(args$assembly_directory)){
  assembly_directory <- paste(prefix, "/assembly/", sep = "")
} else{
  assembly_directory <- args$assembly_directory
}

# blast data
if(is.na(args$blastn_table)){
  blast_blastn_file <- paste(prefix, "/post_assembly/blast/", prefix, ".blastn.filtered.txt", sep="")
} else{
  blast_blastn_file <- args$blastn_table
}

if(is.na(args$blastn_html)){
  blastn_html_link <- paste(prefix, "/post_assembly/blast/", prefix, ".blastn.filtered.html", sep="")
} else{
  blastn_html_link <- args$blastn_html
}

if(is.na(args$megablast_table)){
  blast_megablast_file <- paste(prefix, "/post_assembly/blast/", prefix, ".megablast.filtered.txt", sep="")
} else{
  blast_megablast_file <- args$megablast_table
}

if(is.na(args$megablast_html)){
  megablast_html_link <- paste(prefix, "/post_assembly/blast/", prefix, ".megablast.filtered.html", sep="")
} else{
  megablast_html_link <- args$megablast_html
}

# zoonotic rank data
if(is.na(args$zoonotic_rank_predictions)){
  zoonotic_rank_predictions_file <- paste(prefix, "/post_assembly/zoonotic_rank/", prefix, ".predictions.csv", sep="")
} else{
  zoonotic_rank_predictions_file <- args$zoonotic_rank_predictions
}

rmd_template_file <- args$markdown
render(rmd_template_file,
      output_file=paste(prefix,".virpipe_summary", sep = ""),
      output_dir='.',
      knit_root_dir=getwd(),
      params=list(prefix=prefix,
                  qc_link=qc_link, 
                  map_summary_file=map_summary_file,
                  tax_classify_directory=tax_classify_directory,
                  tax_classify_kraken_link=tax_classify_kraken_link, 
                  assembly_directory=assembly_directory,
                  blast_blastn_file=blast_blastn_file,
		              blastn_html_link=blastn_html_link,
                  blast_megablast_file=blast_megablast_file,
		              megablast_html_link=megablast_html_link,
                  zoonotic_rank_predictions_file=zoonotic_rank_predictions_file))
