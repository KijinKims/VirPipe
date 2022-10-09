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

p <- add_argument(p, "--outdir", help="Output directory")

p <- add_argument(p, "--markdown", help="Rmarkdown template path")

p <- add_argument(p, "--qc_html", help="QC html file path")

p <- add_argument(p, "--map_summary", help="Reference map summary file path")

p <- add_argument(p, "--taxclassify_directory", help="Directory where taxonomic classifcation result is saved")

p <- add_argument(p, "--taxclassify_kraken_html", help="Kraken result html file path")

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
outdir <- args$outdir

# qc data
if(is.na(args$qc_html)){
  qc_link <- paste(outdir, "/qc/", prefix, ".qc.html", sep="")
} else{
  qc_link <- args$qc_html
}

# reference map data

if(is.na(args$map_summary)){
  map_summary_file <- paste(outdir, "/map/", prefix, ".map.filtered.tsv", sep="")
} else{
  map_summary_file <- args$map_summary
}

if(is.na(args$taxclassify_directory)){
  taxclassify_directory <- paste(outdir, "/taxclassify/", sep = "")
} else{
  taxclassify_directory <- args$taxclassify_directory
}

if(is.na(args$taxclassify_kraken_html)){
  taxclassify_kraken_link <- paste(outdir, "/taxclassify/", prefix, ".kraken.html", sep = "")
} else{
  taxclassify_kraken_link <- args$taxclassify_kraken_html
}

if(is.na(args$assembly_directory)){
  assembly_directory <- paste(outdir, "/assembly/", sep = "")
} else{
  assembly_directory <- args$assembly_directory
}

# blast data
if(is.na(args$blastn_table)){
  blast_blastn_file <- paste(outdir, "/postassembly/blast/", prefix, ".blastn.filtered.txt", sep="")
} else{
  blast_blastn_file <- args$blastn_table
}

if(is.na(args$blastn_html)){
  blastn_html_link <- paste(outdir, "/postassembly/blast/", prefix, ".blastn.filtered.html", sep="")
} else{
  blastn_html_link <- args$blastn_html
}

if(is.na(args$megablast_table)){
  blast_megablast_file <- paste(outdir, "/postassembly/blast/", prefix, ".megablast.filtered.txt", sep="")
} else{
  blast_megablast_file <- args$megablast_table
}

if(is.na(args$megablast_html)){
  megablast_html_link <- paste(outdir, "/postassembly/blast/", prefix, ".megablast.filtered.html", sep="")
} else{
  megablast_html_link <- args$megablast_html
}

# zoonotic rank data
if(is.na(args$zoonotic_rank_predictions)){
  zoonotic_rank_predictions_file <- paste(outdir, "/postassembly/zoonotic_rank/", prefix, ".predictions.csv", sep="")
} else{
  zoonotic_rank_predictions_file <- args$zoonotic_rank_predictions
}

rmd_template_file <- args$markdown
render(rmd_template_file,
      output_file=paste(prefix,".virpipe_summary", sep = ""),
      output_dir=outdir,
      knit_root_dir=getwd(),
      quiet=TRUE,
      params=list(prefix=prefix,
                  outdir=outdir,
                  qc_link=qc_link, 
                  map_summary_file=map_summary_file,
                  taxclassify_directory=taxclassify_directory,
                  taxclassify_kraken_link=taxclassify_kraken_link, 
                  assembly_directory=assembly_directory,
                  blast_blastn_file=blast_blastn_file,
		              blastn_html_link=blastn_html_link,
                  blast_megablast_file=blast_megablast_file,
		              megablast_html_link=megablast_html_link,
                  zoonotic_rank_predictions_file=zoonotic_rank_predictions_file))
