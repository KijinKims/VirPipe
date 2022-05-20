nextflow.enable.dsl=2


workflow {
    main:
        Channel.fromPath(params.x).set{map_out}
        filter_map(map_out)
}

workflow filter_map {
    take:
        map_out
    main:
        filtered = filter_map_process(map_out)
}

process filter_map_process {
    tag "${params.prefix}:filter_map"

    publishDir "${params.outdir}/map", mode: 'copy' 

    input:
        path map_out
    output:
        path "${map_out.baseName}.filtered.tsv" optional true
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    
    tb <- read.table("${map_out}", sep="\t", header=TRUE)
    if (nrow(tb)==0){
        quit("no")
    }
    #"SEQ_NAME\tSTART\tEND\tN_READS\tN_COVERED_BASES\tPERCENT_COVERED\tAVG_COV\tAVG_BASEQ\tAVG_MAPQ"
    names(tb)[names(tb) == "END"] <- "LEN"
    tb <- subset(tb, select=-START)
    filtered_tb <- filter(tb, AVG_COV > ${params.min_map_out_avg_cov})
    if (nrow(tb)>0) {
        sorted_filtered_tb <- filtered_tb[order(-filtered_tb\$"PERCENT_COVERED"),]
        write.table(sorted_filtered_tb, file ="${map_out.baseName}.filtered.tsv", row.names=FALSE, sep='\t', quote=FALSE)
    }
    """
}