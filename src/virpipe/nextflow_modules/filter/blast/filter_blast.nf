nextflow.enable.dsl=2

workflow {
    emit:
        filtered

    main:
        Channel.fromPath(params.x).set{blast_out}

        filtered = filter_blast(blast_out)
}

workflow filter_blast {
    take:
        blast_out
    emit:
        filtered
    main:
        filtered = filter_blast_process(blast_out, params.min_blast_aln_len)
}

process filter_blast_process {
    tag "${params.prefix}:filter_blast_process"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy' 

    input:
        path blast_out
        val min_blast_aln_len
    output:
        path "${blast_out.baseName}.filtered.txt" optional true
        
    """
    #!/usr/bin/env Rscript

    library(dplyr)

    tb <- read.table("${blast_out}", sep="\t", header=TRUE)
    filtered_tb <- tb %>%
            filter(
                .data[["ALN_LEN"]] > $min_blast_aln_len
            )
    if (nrow(filtered_tb)>0) {
        sorted_filtered_tb <- tb[order(-tb\$"BITSCORE"),]
        write.table(sorted_filtered_tb, file ="${blast_out.baseName}.filtered.txt", row.names=FALSE, sep='\t', quote=FALSE)
    }
    """
}