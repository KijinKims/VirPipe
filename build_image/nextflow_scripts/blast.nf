nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.fasta).set{contigs}

    blast(contigs)
}

workflow blast {
    take:
        contigs

    emit:
        filtered

    main:

    blast_nucl_outs = Channel.empty()

    Channel.fromPath(params.blast_db_dir, type: 'dir').set{blast_db_dir}
    blastn(contigs, blast_db_dir)
    blast_nucl_outs = blast_nucl_outs.concat(blastn.out)

    Channel.fromPath(params.blast_db_dir, type: 'dir').set{blast_db_dir}
    megablast(contigs, blast_db_dir)
    blast_nucl_outs = blast_nucl_outs.concat(megablast.out)

    filtered = filter_blast_process(blast_nucl_outs, params.min_blast_aln_len)

    Channel.fromPath(params.taxonomizr_db).set{taxonomizr_db}
    match_taxonomy_process(filtered.combine(taxonomizr_db))
}

process blastn {
    tag "${params.prefix}:blastn"

    publishDir "${params.outdir}/postassembly/blast", mode: 'copy'

    input:
        path contigs
        path blast_db_dir
    output:
        path "${params.prefix}.blastn.raw.txt"
    """
    { echo 'QUERY_ID\tREF_ID\tTAX_ID\tREF_TITLE\tPER_IDENT\tQUERY_LEN\tALN_LEN\tMISMATCH\tGAPOPEN\tQUERY_START\tQUERY_END\tREF_START\tREF_END\tEVALUE\tBITSCORE'; \
    blastn -query $contigs -db "${blast_db_dir}/${params.blast_db_name}" -task blastn -evalue ${params.min_evalue} -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 12; \
    } > ${params.prefix}.blastn.txt
    """
}

process megablast {
    tag "${params.prefix}:megablast"

    publishDir "${params.outdir}/postassembly/blast", mode: 'copy'

    input:
        path contigs
        path blast_db_dir
    output:
        path "${params.prefix}.megablast.raw.txt"
    """
    { echo 'QUERY_ID\tREF_ID\tTAX_ID\tREF_TITLE\tPER_IDENT\tQUERY_LEN\tALN_LEN\tMISMATCH\tGAPOPEN\tQUERY_START\tQUERY_END\tREF_START\tREF_END\tEVALUE\tBITSCORE'; \
    blastn -query $contigs -db "${blast_db_dir}/${params.blast_db_name}" -task megablast -evalue ${params.min_evalue} -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 12; \
    } > ${params.prefix}.megablast.txt
    """
}

process filter_blast_process {
    tag "${params.prefix}:filter_blast_process"

    publishDir "${params.outdir}/postassembly/blast", mode: 'copy' 

    input:
        path blast_out
        val min_blast_aln_len
    output:
        path "${blast_out.baseName}.filtered.txt" optional true
        
    """
    #!/usr/bin/env Rscript

    library(dplyr)

    tb <- read.table("${blast_out}", sep="\t", header=TRUE, comment.char = "")
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

process match_taxonomy_process {
    tag "${params.prefix}:match_taxonomy"

    publishDir "${params.outdir}/postassembly/blast", mode: 'copy'

    input:
        tuple path(blast_out), path(taxonomizr_db)
    output:
        path "*"
    """
    Rscript /home/user/custom_scripts/match_taxonomy.R $blast_out ${blast_out.baseName}.html ${taxonomizr_db}
    """
}