nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.x).set{contigs}
    
    zoonotic_rank(contigs)
}

workflow zoonotic_rank {
    take:
        contigs

    main:
    prodigal(contigs)
    prodigal_sco_to_zoonotic_rank_metadata(prodigal.out)
    zoonotic_rank_run(contigs, prodigal_sco_to_zoonotic_rank_metadata.out)
}

process prodigal {
    tag "${params.prefix}:prodigal"

    input:
        path contigs
    output:
        path "${params.prefix}.prodigal.sco"
    """
    prodigal -i $contigs -f sco -o ${params.prefix}.prodigal.sco -p meta
    """
}

process prodigal_sco_to_zoonotic_rank_metadata {
    tag "${params.prefix}:prodigal_sco_to_zoonotic_rank_metadata"

    input:
        path sco
    output:
        path "${params.prefix}.metadata.txt"

    """
    python ~/prodigal_sco_to_zoonotic_rank_metadata.py --input $sco --output ${params.prefix}.metadata.txt
    """
}

process zoonotic_rank_run {
    tag "${params.prefix}:zoonotic_rank_run"

    publishDir "${params.outdir}/post_assembly/zoonotic_rank", mode: 'copy'

    input:
        path contigs
        path metadata
    output:
        path "${params.prefix}*"
    """
    Rscript ~/zoonotic_rank/Scripts/PredictNovel.R fasta $contigs $metadata ${params.prefix}
    """
}