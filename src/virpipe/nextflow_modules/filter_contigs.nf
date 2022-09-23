nextflow.enable.dsl=2

workflow {
    main:
        Channel.fromPath(params.x).set{contigs}
        filter_contigs(contigs)
}

workflow filter_contigs {
    take:
        contigs
    emit:
        filtered_contigs
    main:
        length_filter_contigs(contigs)
        filtered_contigs = rename_contigs(length_filter_contigs.out)
}

process length_filter_contigs {
    tag "${params.prefix}:length_filter_contigs"
    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.temp.fasta")
    """
    reformat.sh in=$contigs out="${params.prefix}.filtered_contigs.temp.fasta" minlength=${params.min_contig_length}
    """
}

process rename_contigs {
    tag "${params.prefix}:rename_contigs"
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.fasta")
    """
    awk '/>/{print ">tig" ++i; next}{print}' < $contigs > "${params.prefix}.filtered_contigs.fasta"
    """
}