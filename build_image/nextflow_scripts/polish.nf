nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.fasta).set{contigs}
    Channel.fromPath(params.reads).set{reads}
    
    polish(contigs,reads)
}

workflow polish {
    take:
        contigs
        reads

    emit:
        contigs

    main:

    paf = ava_map(contigs, reads)
    contigs = racon(contigs, reads, paf)
    contigs = medaka(contigs, reads)
}

process ava_map {
    tag "${params.prefix}:ava_map"

    input:
        path contigs
        path reads
    output:
        path "${params.prefix}.paf"
    """
    minimap2 -x ava-ont $contigs $reads -t 12 > ${params.prefix}.paf
    """
}

process racon {
    tag "${params.prefix}:racon"

    input:
        path contigs
        path reads
        path paf
    output:
        path "${params.prefix}.polished_contigs.fasta"
    """
    racon -t 12 $reads $paf $contigs > ${params.prefix}.polished_contigs.fasta
    """
}

process medaka {
    tag "${params.prefix}:medaka"

    publishDir "$params.outdir/polish", mode: 'copy'
    
    input:
        path contigs
        path reads
    output:
        path "${params.prefix}.polished_contigs.fasta"
    """
    medaka_consensus -i $reads -d $contigs \
         -o ${params.prefix} \
         -t 12 -m ${params.medaka_model}
    """
}