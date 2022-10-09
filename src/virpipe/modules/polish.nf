nextflow.enable.dsl=2

workflow {
    Channel.fromPath(params.x).set{contigs}
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

    tool = params.tool.tokenize()

    if (tool.contains("racon")) {
        paf = ava_map(contigs, reads)
        contigs = racon(contigs, reads, paf)
    }
    
    if (tool.contains("medaka")) {
        contigs = medaka(contigs, reads)
    }
}

process ava_map {
    tag "${params.prefix}:ava_map"

    input:
        path contigs
        path reads
    output:
        path "${params.prefix}.paf"
    """
    minimap2 -x ava-ont $contigs $reads -t ${params.threads} > ${params.prefix}.paf
    """
}

process racon {
    tag "${params.prefix}:racon"

    publishDir "${params.outdir}/polish", mode: 'copy', saveAs: { filename -> params.tool.tokenize().contains("medaka") ? "${params.prefix}.racon_polished_contigs.fasta" : filename }

    input:
        path contigs
        path reads
        path paf
    output:
        path "${params.prefix}.polished_contigs.fasta"
    """
    racon -t ${params.threads} $reads $paf $contigs > ${params.prefix}.polished_contigs.fasta
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
         -t ${params.threads} -m ${params.medaka_model}
    """
}