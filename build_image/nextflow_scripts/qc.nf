nextflow.enable.dsl=2

workflow {
    main:
        if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).set{fastq_pair}
            qc_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
            qc_nanopore(fastq)
        }
}

workflow qc_illumina {
    take:
        fastq_pair
    
    main:
        fastqc(fastq_pair)
        multiqc(fastqc.out.collect())
}

workflow qc_nanopore {
    take:
        fastq
    
    main:
        nanoplot(fastq)
}

process fastqc {
    tag "${params.prefix}:fastqc"

    input:
        path f
    output:
        path "*.zip"
    """
    fastqc $f
    """
}

process multiqc {
    tag "${params.prefix}:multiqc"
    publishDir "${params.outdir}/qc", mode: 'copy', saveAs: { filename -> "${params.prefix}.qc.html"}

    input:
        path f
    output:
        path "multiqc_report.html"
    """
    multiqc ${f.join(' ')}
    """
}

process nanoplot {
    tag "${params.prefix}:nanoplot"
    publishDir "${params.outdir}/qc", mode: 'copy', saveAs: { filename -> "${params.prefix}.qc.html"}

    input:
        path f
    output:
        path "NanoPlot-report.html"
    script:
        """
        #!/usr/bin/env bash
        NanoPlot --fastq $f --minlength ${params.nanopore_min_read_length} --minqual ${params.nanopore_min_read_quality}
        """
}