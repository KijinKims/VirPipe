nextflow.enable.dsl=2

workflow {
    main:
        if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).set{fastq_pair}
            qc_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            qc_nanopore(fastx)
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
        fastx
    
    main:

        nanoplot(fastx)
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

        fullfile=\$(basename -- "$f")
        ext="\${fullfile##*.}"
        filename="\${fullfile%.*}"
        if [[ $f == *gz || $f == *gunzip ]] ; then
            real_ext="\${filename##*.}"
        else
            real_ext=\$ext
        fi

        if [[ \${real_ext} == "fastq" || \${real_ext} == "fq" ]]; then
            NanoPlot --fastq $f --minlength ${params.nanopore_min_read_length} --minqual ${params.nanopore_min_read_quality}
        else
            NanoPlot --fasta $f --minlength ${params.nanopore_min_read_length} --minqual ${params.nanopore_min_read_quality}
        fi
        """
}