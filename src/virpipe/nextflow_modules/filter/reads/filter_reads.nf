nextflow.enable.dsl=2

workflow {
    main:

        if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            filter_reads_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            filter_reads_nanopore(fastx)
        }
}


workflow filter_reads_illumina {
    take:
        fastq_pair
    emit:
        filtered

    main:
        decompress_pair(fastq_pair)
        filter_illumina(decompress_pair.out)
        filtered = trim_illumina(filter_illumina.out)
}

workflow filter_reads_nanopore {
    take:
        fastx
    emit:
        filtered
    main:
        decompress_single(fastx)
        filtered = filter_nanopore(decompress_single.out)
}

process decompress_pair {
    tag "${params.prefix}:decompress_pair"
    stageInMode "copy"

    input:
        tuple path(pe1), path(pe2)
    output:
        tuple path("${pe1.simpleName}.fastq"), path("${pe2.simpleName}.fastq")
    script:

    if (pe1.Extension == "gz" || pe1.Extension == "gunzip")
        """
            gunzip $pe1
            gunzip $pe2
        """
    else
        """
        """

}
process decompress_single {
    tag "${params.prefix}:decompress_single"
    stageInMode "copy"

    input:
        path fastx
    output:
        path "*", includeInputs: true
    script:

    if (fastx.Extension == "gz" || fastx.Extension == "gunzip")
        """
            gunzip -c $fastx > ${fastx.simpleName}
        """
    else
        """
        """
}

process filter_illumina {
    tag "${params.prefix}:filter_illumina"

    input:
        tuple path(pe1), path(pe2)
    output:
        tuple path("${params.prefix}.filtered_1.fastq"), path("${params.prefix}.filtered_2.fastq")
    """
    prinseq-lite.pl -fastq $pe1 -fastq2 $pe2 -min_qual_score ${params.illumina_min_read_quality} -derep 12 -out_good ${params.prefix}.filtered -out_bad null
    """
}

process trim_illumina {
    tag "${params.prefix}:trim_illumina"

    if(params.saveFiltered){
        publishDir path: "$params.outdir/filter", mode: 'copy'
    }
    input:
        tuple path(pe1), path(pe2)
    output:
        tuple path("${params.prefix}.filtered_trimmed_1.fastq"), path("${params.prefix}.filtered_trimmed_2.fastq")
    """
    trimmomatic PE -phred33 $pe1 $pe2 ${params.prefix}.filtered_trimmed_1.fastq unpaired.fastq ${params.prefix}.filtered_trimmed_2.fastq unpaired.fq ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process filter_nanopore {
    tag "${params.prefix}:filter_nanopore"

    if(params.saveFiltered){
        publishDir path: "$params.outdir/filter", mode: 'copy'
    }

    input:
        path single
    output:
        path("${params.prefix}.filtered.${single.extension}")
    """
    NanoFilt --length ${params.nanopore_min_read_length} --readtype 1D --quality ${params.nanopore_min_read_quality} $single > ${params.prefix}.filtered.${single.extension}
    """
}