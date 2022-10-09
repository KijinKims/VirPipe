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
        summary_illumina(decompress_pair.out, filtered)
}

workflow filter_reads_nanopore {
    take:
        fastx
    emit:
        filtered
    main:
        decompress_single(fastx)
        filtered = filter_nanopore(decompress_single.out)
        summary_nanopore(decompress_single.out, filtered)
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

process summary_illumina {
    tag "${params.prefix}:filter_reads_summary_illumina"
    publishDir path: "$params.outdir/filter", mode: 'copy'
    input:
        tuple path(prefilter_reads_1), path(prefilter_reads_2) 
        tuple path(postfilter_reads_1), path(postfilter_reads_2) 
    output:
        path "filter_reads_summary.txt"
    """
    echo "#min_qual=$params.illumina_min_read_quality" >> filter_reads_summary.txt
    echo "file\tpre\tpost"  >> filter_reads_summary.txt

    pre_1=\$(bc <<< \$(cat $prefilter_reads_1 | wc -l)/4)
    post_1=\$(bc <<< \$(cat $postfilter_reads_1 | wc -l)/4)
    line_1="pe1\t"
    line_1+="\$pre_1\t"
    line_1+="\$post_1"
    echo \$line_1 >> filter_reads_summary.txt

    pre_2=\$(bc <<< \$(cat $prefilter_reads_2 | wc -l)/4)
    post_2=\$(bc <<< \$(cat $postfilter_reads_2 | wc -l)/4)
    line_2="pe1\t"
    line_2+="\$pre_2\t"
    line_2+="\$post_2"
    echo \$line_2 >> filter_reads_summary.txt
    """
}

process summary_nanopore {
    tag "${params.prefix}:filter_reads_summary_nanopore"
    publishDir path: "$params.outdir/filter", mode: 'copy'
    input:
        path prefilter_reads
        path postfilter_reads
    output:
        path "filter_reads_summary.txt"
    """
    echo "#min_qual=$params.nanopore_min_read_quality" >> filter_reads_summary.txt
    echo "#min_len=$params.nanopore_min_read_length" >> filter_reads_summary.txt
    echo "file\tpre\tpost"  >> filter_reads_summary.txt

    pre=\$(bc <<< \$(cat $prefilter_reads | wc -l)/4)
    post=\$(bc <<< \$(cat $postfilter_reads | wc -l)/4)
    line="single\t"
    line+="\$pre\t"
    line+="\$post"
    echo \$line >> filter_reads_summary.txt
    """
}