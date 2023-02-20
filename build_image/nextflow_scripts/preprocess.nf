nextflow.enable.dsl=2

workflow {
    main:
        if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).buffer(size:2).set{fastq_pair}
            preprocess_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
            preprocess_nanopore(fastq)
        }
}


workflow preprocess_illumina {
    take:
        fastq_pair
    emit:
        preprocessed
    main:
        decompress_pair(fastq_pair)
        filter_illumina(decompress_pair.out)
        preprocessed = trim_illumina(filter_illumina.out)
        summary_illumina(decompress_pair.out, preprocessed)
}

workflow preprocess_nanopore {
    take:
        fastq
    emit:
        preprocessed
    main:
        decompress_single(fastq)
        preprocessed = filter_nanopore(decompress_single.out)
        summary_nanopore(decompress_single.out, preprocessed)
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
        """
    else
        """
        """

    if (pe2.Extension == "gz" || pe2.Extension == "gunzip")
        """
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
        path fastq
    output:
        path "*", includeInputs: true
    script:

    if (fastq.Extension == "gz" || fastq.Extension == "gunzip")
        """
            gunzip $fastq
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
        tuple path("${params.prefix}.preprocessed_1.fastq"), path("${params.prefix}.preprocessed_2.fastq")
    """
    prinseq-lite.pl -fastq $pe1 -fastq2 $pe2 -min_qual_score ${params.illumina_min_read_quality} -derep 12 -out_good ${params.prefix}.preprocessed -out_bad null
    """
}

process trim_illumina {
    tag "${params.prefix}:trim_illumina"

    if(params.save_preprocessed){
        publishDir path: "$params.outdir/preprocess", mode: 'copy'
    }
    input:
        tuple path(pe1), path(pe2)
    output:
        tuple path("${params.prefix}.preprocessed_trimmed_1.fastq"), path("${params.prefix}.preprocessed_trimmed_2.fastq")
    """
    trimmomatic PE -phred33 $pe1 $pe2 ${params.prefix}.preprocessed_trimmed_1.fastq unpaired.fastq ${params.prefix}.preprocessed_trimmed_2.fastq unpaired.fq ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process filter_nanopore {
    tag "${params.prefix}:filter_nanopore"

    if(params.save_preprocessed){
        publishDir path: "$params.outdir/preprocess", mode: 'copy'
    }

    input:
        path single
    output:
        path("${params.prefix}.preprocessed.${single.extension}")
    """
    NanoFilt --length ${params.nanopore_min_read_length} --readtype 1D --quality ${params.nanopore_min_read_quality} $single > ${params.prefix}.preprocessed.${single.extension}
    """
}

process summary_illumina {
    tag "${params.prefix}:preprocess_summary_illumina"
    publishDir path: "$params.outdir/preprocess", mode: 'copy'
    input:
        tuple path(prepreprocess_1), path(prepreprocess_2) 
        tuple path(postpreprocess_1), path(postpreprocess_2) 
    output:
        path "preprocess_summary.txt"
    """
    echo "#min_qual=$params.illumina_min_read_quality" >> preprocess_summary.txt
    echo "file\tpre\tpost"  >> preprocess_summary.txt

    pre_1=\$(bc <<< \$(cat $prepreprocess_1 | wc -l)/4)
    post_1=\$(bc <<< \$(cat $postpreprocess_1 | wc -l)/4)
    line_1="pe1\t"
    line_1+="\$pre_1\t"
    line_1+="\$post_1"
    echo \$line_1 >> preprocess_summary.txt

    pre_2=\$(bc <<< \$(cat $prepreprocess_2 | wc -l)/4)
    post_2=\$(bc <<< \$(cat $postpreprocess_2 | wc -l)/4)
    line_2="pe1\t"
    line_2+="\$pre_2\t"
    line_2+="\$post_2"
    echo \$line_2 >> preprocess_summary.txt
    """
}

process summary_nanopore {
    tag "${params.prefix}:preprocess_summary_nanopore"
    publishDir path: "$params.outdir/preprocess", mode: 'copy'
    input:
        path prepreprocess
        path postpreprocess
    output:
        path "preprocess_summary.txt"
    """
    echo "#min_qual=$params.nanopore_min_read_quality" >> preprocess_summary.txt
    echo "#min_len=$params.nanopore_min_read_length" >> preprocess_summary.txt
    echo "file\tpre\tpost"  >> preprocess_summary.txt

    pre=\$(bc <<< \$(cat $prepreprocess | wc -l)/4)
    post=\$(bc <<< \$(cat $postpreprocess | wc -l)/4)
    line="single\t"
    line+="\$pre\t"
    line+="\$post"
    echo \$line >> preprocess_summary.txt
    """
}