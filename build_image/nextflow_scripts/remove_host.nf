nextflow.enable.dsl=2

workflow {
    main:
        if (params.host_genome != null) {
            Channel.fromPath(params.host_genome).set{host_genome}
        } else {
            Channel.empty().set{host_genome}
        }

        if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).buffer(size:2).set{fastq_pair}
            remove_host_illumina(fastq_pair, host_genome)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
            remove_host_nanopore(fastq, host_genome)
        }
}

workflow remove_host_illumina {
    take:
        fastq_pair
        host_genome
    emit:
        host_removed
    main:
        host_map_pair(fastq_pair, host_genome)
        host_removed = extract_not_mapped_reads_pair(host_map_pair.out)
        summary_illumina(fastq_pair, host_removed)
}

workflow remove_host_nanopore {
    take:
        fastq
        host_genome
    emit:
        host_removed
    main:
        host_map_single(fastq, host_genome)
        host_removed = extract_not_mapped_reads_single(host_map_single.out)
        summary_nanopore(fastq, host_removed)
}

process host_map_pair {
    tag "${params.prefix}:host_map_pair"

    input:
        tuple path(pe1), path(pe2)
        path host_genome
    output:
        path "aln.sorted.bam"
    """
    minimap2 -ax sr $host_genome $pe1 $pe2 | samtools view -Sb - | samtools view -b -f 12 -F 256 | samtools sort -o aln.sorted.bam
    """
}

process host_map_single {
    tag "${params.prefix}:host_map_single"

    input:
        path single 
        path host_genome
    output:
        path "aln.sam"
        val "$ext"
    
    script:
    ext = single.Extension
    """
    minimap2 -ax map-ont $host_genome $single > aln.sam
    """
}

process extract_not_mapped_reads_pair {
    tag "${params.prefix}:extract_not_mapped_reads_pair"

    if(params.save_host_removed){
        publishDir path: "$params.outdir/host_removed", mode: 'copy'
    }
        
    input:
        path bam
    output:
        tuple path("${params.prefix}.host_removed_1.fastq"), path("${params.prefix}.host_removed_2.fastq")
    """
    bedtools bamtofastq -i $bam -fq ${params.prefix}.host_removed_1.fastq -fq2 ${params.prefix}.host_removed_2.fastq
    """
}

process extract_not_mapped_reads_single {
    tag "${params.prefix}:extract_not_mapped_reads_single"

    if(params.save_host_removed){
        publishDir path: "$params.outdir/host_removed", mode: 'copy'
    }
            
    input:
        path sam
        val ext
    output:
        path("${params.prefix}.host_removed.$ext")
    script:
    if (ext == "fastq" || ext =="fq")
    """
    samtools view -Sb $sam | samtools sort | samtools fastq -f 4 - > "${params.prefix}.host_removed.$ext"
    """
    else
    """
    samtools view -Sb $sam | samtools sort | samtools fasta -f 4 - > "${params.prefix}.host_removed.$ext"
    """

}

process summary_illumina {
    tag "${params.prefix}:remove_host_summary_illumina"
    publishDir path: "$params.outdir/host_removed", mode: 'copy'
    input:
        tuple path(prepreprocess_1), path(prepreprocess_2) 
        tuple path(postpreprocess_1), path(postpreprocess_2) 
    output:
        path "remove_host_summary.txt"
    """
    echo "#host_genome_path=$params.host_genome.Name" >> remove_host_summary.txt
    echo "file\tpre\tpost"  >> remove_host_summary.txt

    pre_1=\$(bc <<< \$(cat $prepreprocess_1 | wc -l)/4)
    post_1=\$(bc <<< \$(cat $postpreprocess_1 | wc -l)/4)
    line_1="pe1\t"
    line_1+="\$pre_1\t"
    line_1+="\$post_1"
    echo \$line_1 >> remove_host_summary.txt

    pre_2=\$(bc <<< \$(cat $prepreprocess_2 | wc -l)/4)
    post_2=\$(bc <<< \$(cat $postpreprocess_2 | wc -l)/4)
    line_2="pe1\t"
    line_2+="\$pre_2\t"
    line_2+="\$post_2"
    echo \$line_2 >> remove_host_summary.txt
    """
}

process summary_nanopore {
    tag "${params.prefix}:remove_host_summary_nanopore"
    publishDir path: "$params.outdir/host_removed", mode: 'copy'
    input:
        path prepreprocess
        path postpreprocess
    output:
        path "remove_host_summary.txt"
    """
    echo "#host_genome_path=$params.host_genome.Name" >> remove_host_summary.txt
    echo "file\tpre\tpost"  >> remove_host_summary.txt

    pre=\$(bc <<< \$(cat $prepreprocess | wc -l)/4)
    post=\$(bc <<< \$(cat $postpreprocess | wc -l)/4)
    line="single\t"
    line+="\$pre\t"
    line+="\$post"
    echo \$line >> remove_host_summary.txt
    """
}