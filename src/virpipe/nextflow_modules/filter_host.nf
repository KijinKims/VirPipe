nextflow.enable.dsl=2

workflow {
    main:

        if (params.host_genome != null) {
            Channel.fromPath(params.host_genome).set{host_genome}
        } else {
            Channel.empty().set{host_genome}
        }

        if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            filter_illumina(fastq_pair, host_genome)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            filter_nanopore(fastx, host_genome)
        }
}

workflow filter_host_illumina {
    take:
        fastq_pair
        host_genome
    emit:
        host_filtered
    main:
        host_map_pair(fastq_pair, host_genome)
        host_filtered = extract_not_mapped_reads_pair(host_map_pair.out)
        summary_illumina(fastq_pair, host_filtered)
}

workflow filter_host_nanopore {
    take:
        fastx
        host_genome
    emit:
        host_filtered
    main:
        host_map_single(fastx, host_genome)
        host_filtered = extract_not_mapped_reads_single(host_map_single.out)
        summary_nanopore(fastx, host_filtered)
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

    if(params.saveHostFiltered){
        publishDir path: "$params.outdir/filter", mode: 'copy'
    }
        
    input:
        path bam
    output:
        tuple path("${params.prefix}.host_filtered_1.fastq"), path("${params.prefix}.host_filtered_2.fastq")
    """
    bedtools bamtofastq -i $bam -fq ${params.prefix}.host_filtered_1.fastq -fq2 ${params.prefix}.host_filtered_2.fastq
    """
}

process extract_not_mapped_reads_single {
    tag "${params.prefix}:extract_not_mapped_reads_single"

    if(params.saveHostFiltered){
        publishDir path: "$params.outdir/filter", mode: 'copy'
    }
            
    input:
        path sam
        val ext
    output:
        path("${params.prefix}.host_filtered.$ext")
    script:
    if (ext == "fastq" || ext =="fq")
    """
    samtools view -Sb $sam | samtools sort | samtools fastq -f 4 - > "${params.prefix}.host_filtered.$ext"
    """
    else
    """
    samtools view -Sb $sam | samtools sort | samtools fasta -f 4 - > "${params.prefix}.host_filtered.$ext"
    """

}

process summary_illumina {
    tag "${params.prefix}:filter_host_summary_illumina"
    publishDir path: "$params.outdir/filter", mode: 'copy'
    input:
        tuple path(prefilter_reads_1), path(prefilter_reads_2) 
        tuple path(postfilter_reads_1), path(postfilter_reads_2) 
    output:
        path "filter_host_summary.txt"
    """
    echo "#host_genome_path=$params.host_genome" >> filter_host_summary.txt
    echo "file\tpre\tpost"  >> filter_host_summary.txt

    pre_1=\$(bc <<< \$(cat $prefilter_reads_1 | wc -l)/4)
    post_1=\$(bc <<< \$(cat $postfilter_reads_1 | wc -l)/4)
    line_1="pe1\t"
    line_1+="\$pre_1\t"
    line_1+="\$post_1"
    echo \$line_1 >> filter_host_summary.txt

    pre_2=\$(bc <<< \$(cat $prefilter_reads_2 | wc -l)/4)
    post_2=\$(bc <<< \$(cat $postfilter_reads_2 | wc -l)/4)
    line_2="pe1\t"
    line_2+="\$pre_2\t"
    line_2+="\$post_2"
    echo \$line_2 >> filter_host_summary.txt
    """
}

process summary_nanopore {
    tag "${params.prefix}:filter_host_summary_nanopore"
    publishDir path: "$params.outdir/filter", mode: 'copy'
    input:
        path prefilter_reads
        path postfilter_reads
    output:
        path "filter_host_summary.txt"
    """
    echo "#host_genome_path=$params.host_genome" >> filter_host_summary.txt
    echo "file\tpre\tpost"  >> filter_host_summary.txt

    pre=\$(bc <<< \$(cat $prefilter_reads | wc -l)/4)
    post=\$(bc <<< \$(cat $postfilter_reads | wc -l)/4)
    line="single\t"
    line+="\$pre\t"
    line+="\$post"
    echo \$line >> filter_host_summary.txt
    """
}