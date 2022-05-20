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
            filter_host_illumina(fastq_pair, host_genome)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            filter_host_nanopore(fastx, host_genome)
        } else { //hybrid
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            Channel.fromPath(params.y).set{fastx}
            filter_host_hybrid(fastq_pair, fastx, host_genome)
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
}

workflow filter_host_hybrid {
    take:
        fastq_pair
        fastx
        host_genome

    emit:
        host_filtered

    main:
        res1 = filter_host_illumina(fastq_pair, host_genome)
        res2 = filter_host_nanopore(fastx, host_genome)
        host_filtered = res1.combine(res2)
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
    tag "${params.prefix}:host_map_pair"

    input:
        path single 
        path host_genome
    output:
        path "aln.sorted.bam"
        val "$ext"
    
    script:
    ext = single.Extension
    """
    minimap2 -ax map-ont $host_genome $single | samtools view -Sb - | samtools sort -o aln.sorted.bam
    """
}

process extract_not_mapped_reads_pair {
    tag "${params.prefix}:extract_not_mapped_reads_pair"

    if(params.saveHostFiltered){
        publishDir path: "$params.outdir/filter_host", mode: 'copy'
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
        publishDir path: "$params.outdir/filter_host", mode: 'copy'
    }
            
    input:
        path bam
        val ext
    output:
        path("${params.prefix}.host_filtered.$ext")
    script:
    if (ext == "fastq" || ext =="fq")
    """
    samtools fastq -f 4 aln.sorted.bam > "${params.prefix}.host_filtered.$ext"
    """
    else
    """
    samtools fasta -f 4 aln.sorted.bam > "${params.prefix}.host_filtered.$ext"
    """

}