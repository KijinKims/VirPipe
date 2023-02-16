nextflow.enable.dsl=2

workflow {
    main:
        if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).buffer(size:2).set{fastq_pair}
            map_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
            map_nanopore(fastq)
        }
}

workflow map_illumina {
    take:
        fastq_pair

    emit:
        filtered_map_out
    
    main:
        refs = map_refs_parse()
        refs.combine(fastq_pair).set{ref_fastq_pair}
        map_pair(ref_fastq_pair)
        qualimap(map_pair.out.flatten())
        bamcov(map_pair.out.flatten())
        Channel.from("SEQ_NAME\tSTART\tEND\tN_READS\tN_COVERED_BASES\tPERCENT_COVERED\tAVG_DEP\tAVG_BASEQ\tAVG_MAPQ").set{bamcov_header}
        bamcov.out.map{ it.text }.set{ bamcov_outputs }
        bamcov_header.concat(bamcov_outputs).collectFile(name: "${params.prefix}.map.tsv", newLine: true, storeDir: "${params.outdir}/map", sort: false).set{map_out}
        filtered_map_out = filter_map_process(map_out)
}

workflow map_nanopore {
    take:
        fastq

    emit:
        filtered_map_out
    
    main:
        refs = map_refs_parse()
        refs.combine(fastq).set{ref_fastq}
        map_single(ref_fastq)
        qualimap(map_single.out.flatten())
        bamcov(map_single.out.flatten())
        Channel.from("SEQ_NAME\tSTART\tEND\tN_READS\tN_COVERED_BASES\tPERCENT_COVERED\tAVG_DEP\tAVG_BASEQ\tAVG_MAPQ").set{bamcov_header}
        bamcov.out.map{ it.text }.set{ bamcov_outputs }
        bamcov_header.concat(bamcov_outputs).collectFile(name: "${params.prefix}.map.tsv", newLine: true, storeDir: "${params.outdir}/map", sort: false).set{map_out}
        filtered_map_out = filter_map_process(map_out)
}

workflow map_refs_parse {
    emit:
        refs
    
    main:
        if (params.ref) {
            Channel.fromPath(params.ref.tokenize()).set{refs}
        } else {
            Channel.empty().set{refs}
        }

        if (params.dir_ref) {
            Channel.fromPath(["${params.dir_ref}/*.fa", "${params.dir_ref}/*.fasta"])
                    .set { dir_ref_ch }
            refs.concat(dir_ref_ch).set{refs}
        }
    
    refs.ifEmpty{ println "WARN: No reference genome is given. Reference mapping will not happen."}
}

process map_pair {
    tag "${params.prefix}:map_pair"

    if(params.save_bam){
        publishDir path: "$params.outdir/map", mode: 'copy'
    }

    input:
        tuple path(ref), path(pe1), path(pe2)
    output:
        path '*bam'
    """
    minimap2 -ax sr $ref $pe1 $pe2 | samtools view -Sb - | samtools sort -o "${params.prefix}_${ref.simpleName}.bam"
    """
}

process map_single {
    tag "${params.prefix}:map_single"

    if(params.save_bam){
        publishDir path: "$params.outdir/map", mode: 'copy'
    }

    input:
        tuple path(ref), path(single)
    output:
        path '*bam'
    """
    minimap2 -ax map-ont $ref $single | samtools view -Sb - | samtools sort -o "${params.prefix}_${ref.simpleName}.bam"
    """
}

process bamcov {
    tag "${params.prefix}:bamcov"

    input:
        path bam
    output:
        path "${bam.simpleName}.bamcov.txt"
    """
    /home/user/bamcov -H $bam > ${bam.simpleName}.bamcov.txt
    """
}

process qualimap {
    tag "${params.prefix}:qualimap"

    publishDir "${params.outdir}/map", mode: 'copy', saveAs: { filename -> "${bam.simpleName}.pdf"}

    input:
        path bam
    output:
        path "*/qualimapReport.pdf"
    """
    qualimap bamqc -bam $bam -outfile qualimapReport.pdf
    """
}

process filter_map_process {
    tag "${params.prefix}:filter_map"

    publishDir "${params.outdir}/map", mode: 'copy' 

    input:
        path map_out
    output:
        path "${map_out.baseName}.filtered.tsv" optional true
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    
    tb <- read.table("${map_out}", sep="\t", header=TRUE, comment.char = "")
    if (nrow(tb)==0){
        quit("no")
    }
    #"SEQ_NAME\tSTART\tEND\tN_READS\tN_COVERED_BASES\tPERCENT_COVERED\tAVG_DEP\tAVG_BASEQ\tAVG_MAPQ"
    names(tb)[names(tb) == "END"] <- "LEN"
    tb <- subset(tb, select=-START)
    filtered_tb <- filter(tb, AVG_DEP > ${params.min_avg_dep})
    if (nrow(tb)>0) {
        sorted_filtered_tb <- filtered_tb[order(-filtered_tb\$"PERCENT_COVERED"),]
        write.table(sorted_filtered_tb, file ="${map_out.baseName}.filtered.tsv", row.names=FALSE, sep='\t', quote=FALSE)
    }
    """
}