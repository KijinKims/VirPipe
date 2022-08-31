nextflow.enable.dsl=2

workflow {
    main:

        if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            consensus_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            consensus_nanopore(fastx)
        } else { //hybrid
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            Channel.fromPath(params.y).set{fastx}
            consensus_hybrid(fastq_pair, fastx)
        }
}

workflow consensus_illumina {
    take:
        fastq_pair
    
    main:
        refs = consensus_refs_parse()
        ref.combine(fastq_pair).flatten().buffer(size:3).set{ref_fastq}
        consensus_illumina_process(ref_fastq)
}

workflow consensus_nanopore {
    take:
        fastx
    
    main:
        refs = consensus_refs_parse()
        refs.combine(fastx).set{ref_fastx}
        consensus_nanopore_process(ref_fastx)
}

workflow consensus_refs_parse {
    emit:
        refs
    
    main:
        if (params.ref  != null) {
            Channel.fromPath(params.ref.tokenize()).set{refs}
        } else {
            Channel.empty().set{refs}
        }
        refs.ifEmpty{ println "WARN: No reference genome is given. No consensus will not be generated."}
}

process consensus_illumina_process {
    tag "${params.prefix}:consensus_illumina_process"

    input:
        tuple path(pe1), path(pe2), path(ref)
    output:
        path "${params.prefix}_${ref.simpleName}.consensus.fasta"
    """
    """
}

process consensus_nanopore_process {
    publishDir "${params.outdir}/consensus", mode: 'copy'
    tag "${params.prefix}:consensus_nanopore_process"

    input:
        tuple path(ref), path(single)
    output:
        path "${params.prefix}_${ref.simpleName}.consensus.fasta"
    """
    # variant call
    mini_align -t ${params.threads} -p ${params.prefix}_${ref.simpleName} -i $single -r $ref -m -f -M 2 -S 4 -O 4,24 -E 2,1
    medaka consensus --model r941_prom_variant_g360 --threads ${params.threads} --batch_size 100 ${params.prefix}_${ref.simpleName}.bam ${params.prefix}_${ref.simpleName}.hdf
    medaka variant $ref ${params.prefix}_${ref.simpleName}.hdf ${params.prefix}_${ref.simpleName}.medaka.vcf
    medaka tools annotate --dpsp ${params.prefix}_${ref.simpleName}.medaka.vcf $ref ${params.prefix}_${ref.simpleName}.bam ${params.prefix}_${ref.simpleName}.medaka.annotated.vcf
    bcftools view -i 'QUAL>${params.variant_quality_threshold} & INFO/DP>${params.variant_depth_threshold}' ${params.prefix}_${ref.simpleName}.medaka.annotated.vcf > ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf
    bgzip -f ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf
    tabix -p vcf ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf.gz

    # apply variants to reference to generate consensus
    cat $ref | bcftools consensus ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf.gz > ${ref.simpleName}.consensus.fasta
    samtools faidx ${ref.simpleName}.consensus.fasta
    awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' ${ref.simpleName}.consensus.fasta.fai > ${ref.simpleName}.consensus.fasta.bed

    # drop low coverage region
    bedtools genomecov -bga -ibam ${params.prefix}_${ref.simpleName}.bam | awk '\$4 < ${params.low_cov_threshold}' | bedtools merge -i - > low_cov_${params.low_cov_threshold}.bed
    bedtools intersect -a low_cov_${params.low_cov_threshold}.bed -b ${ref.simpleName}.consensus.fasta.bed > consensus_low_cov_${params.low_cov_threshold}.bed
    bedtools maskfasta -fi ${ref.simpleName}.consensus.fasta -bed consensus_low_cov_${params.low_cov_threshold}.bed -fo ${params.prefix}_${ref.simpleName}.consensus.fasta

    # change fasta header
    header=">${params.prefix}_${ref.simpleName}_consensus low_cov_thrs=${params.low_cov_threshold} var_qual_thrs=${params.variant_quality_threshold} var_depth_thrs=${params.variant_depth_threshold}"
    sed -i "1s/.*/\$header/" ${params.prefix}_${ref.simpleName}.consensus.fasta
    """
}
