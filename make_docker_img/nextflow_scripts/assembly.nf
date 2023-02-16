nextflow.enable.dsl=2

workflow {
    emit:
        contigs

    main:

        if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).buffer(size:2).set{fastq_pair}
            contigs = assembly_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
            contigs = assembly_nanopore(fastq)
        }
}

workflow assembly_illumina {
    take:
        fastq_pair
    emit:
        filtered_contigs
    main:
        
        contigs = spades(fastq_pair)
        length_filter_contigs(contigs)
        filtered_contigs = rename_contigs(length_filter_contigs.out)
}

workflow assembly_nanopore {
    take:
        fastq
    emit:
        filtered_contigs
    main:
        
        if (params.tool == "canu") {
            contigs = canu(fastq)
        } else if (params.tool == "megahit") {
            contigs = megahit(fastq)
        } else {
            contigs = flye(fastq)
        }

        length_filter_contigs(contigs)
        filtered_contigs = rename_contigs(length_filter_contigs.out)
}

process spades {
    tag "${params.prefix}:spades"

    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}
    
    input:
        tuple path(pe1), path(pe2)
    output:
        path "${params.prefix}/contigs.fasta"
    script:

    """
    spades.py --pe1-1 $pe1 --pe1-2 $pe2 -o ${params.prefix} -t ${params.threads} $params.spades_options
    """
}

process megahit {
    tag "${params.prefix}:megahit"
    
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        path fastq
    output:
        path "${params.prefix}/${params.prefix}.contigs.fa"
    script:

    """
    megahit -r $fastq -o $params.prefix --out-prefix $params.prefix -t ${params.threads} $params.megahit_options
    """
}

process canu {
    tag "${params.prefix}:canu"
    
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        path fastq
    output:
        path "${params.prefix}/${params.prefix}.contigs.fasta"
    script:

    """
    canu -p $params.prefix -d $params.prefix -nanopore $fastq $params.canu_options
    """
}

process flye {
    tag "${params.prefix}:flye"
    
    publishDir "${params.outdir}/assembly", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        path fastq
    output:
        path "${params.prefix}/assembly.fasta" optional true
    script:

    """
    flye --nano-raw $fastq --out-dir $params.prefix --threads ${params.threads} $params.flye_options
    """
}

process length_filter_contigs {
    tag "${params.prefix}:length_filter_contigs"
    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.temp.fasta")
    """
    reformat.sh in=$contigs out="${params.prefix}.filtered_contigs.temp.fasta" minlength=${params.min_contig_length}
    """
}

process rename_contigs {
    tag "${params.prefix}:rename_contigs"
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.fasta")
    """
    awk '/>/{print ">tig" ++i; next}{print}' < $contigs > "${params.prefix}.filtered_contigs.fasta"
    """
}