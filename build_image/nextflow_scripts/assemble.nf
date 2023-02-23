nextflow.enable.dsl=2

workflow {
    emit:
        contigs

    main:

        if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).buffer(size:2).set{fastq_pair}
            contigs = assemble_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
            contigs = assemble_nanopore(fastq)
        }
}

workflow assemble_illumina {
    take:
        fastq_pair
    emit:
        filtered_contigs
    main:
        
        contigs = spades(fastq_pair)
        length_filter_contigs(contigs)
        filtered_contigs = rename_contigs(length_filter_contigs.out)
}

workflow assemble_nanopore {
    take:
        fastq
    emit:
        filtered_contigs
    main:
        
        if (params.assembly_tool == "canu") {
            contigs = canu(fastq)
        } else if (params.assembly_tool == "megahit") {
            contigs = megahit(fastq)
        } else {
            contigs = flye(fastq)
        }

        length_filter_contigs(contigs)
        filtered_contigs = rename_contigs(length_filter_contigs.out)
}

process spades {
    tag "${params.prefix}:spades"

    publishDir "${params.outdir}/assemble", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}
    
    input:
        tuple path(pe1), path(pe2)
    output:
        path "${params.prefix}/contigs.fasta"
    script:

    """
    spades.py --pe1-1 $pe1 --pe1-2 $pe2 -o ${params.prefix} -t 12 $params.spades_options
    """
}

process megahit {
    tag "${params.prefix}:megahit"
    
    publishDir "${params.outdir}/assemble", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        path fastq
    output:
        path "${params.prefix}/${params.prefix}.contigs.fa"
    script:

    """
    megahit -r $fastq -o $params.prefix --out-prefix $params.prefix -t 12 $params.megahit_options
    """
}

process canu {
    tag "${params.prefix}:canu"
    
    publishDir "${params.outdir}/assemble", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

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
    errorStrategy 'ignore'
    
    tag "${params.prefix}:flye"
    
    publishDir "${params.outdir}/assemble", mode: 'copy', saveAs: { filename -> "${params.prefix}.contigs.fasta"}

    input:
        path fastq
    output:
        path "${params.prefix}/assembly.fasta" optional true
    script:

    """
    flye --nano-raw $fastq --out-dir $params.prefix --threads 12 $params.flye_options
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
    publishDir "${params.outdir}/assemble", mode: 'copy'
    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.fasta")
    """
    awk '/>/{print ">tig" ++i; next}{print}' < $contigs > "${params.prefix}.filtered_contigs.fasta"
    """
}