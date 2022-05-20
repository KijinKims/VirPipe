nextflow.enable.dsl=2

workflow {
    main:
        if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
            tax_classify_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
            tax_classify_nanopore(fastx)
        }
}

workflow tax_classify_illumina {
    take:
        fastq_pair

    main:
        tool = params.tool.tokenize()

        krona_inputs = Channel.empty()

        if ( tool.contains('kraken2') ) {
            Channel.fromPath(params.kraken2_db, type: 'dir').set{kraken2_db}
            kraken2_pair(fastq_pair, kraken2_db)
            kraken_krona_input = kraken2krona(kraken2_pair.out)
            krona_inputs = krona_inputs.concat(kraken_krona_input)

            kreport2list(kraken2_pair.out)
            metacomp(kreport2list.out)
        }

        krona(krona_inputs)
}

workflow tax_classify_nanopore {
    take:
        fastx

    main:
        tool = params.tool.tokenize()
        
        krona_inputs = Channel.empty()

        if ( tool.contains('kraken2') ) {
            Channel.fromPath(params.kraken2_db, type: 'dir').set{kraken2_db}
            kraken2_single(fastx, kraken2_db)
            kraken_krona_input = kraken2krona(kraken2_single.out)
            krona_inputs = krona_inputs.concat(kraken_krona_input)

            kreport2list(kraken2_single.out)
            metacomp(kreport2list.out)
        }

        krona(krona_inputs)
}

process kraken2_pair {
    tag "${params.prefix}:kraken2_pair"

    input:
        tuple path(pe1), path(pe2)
        path kraken2_db
    output:
        path "${params.prefix}.kraken_report.csv"
    script:

    """
    kraken2 --db $kraken2_db \
        --report ${params.prefix}.kraken_report.csv \
        --paired --threads ${params.threads} --confidence ${params.kraken2_confidence_threshold}\
        $pe1 $pe2
    """
}

process kraken2_single {
    tag "${params.prefix}:kraken2_single"

    input:
        path single
        path kraken2_db
    output:
        path "${params.prefix}.kraken_report.csv"

    """
    kraken2 --db $kraken2_db \
        --report ${params.prefix}.kraken_report.csv \
        --threads ${params.threads} --confidence ${params.kraken2_confidence_threshold}\
        $single
    """
}

process kraken2krona {
    tag "${params.prefix}:kraken2krona"

    input:
        path kraken_report
    output:
        tuple path("${params.prefix}.kraken"), val("kraken")
    """
    kreport2krona.py -r $kraken_report -o ${params.prefix}.kraken
    """
}

process kreport2list {
    tag "${params.prefix}:kreport2list"

    input:
        path kraken_report
    output:
        path "${params.prefix}.list"

    """
    perl ~/convert_krakenRep2list.pl < $kraken_report > ${params.prefix}.list
    """
}

process metacomp {
    tag "${params.prefix}:metacomp"
    
    publishDir "${params.outdir}/tax_classify", mode: 'copy'

    input:
        path kraken_list
    output:
        path "${params.prefix}_family.svg"
        path "${params.prefix}_genus.svg"
        path "${params.prefix}_species.svg"
    """
    Rscript ~/generate_heatmap.R ${kraken_list} ${params.prefix}
    """
}

process krona {
    tag "${params.prefix}:krona"
    publishDir "${params.outdir}/tax_classify", mode: 'copy'

    input:
        tuple path(krona_input), val(tool)
    output:
        path "${params.prefix}.${tool}.html"
    """
    ktImportText $krona_input -o ${params.prefix}.${tool}.html
    """
}