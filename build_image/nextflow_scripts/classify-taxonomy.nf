nextflow.enable.dsl=2

workflow {
    main:
        if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).buffer(size:2).set{fastq_pair}
            classify_taxonomy_illumina(fastq_pair)
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
            classify_taxonomy_nanopore(fastq)
        }
}

workflow classify_taxonomy_illumina {
    take:
        fastq_pair

    main:
        krona_inputs = Channel.empty()

        Channel.fromPath(params.kraken2_db, type: 'dir').set{kraken2_db}
        kraken2_pair(fastq_pair, kraken2_db)
        kraken_krona_input = kraken2krona(kraken2_pair.out)
        krona_inputs = krona_inputs.concat(kraken_krona_input)

        kreport2list(kraken2_pair.out)
        metacomp(kreport2list.out)

        krona(krona_inputs)
}

workflow classify_taxonomy_nanopore {
    take:
        fastq

    main:
        krona_inputs = Channel.empty()

        Channel.fromPath(params.centrifuge_db_dir, type: 'dir').set{centrifuge_db_dir}
        centrifuge(fastq, centrifuge_db_dir)
        kraken_krona_input = kraken2krona(centrifuge.out)
        krona_inputs = krona_inputs.concat(kraken_krona_input)

        kreport2list(centrifuge.out)
        metacomp(kreport2list.out)

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
        --paired --threads 12 --confidence ${params.kraken2_confidence_threshold}\
        $pe1 $pe2
    """
}

process centrifuge {
    tag "${params.prefix}:centrifuge"

    input:
        path single
        path centrifuge_db_dir
    output:
        path "${params.prefix}.centrifuge.kraken_report.txt"
    """
    centrifuge -x ${centrifuge_db_dir}/${params.centrifuge_db_name} -U $single -S ${params.prefix}.centrifuge.txt
    centrifuge-kreport -x ${centrifuge_db_dir}/${params.centrifuge_db_name} ${params.prefix}.centrifuge.txt > ${params.prefix}.centrifuge.kraken_report.txt
    """
}

process kraken2krona {
    tag "${params.prefix}:kraken2krona"

    input:
        path kraken_report
    output:
        path "${params.prefix}.kraken"
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
    perl /home/user/custom_scripts/convert_krakenRep2list.pl < $kraken_report > ${params.prefix}.list
    """
}

process metacomp {
    tag "${params.prefix}:metacomp"
    
    publishDir "${params.outdir}/classify_taxonomy", mode: 'copy'

    input:
        path kraken_list
    output:
        path "${params.prefix}_family.svg"
        path "${params.prefix}_genus.svg"
        path "${params.prefix}_species.svg"
    """
    Rscript /home/user/custom_scripts/generate_heatmap.R ${kraken_list} ${params.prefix}
    """
}

process krona {
    tag "${params.prefix}:krona"
    publishDir "${params.outdir}/classify_taxonomy", mode: 'copy'

    input:
        path krona_input
    output:
        path "${params.prefix}.krona.html"
    """
    ktImportText $krona_input -o ${params.prefix}.krona.html
    """
}