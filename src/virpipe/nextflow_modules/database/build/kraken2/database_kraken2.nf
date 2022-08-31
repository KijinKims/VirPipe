nextflow.enable.dsl = 1

if (params.sequence) {
    Channel.fromPath(params.sequence.tokenize()).set{seqs}
} else {
    Channel.empty().set{seqs}
}

process build_database {
    publishDir "${params.kraken2dbdir}" , mode: 'move'
    tag "${params.prefix}:build_kraken2_database"

    input:
        path seqList from seqs.collect()
    output:
        path "${params.db}/*"

    script:
    if (params.library)
    """
    kraken2-build --download-taxonomy --db ${params.db} --use-ftp

    kraken2-build --download-library ${params.library} --db ${params.db}

    for seq in ${seqList}
    do
        kraken2-build --add-to-library \$seq --db ${params.db}
    done

    kraken2-build --build --db ${params.db}
    """
    else
    """
    kraken2-build --download-taxonomy --db ${params.db} --use-ftp

    for seq in ${seqList}
    do
        kraken2-build --add-to-library \$seq --db ${params.db}
    done

     kraken2-build --build --db ${params.db}
    """
}