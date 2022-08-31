nextflow.enable.dsl=2

workflow {
    build_database()
}

process build_database {
    publishDir "${params.taxonomizrdb}", mode: 'move'
    tag "${params.prefix}:build_taxonomizr_database"

    output:
        path "accessionTaxa.sql"
    """
    #!/usr/bin/env Rscript

    library(taxonomizr)
    prepareDatabase('accessionTaxa.sql')
    """
}