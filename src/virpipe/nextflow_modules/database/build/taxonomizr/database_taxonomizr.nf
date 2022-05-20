nextflow.enable.dsl=2

workflow {

    main:
        make_db()
}

process make_db {
    publishDir "${params.outdir}", mode: 'copy'

    output:
        path "accessionTaxa.sql"
    """
    #!/usr/bin/env Rscript

    library(taxonomizr)
    prepareDatabase('accessionTaxa.sql')
    """
}