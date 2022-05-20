nextflow.enable.dsl=2

workflow {
    main:
        Channel.fromPath(params.x).set{blast_out}
        match_taxonomy(blast_out)
}

workflow match_taxonomy {
    take:
        blast_out
    main:
        Channel.fromPath(params.taxonomizr_db).set{taxonomizr_db}
        matched = match_taxonomy_process(blast_out.combine(taxonomizr_db))
    
}

process match_taxonomy_process {
    tag "${params.prefix}:match_taxonomy"

    publishDir "${params.outdir}/post_assembly/blast", mode: 'copy'

    input:
        tuple path(blast_out), path(taxonomizr_db)
    output:
        path "*"
    """
    Rscript ~/match_taxonomy.R $blast_out ${blast_out.baseName}.html ${taxonomizr_db}
    """
}