nextflow.enable.dsl=2

workflow {
    main:
        infile = Channel.fromPath(params.in)
        accession2taxid = Channel.fromPath(params.accession2taxid)

        generate_taxid_map(accession2taxid)
        makeblastdb(infile, generate_taxid_map.out)
}

process generate_taxid_map {
    stageInMode 'copy'

    input:
        path accession2taxid
    output:
        path "taxidmapfile"
    """
    gunzip -c $accession2taxid | sed '1d' | awk '{print \$2" "\$3}' > taxidmapfile
    """
}

process makeblastdb {
    publishDir "${params.outdir}" , mode: 'copy'

    input:
        path infile
        path taxidmap
    output:
        path "*"
    script:
    if (params.parse_seqids)
    """
    makeblastdb -in ${infile} -parse_seqids -blastdb_version 5 -taxid_map ${taxidmap} -title ${params.title} -dbtype ${params.dbtype} -out ${params.dbname}
    """
    else
    """
    makeblastdb -in ${infile} -blastdb_version 5 -taxid_map ${taxidmap} -title ${params.title} -dbtype ${params.dbtype} -out ${parmas.dbname}
    """
}