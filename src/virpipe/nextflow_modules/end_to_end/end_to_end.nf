nextflow.enable.dsl=2
include { qc_illumina; qc_nanopore } from '../qc/qc'
include { filter_reads_illumina; filter_reads_nanopore } from '../filter/reads/filter_reads'
include { filter_host_illumina; filter_host_nanopore } from '../filter/host/filter_host'
include { map_illumina; map_nanopore } from '../map/map'
include { filter_map } from '../filter/map/filter_map'
include { tax_classify_illumina; tax_classify_nanopore } from '../tax_classify/tax_classify' addParams(tool: 'kraken2')
include { assembly_illumina; assembly_nanopore } from '../assembly/assembly' addParams(tool: params.assembly_tool)
include { polish } from '../polish/polish' addParams(tool: 'racon medaka')
include { filter_contigs } from '../filter/contigs/filter_contigs'
include { blast } from '../post_assembly/blast/post_assembly_blast' addParams(tool: 'blastn megablast')
include { filter_blast } from '../filter/blast/filter_blast'
include { match_taxonomy } from '../report/blast/report_blast'
include { zoonotic_rank } from '../post_assembly/zoonosis/post_assembly_zoonosis' addParams(tool: 'zoonotic_rank')

workflow {

    main:

    if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
        }

    if (params.platform == 'illumina') {
        qc_illumina(fastq_pair)
        filter_reads_illumina(fastq_pair)

        if (params.host_genome != null) {
            Channel.fromPath(params.host_genome).set{host_genome}
            filter_completed = filter_host_illumina(filter_reads_illumina.out, host_genome)
        } else {
            filter_completed = filter_reads_illumina.out
        }

        map_illumina(filter_completed)
        filter_map(map_illumina.out)

        tax_classify_illumina(filter_completed)

        contigs = assembly_illumina(filter_completed)

    } else if (params.platform == 'nanopore') {
        
        qc_nanopore(fastx)
        filter_reads_nanopore(fastx)

        if (params.host_genome) {
            Channel.fromPath(params.host_genome).set{host_genome}
            filter_completed = filter_host_nanopore(filter_reads_nanopore.out, host_genome)
        } else {
            filter_completed = filter_reads_nanopore.out
        }

        map_nanopore(filter_completed)
        filter_map(map_nanopore.out)

        tax_classify_nanopore(filter_completed)

        contigs = assembly_nanopore(filter_completed)
        }
    
    filter_contigs(contigs)

    blast(filter_contigs.out)
    filter_blast(blast.out)
    match_taxonomy(filter_blast.out)

    zoonotic_rank(filter_contigs.out)
}