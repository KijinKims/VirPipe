nextflow.enable.dsl=2
include { qc_illumina; qc_nanopore } from './qc'
include { filter_reads_illumina; filter_reads_nanopore } from './filter_reads'
include { filter_host_illumina; filter_host_nanopore } from './filter_host'
include { map_illumina; map_nanopore } from './map'
include { filter_map } from './filter_map'
include { taxclassify_illumina; taxclassify_nanopore } from './taxclassify' addParams(tool: 'kraken2')
include { assembly_illumina; assembly_nanopore } from './assembly' addParams(tool: params.assembly_tool)
include { polish } from './polish' addParams(tool: 'racon medaka')
include { filter_contigs } from './filter_contigs'
include { blast } from './postassembly_blast' addParams(tool: 'blastn megablast')
include { filter_blast } from './filter_blast'
include { match_taxonomy } from './report_blast'
include { zoonotic_rank } from './postassembly_zoonosis' addParams(tool: 'zoonotic_rank')

workflow {

    main:

    if (params.platform == 'illumina') {
            Channel.fromPath([params.x, params.x2]).buffer(size:2).set{fastq_pair}
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.x).set{fastx}
        }

    if (params.platform == 'illumina') {

        if(!params.skip_qc){
            qc_illumina(fastq_pair)
        }

        if(!params.skip_filter_read){
            reads = filter_reads_illumina(fastq_pair)
        }
        else{
            reads = fastq_pair
        }

        if (params.host_genome == null | params.skip_host_genome) {
            filter_completed = filter_reads_illumina.out
        } else {
            Channel.fromPath(params.host_genome).set{host_genome}
            filter_completed = filter_host_illumina(filter_reads_illumina.out, host_genome)
        }

        if(!params.skip_map){
            map_illumina(filter_completed)
            filter_map(map_illumina.out)
        }

        if(!params.skip_taxclassify){
            taxclassify_illumina(filter_completed)
        }

        if(!params.skip_assembly){
            contigs = assembly_illumina(filter_completed)
            filter_contigs(contigs)

            if(!params.skip_blast){
                blast(filter_contigs.out)
                filter_blast(blast.out)
                match_taxonomy(filter_blast.out)
            }

            if(!params.skip_zoonosis){
                zoonotic_rank(filter_contigs.out)
            }
        }

    } else if (params.platform == 'nanopore') {
        
        if(!params.skip_qc){
            qc_nanopore(fastx)
        }

        if(!params.skip_filter_read){
            reads = filter_reads_nanopore(fastx)
        }
        else{
            reads = fastq_pair
        }

        if (params.host_genome == null | params.skip_host_genome) {
            filter_completed = filter_reads_nanopore.out
        } else {
            Channel.fromPath(params.host_genome).set{host_genome}
            filter_completed = filter_host_nanopore(filter_reads_nanopore.out, host_genome)
        }

        if(!params.skip_map){
            map_nanopore(filter_completed)
            filter_map(map_nanopore.out)
        }

        if(!params.skip_taxclassify){
            taxclassify_nanopore(filter_completed)
        }

        if(!params.skip_assembly){
            contigs = assembly_nanopore(filter_completed)

            if(!params.skip_polish){
                polished_contigs = polish(contigs)
            }
            else {
                polished_contigs = contigs
            }

            filter_contigs(polished_contigs, fastx)
            if(!params.skip_blast){
                blast(filter_contigs.out)
                filter_blast(blast.out)
                match_taxonomy(filter_blast.out)
            }

            if(!params.skip_zoonosis){
                zoonotic_rank(filter_contigs.out)
            }
        }
    }
    

}