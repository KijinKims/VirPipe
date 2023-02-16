nextflow.enable.dsl=2
include { qc_illumina; qc_nanopore } from './qc'
include { preprocess_illumina; preprocess_nanopore } from './preprocess'
include { remove_host_illumina; remove_host_nanopore } from './remove_host'
include { map_illumina; map_nanopore } from './map'
include { taxclassify_illumina; taxclassify_nanopore } from './taxclassify'
include { assembly_illumina; assembly_nanopore } from './assembly' addParams(tool: params.assembly_tool)
include { polish } from './polish'
include { blast, match_taxonomy } from './blast' addParams(tool: 'blastn megablast')
include { zoonotic_rank } from './zoonosis'

workflow {

    main:

    if (params.platform == 'illumina') {
            Channel.fromPath([params.fastq, params.fastq2]).buffer(size:2).set{fastq_pair}
        } else if (params.platform == 'nanopore') {
            Channel.fromPath(params.fastq).set{fastq}
        }

    if (params.platform == 'illumina') {

        if(!params.skip_qc){
            qc_illumina(fastq_pair)
        }

        if(!params.skip_preprocess){
            reads = preprocess_illumina(fastq_pair)
        }
        else{
            reads = fastq_pair
        }

        if (params.host_genome == null | params.skip_remove_host) {
            preprocess_completed = preprocess_illumina.out
        } else {
            Channel.fromPath(params.host_genome).set{host_genome}
            preprocess_completed = remove_host_illumina(preprocess_illumina.out, host_genome)
        }

        if(!params.skip_map){
            map_illumina(preprocess_completed)
        }

        if(!params.skip_taxclassify){
            taxclassify_illumina(preprocess_completed)
        }

        if(!params.skip_assembly){
            contigs = assembly_illumina(preprocess_completed)

            if(!params.skip_blast){
                blast(contigs.out)
                match_taxonomy(blast.out)
            }

            if(!params.skip_zoonosis){
                zoonotic_rank(contigs.out)
            }
        }

    } else if (params.platform == 'nanopore') {
        
        if(!params.skip_qc){
            qc_nanopore(fastq)
        }

        if(!params.skip_preprocess){
            reads = preprocess_nanopore(fastq)
        }
        else{
            reads = fastq_pair
        }

        if (params.host_genome == null | params.skip_remove_host) {
            preprocess_completed = preprocess_nanopore.out
        } else {
            Channel.fromPath(params.host_genome).set{host_genome}
            preprocess_completed = remove_host_nanopore(preprocess_nanopore.out, host_genome)
        }

        if(!params.skip_map){
            map_nanopore(preprocess_completed)
        }

        if(!params.skip_taxclassify){
            taxclassify_nanopore(preprocess_completed)
        }

        if(!params.skip_assembly){
            contigs = assembly_nanopore(preprocess_completed)

            if(!params.skip_polish){
                polished_contigs = polish(contigs)
            }
            else {
                polished_contigs = contigs
            }

            if(!params.skip_blast){
                blast(polished_contigs.out)
                match_taxonomy(blast.out)
            }

            if(!params.skip_zoonosis){
                zoonotic_rank(polished_contigs.out)
            }
        }
    }
}