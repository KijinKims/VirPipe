nextflow.enable.dsl=2
include { qc_illumina; qc_nanopore } from './qc'
include { preprocess_illumina; preprocess_nanopore } from './preprocess'
include { remove_host_illumina; remove_host_nanopore } from './remove_host'
include { map_illumina; map_nanopore } from './map'
include { classify_taxonomy_illumina; classify_taxonomy_nanopore } from './classify_taxonomy'
include { assemble_illumina; assemble_nanopore } from './assemble' addParams(tool: params.assembly_tool)
include { polish } from './polish'
include { blast } from './blast'
include { zoonotic_rank } from './zoonotic_rank'

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
            preprocess_completed = reads
        } else {
            Channel.fromPath(params.host_genome).set{host_genome}
            preprocess_completed = remove_host_illumina(reads, host_genome)
        }

        if(!params.skip_map){
            map_illumina(preprocess_completed)
        }

        if(!params.skip_classify_taxonomy){
            classify_taxonomy_illumina(preprocess_completed)
        }

        if(!params.skip_assemble){
            contigs = assemble_illumina(preprocess_completed)

            if(!params.skip_blast){
                blast(contigs)
            }

            if(params.do_zoonotic_rank){
                zoonotic_rank(contigs)
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
            reads = fastq
        }

        if (params.host_genome == null | params.skip_remove_host) {
            preprocess_completed = reads
        } else {
            Channel.fromPath(params.host_genome).set{host_genome}
            preprocess_completed = remove_host_nanopore(reads, host_genome)
        }

        if(!params.skip_map){
            map_nanopore(preprocess_completed)
        }

        if(!params.skip_classify_taxonomy){
            classify_taxonomy_nanopore(preprocess_completed)
        }

        if(!params.skip_assemble){
            contigs = assemble_nanopore(preprocess_completed)

            if(!params.skip_polish){
                polished_contigs = polish(contigs, preprocess_completed)
            }
            else {
                polished_contigs = contigs
            }

            if(!params.skip_blast){
                blast(polished_contigs)
            }

            if(params.do_zoonotic_rank){
                zoonotic_rank(polished_contigs)
            }
        }
    }
}