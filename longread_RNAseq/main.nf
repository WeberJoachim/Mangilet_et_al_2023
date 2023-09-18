#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {basecall}                      from "./modules/flep"
include {convert_fasta}                 from "./modules/flep"
include {align_to_genome}               from "./modules/flep"
include {merge_bams as merge_bams1}     from "./modules/flep" 
include {merge_bams as merge_bams2}     from "./modules/flep" 
include {merge_bams as merge_bams3}     from "./modules/flep" 


workflow {

        main:     
            // Input
            read_folder = Channel.fromPath(params.reads_folder_fast5, type: 'dir').map{
	        tuple(it.toString().split('/')[6], it)
    	    }
            genome      = Channel.fromPath(params.genome)
           

           
            /// Workflow
            basecall(read_folder)
	        convert_fasta(basecall.out)
            align_to_genome(convert_fasta.out, genome.collect())
            merge_bams1(align_to_genome.out.filter(~/.*PAG73257.*/).collect())
            merge_bams2(align_to_genome.out.filter(~/.*PAG71173.*/).collect())
            merge_bams3(align_to_genome.out.filter(~/.*PAG72817.*/).collect())
           
}