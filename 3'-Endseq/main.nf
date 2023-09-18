#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {fastqc_before}                 from "./modules/quality_control_and_preprocessing"
include {fastqc_after}                  from "./modules/quality_control_and_preprocessing"
include {multiqc_before}                from "./modules/quality_control_and_preprocessing"
include {multiqc_after}                 from "./modules/quality_control_and_preprocessing"
include {trim_galore}                   from "./modules/quality_control_and_preprocessing"
include {make_tRNA_rRNA_reference}      from "./modules/quality_control_and_preprocessing"
include {hisat2_index}                  from "./modules/quality_control_and_preprocessing"
include {hisat2_align}                  from "./modules/quality_control_and_preprocessing"
include {preprocess_genome}             from "./modules/align"
include {hisat2_index_genome}           from "./modules/align"
include {hisat2_genome_align}           from "./modules/align"
include {sam_to_bam}                    from "./modules/sam_to_bam"
include {sort_and_index}                from "./modules/sam_to_bam"
include {deeptools_multibamsummary}     from "./modules/deeptools"
include {deeptools_plots}               from "./modules/deeptools"



workflow {
        main:     
            // Input
            reads                       = Channel.fromPath(params.reads)
	        latest_transcriptome        = Channel.fromPath(params.latest_transcriptome)
            gene_types_araport11        = Channel.fromPath(params.gene_types_araport11)
	        genome                      = Channel.fromPath(params.genome)
                   


            /// Workflow
            // Initial Qualitycontrol
            fastqc_before(reads)
            multiqc_before(fastqc_before.out.collect())

            /// Trimming
            trim_galore(reads)
    
            /// filter tRNA and rRNA          
            make_tRNA_rRNA_reference(latest_transcriptome, gene_types_araport11)
            hisat2_index(make_tRNA_rRNA_reference.out.ref)
            hisat2_align(trim_galore.out.trimmed_reads, hisat2_index.out.collect())

            // Qualitycontrol after trimming and filtering
            fastqc_after(hisat2_align.out.filtered_reads)
            multiqc_after(fastqc_after.out.collect())

            // align
            preprocess_genome(genome)
            hisat2_index_genome(preprocess_genome.out)
            hisat2_genome_align(hisat2_align.out.filtered_reads, hisat2_index_genome.out.collect())
            
            // sam to bam - sort and index
            sam_to_bam(hisat2_genome_align.out.sam)
            sort_and_index(sam_to_bam.out.bams)
		    
            // deeptools multibamsummary
            deeptools_multibamsummary(sort_and_index.out.sorted_bam.collect(), sort_and_index.out.sorted_bam_index.collect())
            deeptools_plots(deeptools_multibamsummary.out.numpy_array)
	    
            

	
}           