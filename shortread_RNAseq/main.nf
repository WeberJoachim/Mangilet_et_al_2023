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
include {salmon_pre_index}              from "./modules/salmon"
include {salmon_index}                  from "./modules/salmon"
include {salmon_quant}                  from "./modules/salmon"
include {preprocess_genome}		        from "./modules/align"
include {hisat2_index_genome}		    from "./modules/align"
include {hisat2_genome_align}		    from "./modules/align"
include {sam_to_bam}	                from "./modules/sam_to_bam"
include {process_transcriptome_gtf}     from "./modules/rmats"
include {rmats_U1_C}                    from "./modules/rmats"
include {rmats_U1_70K}                  from "./modules/rmats"




workflow {
        main:     
            // Input
            reads                       = Channel.fromFilePairs(params.reads)
	        latest_transcriptome        = Channel.fromPath(params.latest_transcriptome)
            gene_types_araport11        = Channel.fromPath(params.gene_types_araport11)
	        genome                      = Channel.fromPath(params.genome)
            genome_gff3                 = Channel.fromPath(params.genome_gff3)
            samplesheet                 = Channel.fromPath(params.samplesheet)
            latest_transcriptome_gtf    = Channel.fromPath(params.latest_transcriptome_gtf)


            /// Workflow
            // Initial Qualitycontrol
            fastqc_before(reads)
            multiqc_before(fastqc_before.out.collect())

            /// Trimming
            trim_galore(reads)
    
            /// filter tRNA and rRNA          
            make_tRNA_rRNA_reference(latest_transcriptome.collect(), gene_types_araport11)
            hisat2_index(make_tRNA_rRNA_reference.out.ref)
            hisat2_align(trim_galore.out.trimmed_reads_1, trim_galore.out.trimmed_reads_2, hisat2_index.out.collect())

            // Qualitycontrol after trimming and filtering
            fastqc_after(hisat2_align.out.filtered_reads_1, hisat2_align.out.filtered_reads_2)
            multiqc_after(fastqc_after.out.collect())


            //quantify
            salmon_pre_index(preprocess_genome.out.collect(), latest_transcriptome.collect())
            salmon_index(salmon_pre_index.out.gentrome, salmon_pre_index.out.decoys)
            salmon_quant(hisat2_align.out.filtered_reads_1, hisat2_align.out.filtered_reads_2, salmon_index.out)

            // align
            preprocess_genome(genome.collect())
            hisat2_index_genome(preprocess_genome.out)
            hisat2_genome_align(hisat2_align.out.filtered_reads_1, hisat2_align.out.filtered_reads_2, hisat2_index_genome.out.collect())
            sam_to_bam(hisat2_genome_align.out.sam)

          	    
            // rmats
   
            process_transcriptome_gtf(latest_transcriptome_gtf)
            rmats_U1_C(process_transcriptome_gtf.out.collect(), sam_to_bam.out.filter(~/^\[WT_.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect(), sam_to_bam.out.filter(~/^\[U1_C.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect())
            rmats_U1_70K(process_transcriptome_gtf.out.collect(), sam_to_bam.out.filter(~/^\[WT_.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect(), sam_to_bam.out.filter(~/^\[U1_70K.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect())

            
	
	    
}           