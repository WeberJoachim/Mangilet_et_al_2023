#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {fastqc_before}                                     from "./modules/quality_control_and_preprocessing"
include {fastqc_after}                                      from "./modules/quality_control_and_preprocessing"
include {multiqc_before}                                    from "./modules/quality_control_and_preprocessing"
include {multiqc_after}                                     from "./modules/quality_control_and_preprocessing"
include {trim_galore}                                       from "./modules/quality_control_and_preprocessing"
include {make_tRNA_rRNA_reference}                          from "./modules/quality_control_and_preprocessing"
include {hisat2_index}                                      from "./modules/quality_control_and_preprocessing"
include {hisat2_align}                                      from "./modules/quality_control_and_preprocessing"


include {preprocess_genome}		                    from "./modules/align"
include {hisat2_index_genome}		                    from "./modules/align"
include {hisat2_genome_align}		                    from "./modules/align"

include {sam_to_bam}	                                    from "./modules/sam_to_bam"

include {process_transcriptome_gtf}                         from "./modules/rmats"
include {rmats_roots_ago4_1}                                        from "./modules/rmats"
include {rmats_roots_ago1_27}                                      from "./modules/rmats"
include {rmats_shoots_ago4_1}                                      from "./modules/rmats"
include {rmats_shoots_ago1_27}                                      from "./modules/rmats"
include {sort_and_index}					from "./modules/sam_to_bam"


workflow {
        main:     

        // Input

 
        reads 					= Channel.fromFilePairs(params.reads)
	latest_transcriptome                    = Channel.fromPath(params.latest_transcriptome)
        gene_types_araport11                    = Channel.fromPath(params.gene_types_araport11)
        genome                                  = Channel.fromPath(params.genome)
        genome_gff3                             = Channel.fromPath(params.genome_gff3)
        samplesheet                             = Channel.fromPath(params.samplesheet)
        latest_transcriptome_gtf                = Channel.fromPath(params.latest_transcriptome_gtf)


        intron_gtf				= Channel.fromPath(params.intron_gtf)
        



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

        //preprocess genome
        preprocess_genome(genome.collect())

        //quantify

        //align
        hisat2_index_genome(preprocess_genome.out)
        hisat2_genome_align(hisat2_align.out.filtered_reads_1, hisat2_align.out.filtered_reads_2, hisat2_index_genome.out.collect())
        sam_to_bam(hisat2_genome_align.out.sam)

        //rmats
        process_transcriptome_gtf(latest_transcriptome_gtf)
        rmats_shoots_ago4_1(process_transcriptome_gtf.out.collect(), sam_to_bam.out.filter {it -> it[0] in ['SRR10695012', 'SRR10695013', 'SRR10694996']}.flatten().collect(), sam_to_bam.out.filter {it -> it[0] in ['SRR10695002', 'SRR10695003', 'SRR10695004']}.flatten().collect())

        rmats_shoots_ago1_27(process_transcriptome_gtf.out.collect(),  sam_to_bam.out.filter {it -> it[0] in ['SRR10695012', 'SRR10695013', 'SRR10694996']}.flatten().collect(), sam_to_bam.out.filter {it -> it[0] in ['SRR10695007', 'SRR10695008', 'SRR10695009']}.flatten().collect())

	rmats_roots_ago4_1(process_transcriptome_gtf.out.collect(),  sam_to_bam.out.filter {it -> it[0] in ['SRR10694978', 'SRR10694985', 'SRR10695010']}.flatten().collect(), sam_to_bam.out.filter {it -> it[0] in ['SRR10694999', 'SRR10695000', 'SRR10695001']}.flatten().collect())

	rmats_roots_ago1_27(process_transcriptome_gtf.out.collect(),  sam_to_bam.out.filter {it -> it[0] in ['SRR10694978', 'SRR10694985', 'SRR10695010']}.flatten().collect(), sam_to_bam.out.filter {it -> it[0] in ['SRR10695005', 'SRR10695006', 'SRR10695011']}.flatten().collect())
        
	sort_and_index(sam_to_bam.out)        
         
}          
