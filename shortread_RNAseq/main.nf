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

include {salmon_pre_index}                                  from "./modules/salmon"
include {salmon_index}                                      from "./modules/salmon"
include {salmon_quant}                                      from "./modules/salmon"

include {preprocess_genome}		                            from "./modules/align"
include {hisat2_index_genome}		                        from "./modules/align"
include {hisat2_genome_align}		                        from "./modules/align"

include {sam_to_bam}	                                    from "./modules/sam_to_bam"

include {process_transcriptome_gtf}                         from "./modules/rmats"
include {rmats_U1_C}                                        from "./modules/rmats"
include {rmats_U1_70K}                                      from "./modules/rmats"


include {homer_findMotifsGenome_custom_background}          from "./modules/homer"
include {homer_findMotifsGenome_no_background}              from "./modules/homer"

include {preprocess_and_extend_bed}                         from "./modules/quality_control_and_preprocessing"
include {preprocess_and_extend_bed_custom_background}       from "./modules/quality_control_and_preprocessing"

include {bedtools_getfasta}                                 from "./modules/bedtools"
include {bedtools_getfasta as bedtools_getfasta2}           from "./modules/bedtools"

include {bioconvert_fastq}                                  from "./modules/bioconvert"
include {bioconvert_fastq as bioconvert_fastq2}             from "./modules/bioconvert"

include {homer_buildMotif_AAUAAA}                           from "./modules/homer"
include {homer_buildMotif_UUGUUU}                           from "./modules/homer"

include {homer_countMotifs}                                 from "./modules/homer"
include {homer_countMotifs as homer_countMotifs2}           from "./modules/homer"

include {homer_buildMotif_cooccurence}                      from "./modules/homer"

include {homer_co_countMotifs as homer_co_countMotifs0}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs1}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs2}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs3}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs4}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs5}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs6}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs7}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs8}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs9}     from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs10}    from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs11}    from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs12}    from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs13}    from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs14}    from "./modules/homer"
include {homer_co_countMotifs as homer_co_countMotifs15}    from "./modules/homer"




workflow {
        main:     
            // Input
            reads                               = Channel.fromFilePairs(params.reads)
	        latest_transcriptome                = Channel.fromPath(params.latest_transcriptome)
            gene_types_araport11                = Channel.fromPath(params.gene_types_araport11)
	        genome                              = Channel.fromPath(params.genome)
            genome_gff3                         = Channel.fromPath(params.genome_gff3)
            samplesheet                         = Channel.fromPath(params.samplesheet)
            latest_transcriptome_gtf            = Channel.fromPath(params.latest_transcriptome_gtf)

            bed_background_rnd_intron           = Channel.fromPath(params.bed_background_rnd_intron).map{tuple(it.getSimpleName(), it)}
            bed_background_distal_pAs           = Channel.fromPath(params.bed_background_distal_pAs).map{tuple(it.getSimpleName(), it)}
            bed_prox_rep_pA_composite_u170k     = Channel.fromPath(params.bed_prox_rep_pA_composite_u170k).map{tuple(it.getSimpleName(), it)}
            bed_prox_rep_pA_composite_u1c       = Channel.fromPath(params.bed_prox_rep_pA_composite_u1c).map{tuple(it.getSimpleName(), it)}

            bed_all_introns			    	    = Channel.fromPath(params.bed_all_introns).map{tuple(it.getSimpleName(), it)}
            bed_distal_pA_where_comp_rep_u1c	= Channel.fromPath(params.bed_distal_pA_where_comp_rep_u1c).map{tuple(it.getSimpleName(), it)}
            bed_distal_pA_where_comp_rep_u170k	= Channel.fromPath(params.bed_distal_pA_where_comp_rep_u170k).map{tuple(it.getSimpleName(), it)}
            motifsize                           = Channel.of(params.motifsize)
            a_rich_motifs                       = Channel.of(params.A_rich_motifs)
            u_rich_motifs                       = Channel.of(params.U_rich_motifs)


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
            salmon_pre_index(preprocess_genome.out.collect(), latest_transcriptome.collect())
            salmon_index(salmon_pre_index.out.gentrome, salmon_pre_index.out.decoys)
            salmon_quant(hisat2_align.out.filtered_reads_1, hisat2_align.out.filtered_reads_2, salmon_index.out)

            // align
            hisat2_index_genome(preprocess_genome.out)
            hisat2_genome_align(hisat2_align.out.filtered_reads_1, hisat2_align.out.filtered_reads_2, hisat2_index_genome.out.collect())
            sam_to_bam(hisat2_genome_align.out.sam)

          	    
            // rmats
            process_transcriptome_gtf(latest_transcriptome_gtf)
            rmats_U1_C(process_transcriptome_gtf.out.collect(), sam_to_bam.out.filter(~/^\[WT_.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect(), sam_to_bam.out.filter(~/^\[U1_C.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect())
            rmats_U1_70K(process_transcriptome_gtf.out.collect(), sam_to_bam.out.filter(~/^\[WT_.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect(), sam_to_bam.out.filter(~/^\[U1_70K.*/ ).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect())

            
            //motifanalysis - You can only run this, if you ran the Deseq2 analysis in the .pynb file. You also need to create the BED-files and Fastafiles as described in the .pynb file to be able to run this
            
            background = bed_background_rnd_intron.concat(bed_background_distal_pAs, bed_all_introns, bed_distal_pA_where_comp_rep_u1c, bed_distal_pA_where_comp_rep_u170k).ifEmpty(false)
            regions = bed_prox_rep_pA_composite_u170k.concat(bed_prox_rep_pA_composite_u1c).ifEmpty(false)

            if(regions != false){
                

                homer_buildMotif_AAUAAA()
                homer_buildMotif_UUGUUU()


                if (background == false){

                    preprocess_and_extend_bed(regions)
                    homer_findMotifsGenome_no_background(preprocess_and_extend_bed.out, preprocess_genome.out.collect(), motifsize)
                    bedtools_getfasta(preprocess_and_extend_bed.out, preprocess_genome.out)
                    bioconvert_fastq(bedtools_getfasta.out)

                }else{

                    preprocess_and_extend_bed(regions)
                    preprocess_and_extend_bed_custom_background(background)

                    //cartesian product
                    homer_findMotifsGenome_custom_background(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), preprocess_genome.out.collect(), motifsize.collect())
                    bedtools_getfasta(preprocess_and_extend_bed.out, preprocess_genome.out.collect())
                    bedtools_getfasta2(preprocess_and_extend_bed_custom_background.out, preprocess_genome.out.collect())
                    bioconvert_fastq(bedtools_getfasta.out)
                    bioconvert_fastq2(bedtools_getfasta2.out)

                }

                
                homer_countMotifs( preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), preprocess_genome.out.collect(), homer_buildMotif_AAUAAA.out.motif, homer_buildMotif_AAUAAA.out.mask)
                homer_countMotifs2(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), preprocess_genome.out.collect(), homer_buildMotif_UUGUUU.out.motif, homer_buildMotif_UUGUUU.out.mask)

                homer_buildMotif_cooccurence(a_rich_motifs.flatten().combine(u_rich_motifs.flatten()))
		
		

                homer_co_countMotifs0(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.zero.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs1(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.one.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs2(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.two.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs3(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.three.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs4(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.four.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs5(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.five.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs6(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.six.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs7(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.seven.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs8(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.eight.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs9(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out), homer_buildMotif_cooccurence.out.nine.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs10(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out),homer_buildMotif_cooccurence.out.ten.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs11(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out),homer_buildMotif_cooccurence.out.eleven.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs12(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out),homer_buildMotif_cooccurence.out.twelve.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs13(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out),homer_buildMotif_cooccurence.out.thirteen.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs14(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out),homer_buildMotif_cooccurence.out.fourteen.collect(), preprocess_genome.out.collect())
		homer_co_countMotifs15(preprocess_and_extend_bed.out.combine(preprocess_and_extend_bed_custom_background.out),homer_buildMotif_cooccurence.out.fifteen.collect(), preprocess_genome.out.collect())

            } 
}          
