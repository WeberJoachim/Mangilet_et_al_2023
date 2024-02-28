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



include {preprocess_and_extend_bed}                         from "./modules/quality_control_and_preprocessing"
include {preprocess_and_find_intersecting_feature_bed}      from "./modules/quality_control_and_preprocessing"

include {getnonpolyAintrons}				    from "./modules/quality_control_and_preprocessing"



include {bedtools_getfasta}                                 from "./modules/bedtools"

include {bioconvert_fastq}                                  from "./modules/bioconvert"

include {homer_buildMotif_AAUAAA}                           from "./modules/homer"
include {homer_buildMotif_UUGUUU}                           from "./modules/homer"
include {homer_buildMotif_UGUA}                             from "./modules/homer"
include {homer_buildMotif_YA}                               from "./modules/homer"

include {homer_count_coMotifs}                              from "./modules/homer"
include {R_extract_infos_motif}				    from "./modules/homer"




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
            intron_gtf				= Channel.fromPath(params.intron_gtf)
	    bed_all_pA_plantapadb		= Channel.fromPath(params.pAs_plantapadb)
            denovo                              = Channel.of(params.denovo)
	    rscript				= Channel.fromPath(params.rscript)

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

            
            //motifanalysis - You can only run this, if you ran the Deseq2 analysis in the .pynb file. You also need to create the BED-files and Fastafiles as described in the .pynb file from the 3' end seq part of this study to be able to run this
            
            regions = bed_prox_rep_pA_composite_u170k.concat(bed_prox_rep_pA_composite_u1c, bed_background_distal_pAs,  bed_distal_pA_where_comp_rep_u1c, bed_distal_pA_where_comp_rep_u170k).ifEmpty(false)
 

            if(regions != false){
                

                preprocess_and_extend_bed(regions)
                bedtools_getfasta(preprocess_and_extend_bed.out, preprocess_genome.out.collect())
                bioconvert_fastq(bedtools_getfasta.out)
               

                preprocess_and_find_intersecting_feature_bed(regions, intron_gtf.collect())

		getnonpolyAintrons(intron_gtf.collect(), bed_all_pA_plantapadb.collect(), preprocess_and_find_intersecting_feature_bed.out.filter(~/.*u1c.*/).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect(), preprocess_and_find_intersecting_feature_bed.out.filter(~/.*u170k.*/).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect(), preprocess_and_find_intersecting_feature_bed.out.filter(~/.*WT_consensus.*/).toSortedList({a,b -> a[0] <=> b[0] }).flatten().collect())
		

                homer_buildMotif_AAUAAA()
                homer_buildMotif_UUGUUU()
                homer_buildMotif_UGUA()
                homer_buildMotif_YA()

                homer_count_coMotifs(preprocess_and_find_intersecting_feature_bed.out.concat(getnonpolyAintrons.out), homer_buildMotif_AAUAAA.out.motif, homer_buildMotif_UUGUUU.out.motif, homer_buildMotif_UGUA.out, homer_buildMotif_YA.out, preprocess_genome.out.collect())
               
		R_extract_infos_motif(homer_count_coMotifs.out, rscript.collect())
            } 
}          
