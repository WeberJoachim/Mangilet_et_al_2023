#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {fastqc_before}                                                         from "./modules/quality_control_and_preprocessing"
include {fastqc_after}                                                          from "./modules/quality_control_and_preprocessing"
include {multiqc_before}                                                        from "./modules/quality_control_and_preprocessing"
include {multiqc_after}                                                         from "./modules/quality_control_and_preprocessing"
include {trim_galore}                                                           from "./modules/quality_control_and_preprocessing"

include {preprocess_genome}		                                                from "./modules/align"
include {hisat2_index_genome}		                                            from "./modules/align"
include {hisat2_genome_align}		                                            from "./modules/align"

include {sam_to_bam}	                                                        from "./modules/sam_to_bam"
include {sort_and_index}                                                        from "./modules/sam_to_bam"
include {extractChromsizes}                                                     from "./modules/sam_to_bam"


include {macs2_peakcall as macs2_peakcall_wt_input_rep1}                        from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_wt_input_rep2}                        from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_wt_input_rep3}                        from "./modules/peakcall"

include {macs2_peakcall as macs2_peakcall_u1c_input_rep1}                       from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u1c_input_rep2}                       from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u1c_input_rep3}                       from "./modules/peakcall"

include {macs2_peakcall as macs2_peakcall_u170k_input_rep1}                     from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u170k_input_rep2}                     from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u170k_input_rep3}                     from "./modules/peakcall"

include {macs2_peakcall as macs2_peakcall_wt_igg_input_rep1}                    from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_wt_igg_input_rep2}                    from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_wt_igg_input_rep3}                    from "./modules/peakcall"

include {macs2_peakcall as macs2_peakcall_u1c_igg_input_rep1}                   from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u1c_igg_input_rep2}                   from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u1c_igg_input_rep3}                   from "./modules/peakcall"

include {macs2_peakcall as macs2_peakcall_u170k_igg_input_rep1}                 from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u170k_igg_input_rep2}                 from "./modules/peakcall"
include {macs2_peakcall as macs2_peakcall_u170k_igg_input_rep3}                 from "./modules/peakcall"


include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_wt_rep1}                  from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_wt_rep2}                  from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_wt_rep3}                  from "./modules/peakcall"

include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_u1c_rep1}                 from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_u1c_rep2}                 from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_u1c_rep3}                 from "./modules/peakcall"

include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_u170k_rep1}               from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_u170k_rep2}               from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_treat_input_u170k_rep3}               from "./modules/peakcall"

include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_wt_rep1}                    from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_wt_rep2}                    from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_wt_rep3}                    from "./modules/peakcall"

include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_u1c_rep1}                   from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_u1c_rep2}                   from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_u1c_rep3}                   from "./modules/peakcall"

include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_u170k_rep1}                 from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_u170k_rep2}                 from "./modules/peakcall"
include {macs2_bdgcmp_pre as macs2_bdgcmp_igg_input_u170k_rep3}                 from "./modules/peakcall"

include {macs2_bdgcmp_post as macs2_bdgcmp_post_wt_rep1}                        from "./modules/peakcall"
include {macs2_bdgcmp_post as macs2_bdgcmp_post_wt_rep2}                        from "./modules/peakcall"
include {macs2_bdgcmp_post as macs2_bdgcmp_post_wt_rep3}                        from "./modules/peakcall"

include {macs2_bdgcmp_post as macs2_bdgcmp_post_u1c_rep1}                       from "./modules/peakcall"
include {macs2_bdgcmp_post as macs2_bdgcmp_post_u1c_rep2}                       from "./modules/peakcall"
include {macs2_bdgcmp_post as macs2_bdgcmp_post_u1c_rep3}                       from "./modules/peakcall"

include {macs2_bdgcmp_post as macs2_bdgcmp_post_u170k_rep1}                     from "./modules/peakcall"
include {macs2_bdgcmp_post as macs2_bdgcmp_post_u170k_rep2}                     from "./modules/peakcall"
include {macs2_bdgcmp_post as macs2_bdgcmp_post_u170k_rep3}                     from "./modules/peakcall"


include {bedgraphToBigwig}                                                      from "./modules/deeptools"
include {multibigwigsummary}                                                    from "./modules/deeptools"
include {plotCorrelation}                                                       from "./modules/deeptools"


include {wiggletools_mean_wt}           					                    from "./modules/ucsc"
include {wiggletools_mean_u1c}          					                    from "./modules/ucsc"
include {wiggletools_mean_u170k}        					                    from "./modules/ucsc"


include {bedgraphToBigwig as bdgToBw_WT}                                        from "./modules/deeptools"
include {bedgraphToBigwig as bdgToBw_U1C}                                       from "./modules/deeptools"
include {bedgraphToBigwig as bdgToBw_U170K}                                     from "./modules/deeptools"

include {preprocess_bed}                                                        from "./modules/deeptools"
include {computeMatrix}								                            from "./modules/deeptools"
include {plotProfile}								                            from "./modules/deeptools"


include {sub_bed_apa_repressed_composite}                                       from "./modules/deeptools"
include {computeMatrix as compute_matrix_sub}					                from "./modules/deeptools"
include {plotProfile as plotProfile_sub}					                    from "./modules/deeptools"



workflow {
        main:     
            // Input
            reads       = Channel.fromFilePairs(params.reads)
	        genome      = Channel.fromPath(params.genome)
            gff         = Channel.fromPath(params.gff)
            gene_list   = Channel.fromPath(params.gene_list)


            /// Workflow
            // Initial Qualitycontrol
            fastqc_before(reads)
            multiqc_before(fastqc_before.out.collect())

            /// Trimming
            trim_galore(reads)
    
                 
       
            // Qualitycontrol after trimming and filtering
            fastqc_after(trim_galore.out.trimmed_reads_1,trim_galore.out.trimmed_reads_2)
            multiqc_after(fastqc_after.out.collect())


            // align
            preprocess_genome(genome.collect())
            hisat2_index_genome(preprocess_genome.out)
            hisat2_genome_align(trim_galore.out.trimmed_reads_1, trim_galore.out.trimmed_reads_2, hisat2_index_genome.out.collect())
            
            // sort and index
            sam_to_bam(hisat2_genome_align.out.sam)
            sort_and_index(sam_to_bam.out.bams)

          
            //peakcall for pileups and peaks Treatment vs Input
	        macs2_peakcall_wt_input_rep1(sort_and_index.out.filter(~/^\[WT_Pol_R1.*/), sort_and_index.out.filter(~/^\[WT_Input_R1.*/))
            macs2_peakcall_wt_input_rep2(sort_and_index.out.filter(~/^\[WT_Pol_R2.*/), sort_and_index.out.filter(~/^\[WT_Input_R2.*/))
            macs2_peakcall_wt_input_rep3(sort_and_index.out.filter(~/^\[WT_Pol_R3.*/), sort_and_index.out.filter(~/^\[WT_Input_R3.*/))

            macs2_peakcall_u1c_input_rep1(sort_and_index.out.filter(~/^\[U1C_Pol_R1.*/), sort_and_index.out.filter(~/^\[WT_Input_R1.*/))
            macs2_peakcall_u1c_input_rep2(sort_and_index.out.filter(~/^\[U1C_Pol_R2.*/), sort_and_index.out.filter(~/^\[WT_Input_R2.*/))
            macs2_peakcall_u1c_input_rep3(sort_and_index.out.filter(~/^\[U1C_Pol_R3.*/), sort_and_index.out.filter(~/^\[WT_Input_R3.*/))

            macs2_peakcall_u170k_input_rep1(sort_and_index.out.filter(~/^\[70K_Pol_R1.*/), sort_and_index.out.filter(~/^\[WT_Input_R1.*/))
            macs2_peakcall_u170k_input_rep2(sort_and_index.out.filter(~/^\[70K_Pol_R2.*/), sort_and_index.out.filter(~/^\[WT_Input_R2.*/))
            macs2_peakcall_u170k_input_rep3(sort_and_index.out.filter(~/^\[70K_Pol_R3.*/), sort_and_index.out.filter(~/^\[WT_Input_R3.*/))


            //peakcall for pileups and peaks IGG vs Input
            macs2_peakcall_wt_igg_input_rep1(sort_and_index.out.filter(~/^\[WT_IGG_R1.*/), sort_and_index.out.filter(~/^\[WT_Input_R1.*/))
            macs2_peakcall_wt_igg_input_rep2(sort_and_index.out.filter(~/^\[WT_IGG_R2.*/), sort_and_index.out.filter(~/^\[WT_Input_R2.*/))
            macs2_peakcall_wt_igg_input_rep3(sort_and_index.out.filter(~/^\[WT_IGG_R3.*/), sort_and_index.out.filter(~/^\[WT_Input_R3.*/))

            macs2_peakcall_u1c_igg_input_rep1(sort_and_index.out.filter(~/^\[U1C_IGG_R1.*/), sort_and_index.out.filter(~/^\[WT_Input_R1.*/))
            macs2_peakcall_u1c_igg_input_rep2(sort_and_index.out.filter(~/^\[U1C_IGG_R2.*/), sort_and_index.out.filter(~/^\[WT_Input_R2.*/))
            macs2_peakcall_u1c_igg_input_rep3(sort_and_index.out.filter(~/^\[U1C_IGG_R3.*/), sort_and_index.out.filter(~/^\[WT_Input_R3.*/))

            macs2_peakcall_u170k_igg_input_rep1(sort_and_index.out.filter(~/^\[70K_IGG_R1.*/), sort_and_index.out.filter(~/^\[WT_Input_R1.*/))
            macs2_peakcall_u170k_igg_input_rep2(sort_and_index.out.filter(~/^\[70K_IGG_R2.*/), sort_and_index.out.filter(~/^\[WT_Input_R2.*/))
            macs2_peakcall_u170k_igg_input_rep3(sort_and_index.out.filter(~/^\[70K_IGG_R3.*/), sort_and_index.out.filter(~/^\[WT_Input_R3.*/))



            //bdgcmp Treatment vs Input
            macs2_bdgcmp_treat_input_wt_rep1(macs2_peakcall_wt_input_rep1.out.bedgraph_treatment, macs2_peakcall_wt_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_treat_input_wt_rep2(macs2_peakcall_wt_input_rep2.out.bedgraph_treatment, macs2_peakcall_wt_input_rep2.out.bedgraph_control)
            macs2_bdgcmp_treat_input_wt_rep3(macs2_peakcall_wt_input_rep3.out.bedgraph_treatment, macs2_peakcall_wt_input_rep3.out.bedgraph_control)

            macs2_bdgcmp_treat_input_u1c_rep1(macs2_peakcall_u1c_input_rep1.out.bedgraph_treatment, macs2_peakcall_u1c_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_treat_input_u1c_rep2(macs2_peakcall_u1c_input_rep2.out.bedgraph_treatment, macs2_peakcall_u1c_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_treat_input_u1c_rep3(macs2_peakcall_u1c_input_rep3.out.bedgraph_treatment, macs2_peakcall_u1c_input_rep1.out.bedgraph_control)

            macs2_bdgcmp_treat_input_u170k_rep1(macs2_peakcall_u170k_input_rep1.out.bedgraph_treatment, macs2_peakcall_u170k_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_treat_input_u170k_rep2(macs2_peakcall_u170k_input_rep2.out.bedgraph_treatment, macs2_peakcall_u170k_input_rep2.out.bedgraph_control)
            macs2_bdgcmp_treat_input_u170k_rep3(macs2_peakcall_u170k_input_rep3.out.bedgraph_treatment, macs2_peakcall_u170k_input_rep3.out.bedgraph_control)



            //bdgcmp IGG vs Input
            macs2_bdgcmp_igg_input_wt_rep1(macs2_peakcall_wt_igg_input_rep1.out.bedgraph_treatment, macs2_peakcall_wt_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_igg_input_wt_rep2(macs2_peakcall_wt_igg_input_rep2.out.bedgraph_treatment, macs2_peakcall_wt_input_rep2.out.bedgraph_control)
            macs2_bdgcmp_igg_input_wt_rep3(macs2_peakcall_wt_igg_input_rep3.out.bedgraph_treatment, macs2_peakcall_wt_input_rep3.out.bedgraph_control)

            macs2_bdgcmp_igg_input_u1c_rep1(macs2_peakcall_u1c_igg_input_rep1.out.bedgraph_treatment, macs2_peakcall_u1c_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_igg_input_u1c_rep2(macs2_peakcall_u1c_igg_input_rep2.out.bedgraph_treatment, macs2_peakcall_u1c_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_igg_input_u1c_rep3(macs2_peakcall_u1c_igg_input_rep3.out.bedgraph_treatment, macs2_peakcall_u1c_input_rep1.out.bedgraph_control)

            macs2_bdgcmp_igg_input_u170k_rep1(macs2_peakcall_u170k_igg_input_rep1.out.bedgraph_treatment, macs2_peakcall_u170k_input_rep1.out.bedgraph_control)
            macs2_bdgcmp_igg_input_u170k_rep2(macs2_peakcall_u170k_igg_input_rep2.out.bedgraph_treatment, macs2_peakcall_u170k_input_rep2.out.bedgraph_control)
            macs2_bdgcmp_igg_input_u170k_rep3(macs2_peakcall_u170k_igg_input_rep3.out.bedgraph_treatment, macs2_peakcall_u170k_input_rep3.out.bedgraph_control)



            //second bdgcmp 
            macs2_bdgcmp_post_wt_rep1(macs2_bdgcmp_treat_input_wt_rep1.out.bedgraph, macs2_bdgcmp_igg_input_wt_rep1.out.bedgraph)
            macs2_bdgcmp_post_wt_rep2(macs2_bdgcmp_treat_input_wt_rep2.out.bedgraph, macs2_bdgcmp_igg_input_wt_rep2.out.bedgraph)
            macs2_bdgcmp_post_wt_rep3(macs2_bdgcmp_treat_input_wt_rep3.out.bedgraph, macs2_bdgcmp_igg_input_wt_rep3.out.bedgraph)

            macs2_bdgcmp_post_u1c_rep1(macs2_bdgcmp_treat_input_u1c_rep1.out.bedgraph, macs2_bdgcmp_igg_input_u1c_rep1.out.bedgraph)
            macs2_bdgcmp_post_u1c_rep2(macs2_bdgcmp_treat_input_u1c_rep2.out.bedgraph, macs2_bdgcmp_igg_input_u1c_rep2.out.bedgraph)
            macs2_bdgcmp_post_u1c_rep3(macs2_bdgcmp_treat_input_u1c_rep3.out.bedgraph, macs2_bdgcmp_igg_input_u1c_rep3.out.bedgraph)

            macs2_bdgcmp_post_u170k_rep1(macs2_bdgcmp_treat_input_u170k_rep1.out.bedgraph, macs2_bdgcmp_igg_input_u170k_rep1.out.bedgraph)
            macs2_bdgcmp_post_u170k_rep2(macs2_bdgcmp_treat_input_u170k_rep2.out.bedgraph, macs2_bdgcmp_igg_input_u170k_rep2.out.bedgraph)
            macs2_bdgcmp_post_u170k_rep3(macs2_bdgcmp_treat_input_u170k_rep3.out.bedgraph, macs2_bdgcmp_igg_input_u170k_rep3.out.bedgraph)

           

            //bedgraphs to bigwigs           
	        extractChromsizes(preprocess_genome.out.collect())
	        bedgraphToBigwig(macs2_bdgcmp_post_wt_rep1.out.bedgraph | combine(macs2_bdgcmp_post_wt_rep2.out.bedgraph) | combine(macs2_bdgcmp_post_wt_rep3.out.bedgraph) | combine(macs2_bdgcmp_post_u1c_rep1.out.bedgraph) | combine(macs2_bdgcmp_post_u1c_rep2.out.bedgraph) | combine(macs2_bdgcmp_post_u1c_rep3.out.bedgraph) | combine(macs2_bdgcmp_post_u170k_rep1.out.bedgraph) | combine(macs2_bdgcmp_post_u170k_rep2.out.bedgraph) | combine(macs2_bdgcmp_post_u170k_rep3.out.bedgraph) |flatten |collate(2) | toSortedList({a,b -> a[0] <=> b[0]}) |flatten| collate(2), extractChromsizes.out.collect())

            //qualitycontrol
            multibigwigsummary(bedgraphToBigwig.out.collect())
            plotCorrelation(multibigwigsummary.out.npz)
	    
            //merge bigwigs (has bedgraph as output)
            wiggletools_mean_wt(bedgraphToBigwig.out.filter(~/.*WT_Pol_R.*/).collect())
            wiggletools_mean_u1c(bedgraphToBigwig.out.filter(~/.*U1C_Pol_R[23].*/).collect())
            wiggletools_mean_u170k(bedgraphToBigwig.out.filter(~/.*70K_Pol_R.*/).collect())

            //convert begraphs from merging to bigwigs
            bdgToBw_WT(wiggletools_mean_wt.out, extractChromsizes.out.collect())
            bdgToBw_U1C(wiggletools_mean_u1c.out, extractChromsizes.out.collect())
            bdgToBw_U170K(wiggletools_mean_u170k.out, extractChromsizes.out.collect())

            //plotProfile over all genes
            preprocess_bed(gff)
            computeMatrix(bdgToBw_WT.out, bdgToBw_U1C.out, bdgToBw_U170K.out, preprocess_bed.out.collect())
            plotProfile(computeMatrix.out.matrix)

            //plotProfile over subset
            sub_bed_apa_repressed_composite(preprocess_bed.out, gene_list.collect())
            compute_matrix_sub(bdgToBw_WT.out, bdgToBw_U1C.out, bdgToBw_U170K.out, sub_bed_apa_repressed_composite.out.collect())
            plotProfile_sub(compute_matrix_sub.out.matrix)
	    
}           