#!/bin/bash


#run this after the nextflow pipeline finished
cat results/extracted_motif_info/collapsed_introns_without_pA*.bed              >> collapsed_introns_without_pA
cat results/extracted_motif_info/distal_pA_where_comp_rep_u1c_new*.bed          >> distal_pA_where_comp_rep_u1c_new
cat results/extracted_motif_info/distal_pA_where_comp_rep_u170k_new*.bed        >> distal_pA_where_comp_rep_u170k_new
cat results/extracted_motif_info/non_affected_WT_consensus_pA_sites_new*.bed    >> non_affected_WT_consensus_pA_sites_new
cat results/extracted_motif_info/proximal_pA_comp_rep_u1c_new*.bed              >> proximal_pA_comp_rep_u1c_new
cat results/extracted_motif_info/proximal_pA_comp_rep_u170k_new*.bed            >> proximal_pA_comp_rep_u170k_new




