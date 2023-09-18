#!/usr/bin/env nextflow

nextflow.enable.dsl=2




process macs2_peakcall {
    label "peakcall"

    input:
        tuple val(name_treatment), path(sorted_indexed_bam_treatment), path(sorted_indexed_bam_treatment_index)
        tuple val(name_control), path(sorted_indexed_bam_control), path(sorted_indexed_bam_control_index) 
       
    output:
        tuple val(name_control), path("*_control_lambda.bdg"), emit: bedgraph_control
        tuple val(name_treatment), path("*_treat_pileup.bdg"), emit: bedgraph_treatment
        tuple val(name_treatment), path("*.narrowPeak"), emit: peaks
        tuple val(name_treatment), path("*summits.bed"), emit: summits


    script:
        """     
        macs2 callpeak -t ${sorted_indexed_bam_treatment} -c ${sorted_indexed_bam_control} -f BAMPE -n scaled_${name_treatment} --gsize 1.25e8 -B --SPMR --outdir .
        """
}


process macs2_bdgcmp_pre {

label "bdgcmp"

    input:
        tuple val(name_treatment), path(pileup_treatment) 
        tuple val(name_control), path(pileup_control) 
       
    output:
        tuple val(name_treatment), path("*.bdg"), emit: bedgraph 
        path("*.log")

    script:
        """     
        macs2 bdgcmp -t ${pileup_treatment} -c ${pileup_control} -o ./${name_treatment}_scaled_FE_input.bdg -m FE 2> ${name_treatment}_FE_Input.log
        """




}



process macs2_bdgcmp_post {

label "bdgcmp"

    input:
        tuple val(name_treatment), path(pileup_treatment) 
        tuple val(name_control), path(pileup_control) 
       
    output:
        tuple val(name_treatment), path("*.bdg"), emit: bedgraph 
        path("*.log")

    script:
        """     
        macs2 bdgcmp -t ${pileup_treatment} -c ${pileup_control} -o ./${name_treatment}_subtract_input_scaled_IGG.bdg -m subtract  2> ${name_treatment}_subtract_input_scaled_IGG.log
        """




}
