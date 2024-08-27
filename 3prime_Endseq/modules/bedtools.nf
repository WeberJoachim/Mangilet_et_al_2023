#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bedtools_closest {

    publishDir 'results/distances', mode: 'copy'
    label 'bedtools'

    
    input:
        path(intron_gtf)
        path(bed)

    output:
        path "*.out"

    script:
        """
        bedtools intersect -wa -wb -a ${intron_gtf} -b ${bed} > distance_intronboundaries_to_${bed.getSimpleName()}.out
        """
}
