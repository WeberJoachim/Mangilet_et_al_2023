#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process python_add_introns {

    publishDir 'results/'
    label 'python'

    input:
        path(gtf)
    
    output:
        "*"

    script:
        """

        python3 add_introns.py ${gtf} ${gtf.getSimpleName()}_introns.gtf

        """
}