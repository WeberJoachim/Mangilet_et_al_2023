#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process python_add_introns {

    publishDir 'results/', mode: 'copy', overwrite: false
    label 'python'

    input:
        path(gtf)
    
    output:
        path("*_introns.gtf"), emit: intron_gtf

    script:
        """

        add_introns.py ${gtf} ${gtf.getSimpleName()}_introns.gtf
        

        """
}
