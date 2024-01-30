#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bedtools_getfasta {

    input:
        tuple val(name), path(regions)
        path(genome)

    output:
        tuple val(name), path("*fa")

    script:
        """

        bedtools getfasta -fi ${genome} -bed ${regions} -s -fo ${name}.fa

        """
}
