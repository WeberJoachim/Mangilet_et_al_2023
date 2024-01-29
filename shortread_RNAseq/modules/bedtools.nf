#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process bedtools_getfasta {

    input:
        path(regions)
        path(genome)

    output:
        path("*fa")

    script:
        """

        bedtools getfasta -fi ${genome} -bed ${regions} -s -fo ${regions.getSimpleName()}.fa

        """
}
