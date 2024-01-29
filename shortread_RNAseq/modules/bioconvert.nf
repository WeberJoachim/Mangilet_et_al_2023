#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process bioconvert_fastq {

    input:
        path(fasta)

    output:
        "*"

    script:
        """

        bioconvert fasta2fastq ${fasta} ${fasta.getSimpleName()}.fq
        
        """
}
