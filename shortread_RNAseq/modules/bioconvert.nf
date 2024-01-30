#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process bioconvert_fastq {

    input:
        tuple val(name), path(fasta)

    output:
        "*.fq"

    script:
        """

        bioconvert fasta2fastq ${fasta} ${name}.fq
        
        """
}
