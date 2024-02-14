#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process bioconvert_fastq {
	

    publishDir 'results/fastq_for_nucfreq/'

    input:
        tuple val(name), path(fasta)

    output:
        path("*.fq", emit: fastq)

    script:
        """

        bioconvert fasta2fastq ${fasta} ${name}.fq
        
        """
}
