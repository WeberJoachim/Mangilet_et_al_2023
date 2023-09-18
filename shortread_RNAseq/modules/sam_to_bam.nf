#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process sam_to_bam {
    
    input:
        tuple val(name), path(alignment)
        
    output:
        tuple val(name), path("*.bam"), emit: bams
	     
    script:
    """
    samtools view -b -@ 3 -o ${name}.bam ${alignment}
    """
}


