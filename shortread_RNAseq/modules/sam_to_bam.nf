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



process sort_and_index {
 
    publishDir "results/alignments/sorted_bams/"

    input:
	    tuple val(name), path(bams)
	
    output:
        tuple val(name), path("*sorted.bam"), path("*sorted.bam.bai"), emit: sorted_bam_and_index
	    

    script:
	"""
	samtools sort -@ 4 -m 4G ${bams} > ${name}_sorted.bam
	samtools index -@ 3 ${name}_sorted.bam
	"""
}

