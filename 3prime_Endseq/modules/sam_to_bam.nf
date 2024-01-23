#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process sam_to_bam {
    
    input:
        path(alignment)
        
    output:
        path("*.bam"), emit: bams
	     
    script:
    """
    samtools view -b -@ 3 -o ${alignment.getSimpleName()}.bam ${alignment}
    """
}



process sort_and_index {
 
    publishDir "results/alignments/sorted_bams/"

    input:
	    path(bams)
	
    output:
        path("*sorted.bam"), emit: sorted_bam
	    path("*sorted.bam.bai"), emit: sorted_bam_index

    script:
	"""
	samtools sort -@ 4 -m 4G ${bams.getSimpleName()}.bam > ${bams.getSimpleName()}_sorted.bam
    samtools index -@ 3 ${bams.getSimpleName()}_sorted.bam
	"""
}
