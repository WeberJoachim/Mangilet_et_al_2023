#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process fastqc_before {
    label "fastqc"

    input:
        tuple val(name), path(reads) 
       
    output:
        path "*_fastqc.zip" 

    script:
        """     
        fastqc -t 4 -q ${reads} 
        """
}



process multiqc_before {

    publishDir = "$baseDir/results/qc"
    label "multiqc"

    input:
        path "*_fastqc.zip" 
       
    output:
        path "*multiqc_report.html"
        
    script:
        """              
        multiqc . -n "raw_multiqc_report.html"
        """
}

process fastqc_after {
    label "fastqc"

    input:
        tuple val(name), path(trimmed_reads_1)
	    tuple val(name), path(trimmed_reads_2) 
       
    output:
        path "*_fastqc.zip" 

    script:
        """     
        fastqc -t 4 -q ${trimmed_reads_1} ${trimmed_reads_2} 
        """
}

process multiqc_after {

    publishDir = "$baseDir/results/qc"
    label "multiqc"
    
    input:
        path "*_fastqc.zip" 
       
    output:
        path "*multiqc_report.html"
        
    script:
        """              
        multiqc . -n "after_processing_multiqc_report.html"
        """
}




process trim_galore {
    
    publishDir "$baseDir/results/trimming", pattern: "*.txt", saveAs: {filename -> "${name}_trimming.log"}
 
    input:
        tuple val(name), path(reads) 
       
    output:
        tuple val(name), path("*_val_1.fq.gz"), emit: trimmed_reads_1
    	tuple val(name), path("*_val_2.fq.gz"), emit: trimmed_reads_2
        path "*.txt", emit: logs
      
    script:
        """
        trim_galore --quality 20 --phred33 ${reads} --cores 4 --paired --basename ${name}
        """
}


