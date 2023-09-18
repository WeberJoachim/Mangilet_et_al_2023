#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process fastqc_before {
    label "fastqc"

    input:
        path(reads) 
       
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
        path(trimmed_reads)
	
       
    output:
        path "*_fastqc.zip" 

    script:
        """     
        fastqc -t 4 -q ${trimmed_reads}
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
        path(reads) 
       
    output:
        path("*.fq.gz"), emit: trimmed_reads
    	path "*.txt", emit: logs
      
    script:
        """
        trim_galore --quality 20 --phred33 ${reads} --cores 4 --basename ${reads.getSimpleName()}
        """
}


process make_tRNA_rRNA_reference {

    input:
        path(latest_transcriptome)
        path(gene_types_araport11)

    output:
        path "tRNA_rRNA_reference_Atrtd3.fa", emit: ref

    script:
        """
        grep pre_trna ${gene_types_araport11} | cut -f 1 > filter_these_rnas.txt
        grep ribosomal_rna ${gene_types_araport11} | cut -f 1 >> filter_these_rnas.txt


        seqkit grep -r -f filter_these_rnas.txt ${latest_transcriptome} -o tRNA_rRNA_reference_Atrtd3.fa
        """    
}


process hisat2_index {
    
    input:     
        path(ref)

    output:
        path "*" 

    script:
        """       
        hisat2-build -p 16 ${ref} index_hisat_${ref.getSimpleName()}
        """
}

process hisat2_align {

	publishDir "results/filtering/", pattern: "*.log"
   	 
    input:
        path(trimmed_reads)
	    path(index) 
       
    output:
        path("*fq.gz"), emit: filtered_reads
	    path("*.log"), emit: logs
          
    script:
        """            
        hisat2 -p 16 --min-intronlen 20 --max-intronlen 1000 -x ${index[0].getSimpleName()} -U ${trimmed_reads} --un-gz ${trimmed_reads.getSimpleName()}_un.fq.gz 2> ${trimmed_reads.getSimpleName()}.log
        """
}
