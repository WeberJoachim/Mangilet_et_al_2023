#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process preprocess_genome {

    input:
        path(genome)

    output:
        path "*.fas"

    script:
        """
        sed -i s/Chr1/1/g ${genome}
        sed -i s/Chr2/2/g ${genome}
        sed -i s/Chr3/3/g ${genome}
        sed -i s/Chr4/4/g ${genome} 
        sed -i s/Chr5/5/g ${genome}
        sed -i s/ChrC/chloroplast/g ${genome}
        sed -i s/ChrM/mitochondria/g ${genome}

       	seqkit replace -p " TPA_asm: Arabidopsis thaliana ecotype Col-0 mitochondrion, complete genome, " -r " " ${genome} | seqkit seq -w 80 > ${genome.getSimpleName()}_mod.fas
        """
}

process hisat2_index_genome {
    
    input:     
        path(ref)

    output:
        path "*" 

    script:
        """       
        hisat2-build -p 16 ${ref} index_hisat_${ref.getSimpleName()}
        """
}

process hisat2_genome_align {

	publishDir "results/align/", pattern: "*.log"
   	 
    input:
        path(trimmed_reads)
	    path(index) 
       
    output:
        path("*.sam"), emit: sam
        path("*.log"), emit: logs
          
    script:
        """            
        hisat2 -p 16 --min-intronlen 20 --max-intronlen 1000 -x ${index[0].getSimpleName()} -U ${trimmed_reads} -S ${trimmed_reads.getSimpleName()}.sam 2> ${trimmed_reads.getSimpleName()}.log
        """
}





