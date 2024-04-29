#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process salmon_pre_index {
    
    input:
        path(genome)
        path(transcriptome)
        
    output:
        path "gentrome_TAIR10_Atrtd3.fa", emit: gentrome
        path "decoys.txt", emit: decoys

    script:
        """
        grep "^>" <${genome} | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt
        cat ${transcriptome} ${genome} > gentrome_TAIR10_Atrtd3.fa
        """      
}

process salmon_index {

    input:
        path(gentrome)
        path(decoys)
        
    output:
        path "*"
        
    script:
        """
        salmon index -t ${gentrome} -d ${decoys} -p 16 -i salmon_index_k31 -k 31
        """
}


process salmon_quant {
    
    publishDir "$baseDir/results/quant"

    input:
        tuple val(name), path(filtered_reads_1)
        tuple val(name), path(filtered_reads_2)
        path(index)
            
    output:
        path "*"
		
    script:
        """
        salmon quant -l A -i ${index} -1 ${filtered_reads_1} -2 ${filtered_reads_2} --seqBias --gcBias -o ${name}_quant
        """
}
