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
        tuple val(name), path(trimmed_reads_1)
	    tuple val(name), path(trimmed_reads_2)
        path(index) 
       
    output:
        tuple val(name), path("*1.fq.gz"), emit: filtered_reads_1
	tuple val(name), path("*2.fq.gz"), emit: filtered_reads_2
        path("*.log"), emit: logs
          
    script:
        """            
        hisat2 -p 16 --min-intronlen 20 --max-intronlen 1000 -x ${index[0].getSimpleName()} -1 ${trimmed_reads_1} -2 ${trimmed_reads_2} --un-conc-gz ${name}_R%.fq.gz 2> ${name}.log
        """
}




process preprocess_and_find_intersecting_feature_bed {

    label "preprocess_bed"

    input:
        tuple val(name), path(regions)
        path(intron_gtf)

    output:
        tuple val(name), path("*_new.bed"), emit: new_bed

    script:
        """
        bedtools intersect -s -wa -a ${intron_gtf} -b ${background} > ${name}_new.bed

        sed -i s/Chr1/1/g ${name}_new.bed
        sed -i s/Chr2/2/g ${name}_new.bed
        sed -i s/Chr3/3/g ${name}_new.bed
        sed -i s/Chr4/4/g ${name}_new.bed
        sed -i s/Chr5/5/g ${name}_new.bed
        sed -i s/ChrC/chloroplast/g ${name}_new.bed
        sed -i s/ChrM/mitochondria/g ${name}_new.bed

        """
}


process preprocess_and_extend_bed {

    label "preprocess_bed"

    input:
        tuple val(name), path(positions)

    output:
        tuple val(name), path("*_ext.bed"), emit: extended_bed

    script:
        """

        awk -v OFS='\\t' '
        {
            if (\$6 == "-") {
                \$2 = \$2 - 49
                \$3 = \$3 + 50
            } else {
                \$2 = \$2 - 50
                \$3 = \$3 + 49
            }
            print \$0

        }' "${positions}" > "${name}_ext.bed"

        sed -i s/Chr1/1/g ${name}_ext.bed
        sed -i s/Chr2/2/g ${name}_ext.bed
        sed -i s/Chr3/3/g ${name}_ext.bed
        sed -i s/Chr4/4/g ${name}_ext.bed
        sed -i s/Chr5/5/g ${name}_ext.bed
        sed -i s/ChrC/chloroplast/g ${name}_ext.bed
        sed -i s/ChrM/mitochondria/g ${name}_ext.bed

        """
}