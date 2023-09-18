#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process process_transcriptome_gtf {

    input:
        path(gtf)
        
        
    output:
        path "*.gtf"
        
	     
    script:
        """
        sed -i s/^Chr1/1/g ${gtf}
        sed -i s/^Chr2/2/g ${gtf}
        sed -i s/^Chr3/3/g ${gtf}
        sed -i s/^Chr4/4/g ${gtf}
        sed -i s/^Chr5/5/g ${gtf}
        sed -i s/^ChrC/chloroplast/g ${gtf}
        sed -i s/^ChrM/mitochondria/g ${gtf}

	
        """

}


process rmats_U1_C {

    label "rmats"
    publishDir "results/rmats/U1C"

    input:
        path(tair10_gtf)
        tuple val(name1), path(bam1, stageAs: "WT_1.bam"), val(name2), path(bam2, stageAs: "WT_2.bam"), val(name3), path(bam3, stageAs: "WT_3.bam")

        tuple val(name4), path(bam4, stageAs: "U1C_1.bam"), val(name5), path(bam5, stageAs: "U1C_2.bam"), val(name6), path(bam6, stageAs: "U1C_3.bam")


        
    output:
        path "*"
        
	     
    script:
        """
	
	echo "WT_1.bam,WT_2.bam,WT_3.bam" > wt.txt
	echo "U1C_1.bam,U1C_2.bam,U1C_3.bam" > u1c.txt	


        rmats.py --b1 wt.txt --b2 u1c.txt --gtf ${tair10_gtf} -t paired --readLength 133 --variable-read-length --nthread 8 --novelSS --mil 20 --mel 15000 --od WT_vs_U1C --tmp temp
        """
}


process rmats_U1_70K {

    label "rmats"
    publishDir "results/rmats/U170K"

    input:
        path(tair10_gtf)
	tuple val(name1), path(bam1, stageAs: "WT_1.bam"), val(name2), path(bam2, stageAs: "WT_2.bam"), val(name3), path(bam3, stageAs: "WT_3.bam") 
        tuple val(name4), path(bam4, stageAs: "U170K_1.bam"), val(name5), path(bam5, stageAs: "U170K_2.bam"), val(name6), path(bam6, stageAs: "U170K_3.bam") 
        
    output:
        path "*"
        
	     
    script:
        """

	echo "WT_1.bam,WT_2.bam,WT_3.bam" > wt.txt
        echo "U170K_1.bam,U170K_2.bam,U170K_3.bam" > u170k.txt


        rmats.py --b1 wt.txt --b2 u170k.txt --gtf ${tair10_gtf} -t paired --readLength 133 --variable-read-length --nthread 8 --novelSS --mil 20 --mel 15000 --od WT_vs_U170K --tmp temp
        """
}
