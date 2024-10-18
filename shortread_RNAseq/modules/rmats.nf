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
	
	    mv ${gtf} ${gtf.getSimpleName()}_noChr.gtf
	
        """

}


process rmats_shoots_ago4_1 {

    label "rmats"
    publishDir "results/rmats/shoots_WT_vs_ago4-1"

    input:
        path(tair10_gtf)
        tuple val(name1), path(bam1, stageAs: "WT_shoots_1.bam"), val(name2), path(bam2, stageAs: "WT_shoots_2.bam"), val(name3), path(bam3, stageAs: "WT_shoots_3.bam")
        tuple val(name4), path(bam4, stageAs: "ago4-1_shoots_1.bam"), val(name5), path(bam5, stageAs: "ago4-1_shoots_2.bam"), val(name6), path(bam6, stageAs: "ago4-1_shoots_3.bam")
 
    output:
        path "*"
        
	     
    script:
        """
	    echo "WT_shoots_1.bam,WT_shoots_2.bam,WT_shoots_3.bam" > wt.txt
	    echo "ago4-1_shoots_1.bam,ago4-1_shoots_2.bam,ago4-1_shoots_3.bam" > ago4-1.txt	

        rmats.py --b1 wt.txt --b2 ago4-1.txt --gtf ${tair10_gtf} -t paired --readLength 150 --variable-read-length --nthread 8 --novelSS --mil 20 --mel 15000 --od shoots_WT_vs_ago4-1 --tmp temp
        """
}





process rmats_roots_ago4_1 {

    label "rmats"
    publishDir "results/rmats/roots_WT_vs_ago4-1"

    input:
        path(tair10_gtf)
        tuple val(name1), path(bam1, stageAs: "WT_roots_1.bam"), val(name2), path(bam2, stageAs: "WT_roots_2.bam"), val(name3), path(bam3, stageAs: "WT_roots_3.bam")
        tuple val(name4), path(bam4, stageAs: "ago4-1_roots_1.bam"), val(name5), path(bam5, stageAs: "ago4-1_roots_2.bam"), val(name6), path(bam6, stageAs: "ago4-1_roots_3.bam")

    output:
        path "*"


    script:
        """
            echo "WT_roots_1.bam,WT_roots_2.bam,WT_roots_3.bam" > wt.txt
            echo "ago4-1_roots_1.bam,ago4-1_roots_2.bam,ago4-1_roots_3.bam" > ago4-1.txt

        rmats.py --b1 wt.txt --b2 ago4-1.txt --gtf ${tair10_gtf} -t paired --readLength 150 --variable-read-length --nthread 8 --novelSS --mil 20 --mel 15000 --od roots_WT_vs_ago4-1 --tmp temp
        """
}

process rmats_shoots_ago1_27 {

    label "rmats"
    publishDir "results/rmats/shoots_WT_vs_ago1-27"

    input:
        path(tair10_gtf)
        tuple val(name1), path(bam1, stageAs: "WT_shoots_1.bam"), val(name2), path(bam2, stageAs: "WT_shoots_2.bam"), val(name3), path(bam3, stageAs: "WT_shoots_3.bam")
        tuple val(name4), path(bam4, stageAs: "ago1-27_shoots_1.bam"), val(name5), path(bam5, stageAs: "ago1-27_shoots_2.bam"), val(name6), path(bam6, stageAs: "ago1-27_shoots_3.bam")

    output:
        path "*"


    script:
        """
            echo "WT_shoots_1.bam,WT_shoots_2.bam,WT_shoots_3.bam" > wt.txt
            echo "ago1-27_shoots_1.bam,ago1-27_shoots_2.bam,ago1-27_shoots_3.bam" > ago1-27.txt

        rmats.py --b1 wt.txt --b2 ago1-27.txt --gtf ${tair10_gtf} -t paired --readLength 150 --variable-read-length --nthread 8 --novelSS --mil 20 --mel 15000 --od shoots_WT_vs_ago1-27 --tmp temp
        """
}





process rmats_roots_ago1_27 {

    label "rmats"
    publishDir "results/rmats/roots_WT_vs_ago1-27"

    input:
        path(tair10_gtf)
        tuple val(name1), path(bam1, stageAs: "WT_roots_1.bam"), val(name2), path(bam2, stageAs: "WT_roots_2.bam"), val(name3), path(bam3, stageAs: "WT_roots_3.bam")
        tuple val(name4), path(bam4, stageAs: "ago1-27_roots_1.bam"), val(name5), path(bam5, stageAs: "ago1-27_roots_2.bam"), val(name6), path(bam6, stageAs: "ago1-27_roots_3.bam")

    output:
        path "*"


    script:
        """
            echo "WT_shoots_1.bam,WT_shoots_2.bam,WT_shoots_3.bam" > wt.txt
            echo "ago1-27_roots_1.bam,ago1-27_roots_2.bam,ago1-27_roots_3.bam" > ago1-27.txt

        rmats.py --b1 wt.txt --b2 ago1-27.txt --gtf ${tair10_gtf} -t paired --readLength 150 --variable-read-length --nthread 8 --novelSS --mil 20 --mel 15000 --od roots_WT_vs_ago1-27 --tmp temp
        """

}
