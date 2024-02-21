#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process homer_findMotifsGenome_custom_background {

    publishDir 'results/motif_finding/'
    label "homer_findMotifsGenome"

    input:
        tuple val(name_regions), path(regions), val(name_background), path(background)
        path(genome)
        val(motifsize)

    output:
        path('*', type: 'dir', emit: results_folder)

    script:
        """

        findMotifsGenome.pl ${regions} ${genome} ./${name_regions}_vs_${name_background} -size given -len ${motifsize} -bg ${background}

        """
}



process homer_findMotifsGenome_no_background {
    
    publishDir 'results/motif_finding/'
    label "homer_findMotifsGenome"

    input:
        tuple val(name_regions), path(regions)
        path(genome)
        val(motifsize)

    output:
        "*"

    script:
        """

        findMotifsGenome.pl ${regions} ${genome} ./${name}_vs_random_bg -len ${motifsize} -size given

        """
}



process homer_buildMotif_AAUAAA {

    label "homer_buildMotif"

    output:
        path("*AAUAAA.motif"), emit: motif
        path("*AAAAAA.motif"), emit: mask

    script:
        """

        seq2profile.pl AATAAA 1 AAUAAA > polyA_motif_AAUAAA.motif
        seq2profile.pl AAAAAA 0 AAAAAA > polyA_motif_AAAAAA.motif

        """

}


process homer_buildMotif_UUGUUU {

    label "homer_buildMotif"

    output:
        path("*UUGUUU.motif"), emit: motif
        path("*UUUUUU.motif"), emit: mask


    script:
        """

        seq2profile.pl TTGTTT 1 UUGUUU > polyA_motif_UUGUUU.motif
        seq2profile.pl TTTTTT 0 UUUUUU > polyA_motif_UUUUUU.motif
        """

}


process homer_buildMotif_UGUA {

    label "homer_buildMotif"

    output:
        path("*UGUA.motif"), emit: motif


    script:
        """

        seq2profile.pl TGTA 0 TGTA > polyA_motif_UGUA.motif
        
        """

}



process homer_countMotifs {


    label "homer_countMotifs"

    input:
        tuple val(name_regions), path(regions), val(name_background), path(background)
        path(genome)
        path(motif)
        path(mask)

    output:
        path("*"), emit: outputfolder

    script:
        """

        findMotifsGenome.pl ${regions} ${genome} ./${name_regions}_vs_${name_background} -norevopp -mis 0 -size given -len 7 -mknown ${motif} -maskMotif ${mask} -bg ${background}

        """

}


process homer_co_countMotifs {

    publishDir 'results/co_count_motifs/'
    label "homer_countMotifs"

    input:
        tuple val(name_regions), path(regions), val(name_background), path(background)
	path '*.motif'
        path(genome)


    output:
        path("co_occurence_*"), emit: outputfolder
	

    script:
        """
	
	cat *.motif >> co_motif_Urich_Arich.motif
	spacer=\$(head co_motif_Urich_Arich.motif -n 1 | awk '{ print \$2 }' | cut -c 7- | rev | cut -c 7- | rev | tr -d "\n"  | wc -c) 

        findMotifsGenome.pl ${regions} ${genome} ./co_occurence_\${spacer}_${name_regions}_${name_background} -norevopp -mknown co_motif_Urich_Arich.motif -bg ${background} -nomotif
        """

}



process homer_buildMotif_cooccurence{

    label "homer_buildMotif"


    input:
        tuple val(a_rich_motif), val(u_rich_motif)


    output:

	path("*_0N_*"),  emit: zero
        path("*_1N_*"),  emit: one
        path("*_2N_*"),  emit: two
        path("*_3N_*"),  emit: three
        path("*_4N_*"),  emit: four
        path("*_5N_*"),  emit: five
        path("*_6N_*"),  emit: six
        path("*_7N_*"),  emit: seven
        path("*_8N_*"),  emit: eight
        path("*_9N_*"),  emit: nine
        path("*_10N_*"), emit: ten
        path("*_11N_*"), emit: eleven
        path("*_12N_*"), emit: twelve
        path("*_13N_*"), emit: thirteen
        path("*_14N_*"), emit: fourteen
        path("*_15N_*"), emit: fifteen


    script:

        """
    
        stretch=""

        for ((i=0; i<=15;i++));do

            seq2profile.pl ${a_rich_motif}\${stretch}${u_rich_motif} 0 > coMotif_${a_rich_motif}_\${i}N_${u_rich_motif}.motif
            stretch+="N"

        done
        """
}



process homer_count_coMotifs {

    publishDir 'results/count_coMotifs/'
    label "homer_countMotifs"

    input:
        tuple val(name_regions), path(regions), val(name_background), path(background)
	    path(arich_motif)
        path(arich_mask)
        path(urich_motif)
        path(urich_mask)
        path(use_motif)
        path(genome)


    output:
        path("*countRegions_with_Motif_*.txt")
	

    script:
        """
		    
        annotatePeaks.pl ${regions} ${genome} -norevopp -m ${arich_motif} ${urich_motif} ${use_motif} -mask > ${name_regions}_countRegions_with_Motif.txt
        annotatePeaks.pl ${background} ${genome} -norevopp -m ${arich_motif} ${urich_motif} ${use_motif} -mask > ${name_background}_countRegions_with_Motif.txt

        """

}

