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
        path(genome)
        tuple val(name), path(motif), path(maskA), path(maskAT), path(maskT)


    output:
        path("co_occurence_*"), emit: outputfolder

    script:
        """

        findMotifsGenome.pl ${regions} ${genome} ./co_occurence_${name}_${name_regions}_${name_background} -norevopp -size given -mknown ${motif} -maskMotif ${maskA} ${maskT} ${maskAT} -bg ${background}

        """

}



process homer_buildMotif_cooccurence{

    label "homer_buildMotif"

    output:
        tuple val("1N"),  path("*_1N.motif"),  path("*_1N_mask_A.motif"),  path("*_1N_mask_AT.motif"),  path("*_1N_mask_T.motif"),  emit: one
        tuple val("2N"),  path("*_2N.motif"),  path("*_2N_mask_A.motif"),  path("*_2N_mask_AT.motif"),  path("*_2N_mask_T.motif"),  emit: two
        tuple val("3N"),  path("*_3N.motif"),  path("*_3N_mask_A.motif"),  path("*_3N_mask_AT.motif"),  path("*_3N_mask_T.motif"),  emit: three
        tuple val("4N"),  path("*_4N.motif"),  path("*_4N_mask_A.motif"),  path("*_4N_mask_AT.motif"),  path("*_4N_mask_T.motif"),  emit: four
        tuple val("5N"),  path("*_5N.motif"),  path("*_5N_mask_A.motif"),  path("*_5N_mask_AT.motif"),  path("*_5N_mask_T.motif"),  emit: five
        tuple val("6N"),  path("*_6N.motif"),  path("*_6N_mask_A.motif"),  path("*_6N_mask_AT.motif"),  path("*_6N_mask_T.motif"),  emit: six
        tuple val("7N"),  path("*_7N.motif"),  path("*_7N_mask_A.motif"),  path("*_7N_mask_AT.motif"),  path("*_7N_mask_T.motif"),  emit: seven
        tuple val("8N"),  path("*_8N.motif"),  path("*_8N_mask_A.motif"),  path("*_8N_mask_AT.motif"),  path("*_8N_mask_T.motif"),  emit: eight
        tuple val("9N"),  path("*_9N.motif"),  path("*_9N_mask_A.motif"),  path("*_9N_mask_AT.motif"),  path("*_9N_mask_T.motif"),  emit: nine
        tuple val("10N"), path("*_10N.motif"), path("*_10N_mask_A.motif"), path("*_10N_mask_AT.motif"), path("*_10N_mask_T.motif"), emit: ten
        tuple val("11N"), path("*_11N.motif"), path("*_11N_mask_A.motif"), path("*_11N_mask_AT.motif"), path("*_11N_mask_T.motif"), emit: eleven
        tuple val("12N"), path("*_12N.motif"), path("*_12N_mask_A.motif"), path("*_12N_mask_AT.motif"), path("*_12N_mask_T.motif"), emit: twelve
        tuple val("13N"), path("*_13N.motif"), path("*_13N_mask_A.motif"), path("*_13N_mask_AT.motif"), path("*_13N_mask_T.motif"), emit: thirteen
        tuple val("14N"), path("*_14N.motif"), path("*_14N_mask_A.motif"), path("*_14N_mask_AT.motif"), path("*_14N_mask_T.motif"), emit: fourteen
        tuple val("15N"), path("*_15N.motif"), path("*_15N_mask_A.motif"), path("*_15N_mask_AT.motif"), path("*_15N_mask_T.motif"), emit: fifteen


    script:

        """
        Amotif="AATAAA"
        Amask="AAAAAA"
        Umotif="TTGTTT"
        Umask="TTTTTT"
        stretch=""


        for ((i=1; i<=15;i++));do

            stretch+="N"
            seq2profile.pl \$Amotif\${stretch}\${Umotif} 2 > coMotif_\${i}N.motif
            seq2profile.pl \$Amask\${stretch}\${Umotif} 0 > coMotif_\${i}N_mask_A.motif
            seq2profile.pl \$Amotif\${stretch}\${Umask} 0 > coMotif_\${i}N_mask_T.motif
            seq2profile.pl \$Amask\${stretch}\${Umask} 0 > coMotif_\${i}N_mask_AT.motif

        done
        """
}

