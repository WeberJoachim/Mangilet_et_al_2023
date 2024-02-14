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


    input:
        tuple val(a_rich_motif), val(u_rich_motif)


    output:

        tuple val("0N"),  path("*_0N_*"),  emit: zero
        tuple val("1N"),  path("*_1N_*"),  emit: one
        tuple val("2N"),  path("*_2N_*"),  emit: two
        tuple val("3N"),  path("*_3N_*"),  emit: three
        tuple val("4N"),  path("*_4N_*"),  emit: four
        tuple val("5N"),  path("*_5N_*"),  emit: five
        tuple val("6N"),  path("*_6N_*"),  emit: six
        tuple val("7N"),  path("*_7N_*"),  emit: seven
        tuple val("8N"),  path("*_8N_*"),  emit: eight
        tuple val("9N"),  path("*_9N_*"),  emit: nine
        tuple val("10N"), path("*_10N_*"), emit: ten
        tuple val("11N"), path("*_11N_*"), emit: eleven
        tuple val("12N"), path("*_12N_*"), emit: twelve
        tuple val("13N"), path("*_13N_*"), emit: thirteen
        tuple val("14N"), path("*_14N_*"), emit: fourteen
        tuple val("15N"), path("*_15N_*"), emit: fifteen


    script:

        """
    
        stretch=""

        for ((i=0; i<=15;i++));do

            stretch+="N"
            seq2profile.pl ${a_rich_motif}\${stretch}${u_rich_motif} 0 > coMotif_${a_rich_motif}_\${i}N_${u_rich_motif}.motif
        
        done
        """
}
