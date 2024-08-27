#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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

process homer_buildMotif_YA {

    label "homer_buildMotif"

    output:
        path("*YA.motif"), emit: motif


    script:
        """
        seq2profile.pl YA 0 YA > polyA_motif_YA.motif
        """

}


process homer_count_coMotifs {

    publishDir 'results/count_coMotifs/'
    label "homer_countMotifs"

    input:
        tuple val(name_regions), path(regions)
        path(arich_motif)
        path(urich_motif)
        path(use_motif)
	    path(cleavage_motif)
        path(genome)


    output:
        path("*countRegions_with_Motif.txt"), emit: info
        val(name_regions), emit: name
        

    script:
        """   
        annotatePeaks.pl ${regions} ${genome} -m ${arich_motif} ${urich_motif} ${use_motif} ${cleavage_motif} -mask > ${name_regions}_countRegions_with_Motif.txt
        """

}


process py_extract_infos_motif {

    publishDir 'results/extracted_motif_info', saveAs: {filename -> "${filename.substring(0,filename.length()-4)}_${task.index}.bed"}
    
    input: 
        file x
        path(pyscript)

    output:
	    path("*.bed")
	

    script:
	"""
	name=\$(head ${x} -n 1 | cut -d " " -f 3 | cut -d . -f 1)
		
	python3 ${pyscript} ${x}
	mv *.bed \${name}.bed
	"""

}
