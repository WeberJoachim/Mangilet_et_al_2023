#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process wiggletools_mean_wt {
    
    label "wiggletools_mean"
    input:
        tuple path(bw1), path(bw2), path(bw3)
	       

    output:
        path("*.bg")

        	     
    script:
        """
	
        wiggletools write_bg WT_mean.bg mean ${bw1} ${bw2} ${bw3} 

        """
}


process wiggletools_mean_u170k {
    
    label "wiggletools_mean"
   
    input:
        tuple path(bw1), path(bw2), path(bw3)
	       

    output:
        path("*.bg")

        	     
    script:
        """
	
        wiggletools write_bg U170K_mean.bg mean ${bw1} ${bw2} ${bw3} 

        """
}

process wiggletools_mean_u1c {

     label "wiggletools_mean"
    
     input:
	tuple path(bw1), path(bw2)	
	
     output:
	path("*.bg")


     script:
	"""
	wiggletools write_bg U1_C_mean_rep1_2.bg mean ${bw1} ${bw2}
	"""