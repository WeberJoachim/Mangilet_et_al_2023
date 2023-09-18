#!/usr/bin/env nextflow

nextflow.enable.dsl=2



process deeptools_multibamsummary {

    publishDir "results/deeptools"
    label "deeptools"

    input:
        path(sorted_bams)
	    path(bambais) 
               
    output:
        path("*.npz"), emit: numpy_array
	    path("*.log")

    script:
    """

    deeptools --version > deeptools.log
	multiBamSummary bins --smartLabels -p 8 --bamfiles ${sorted_bams} -o results.npz
  
    """
}



process deeptools_plots{

    publishDir "results/deeptools"
    label "deeptools"

    input:
        path(numpy_array)
 
               
    output:
        path("heatmap.svg")
	    path("scatterplot.png")

    script:
        """

        plotCorrelation --corData ${numpy_array} --whatToPlot heatmap --corMethod spearman --plotFileFormat svg --colorMap Dark2 -o heatmap.svg
	    plotCorrelation --corData ${numpy_array} --whatToPlot scatterplot --corMethod spearman --plotFileFormat png -o scatterplot.png
        plotPCA -in ${numpy_array} -o PCA_1_2.png -T "PCA of PC1 and PC2"
        plotPCA -in ${numpy_array} -o PCA_2_3.png -T "PCA of PC2 and PC3" --PCs 2 3
  
        """
}
