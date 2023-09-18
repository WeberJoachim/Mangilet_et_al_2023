#!/usr/bin/env nextflow

nextflow.enable.dsl=2




process bedgraphToBigwig {
	label "bedgraphToBigwig"
	publishDir 'results/bigwigs/', mode: 'copy'

    input:
        tuple val(name), path(bedgraph)
        path(chromsizes)

    output:
        tuple val (name), path("*.bw")

    script:
        """
	sort -k1,1 -k2,2n ${bedgraph} > ${name}_sorted.bdg
        bedGraphToBigWig ${name}_sorted.bdg ${chromsizes} ${name}.bw
	rm ${name}_sorted.bdg
        """
}



process multibigwigsummary {
    publishDir "results/deeptools/"

    input:
        tuple val(name1), path(bedgraph1), val(name2), path(bedgraph2), val(name3), path(bedgraph3), val(name4), path(bedgraph4), val(name5), path(bedgraph5), val(name6), path(bedgraph6), val(name7), path(bedgraph7), val(name8), path(bedgraph8), val(name9), path(bedgraph9)

    output:
        path("*.npz"), emit: npz

    script:
   	"""
    	multiBigwigSummary bins -b ${bedgraph1} ${bedgraph2} ${bedgraph3} ${bedgraph4} ${bedgraph5} ${bedgraph6} ${bedgraph7} ${bedgraph8} ${bedgraph9} -o wigcorrelate.npz
 
	"""
}


process plotCorrelation {
    publishDir "results/deeptools"

    input:
        path(npz)

    output:
        path("heatmap.svg")
        path("scatterplot.png")
        path("PCA_1_2.png")
        path("PCA_2_3.png")

        
    script:
    """
        plotCorrelation --corData ${npz} --whatToPlot heatmap --corMethod spearman --plotFileFormat svg --colorMap Dark2 -o heatmap.svg
        plotCorrelation --corData ${npz} --whatToPlot scatterplot --corMethod spearman --plotFileFormat png -o scatterplot.png
        plotPCA -in ${npz} -o PCA_1_2.png -T "PCA of PC1 and PC2"
        plotPCA -in ${npz} -o PCA_2_3.png -T "PCA of PC2 and PC3" --PCs 2 3

 
    """
}


process preprocess_bed{

    
    input:
        path(gff)

    output:
        path("TAIR10_genes.bed"), emit: bed_all

        
    script:
    """
        

        gff2bed < ${gff} > pre_TAIR10_genes.bed
	grep gene pre_TAIR10_genes.bed > TAIR10_genes.bed
        sed -i s/Chr1/1/g TAIR10_genes.bed
        sed -i s/Chr2/2/g TAIR10_genes.bed
        sed -i s/Chr3/3/g TAIR10_genes.bed
        sed -i s/Chr4/4/g TAIR10_genes.bed
        sed -i s/Chr5/5/g TAIR10_genes.bed
        sed -i s/ChrC/chloroplast/g TAIR10_genes.bed
        sed -i s/ChrM/mitochondria/g TAIR10_genes.bed
 
    """
}

process sub_bed_apa_repressed_composite{

    
    input:
        path(bed)
        path(gene_list)

    output:
        path("TAIR10_genes_apa_enhanced_same.bed"), emit: bed_apa_repressed_composite

        
    shell:
    '''
        for i in $(cat !{gene_list}); do grep ${i} !{bed} >> TAIR10_genes_apa_enhanced_same.bed; done

        
    '''


}




process computeMatrix {
    publishDir "results/deeptools"

    input:
        tuple val(name1), path(bw1)
        tuple val(name2), path(bw2)
        tuple val(name3), path(bw3)
        path(bed)

    output:
        tuple val("TAIR10_gene_all"), path("*.gz"), emit: matrix
        tuple val("TAIR10_gene_all"), path("*.bed"), emit: sorted_regions_bed
       
        
    script:
    """
    	computeMatrix scale-regions -S ${bw1} ${bw2} ${bw3} -R ${bed} -p 8 -m 1000 -a 1000 -b 1000 -bs 1 --skipZeros -out matrix_${bed.getSimpleName()}.gz --outFileSortedRegions matrix_${bed.getSimpleName()}_sortedregions.bed
 
    """
}


process plotProfile {
    publishDir "results/deeptools"

    input:
        tuple val(name), path(npz)
       

    output:
        path("*")
       
        
    script:
    """
    plotProfile -m ${npz} --perGroup --regionsLabel " " -y "PolII occupancy"  -o ${name}.pdf --yMin -0.5 
 
    """
}
