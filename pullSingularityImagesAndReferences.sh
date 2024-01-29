#!/bin/bash


#singularity images
singularity pull fastqc-0.11.9-hdfd78af_1.sif https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1 
singularity pull multiqc-1.13-pyhdfd78af_0.sif https://depot.galaxyproject.org/singularity/multiqc%3A1.13--pyhdfd78af_0
singularity pull trim-galore-0.6.7-hdfd78af_0.sif https://depot.galaxyproject.org/singularity/trim-galore%3A0.6.7--hdfd78af_0
singularity pull seqkit-2.3.1-h9ee0642_0.sif https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0
singularity pull hisat2-2.2.1-h87f3376_4.sif https://depot.galaxyproject.org/singularity/hisat2%3A2.2.1--h87f3376_4
singularity pull rmats-4.1.2-py39hd863d43_3.sif https://depot.galaxyproject.org/singularity/rmats%3A4.1.2--py39hd863d43_3
singularity pull wiggletools_1.2.8-hedeb25_0.sif https://depot.galaxyproject.org/singularity/wiggletools%3A1.2.8--hebded25_0
singularity pull bedgraphtobigwig-377--ha8a8165_3.sif https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig%3A377--ha8a8165_3
singularity pull pyfaidx-0.7.2.1-pyh7cba7a3_1.sif https://depot.galaxyproject.org/singularity/pyfaidx%3A0.7.2.1--pyh7cba7a3_1
singularity pull samtools-1.9-h91753b0_8.sif https://depot.galaxyproject.org/singularity/samtools%3A1.9--h91753b0_8
singularity pull deeptools-3.5.1--py0.sif https://depot.galaxyproject.org/singularity/deeptools%3A3.5.1--py_0
singularity pull deeptools-3.5.2--pyhdfd78af_1.sif https://depot.galaxyproject.org/singularity/deeptools%3A3.5.2--pyhdfd78af_1
singularity pull macs2-2.9.1--py39hf95cd2a_0.sif https://depot.galaxyproject.org/singularity/macs2%3A2.2.9.1--py39hf95cd2a_0
singularity pull bedgraphtobigwig-377--ha8a8165_3.sif https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig%3A377--ha8a8165_3
singularity pull wiggletools_1.2.8-hedeb25_0.sif https://depot.galaxyproject.org/singularity/wiggletools%3A1.2.8--hebded25_0
singularity pull bedops-2.4.41-h9f5acd7_0.sif https://depot.galaxyproject.org/singularity/bedops%3A2.4.41--h9f5acd7_0
singularity pull salmon-1.9.0-h7e5ed60_1.sif https://depot.galaxyproject.org/singularity/salmon%3A1.9.0--h7e5ed60_1
singularity pull homer-4.11--pl5321h9f5acd7_7.sif https://depot.galaxyproject.org/singularity/homer%3A4.11--pl5321h9f5acd7_7
singularity pull bedtools_2.31.1--hf5e1c6e_0.sif https://depot.galaxyproject.org/singularity/bedtools%3A2.31.1--hf5e1c6e_0
singularity pull bioconvert-1.1.1--pyhdfd78af_0.sif https://depot.galaxyproject.org/singularity/bioconvert%3A1.1.1--pyhdfd78af_0



#references
wget --no-check-certificate https://www.arabidopsis.org/download_files/Public_Data_Releases/TAIR_Data_20211231/Araport11_functional_descriptions_20211231.txt.gz
gunzip Araport11_functional_descriptions_20211231.txt.gz
wget --no-check-certificate https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_29122021.fa
wget --no-check-certificate https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
gunzip TAIR10_chr_all.fas.gz
wget --no-check-certificate https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
wget --no-check-certificate https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_TS_21Feb22_transfix.gtf

