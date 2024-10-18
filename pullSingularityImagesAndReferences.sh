#!/bin/bash


#singularity images
singularity pull fastqc-0.11.9-hdfd78af_1.sif https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1 
singularity pull multiqc-1.13-pyhdfd78af_0.sif https://depot.galaxyproject.org/singularity/multiqc%3A1.13--pyhdfd78af_0
singularity pull trim-galore-0.6.7-hdfd78af_0.sif https://depot.galaxyproject.org/singularity/trim-galore%3A0.6.7--hdfd78af_0
singularity pull seqkit-2.3.1-h9ee0642_0.sif https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0
singularity pull hisat2-2.2.1-h87f3376_4.sif https://depot.galaxyproject.org/singularity/hisat2%3A2.2.1--h87f3376_4
singularity pull rmats-4.1.2-py39hd863d43_3.sif https://depot.galaxyproject.org/singularity/rmats%3A4.1.2--py39hd863d43_3
singularity pull samtools-1.9-h91753b0_8.sif https://depot.galaxyproject.org/singularity/samtools%3A1.9--h91753b0_8


#references
wget --user-agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3" https://www.arabidopsis.org/download_files/Public_Data_Releases/TAIR_Data_20211231/Araport11_functional_descriptions_20211231.txt.gz
gunzip Araport11_functional_descriptions_20211231.txt.gz
wget --no-check-certificate https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_29122021.fa
wget --user-agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3" https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
gunzip TAIR10_chr_all.fas.gz
wget --no-check-certificate https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
wget --no-check-certificate https://ics.hutton.ac.uk/atRTD/RTD3/atRTD3_TS_21Feb22_transfix.gtf

