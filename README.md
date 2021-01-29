# RNAseq_analysis
Script for analysis of RNAseq data obtained from Candida albicans
This script is used in the Kuchler lab (http://cdl.univie.ac.at/) at MFPL (https://www.mfpl.ac.at/de.html).

This repository conatains a pipeline for the primary analysis of Illumina short read sequencing RNAseq data (single-end) obtained from the fungal pathogen Candida albicans. It includes the retrieval of genomic data required for analysis from the Candida Genome Database (CGD, http://www.candidagenome.org/), quality control of the raw data, trimming and mapping of reads, removal of reads mapping to ribosomal RNAs and counting the number of reads mapping to genes. Only the RNAseq raw data in .bam or .fastq format (compressed or uncompressed) have to be provided by the user.

# Tools required for analysis:

samtools (http://www.htslib.org/)

bedtools (https://bedtools.readthedocs.io/en/latest/)

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

MultiQC (https://multiqc.info/)

cutadapt (https://cutadapt.readthedocs.io/en/stable/)

NextGenMap (https://github.com/Cibiv/NextGenMap/wiki)

DeepTools (https://deeptools.readthedocs.io/en/develop/)

HTSeq (https://htseq.readthedocs.io/en/release_0.11.1/#)

All the above-mentioned tools have to be included in yout PATH environment.

# Usage:

Clone the repository by typing "git clone https://github.com/tschemic/RNAseq_analysis.git" and copy the raw data into the RNAseq_analysis directory.

Change the adapter sequence for read trimming in the config_file.txt if necessary. By default it contains the Illumina TrueSeq adapter.
  For adapter sequences see: https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-11.pdf
  
Change into the required_files directory and run the analysis script (by typing: bash analysis_script.sh).

After the pipeline has finished change into the diff_expr_analysis directory and use the edgeR_analysis.R script as a basis for differential expression analysis in R.
