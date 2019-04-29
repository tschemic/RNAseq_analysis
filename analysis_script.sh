#!/bin/bash

# Preparation and setup of required files

WKDIR=$(grep "working directory" ./Config_file.txt | cut -d ":" -f 2)
FILES=$WKDIR/required_files

read -p 'Do you want to retrieve genomic data from the CGD? (yes or no): ' GENEDATA

# Ask for raw data file format

echo 'Are the data in bam or fastq format?'
read -p 'Specify file format (bam or fastq): ' FORMAT

read -p 'Do you want to do a quality control of the raw data (yes or no): ' QCRAW

read -p 'Are the libraries prepared in a strand-specific way? (yes or no): ' STRANDED

if [ $GENEDATA == 'yes' ]
then
	echo 'Retrieving genomic data from CGD.'

	wget http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_chromosomes.fasta.gz  ## include WKDIR/required_files folder in wget!!!
	wget http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff
	cat C_albicans_SC5314_A22_current_features.gff | egrep -v "Ca22chr[1-7R]B" > C_albicans_SC5314_A22_current_features_haploid.gff
	gunzip C_albicans_SC5314_A22_current_chromosomes.fasta.gz
	cat C_albicans_SC5314_A22_current_chromosomes.fasta | egrep ">Ca22chr[1-7RM][A_]" | sed 's/>//g' | sed 's/(/1	/g' | sed 's/ nucleotides)//g' | sed 's/ /	/g' > chrMA.bed
	bedtools getfasta -fi C_albicans_SC5314_A22_current_chromosomes.fasta -bed chrMA.bed | fold -w 60 | sed 's/:1-[0-9]*//g' > C_albicans_SC5314_A22_current_chromosomesAM.fasta
else
	echo 'No genomic data are retrieved.'
fi


GENOME=$WKDIR/required_files/C_albicans_SC5314_A22_current_chromosomesAM.fasta
FEATURES=$WKDIR/required_files/C_albicans_SC5314_A22_current_features_haploid.gff
ADAPT1=$(cat $WKDIR/required_files/Config_file.txt | grep Read1: | cut -d ":" -f 2)
rRNA=$WKDIR/required_files/Ca_A22chrAM_rRNAloci.bed
mkdir $WKDIR/QC
PICARD=$(cat $WKDIR/required_files/Config_file.txt | grep "picard file path:" | cut -d ":" -f 2)


if [ $FORMAT == 'bam' ]
then
	echo 'File format is bam.'
	for i in $WKDIR/*.bam
	do
		bamToFastq -i $i -fq $i.fq
	done
elif [ $FORMAT == 'fastq' ]
then
	echo 'File format is fastq.'
else
	echo 'Invalid file format! Options are "bam" or "fastq".'
	exit
fi

# QC of raw data

if [ $QCRAW == 'yes' ]
then
	echo 'Quality control of raw data:'
	for i in $WKDIR/*q.gz
	do
		fastqc -o $WKDIR/QC $i
	done
else
	echo 'No QC of raw data done.'
fi

# Adapter removal with cutadapt and mapping of all files with NGM

for i in $WKDIR/*q.gz
do
	SNAME=$(echo $i | sed 's:/.*/::g')
	cutadapt -j 5 -q 30 -a $ADAPT1 $i > $i.trimmed.fq.gz 2>$WKDIR/QC/$SNAME.cutadapt.report.txt   # removes Illumina TrueSeq adapters from reads (change -a for different adapters); -j specifies number of cores to use, remove if not sure
	rm $i

	ngm -q $i.trimmed.fq.gz -r $GENOME -o $i.trimmed.fq.bam -b -t 5  # add -p for paired-end data; -t 6 is optional - means 6 threads of the processor are used, if you don't know what to do, remove it; --topn 1 --strata causes ngm to write only uniquely mapping reads to the output
	rm $i.trimmed.fq.gz

	samtools sort $i.trimmed.fq.bam -o $i.trimmed.fq.bam.sort.bam   # sort .bam files using samtools
	rm $i.trimmed.fq.bam

	bedtools intersect -a $i.trimmed.fq.bam.sort.bam -b $rRNA -v > $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam  # removal of reads mapping to rRNA loci
	rm $i.trimmed.fq.bam.sort.bam

	#samtools rmdup -s $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam.rmdup.bam  # removal of duplicated reads

	# Labelling of duplicated reads and removal of optical duplicates
	java -jar $PICARD MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true I=$i.trimmed.fq.bam.sort.bam.rRNAfilt.bam O=$i.trimmed.fq.bam.sort.bam.rRNAfilt.bam.markdup.bam M=$WKDIR/QC/$SNAME.markdup.metrics.txt

	#echo $i >> $WKDIR/QC/flagstat_analysis.txt
	samtools flagstat $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam.markdup.bam >> $WKDIR/QC/$SNAME.markdup.flagstat_analysis.txt   # flagstat analysis

	fastqc -o $WKDIR/QC $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam.markdup.bam

done

multiqc -s -o $WKDIR/QC $WKDIR/QC

# Preparation of coverage files for visualization in IGV

mkdir $WKDIR/IGV_files

for i in $WKDIR/*.markdup.bam
do
	samtools index $i
	SNAME=$(echo $i | sed 's:/.*/::g')
	bamCoverage -b $i -o $WKDIR/IGV_files/$SNAME.bw -p 5 --normalizeUsing CPM
done



mkdir $WKDIR/count

for i in $WKDIR/*.markdup.bam
do
	htseq-count -f bam -s $STRANDED -t gene -i ID $i $FEATURES > $i.count.txt  # read count  for each gene with htseq-count
	mv $i.count.txt $WKDIR/count
done

for i in $WKDIR/count/*.count.txt
do
	head -n -5 $i > $i.crop.txt  # clear count files for flags
done
