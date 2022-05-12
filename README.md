# JoinASM

JoinASM is a scaffolding tool contigs using linked reads.

## Dependancies
- gcc-7.2.0 or higher

## Installation

Clone this repository and run make in the folder

    git clone https://forgemia.inra.fr/seqoccin/axis1-assembly-methodo.git
    cd joinASM
    make

## Creating input files

1st step: Align linked reads to contigs to create a bam file (example with LongRanger on genotoul clster)


    module load bioinfo/bwa-0.7.17
    module load bioinfo/longranger-2.2.2

    longranger mkref contigFileName.fasta
    longranger align --id=outputFolderName --fastqs=fastqPATH --sample=sampleNAME --reference=refdata-contigFileName --localcores=20 --localmem=100

The bam file would be in the outputFolderName/out folder with the name possorted_bam.bam

2nd step: Create molecule tsv file

Using the modified LongRanger version from this repository:
    module load bioinfo/bwa-0.7.17
    module load bioinfo/samtools-1.10

    samtools sort -t BX -@ 20 -m 2G -o bcsorted_bam.bam possorted_bam.bam
    longranger-2.2.2_reportMolecules/longranger reportMolecules --id=outputFolderName \
     --fastqs=fastqPATH \
     --sample=sampleNAME \
     --reference=refdata-contigFileName \
     --vcmode=disable \
     --bam_bc=bcsorted_bam.bam \
     --bam_pos=possorted_bam.bam  \
     --localcores=20 --localmem=100

The molecule file would be in the outputFolderName/out folder with the name
fragments_tsv.tsv

3rd step: Create contig size bed file

    module load bioinfo/samtools-1.10

    samtools faidx contigs.fasta
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' contigs.fasta.fai  > contigs.bed

4th step: Sort both molecules and contig size files with

    sort -k1,1 -s contigs.bed > contigs_size_sort_by_name.bed
    sort -k1,1 -s fragments_tsv.tsv > molecules_sort_by_name.tsv


## Usage

JoinASM requires two files containing the contigs sizes and the molecules information, both sorted by contig name. See **Creating input file** for how to obtain them.

    module load compiler/gcc-7.2.0

    joinASM -w 10000 -a 8 -c contigs.bed -g scaffoldingGraphFile.gfa  -s scaffoldsStructure.txt molecules_sort_by_name.tsv > logFile.txt

The output files are:
  - scaffoldingGraphFile.gfa contains the complete scaffolding graph
  - scaffoldsStructure.txt the order and the sens of contigs in each scaffold

To obtain the fasta file you can use scaffolds_to_fasta.cpp script from this repository:

./scaffolds_to_fasta.o contigs.fasta scaffoldsStructure.txt output_scaffolds.fasta


## SplitASM parameters

    Usage: ./joinASM <option(s)> MOLECULE FILE
    Options:
        -h, --help		Show this help message
        -w, --window INT	Window size for barcode consideration (default 10kb)
        -a, --arcsCondition INT	condition used for connecting two contigs; values{1..8}
        -c, --contigs FILE	Contig bed file name
        -g, --graph FILE	Output gfa file name
