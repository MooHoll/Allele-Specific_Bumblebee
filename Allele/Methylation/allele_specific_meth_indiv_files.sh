### Allele specific methylation pipeline for repro and non-repro bumblebee workers 
### Input files are .bam files from bismark alignment from trimmed raw data

### The following chunks of code can be copy and pasted to make a .PBS script

###------------------------------------------------------------------------------------------

### Convert .bam files from a bismark aligner so methpipe can use the output (.mr = mapped reads)

#!/bin/bash

#PBS -N mr_file_creation
#PBS -l walltime=60:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

cd $PBS_O_WORKDIR

module load methpipe/3.4.2

for file in $(ls *.bam)
do 
    base=$(basename ${file} "_deduplicated_sorted.bam")
    to-mr -o ${base}.mr -m bismark ${file}
done

###------------------------------------------------------------------------------------------

### Remove all NW_ reads from the .mr files otherwise the epistates command will take years to run 
### (even on loadssss of resources). This is because the bumblbee genome contains 5000+ scaffolds
### Attempted and it didn't finish 1 line in 8hrs with 2 cores and 16ppn.

for file in $(ls *.mr)
do
    base=$(basename ${file} ".mr")
    grep "NC_" ${file} > ${base}_NC_only.mr
done

###------------------------------------------------------------------------------------------

### Need to prepare a genome/chromosome folder for epireads command, remove additional text in headers as
### won't 'match' with chromosomes name in .mr file and make directory containing individual chromosomes

sed -i 's/B.*//' GCF_000214255.1_Bter_1.0_genomic.fa
sed -i 's/ //g' GCF_000214255.1_Bter_1.0_genomic.fa

module load perl/5.24.0
perl genome_split.pl GCF_000214255.1_Bter_1.0_genomic.fa
mkdir genome_chromosomes
mv NW* genome_chromosomes
mv NC* genome_chromosomes

### Perl script (taken from: http://seqanswers.com/forums/archive/index.php/t-32162.html):

#!/usr/bin/perl
 
$f = $ARGV[0]; #get the file name
 
open (INFILE, "<$f")
or die "Can't open: $f $!";
 
while (<INFILE>) {
$line = $_;
chomp $line;
if ($line =~ /\>/) { #if has fasta >
close OUTFILE;
$new_file = substr($line,1);
$new_file .= ".fa";
open (OUTFILE, ">$new_file")
or die "Can't open: $new_file $!";
}
print OUTFILE "$line\n";
}
close OUTFILE;


###------------------------------------------------------------------------------------------

### Making .epiread format from .mr files (more readable version for allele-specific analysis)

#!/bin/bash

#PBS -N epiread_creation
#PBS -l walltime=00:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR

module load methpipe/3.4.2

ref=/scratch/monoallelic/hm257/genome/genome_chromosomes

for file in $(ls *only.mr)
do 
    base=$(basename ${file} "_NC_only.mr")
    methstates -c ${ref} -o ${base}.epiread ${file}
done

###------------------------------------------------------------------------------------------

### Identify allelically methylated CpGs

### In the output file, each row represents a CpG pair made by any CpG and its previous CpG, 
### the first three columns indicate the positions of the CpG site, 
### the fourth column is the name including the number of reads covering the CpG pair, 
### the fifth column is the score for ASM, and the last four columns record the number of 
### reads of four different methylation combinations of the CpG pair: 
### methylated methylated (mm), methylated unmethylated (mu), unmethylated methylated (um), 
### or unmethylated unmethylated (uu). 


#!/bin/bash

#PBS -N allele_cpgs
#PBS -l walltime=00:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR

module load methpipe/3.4.2

ref=/scratch/monoallelic/hm257/genome/genome_chromosomes

for file in $(ls *.epiread)
do 
    base=$(basename ${file} ".epiread")
    allelicmeth -c ${ref} -o ${base}.allelicCpG ${file}
done


### Then pull out significant CpGs

for file in $(ls *r.allelicCpG)
do 
    base=$(basename ${file} ".allelicCpG")
    awk '$5<0.05 {print $0}'${file} > ${base}_significant.allelicCpG
done


###------------------------------------------------------------------------------------------

### Identify allelically methylated regions, using only 3 CpGs per window (default is 10), 
### minimum of 10 reads per CpG, amrfinder also corrects for multiple testing using FDR

#!/bin/bash

#PBS -N amr_finder
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

cd $PBS_O_WORKDIR

module load methpipe/3.4.2

ref=/scratch/monoallelic/hm257/genome/genome_chromosomes

for file in $(ls *.epiread)
do 
    base=$(basename ${file} ".epiread")
    amrfinder -m 10 -w 3 -o ${base}.amr -c ${ref} ${file}
done

### Repeated the above also with standard parameters (10CpG per window)




