### Allele specific methylation pipeline for repro and non-repro bumblebee workers 
### Input files are .bam files from bismark alignment from trimmed raw data

### The following chunks of code can be copy and pasted to make a .PBS script

###------------------------------------------------------------------------------------------
###----------------------- FOR MERGED BAM FILES, ONE REPRO ONE NON-REPRO --------------------

# Using deduplicated and sorted .bam files 

module load samtools/1.3.2

samtools merge repro_deduplicated_sorted.bam j1r_deduplicated_sorted.bam j5r_deduplicated_sorted.bam \
j8r_deduplicated_sorted.bam

samtools merge nonrepro_deduplicated_sorted.bam j1nr_deduplicated_sorted.bam j5nr_deduplicated_sorted.bam \
j8nr_deduplicated_sorted.bam

samtools index repro_deduplicated_sorted.bam
samtools index nonrepro_deduplicated_sorted.bam

# Needed to take only the assembled scaffolds otherwise it didn't finish running for one file after 5 days 
# and many resources 
samtools view nonrepro_dedup_sorted.bam NC_015762.1 NC_015763.1 NC_015764.1 NC_015765.1 NC_015766.1 \
NC_015767.1 NC_015768.1 NC_015769.1 NC_015770.1 NC_015771.1 NC_015772.1 NC_015773.1 NC_015774.1 \
NC_015775.1 N_015776.1 NC_015777.1 NC_015778.1 NC_015779.1 -b > nonrepro_dedup_sorted_NConly.bam

samtools view repro_dedup_sorted.bam NC_015762.1 NC_015763.1 NC_015764.1 NC_015765.1 NC_015766.1 \
NC_015767.1 NC_015768.1 NC_015769.1 NC_015770.1 NC_015771.1 NC_015772.1 NC_015773.1 NC_015774.1 \
NC_015775.1 N_015776.1 NC_015777.1 NC_015778.1 NC_015779.1 -b > repro_dedup_sorted_NConly.bam

###------------------------------------------------------------------------------------------

### Convert .bam files from a bismark aligner so methpipe can use the output (.mr = mapped reads)

#!/bin/bash

#PBS -N mr_file_creation
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR

module load methpipe/3.4.2

for file in $(ls *.bam)
do 
    base=$(basename ${file} "_deduplicated_sorted.bam")
    to-mr -o ${base}.mr -m bismark ${file}
done

###------------------------------------------------------------------------------------------

# The rest of the script is standard, see allele_specific_meth_indiv_files.sh