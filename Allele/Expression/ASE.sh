# Data were trimmed and aligned using STAR, then sorted (pre-done for diff exp analysis)

# Following GATK best practices for SNP calling from RNA-seq
#https://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail

# ---------------------------------------------------------
# Add read group
# ---------------------------------------------------------

#!/bin/bash

#PBS -N read_group
#PBS -l walltime=01:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.3.2
module load java/1.8
module load picard/2.6.0

# Add read-group information to the alignments (needed for SNP calling)
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_Aligned.sortedByCoord.out.bam")
    java -jar /cm/shared/apps/picard/2.6.0/picard.jar \
    AddOrReplaceReadGroups \
    I=${file} \
    O=${base}_RG.bam \
    RGID=0001${base} \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=NA \
    RGSM=${base}
done



# ---------------------------------------------------------
# Mark duplicates 
# ---------------------------------------------------------

#!/bin/bash

#PBS -N mark_duplicates
#PBS -l walltime=01:40:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load picard/2.6.0

for file in $(ls *.bam)
do
  	base=$(basename ${file} "_RG.bam")
    java -jar /cm/shared/apps/picard/2.6.0/picard.jar MarkDuplicates \
        I=${file} \
        O=${base}_marked.bam \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT \
        M=${base}_metrics.txt
done


# ---------------------------------------------------------
# Split'n'trim to remove reads overlapping introns and keep only exonic data
# ---------------------------------------------------------

#!/bin/bash

#PBS -N split_n_trim
#PBS -l walltime=03:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load gatk/3.8

REF_FILE=/scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa    

for file in $(ls *.bam)
do
  	base=$(basename ${file} "_marked.bam")
    java -jar /cm/shared/apps/gatk/3.6/GenomeAnalysisTK.jar \
        -T SplitNCigarReads \
        -R ${REF_FILE} \
        -I ${file} \
        -o ${base}_split.bam \
        -rf ReassignOneMappingQuality \
        -RMQF 255 \
        -RMQT 60 \
        -U ALLOW_N_CIGAR_READS
done

# ---------------------------------------------------------
# SNP calling
# ---------------------------------------------------------

#!/bin/bash

#PBS -N snp_calling
#PBS -l walltime=17:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load gatk/3.8

REF_FILE=/scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa    

for file in $(ls *split.bam)
do
  	base=$(basename ${file} "_split.bam")
    java -jar /cm/shared/apps/gatk/3.6/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R ${REF_FILE} \
        -I ${file} \
        -dontUseSoftClippedBases \
        -stand_call_conf 20.0 \
        -o ${base}.vcf
done

# ---------------------------------------------------------
# SNP filtering
# ---------------------------------------------------------

#!/bin/bash

#PBS -N snp_filtering
#PBS -l walltime=00:10:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load gatk/3.8

REF_FILE=/scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa    

for file in $(ls *.vcf)
do
  	base=$(basename ${file} ".vcf")
    java -jar /cm/shared/apps/gatk/3.6/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R ${REF_FILE} \
        -V ${file} \
        -window 35 \
        -cluster 3 \
        -filterName FS \
        -filter "FS > 30.0" \
        -filterName QD \
        -filter "QD < 2.0" \
        -o ${base}_filtered.vcf
done

# ---------------------------------------------------------
# Select biallelic variants to determine alleles 
# ---------------------------------------------------------

# Tried the GATK selection programme but still kept 1/1 SNPs :S 
# Doing own filtering instead

for file in $(ls *filtered.vcf)
do
  	base=$(basename ${file} "_filtered.vcf")
    grep -v "1[/|]1" ${file} > ${base}_filtered_biallelic.vcf
done



# ---------------------------------------------------------
# Count reads to each allele
# ---------------------------------------------------------

#!/bin/bash

#PBS -N counting_reads
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load gatk/3.8

REF_FILE=/scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa    

# NOTE: this ignores 1/2 SNPs
for file in $(ls *.bam)
do
  	base=$(basename ${file} "_RG_split.bam")
    java -jar /cm/shared/apps/gatk/3.6/GenomeAnalysisTK.jar \
        -R ${REF_FILE} \
        -T ASEReadCounter \
        -o ${base}_allele_specific_counts_per_SNP.csv \
        -I ${file} \
        -sites ${base}_filtered_biallelic.vcf  \
        -U ALLOW_N_CIGAR_READS
done