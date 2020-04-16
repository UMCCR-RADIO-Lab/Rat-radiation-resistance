#!/bin/bash
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --nodes=1
#SBATCH -J Cellranger
#SBATCH -p vccc
#SBATCH --mem=200G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=48:00:00
#SBATCH --mail-user=andrew.pattison@unimelb.edu.au
#SBATCH --mail-type=ALL

export PATH=/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/cellranger-3.1.0:$PATH

refdir=/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/reference

# Make the cellranger gtf
# cellranger mkgtf \
# $refdir/Rattus_norvegicus.Rnor_6.0.99.gtf \
# $refdir/Rattus_norvegicus.Rnor_6.0.99_cellranger_filtered.gtf \
# --attribute=gene_biotype:protein_coding

# cellranger mkref \
# --genome=Rnor_6 \
# --fasta=$refdir/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
# --genes=$refdir/Rattus_norvegicus.Rnor_6.0.99_cellranger_filtered.gtf

# cellranger count --id=KWR \
# --fastqs=/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/raw \
# --sample=LPRJ200257 \
# --transcriptome=/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/Rnor_6

cellranger count --id=KWP \
--fastqs=/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/raw \
--sample=LPRJ200258 \
--transcriptome=/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/Rnor_6




