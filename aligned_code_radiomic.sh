
/data/cath/users/Yang/software/plink2 \
  --bfile /data/circ_ukb/AHua/MGB_array/GSA_53k/GSA_53K.hg38 \
  --keep /data/circ_ukb/AHua/MGB_array/GSA_53k/keep_IDs.txt \
  \
  --maf 0.01 \
  --mac 20 \
  --geno 0.1 \
  --hwe 1e-15 \
  --mind 0.1 \
  \
  --write-snplist \
  \
  --out snps_pass_meghana

# step 1
#!/bin/bash
#SBATCH --job-name=regenie_step1
#SBATCH --output=logs/regenie_%x_%j.out
#SBATCH --error=logs/regenie_%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal

# Load modules or activate environment
cd /data/circ_ukb/AHua
wdir=/data/circ_ukb/AHua/Meghana/project/GWAS/regenie/step1/refined
mkdir -p ${wdir}

ARRAY=GSA_53K
phenoCol=IDN
bed=/data/circ_ukb/AHua/MGB_array/GSA_53k/GSA_53K.hg38

# Activate conda environment
export PATH=$HOME/.conda/envs/regenie_env/bin:$PATH

# Run regenie
regenie \
    --step 1 \
    --bed ${bed} \
    --phenoFile /data/circ_ukb/AHua/rad_phenotype/split_features/${phenoCol}.tsv \
    --phenoCol ${phenoCol} \
    \
    --covarFile /data/circ_ukb/AHua/MGB_array/GSA_53k/GSA_53K.covar_numeric.tsv \
    --covarColList batch,gender,age_dnaextraction,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    \
    --out ${wdir}/${ARRAY}_${phenoCol} \
    \
    --bsize 1000 \
    --lowmem \
    --extract /data/circ_ukb/AHua/Meghana/snps_pass_meghana.snplist




#step 2
#!/bin/bash
#SBATCH --job-name=regenie
#SBATCH --partition=bigmem
#SBATCH --time=140:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=3
#SBATCH --array=1-22
#SBATCH --output=logs/regenie_%A_%a.out
#SBATCH --error=logs/regenie_%A_%a.err

export PATH=$HOME/.conda/envs/regenie_env/bin:$PATH

wdir=/data/circ_ukb/AHua/Meghana/project/GWAS/regenie/step2/refined
mkdir -p ${wdir}
ARRAY=GSA_53K
chr=${SLURM_ARRAY_TASK_ID:-22}


PHENOS=(GLRLM.run.length.non.uniformity GLCM.correlation)

for phenoCol in "${PHENOS[@]}"; do
  echo "Running REGENIE step 2 for ${phenoCol} on chr${chr}..."

  regenie \
    --step 2 \
    --bgen /data/circ_ukb/AHua/imputated_MGBB/bgen/GSA_53K.merged.chr${chr}.bgen \
    --sample /data/circ_ukb/AHua/imputated_MGBB/bgen/GSA_53K.merged.chr${chr}.sample \
    --phenoFile /data/circ_ukb/AHua/rad_phenotype/split_features/${phenoCol}.tsv  \
    --phenoCol ${phenoCol} \
    --bsize 200 --lowmem \
    --approx \
    --covarFile /data/circ_ukb/AHua/MGB_array/GSA_53k/GSA_53K.covar_numeric.tsv \
    --covarColList batch,gender,age_dnaextraction,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList "" \
    --pThresh 0.01 \
    --maxstep-null 2 \
    --pred /data/circ_ukb/AHua/Meghana/project/GWAS/regenie/step1/refined/${ARRAY}_${phenoCol}_pred.list \
    --out ${wdir}/${ARRAY}_${phenoCol}_chr${chr}

done

