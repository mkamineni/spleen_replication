# Spleen Replication Analysis

## Description

Replication and cross-cohort validation of genetic associations between spleen imaging-derived phenotypes and coronary artery disease (CAD) in Massachusetts General Brigham Biobank (MGBB), including evaluation of directional consistency across datasets.

## Overview

This repository contains code for replication analyses of spleen radiomic trait–associated genetic variants across independent cohorts. The goal is to evaluate robustness and consistency of genome-wide associations, identify shared loci, and assess directional concordance of effect estimates across studies.

Analyses include genome-wide summary statistic generation, SNP overlap comparisons, directional association testing, and visualization of replicated loci.

## Repository Structure

* **GWAS / cohort analysis**

  * `aligned_code_radiomic.sh`: scripts to run GWAS in MGBB cohort

* **Replication and SNP comparison**

  * `lead_snp_replication.py`: assesses replication of lead SNP associations across UKB and MGBB
  * `overlap_snp_replication_new.py`: assesses replication of independent and significant SNP associations across UKB and MGBB
  * `replicate_snps_directional.py`: assesses directional association of independent and significant SNP associations across UKB and MGBB

* **Genetic data processing**

  * `make_bed_for_liftover.py`: preparation of BED files for genomic liftover
  * `run_liftover.sh`: execution of liftover pipeline for genomic coordinate harmonization

* **Imaging feature extraction**

  * `get_spleen_radiomics_features_water_img.py`: extraction of spleen radiomic features from MRI-derived images with phase water in MGBB cohort
  * `run_spleen_ts_mri.py`: MRI preprocessing and feature generation pipeline for MGBB cohort

## Related Repositories

This repository is part of the *Spleen Radiomics and Coronary Artery Disease* project.

* spleenCADImaging
* spleenCADgenetics
* spleenResultsAnalysis

## Publication

If you use this repository, please cite:

Kamineni M, Raghu V, Truong B, Alaa A, Schuermans A, Friedman S, Reeder C, Bhattacharya R, Libby P, Ellinor PT, Maddah M, Philippakis A, Hornsby W, Yu Z, Natarajan P. *Deep learning-derived splenic radiomics, genomics, and coronary artery disease*. medRxiv. 2024.

## Author

Meghana Kamineni, MD
Harvard Medical School
Massachusetts General Hospital
