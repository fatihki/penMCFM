# penMCFM: Weibull Mixture Cure Frailty Model with High-dimensional Covariates  

 In this work, we introduce a novel extension of the Weibull mixture cure model that incorporates a frailty component, employed to model an underlying latent population heterogeneity with respect to the outcome risk. Additionally, high-dimensional covariates are integrated into both the cure rate and survival part of the model, providing a comprehensive approach to employ the model in the context of high-dimensional omics data. 
 
 We perform variable selection via an **adaptive elastic-net penalization**, and propose a novel approach to inference using the **expectation–maximization (EM) algorithm**. 
 
We apply the novel approach to analyze **RNAseq gene expression data from bulk breast cancer patients included in The Cancer Genome Atlas (TCGA) database**. 

A set of prognostic biomarkers is then derived from selected genes, and subsequently validated via **functional enrichment analysis**. 

**A prognostic risk score index** based on the identified biomarkers is proposed and validated by exploring the patients' survival.
 


## Citation

Kızılaslan, Fatih, Swanson, David M., and Vitelli, Valeria. **"A Weibull Mixture Cure Frailty Model for High-dimensional Covariates"** 

## Data


## R

``functions.R`` includes the related functions in the algorithms.

``plots_sim.R`` draw the figures based on the simulation results. Results can be loaded as data framed from _df_bp_R1.RData_, _df_betap_R1.RData_ and _df.penCox.1se.train.Cstat.RData_.

## Supplementary Material

We present all the files in the Supplementary Material.

## TCGA-BRCA Application

**Results** includes all the selected genes and complete results of GO and KEGG enrichment analyses as in excel files.

**Figures** includes all the figures for the application part of the study.

**Data** includes genes of prior knowledge of BRCA from [CoxReg study of Li and Liu (2021)](https://github.com/zpliulab/CoxReg).

