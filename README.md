# penMCFM: Weibull Mixture Cure Frailty Model with High-dimensional Covariates  

In this work, we introduce a novel extension of the Weibull mixture cure model that incorporates a frailty component, employed to model an underlying latent population heterogeneity with respect to the outcome risk. Additionally, high-dimensional covariates are integrated into both the cure rate and survival part of the model, providing a comprehensive approach to employ the model in the context of high-dimensional omics data. 
 
We perform variable selection via an **adaptive elastic-net penalization**, and propose a novel approach to inference using the **expectation–maximization (EM) algorithm**. 
 
We apply the novel approach to analyze **RNAseq gene expression data from bulk breast cancer patients included in The Cancer Genome Atlas (TCGA) database**. 

A set of prognostic biomarkers is then derived from selected genes, and subsequently validated via **functional enrichment analysis**. 

**A prognostic risk score index** based on the identified biomarkers is proposed and validated by exploring the patients' survival.
 

## Citation

 > Kızılaslan, Fatih, Swanson, David M., and Vitelli, Valeria (2024). **"A Weibull Mixture Cure Frailty Model for High-dimensional Covariates"**. arXiv DOI: [10.48550/arXiv.2401.06575](https://arxiv.org/abs/2401.06575) 


## R

It includes all the source code necessary to run the introduced and utilized methods in Kızılaslan, Swanson and Vitelli (2024).

``application.R`` includes the code for applying the methods to the TCGA-BRCA RNA-seq data.

``EM_penMCFM.R`` includes EM algorithm codes based on the methods employed in the study.

``functions.R`` includes the related functions in the algorithms.

``GMIFS_MCM.R`` includes the related functions for the GMIFS algorithm for the MCM based on the Fu et al. (2022).

``GMIFS_penCox_penMCFM.R`` includes the related functions for the GMIFS and penCox algorithms for the penMCFM.

``plots_sim.R`` draw the figures based on the simulation results. Results can be loaded from data frames _df_bp_R1.RData_, _df_betap_R1.RData_ and _df.penCox.1se.train.Cstat.RData_ in [Simulation>Data](https://github.com/fatihki/penMCFM/tree/main/Simulation/Data).

``sim_penMCFM.R`` is a function designed for running simulations for the MCFM, utilizing various methods.


## Simulation

**Data** includes data frames obtained from simulation results, enabling the generation of plots

**Figures** includes all the plots from the simulation section of the study.


## Supplementary Material

We present all the files in the Supplementary Materials of the study.


## TCGA-BRCA Application

**Data** includes prior knowledge of BRCA genes from [CoxReg study of Li and Liu (2021)](https://github.com/zpliulab/CoxReg).

**Figures** includes all the plots relevant to the application section of the study.

**Results** contains all the selected genes and comprehensive results of Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analyses in excel files.


## References

 > Fu, H., Nicolet, D., Mrózek, K., Stone, R.M., Eisfeld, A.K., Byrd, J.C., Archer, K.J. (2022). Controlled variable selection in Weibull mixture cure models for high‐dimensional data. Statistics in Medicine, 41 (22), 4340-4366. DOI: [10.1002/sim.9513](https://doi.org/10.1002/sim.9513)

 > Li, L., Liu, ZP. (2021). Detecting prognostic biomarkers of breast cancer by regularized Cox proportional hazards models. Journal of Translational Medicine, 19, 514 (2021). DOI: [10.1186/s12967-021-03180-y](https://doi.org/10.1186/s12967-021-03180-y)

