#================================================================================================================
# Application of the methods to RNA-Seq data from TCGA-BRCA 
#
# author: Fatih Kızılaslan (fatih.kizilaslan@medisin.uio.no)
# date: 19-January-2024
#================================================================================================================


# We follow the tutorial for by Zhao et al. (2023) for retrieving TCGA-BRCA data from  https://ocbe-uio.github.io/survomics/survomics.html
# Zhao, Z., Zobolas, J., Zucknick, M., Aittokallio, T. (2023). Tutorial on survival modelling with omics data. arXiv preprint arXiv:2302.12542.


library("TCGAbiolinks")
library("SummarizedExperiment")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("survival")
library("survminer")
library("M3C")
library("glmnet")

# Download TCGA breast cancer (BRCA) mRNA-Seq data using GDC api method
query <- TCGAbiolinks::GDCquery(project = "TCGA-BRCA",
                               data.category = "Transcriptome Profiling",
                               data.type = "Gene Expression Quantification",
                               workflow.type = "STAR - Counts",
                               experimental.strategy = "RNA-Seq",
                               sample.type = c("Primary Tumor"))
TCGAbiolinks::GDCdownload(query = query, method = "api")
dat <- TCGAbiolinks::GDCprepare(query = query)
SummarizedExperiment::assays(dat)$unstranded[1:5, 1:2]

?
?
?
?


# Duplicated genes determination: we use "tcga.brca.data" for finding and removing of duplicated genes from original "dat" file

dim(dat)  # gene expression matrices.
tcga.brca.data = rowData(dat) #includes "gene_name" without other variables! it is only 10 columns
head(tcga.brca.data) # ensembl id and gene id of the first 6 genes.
dim(tcga.brca.data)


rownames(tcga.brca.data) = NULL #remove row names
head(tcga.brca.data)


gene.ids.names = cbind(tcga.brca.data$gene_id, tcga.brca.data$gene_name)
dim(gene.ids.names)
save(gene.ids.names, file = "all.gene.ids.names.TCGA.BRCA.RData")

# Gene types for this data
as.data.frame(tcga.brca.data) %>% group_by(gene_type) %>% summarise(n=n())

# We will use "protein_coding" genes in our study, so we filter this subset of data as in original "dat" file
protein.coding.genes = which(tcga.brca.data$gene_type=="protein_coding") 
data.protein.coding = dat[protein.coding.genes,]
dim(data.protein.coding)    #  19962 x 1111

# Duplicated gene names in data.protein.coding, we have 24 duplicated genes as:
gene.names = as.data.frame(rowData(data.protein.coding )) %>% group_by(gene_name) %>% summarise(n=n()) %>% dplyr::filter(n >1)
gene.names

# gene_ids for these duplicated gene.names of data.protein.coding. 
AA = list()
for (i in 1:24) {
  aa = as.data.frame(rowData(data.protein.coding )) %>% dplyr::filter(gene_name == gene.names$gene_name[i]) 
  AA[[i]] = aa$gene_id
}
head(AA)

# Since the second values of these 24 genes are completely 0 ((except 5 of them), we remove these duplicated values
remove.gene.ids = unlist(AA)[seq(2,48,2)]         # for only 2nd values of each gene
remove.gene.ids 


# Removing duplicated gene ids, namely, "remove.gene.ids" from "data.protein.coding", then we call it "data.protein.coding.rev"
rm.col = c()
for (i in 1:24) {
  
  rm.col[i] = which(rownames(data.protein.coding) == remove.gene.ids[i] )
  
}
data.protein.coding.rev = data.protein.coding[-rm.col,]
dim(data.protein.coding)
dim(data.protein.coding.rev)               # we remove 24 duplicated genes: 19962-24 = 19938, so dimension 19938 x 1111

# Normalizing RNA count data by DESeq2 package
dds.rev <- DESeq2::DESeqDataSetFromMatrix(assays(data.protein.coding.rev)$unstranded, colData = meta.rev, design = ~ 1) 
dds2.rev <- DESeq2::estimateSizeFactors(dds.rev)
RNA_count.rev <- DESeq2::counts(dds2.rev, normalized=TRUE)
RNA_count.rev[1:5, 1:6]
dim(RNA_count.rev)        # 19938 x 1111


# Arrangements of clinical variables for the data:  "data.protein.coding.rev"

# Some gene values of this data 
SummarizedExperiment::assays(data.protein.coding.rev)$unstranded[1:30, 1:2]
# Column names of the data "data.protein.coding.rev"
colnames(colData(data.protein.coding.rev))


# We will use "submitter_id" variable for each subject for combining different the clinical data structures.
  
meta.rev = colData(data.protein.coding.rev)[, c("project_id", "submitter_id","barcode",
                                                "age_at_diagnosis", "age_at_index", "days_to_birth", "days_to_death",
                                                "days_to_last_follow_up", "year_of_birth","year_of_diagnosis","year_of_death", 
                                                "ethnicity", "gender", "vital_status", "paper_BRCA_Subtype_PAM50", "treatments", "race", 
                                                "ajcc_pathologic_m", "ajcc_pathologic_n", "ajcc_pathologic_stage","ajcc_pathologic_t",
                                                "icd_10_code","prior_malignancy", "synchronous_malignancy","prior_treatment",
                                                "initial_weight","paper_pathologic_stage","paper_BRCA_Pathology","morphology",
                                                'paper_DNA.Methylation Clusters','paper_Mutation Clusters','paper_Protein Clusters',
                                                'paper_miRNA Clusters')]
meta.rev$treatments = unlist(lapply(meta.rev$treatments, function(xx){any(xx$treatment_or_therapy == "yes")}))
meta.rev$pharmaT = pharmaT ##defining Pharmaceutical Treatment/Therapy indicator
meta.rev$radiaT = radiaT ##defining Radiation Treatment/Therapy indicator
dim(meta.rev) ## 1111 x  35    


# Survival time variable arrangement:
# Notice that 733th subject' time value is -Inf, and also we have some negative values and 0 values as a time variable for some subjects.
# We need to give a kind of threshold value for the time variable.
# Hence, we use 30 days as a threshold value for the survival time (namely we remove the sample if its time values is less than 30/365.25=0.08213552 or
# use time values are greater than 26/365.25 because after 26 days in the data is 30 days) 

meta.rev$time = apply(meta.rev[, c("days_to_death", "days_to_last_follow_up")], 1, max, na.rm = TRUE) / 365.25  
summary(meta.rev$time[-733])

meta.rev$status = meta.rev$vital_status
meta.rev$age = meta.rev$age_at_diagnosis / 365.25
dim(meta.rev)          # 1111 x 38

# We also remove the "male" subjects and time value of subjects is less than 30 days. ( hence, we remove 1111-1017=94 observations)
clinical.meta.rev = subset(meta.rev, gender == "female" & !duplicated(submitter_id) & time > 26/365.25 & !is.na(age)) ##in the original George use 0 for time threshold!
dim(clinical.meta.rev) # after that dimension is 1017 x 38


# Clinical data using BCR-Biotab
# In this part, we will get the clinical data for TCGA-BRC using BCR-Biotab, then we will combine both clinical data with respect to subjects sample ids.


                                                                 
# Variance filter:
# We apply variance filter to "RNA_count.rev" data, after that we decide how many genomic variables will be used in the analysis.



