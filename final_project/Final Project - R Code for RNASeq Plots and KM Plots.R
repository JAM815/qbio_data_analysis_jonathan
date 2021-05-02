if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
if(!requireNamespace("arsenal"))install.packages(c("arsenal"))
if(!requireNamespace("survival"))install.packages(c("survival"))
if(!requireNamespace("survminer"))install.packages(c("survminer"))
```

Always use library() to load your packages into the current R workspace / environment
```{r}
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(arsenal)
library(survival)
library(survminer)


Actual Analysis
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

install.packages("BiocManager", lib = "C:/Users/Jonathan's Laptop/Documents/R/win-library/4.0")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinksGUI.data")

install.packages("TCGAbiolinks", lib = "C:/Users/Jonathan's\\ Laptop/Documents/R/win-library/4.0")

library(devtools)

library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query) #only need this line of code once to download the data
sum_exp <- GDCprepare(query)
# Create a tutorial on SummarizedExperiment


patient_data <- colData(sum_exp)
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
patient_data$age_category = ifelse(patient_ages < 50, "Young", "Old")

htseq_counts <- assays(sum_exp)$"HTSeq - Counts"

patient_data$ESR1_counts = htseq_counts["ENSG00000091831",]

patient_data$ESR1_counts_log = sapply(htseq_counts["ENSG00000091831",], log10)

#genes <- c("ENSG00000141510", "ENSG00000141736", "ENSG00000121879")
#htseq_counts["ENSG00000141510",]

# Boxplots by age
boxplot(ESR1_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for ESR1 by Age Category")

plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$ESR1_counts)

ERPos_mask <- rowData(sum_exp)$external_gene_name == "TP53"