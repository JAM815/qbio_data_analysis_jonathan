#QBio - Ellison Data Analysis Group - Project 1#

#Before running the code, all of the following need to be downloaded.
#If each of the following are not downloaded, parts of the the code will NOT run

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("devtools")) BiocManager::install(c("devtools"))
if(!requireNamespace("robustbase"))BiocManager::install(c("robustbase"))
library(devtools)
library(robustbase)
if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
if(!requireNamespace("survival"))install.packages(c("survival"))

#Below are the required libraries that need to be loaded in for the Code to function
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(survival)

# In order to generate two of the main plots needed for this analysis, genome
# mutations data is needed. And so the mutect2 pipline needs to accessed to 
# recieve all of the said mutation data. By using the "query" function, this
# information is downloaded and accessible for us to see and analyze. 

mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="mutect2")
maf_dataframe = read.maf(mutation)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query) # Note: this line of code only needs to be ran ONCE
sum_exp <- GDCprepare(query)

# From here, we need to access and find specific information regarding age.
# More specifically, we want to classify patients as either "old", "young" or "mid."

# In order to start this, we need to access the column containing the age at which
# patients first got sick. We can do this by accessing the Column Data of a Summarized
# Experiment and assigning to a variable. And using this variable, we can specifically
# target a column's data and use that to create another column that tells us whether a
# patient is "young", "mid" or "old."

patient_data <- colData(sum_exp) #access the dataframe contain all of the column info
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
# Here we create a New Column in patient_data called "age_category"
# Note: This column NOT be added to colData(sum_exp). 
# Instead it will only be added to the patient_data data table.

patient_data$age_category = ifelse(patient_ages < 40, "Young", ifelse(patient_ages >= 60, "Old", "Mid"))
# The ifelse() form is: ifelse( condition, action when condition is true, action when condition is false ). 
# In this case, Here we have two ifelse() embedded together to create the conditions we want
# for the information we want in this new column created.

patient_data_ <- patient_data[, c("barcode", "age_category")]

patient_data_

# From here, we need to obtain the shortened patient barcodes so that we can compare 
# it with all of the other patients.
short_maf <- substr(maf_dataframe@clinical.data$Tumor_Sample_Barcode, 1,12)

# From here, we want to create a new column in maf_dataframe that is designated for the
# short barcodes.
maf_dataframe@clinical.data$short_barcodes <- short_maf

# From here, we will create a new variable that extracts age category information for 
# each barcode that we have. In essence, we will have each patient's age along with their
# barcode. 

maf_ages <- patient_data[short_maf, "age_category"]

maf_dataframe@clinical.data$Ages <- maf_ages

# From here, we extract barcodes for each age group so that all of the "Young" barcodes
# are with the the "Young" patients and so on for each of the age groups.

young_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Young",]

mid_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Mid",]

old_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Old",]

NA_codes <- maf_dataframe@clinical.data[is.na(maf_dataframe@clinical.data$Ages),]
length(NA_codes$Ages)

# Following this, we create maf subsets for each age group which will be used to help
# generate the first graph needed for this analysis. 

young_maf <- subsetMaf(maf_dataframe, tsb = young_codes$Tumor_Sample_Barcode)

mid_maf <- subsetMaf(maf_dataframe, tsb = mid_codes$Tumor_Sample_Barcode)

old_maf <- subsetMaf(maf_dataframe, tsb = old_codes$Tumor_Sample_Barcode)


# Using all of this information, we can use the oncoplot command to obtain one of the 
# plots of the experiment, oncoplots.

oncoplot(maf = maf_dataframe, draw_titv = TRUE, top = 3)
oncoplot(maf = young_maf, draw_titv = TRUE, top = 3)
oncoplot(maf = mid_maf, draw_titv = TRUE, top = 3)
oncoplot(maf = old_maf, draw_titv = TRUE, top = 3)

# From here, we wanted to create a plot that could tell us the counts of the top three genes that 
# mutated the most in each age group. So then boxplots for top 3 genes for each age group
# were needed to find these results. 
# First we needed to access the counts of all genes and then afterwards we can seperate it
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"

# From the previously generated Oncoplots, we have four Genes of interest: 
# TP53, PIK3CA, TTN, GATA3. Notice that we have four genes instead of three. This because
# "Young" patients have the gene "GATA3" as their third top gene while "old" and "mid" 
# have the gene "TTN" as their third top gene. 

# To access each of these genes, we create a mask that allows us to target specific names
# and tell us whether a patinet had that mutation or did not.

TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TP53_mask]
# 
PIK3CA_mask <- rowData(sum_exp)$external_gene_name == "PIK3CA"
PIK3CA_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[PIK3CA_mask]
# 
TTN_mask <- rowData(sum_exp)$external_gene_name == "TTN"
TTN_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TTN_mask]
# 
GATA3_mask <- rowData(sum_exp)$external_gene_name == "GATA3"
GATA3_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[GATA3_mask]
# Note: external_gene_name does NOT need quotation marks because there are no spaces in 
# the name of the column. 
 
# From here, we extract the data that we want and get the HTseq counts for desired genes
# Note: We decided to apply a log function to each of these so that the graph looks nicer
# visually and so that the outliers are not so obviously present. 

patient_data$TP53_counts_log = sapply(htseq_counts[TP53_ENSG_name,], log10)
patient_data$PIK3CA_counts_log = sapply(htseq_counts[PIK3CA_ENSG_name,], log10)
patient_data$TTN_counts_log = sapply(htseq_counts[TTN_ENSG_name,], log10)
patient_data$GATA3_counts_log = sapply(htseq_counts[GATA3_ENSG_name,], log10)
 
 
# From here, we generate each of the Boxplots by age
jpeg("TP53_log_counts_by_age.jpg")
boxplot(TP53_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for TP53 by Age Category")
dev.off()

jpeg("PIK3CA_log_counts_by_age.jpg")
boxplot(PIK3CA_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for PIK3CA by Age Category")
dev.off()

jpeg("TTN_log_counts_by_age.jpg")
boxplot(TTN_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for TTN by Age Category")
dev.off()

jpeg("GATA3_log_counts_by_age.jpg")
boxplot(GATA3_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for GATA3 by Age Category")
dev.off()
# If running the data locally, comment the "dev.off()" and "jpeg("gene_log_counts_by_age_jpg)" lines


clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", barcode=barcodes_clinic)
GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up" #fixes an error in the name of the column

age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))


subtypes <- TCGAquery_subtype(tumor = "BRCA")
age_subs = subtypes$age_at_initial_pathologic_diagnosis
#also add the age category to the subtypes information
subtypes$age_category = ifelse(age_subs < 40, "Young", ifelse(age_subs >= 60, "Old", "Mid"))


#modify to your desired file path and name to save the figure. 
#On the cluster use ls and rsync to find the file and copy to your local

#TCGAanalyze_survival(clinic, "ethnicity")
#TCGAanalyze_survival(clinic, "menopause_status")
#TCGAanalyze_survival(clinic, "breast_carcinoma_progesterone_receptor_status")

TCGAanalyze_survival(clinic, "ethnicity", filename="./survival_curves/survival_ethnicity.pdf")
TCGAanalyze_survival(clinic, "menopause_status", filename="./survival_curves/survival_menopause.pdf")
TCGAanalyze_survival(clinic, "breast_carcinoma_progesterone_receptor_status",   filename="./survival_curves/survival_progesterone.pdf")
```