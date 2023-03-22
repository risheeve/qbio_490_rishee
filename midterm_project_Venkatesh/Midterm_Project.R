#set working directory & make outputs folder to store plots
setwd("/PATH/TO/qbio_490_name/midterm_project_Venkatesh")
dir.create("outputs")

#loading libraries - Bioconductor, TCGA, maftools, survival, survminer, ggplot2, DESeq2, EnhancedVolcano
install.packages("BiocManager")
BiocManager::install(version = "3.16")

library(BiocManager)
library(TCGAbiolinks)
library(maftools)

if (!require(survival)){
  install.packages("survival")
}
library(survival)

if(!require(survminer)){
  install.packages("survminer")
}
library(survminer)

if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#clinical data query
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
GDCdownload(clin_query)
clinical <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")

clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug") #clinical - drug data query
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation") #clinical - radiation data query

#Mutation MAF data query
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)

maf <- GDCprepare(maf_query)

#if column is named bcr_patient_barcode and tumor_sample_barcode isn't found:
colnames(clinical)[1] = "Tumor_Sample_Barcode"
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

#RNASeq Data query
rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)
rna_clinical <- rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)
rna_counts <- as.data.frame(rna_se@assays@data$unstranded) 
rna_genes <- as.data.frame(rna_se@rowRanges) 


#Question: What is the correlation between age and mutation -  particular types of mutation/gene that gets mutated?

#making a merged data frame between clinical data and maf data
maf_rad_merge <- merge(clinical, maf_object@data, by = "Tumor_Sample_Barcode")

#removing na values from age and applying mask
rad_anatsite_na_mask <- !is.na(maf_rad_merge$age_at_initial_pathologic_diagnosis) 
maf_rad_merge <- maf_rad_merge[rad_anatsite_na_mask, ]

#creating a boxplot & saving file
setwd("outputs")
par(mar=c(1,1,1,1))
jpeg(file="BoxPlot.jpg")
boxplot(maf_rad_merge$age_at_initial_pathologic_diagnosis ~ maf_rad_merge$Variant_Classification,
        data = maf_rad_merge,
        xlab = "Mutation Variant Classification",
        ylab = "Age at Iniital Pathologic Diagnosis",
        main = "Age at Diagnosis vs Mutation Variant Classification")
dev.off()


#creating survival time category for Kaplan Meier
maf_rad_merge$survival_time <- ifelse(is.na(maf_rad_merge$days_to_death),
                                                maf_rad_merge$survival_time <- maf_rad_merge$days_to_last_followup,
                                                maf_rad_merge$survival_time <- maf_rad_merge$days_to_death)

#removing -infinity values from data frame& applying mask
inf_mask <- ifelse(maf_rad_merge$survival_time == "-Inf", F, T)
maf_rad_merge <- maf_rad_merge[inf_mask, ]

#creating death event category in data frame to categorize patients by vital status
maf_rad_merge$death_event <- ifelse(maf_rad_merge$vital_status == "Alive", 
                                              maf_rad_merge$death_event <- FALSE,
                                             maf_rad_merge$death_event <- TRUE)

#creating survival/fit objects for Kaplan Meier plot, then making the plot/saving file
survival_object <- Surv(time = maf_rad_merge$survival_time, event = maf_rad_merge$death_event)
fit_object <- survfit(survival_object~maf_rad_merge$Variant_Classification, maf_rad_merge)
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

jpeg(file = "KMPlot.jpg")
par(mar=c(1,1,1,1))
KM_plot <- survplot$plot + theme_bw() + theme(axis.title=element_text(size=20), axis.text=element_text(size=16),
                                              legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot
dev.off()

#sorting patients as Young or Old based on median age in dataframe; creating median variable, creating category as age category and applying mask to sort patients
median_age <- median(as.numeric(maf_object@clinical.data$age_at_initial_pathologic_diagnosis))
maf_object@clinical.data$age_category <- ifelse(as.numeric(maf_object@clinical.data$age_at_initial_pathologic_diagnosis) > median_age, 'Old', 'Young')

#creating boolean masks to subset data for young / old patients
young_mask <- ifelse(maf_object@clinical.data$age_category == 'Young', T, F)
young_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[young_mask]
young_maf <- subsetMaf(maf = maf_object,
                       tsb = young_patient_barcodes)

old_mask <- ifelse(maf_object@clinical.data$age_category == 'Old', T, F)
old_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[old_mask]
old_maf <- subsetMaf(maf = maf_object, 
                     tsb = old_patient_barcodes)

#creeating a cooncoplot/saving file
jpeg(file = "coOncoplot.jpg")
coOncoplot(m1 = young_maf,
           m2 = old_maf,
           m1Name = 'Young Patients',
           m2Name = 'Old Patients',
           borderCol = NA)
dev.off()

#prepping data - removing NAs, converting datatype to factor
rnaCountAge <- rna_se@assays@data$unstranded[, is.na(rna_se@colData$age_at_index)]
rnaCountAge <- as.data.frame(rnaCountAge)
rna_clinical$age_category <- factor(ifelse(rna_clinical$age_at_index <= 58, "young", "old"))
#pre-filtering - removing genes from rna counts where total sum across all patients less than 10
row_sums <- rowSums(rna_counts)

low_counts_mask <- ifelse(row_sums < 10, F, T)

rna_counts <- rna_counts[low_counts_mask,]

rna_genes <- rna_genes[low_counts_mask,]

#changing dataframe row / column names to match for gene id
rownames(rna_genes) <- rna_genes$gene_id
rownames(rna_counts) <- rownames(rna_genes)
colnames(rna_counts) <- rownames(rna_clinical)

#removing na's
age_na_mask <- !is.na(rna_clinical$age_category)
rna_clinical <- rna_clinical[age_na_mask, ]
rna_counts <- rna_counts[, age_na_mask]

#creating deseq data set dds
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~age_category)

#creating dds object and saving results to a data frame
dds_obj <- DESeq(dds)
resultsNames(dds_obj)
results <- results(dds_obj, format = "DataFrame", contrast = c("age_category", "old", "young"))

results <- data.frame(rna_genes$gene_name, rownames(results), results$log2FoldChange, 
                      results$pvalue, results$padj, -log10(results$padj))

colnames(results) <- c("gene_name", "gene_id", "log2foldchange", "pvalue", "padj", "-log10(padj)")

#analyzing upregulated and downregulated genes in individual dataframes
sig_results_mask <- ifelse(results$padj < 0.05, T, F)
sig_results <- results[sig_results_mask,]
ordered_up <- order(results$log2foldchange, decreasing = TRUE)
up_reg_results <- results[ordered_up,]
up_sig_results <- up_reg_results[up_reg_results$log2foldchange > 1, ]
order_down <- order(results$log2foldchange)
down_reg_results <- results[order_down,]
down_sig_results <- down_reg_results[down_reg_results$log2foldchange < -1,]

par(mar=c(1,1,1,1))

rownames(results) <- results$gene_id
#creating Enhanced Volcano Plot and saving file
jpeg(file = "EnhancedVolcanoPlot.jpg")
EnhancedVolcano(results, lab = rownames(results), x = 'log2foldchange', y = 'pvalue',
                title = 'Age Category Young vs Old',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 4.0)
dev.off()

