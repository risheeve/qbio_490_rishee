setwd("/Users/risheev/Documents/USC/qbio_490_rishee/analysis_data")


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

clinical <- read.csv("/Users/risheev/Documents/USC/qbio_490_rishee/analysis_data/brca_clinical_data.csv")

clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation")

#Q1/2. person_neoplasm_cancer_status - categorical
#na_mask <- is.na(clinic$person_neoplasm_cancer_status)
#sum(na_mask)
#na_mask <- !is.na(clinical$person_neoplasm_cancer_status)
#no_na_neoplasm_clinic <- clinical[na_mask, ]
#Q3/4. radiation dosage - discrete
#na_mask2 <- is.na(clinical_rad$radiation_dosage)
#sum(na_mask2)


clinic_rad_merge <- merge(clinic, clinical_rad, by = "bcr_patient_barcode")

na_mask3 <- !is.na(clinic_rad_merge$person_neoplasm_cancer_status)
clinic_rad_merge <- clinic_rad_merge[na_mask3, ]
na_mask2 <- !is.na(clinical_rad$radiation_dosage)
clinic_rad_merge <- clinic_rad_merge[na_mask2, ]


#Q5. Patients with tumors require a higher radiation dosage; Patients with tumors have lower survival rates of BRCA; Patients with higher radiation dosages have higher survival rates of cancer

clinic_rad_merge$person_neoplasm_cancer_status <- factor(clinic_rad_merge$person_neoplasm_cancer_status)

boxplot(clinic_rad_merge$radiation_dosage~clinic_rad_merge$person_neoplasm_cancer_status,
     main = "Neoplasm Cancer Status vs Radiation Dosage",
     xlab = "Neoplasm Cancer Status", 
     ylab = "Radiation Dosage")


neoplasm_na_mask <- ifelse(is.na(clinic_rad_merge$person_neoplasm_cancer_status), F, T)
neoplasm_cleaned_clinic <- clinic_rad_merge[neoplasm_na_mask, ]

neoplasm_cleaned_clinic$survival_time <- ifelse(is.na(neoplasm_cleaned_clinic$days_to_death),
                                                neoplasm_cleaned_clinic$survival_time <- neoplasm_cleaned_clinic$days_to_last_followup,
                                                neoplasm_cleaned_clinic$survival_time <- neoplasm_cleaned_clinic$days_to_death)

inf_mask <- ifelse(neoplasm_cleaned_clinic$survival_time == "-Inf", F, T)
neoplasm_cleaned_clinic <- neoplasm_cleaned_clinic[inf_mask, ]

na_mask_neoplasm_x <- !is.na(neoplasm_cleaned_clinic$survival_time)
neoplasm_cleaned_clinic <- neoplasm_cleaned_clinic[na_mask_neoplasm_x, ]

neoplasm_cleaned_clinic$death_event <- ifelse(neoplasm_cleaned_clinic$vital_status == "Alive", 
                                              neoplasm_cleaned_clinic$death_event <- FALSE,
                                              neoplasm_cleaned_clinic$death_event <- TRUE)

survival_object <- Surv(time = neoplasm_cleaned_clinic$survival_time, event = neoplasm_cleaned_clinic$death_event)
fit_object <- survfit(survival_object~neoplasm_cleaned_clinic$person_neoplasm_cancer_status, neoplasm_cleaned_clinic)
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
KM_plot <- survplot$plot + theme_bw() + theme(axis.title=element_text(size=20), axis.text=element_text(size=16),
                                          legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot

rad_na_mask <- ifelse(is.na(clinic_rad_merge$radiation_dosage), F, T)
rad_cleaned_clinic <- clinic_rad_merge[rad_na_mask, ]

low <- (as.numeric(as.character(rad_cleaned_clinic$radiation_dosage)))
low_na_mask <- !is.na(rad_cleaned)
rad_cleaned_clinic <- rad_cleaned_clinic[low_na_mask, ]

low_rad_mask <- ifelse(low <= 6000, T, F)
middle_rad_mask <- ifelse((low > 6000) & (low <= 12000), T, F)
high_rad_mask <- ifelse(rad_cleaned_clinic$radiation_dosage > 12000, T, F)
rad_statuss <- ifelse(low_rad_mask, "0-6000", ifelse(middle_rad_mask, "6000-12000", "12000+"))
rad_cleaned_clinic_merge <- merge(rad_cleaned_clinic, rad_statuss)
  
rad_cleaned_clinic$survival_time <- ifelse(is.na(rad_cleaned_clinic$days_to_death),
                                           rad_cleaned_clinic$survival_time <- rad_cleaned_clinic$days_to_last_followup,
                                           rad_cleaned_clinic$survival_time <- rad_cleaned_clinic$days_to_death)

inf_mask <- ifelse(rad_cleaned_clinic$survival_time == "-Inf", F, T)
rad_cleaned_clinic <- rad_cleaned_clinic[inf_mask, ]

rad_cleaned_clinic$death_event <- ifelse(rad_cleaned_clinic$vital_status == "Alive", rad_cleaned_clinic$death_event <- FALSE,
                                         rad_cleaned_clinic$death_event <- TRUE)

survival_object2 <- Surv(time = rad_cleaned_clinic$survival_time, event = rad_cleaned_clinic$death_event)
fit_object <- survfit(survival_object2 ~ rad_cleaned_clinic$rad_status, data = rad_cleaned_clinic)

