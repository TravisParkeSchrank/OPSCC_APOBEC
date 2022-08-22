
library(survival);
## Needed For KAP MEIRE analysis stuff
library(plotrix); 
#Needed for Multi Hist
library(wesanderson)


#load("opscc_tcga_composite_clinical.RData")

#load clinical data
new_df <- read.table(file ="hnsc_firehose_clinical.txt", sep = "\t", header = TRUE)
new_df <- subset(new_df, new_df$hpv_status == "negative")
new_df$sample_id <- toupper(as.character(new_df$bcr_patient_barcode))
hpv_neg <- new_df
hpv_neg_hnscc <- subset(hpv_neg, !(hpv_neg$anatomic_neoplasm_subdivision %in%  c("base of tongue", "oropharynx", "tonsil")))


load("hnsc_tcga_composite_clinical.RData")
hnsc_tcga_composite_clinical$sample_id <- hnsc_tcga_composite_clinical$study_id
hnsc_tcga_composite_clinical <- hnsc_tcga_composite_clinical[, !(names(hnsc_tcga_composite_clinical) %in% c("sample"))]
hnsc_tcga_composite_clinical <- unique(hnsc_tcga_composite_clinical)

merged <- subset(hnsc_tcga_composite_clinical, hnsc_tcga_composite_clinical$study_id %in% hpv_neg_hnscc$sample_id)




load("APOBEC_hpv_neg_data.RData")

APOBEC_hpv_neg_data$study_id <- gsub(APOBEC_hpv_neg_data$sample, pattern = "[.]", replacement = "-")
merged <- merge(merged, APOBEC_hpv_neg_data, by = "study_id")
#merged$number_pack_years_smoked[is.na(merged$number_pack_years_smoked)] <- 0
#merged <- subset(merged, merged$number_pack_years_smoked <  30)
#merged$A3A_or_NFkB <- merged$lt_n05 | merged$A3A

library("survminer")


fit <- survfit(Surv(PFI.time, PFI) ~ A3A, data = merged)
CoxFit <- coxph(Surv(PFI.time, PFI) ~ A3A, data = merged);
print(summary(CoxFit));
my_surv_plot_PFI_A3A <- ggsurvplot(
  fit, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = c("A3B", "A3A"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

#my_surv_plot_PFI
pdf(file = "HPVneg_surv_plot_A3A_PFI.pdf", width = 3, height = 6)
my_surv_plot_PFI_A3A
dev.off()








































fit <- survfit(Surv(RFS.time, RFS) ~ A3A, data = merged)
CoxFit <- coxph(Surv(RFS.time, RFS) ~ A3A, data = merged);
print(summary(CoxFit));
my_surv_plot_rfs_A3A <- ggsurvplot(
  fit, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  #legend.labs = c("WT", "dNFKB"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

#my_surv_plot_rfs
### run from command line
#pdf(file = "nice_surv_plot_RFS.pdf", width = 3, height = 4)
#my_surv_plot_rfs
#dev.off()

fit <- survfit(Surv(DSS.time, DSS) ~ A3A, data = merged)
CoxFit <- coxph(Surv(DSS.time, DSS) ~ A3A, data = merged);
print(summary(CoxFit));
my_surv_plot_DSS_A3A <- ggsurvplot(
  fit, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
#  legend.labs = c("WT", "dNFKB"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)




fit <- survfit(Surv(OS.time.x, OS.x) ~ A3A, data = merged)
CoxFit <- coxph(Surv(OS.time.x, OS.x) ~ A3A, data = merged);
print(summary(CoxFit));
my_surv_plot_OS_A3A <- ggsurvplot(
  fit, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  #  legend.labs = c("WT", "dNFKB"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)


























traf3_mut <- subset(merged, merged$TRAF3 != "")
traf3_del <- subset(merged, merged$TRAF3_value < -1)
cyld_mut <- subset(merged, merged$CYLD != "")
cyld_del <- subset(merged, merged$CYLD_value < -1)

mut_tumors <- c(as.character(traf3_mut$sample_id), as.character(cyld_mut$sample_id),as.character(traf3_del$sample_id),as.character(traf3_del$sample_id))

merged$pway_mut <- merged$sample_id %in% mut_tumors

merged <- subset(merged, merged$TP53 %in% c(NA))

fit_muts <- survfit(Surv(OS.time.x, OS.x) ~ pway_mut, data = merged)
CoxFit_muts <- coxph(Surv(OS.time.x, OS.x) ~ pway_mut, data = merged);
print(summary(CoxFit_muts));
#print(summary(fit_muts))
my_surv_plot_OS_muts <- ggsurvplot(
  fit_muts, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("WT", "dNFKB"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
#my_surv_plot_OS_muts
### run from command line
#pdf(file = "nice_surv_plot_RFS_muts.pdf", width = 3, height = 4)
#my_surv_plot_rfs_muts
#dev.off()

fit_muts <- survfit(Surv(PFI.time, PFI) ~ pway_mut, data = merged)
CoxFit_muts <- coxph(Surv(PFI.time, PFI) ~ pway_mut, data = merged);
print(summary(CoxFit_muts));
#print(summary(fit_muts))
my_surv_plot_pfi_muts <- ggsurvplot(
  fit_muts, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("WT", "dNFKB"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
#my_surv_plot_pfi_muts
### run from command line
#pdf(file = "nice_surv_plot_RFS_muts.pdf", width = 3, height = 4)
#my_surv_plot_rfs_muts
#dev.off()


fit_muts <- survfit(Surv(RFS.time, RFS) ~ pway_mut, data = merged)
CoxFit_muts <- coxph(Surv(RFS.time, RFS) ~ pway_mut, data = merged);
print(summary(CoxFit_muts));
#print(summary(fit_muts))
my_surv_plot_rfs_muts <- ggsurvplot(
  fit_muts, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("WT", "dNFKB"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
#my_surv_plot_rfs_muts
### run from command line
#pdf(file = "nice_surv_plot_RFS_muts.pdf", width = 3, height = 4)
#my_surv_plot_rfs_muts
#dev.off()


fit_muts <- survfit(Surv(DSS.time, DSS) ~ pway_mut, data = merged)
CoxFit_muts <- coxph(Surv(DSS.time, DSS) ~ pway_mut, data = merged);
print(summary(CoxFit_muts));
#print(summary(fit_muts))
my_surv_plot_DSS_muts <- ggsurvplot(
  fit_muts, 
  data = merged, 
  size = 1,                 # change line size
  palette = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]),
  conf.int = FALSE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("WT", "dNFKB"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
#my_surv_plot_DSS_muts
### run from command line
#pdf(file = "nice_surv_plot_RFS_muts.pdf", width = 3, height = 4)
#my_surv_plot_rfs_muts
#dev.off()



clincical_hpv_opscc_alts_preds_merged <- merged
save(clincical_hpv_opscc_alts_preds_merged, file = "clincical_hpv_opscc_alts_preds_merged.RData")
