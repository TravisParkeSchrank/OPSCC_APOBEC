# $sudo mkdir /mnt/ramdisk
# $sudo mount -t tmpfs -o rw,size=64G tmpfs /mnt/ramdisk
# $cp ./ mc3.v0.2.8.PUBLIC.maf /mnt/ramdisk/mc3.v0.2.8.PUBLIC.maf



library(maftools)
library(stringr)
library(PerformanceAnalytics)

##### truncate a sample name to make it in normal TCGA format
truncate_new <- function(x){
  return(str_sub(x, start = 1L, end = 12L))
}

tail_new <- function(x){
  return(str_sub(x, start = 14L, end = 15L))
}

mut_data <- read.maf(maf = "TCGA.HNSC.varscan.5296cf00-4d8c-4db3-80d7-930a4b44f90d.DR-10.0.somatic.maf")

my_data <- data.frame("hugo" = mut_data@data$Hugo_Symbol, "study_id" = mut_data@data$Tumor_Sample_Barcode, "context" = mut_data@data$CONTEXT)
my_data$REF <- mut_data@data$Reference_Allele
my_data$ALT <- mut_data@data$Allele
my_data$VAF <- mut_data@data$t_alt_count / mut_data@data$t_depth
my_data$depth <- mut_data@data$t_depth
my_data$tri_nuc_context <- str_sub(my_data$context, start = 5L, end = 7L)
my_data$study_id <- str_replace_all(my_data$study_id, "[-]", ".")
my_data$study_id <- str_sub(my_data$study_id, start = 1L, end = 15L)

exp_data <- read.table(file = "apobec_transcript_data.txt", sep = "\t", header = TRUE)
row.names(exp_data) <- exp_data$sample
exp_data <- exp_data[, !(names(exp_data) %in% c("sample"))]

#get all ids with viariants 
tumors_with_variants <- names(table(my_data$study_id))

exp_data <- exp_data[, names(exp_data) %in% tumors_with_variants]
var_data <- subset(my_data, my_data$study_id %in% names(exp_data))
var_data$minus2pos <- str_sub(var_data$context, start = 4L, end = 4L)
exp_data <- data.frame(t(exp_data))
names(exp_data) <- c("cont1", "cont2", "cont3", "cont4", "A3A202", "A3A201", "A3A205", "A3B201", "A3B204", "A3B203")

exp_data$all_A3A <- log((2^exp_data$A3A201 + 2^exp_data$A3A202 + 2^exp_data$A3A205), base =2)
exp_data$all_A3B <- log((2^exp_data$A3B201 + 2^exp_data$A3B204 + 2^exp_data$A3B203), base =2)


apobec_expr_data <- exp_data[ , !(names(exp_data) %in% c("cont1", "cont2", "cont3", "cont4"))]
apobec_expr_data <- apobec_expr_data[, c("A3A202", "A3A201", "A3A205","all_A3A", "A3B201", "A3B204", "A3B203", "all_A3B")]

pdf(file = "performance_paris.pdf", height = 12, width = 12)
print(chart.Correlation(apobec_expr_data, histogram = TRUE, method = "pearson" ))
dev.off()

#### GET CLASSIC APOBEC VARS
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G", "T"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(apobec_expr_data, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0


#### GE Classic  APOBEC VARS C>GT with n2pos A
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G", "T"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("A"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_m2p_A")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0

#### GET Classic APOBEC VARS C>GT with n2pos T
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G", "T"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("T"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_m2p_T")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0

#### GET Classic APOBEC VARS C>G with n2pos G
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G", "T"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("G"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_m2p_G")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0


#### GET Classic APOBEC VARS C>G with n2pos C
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G", "T"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("C"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_m2p_C")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0
































#### GET CLASSIC APOBEC VARS C>G
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_CG")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0

#### GET  APOBEC VARS C>G with n2pos A
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("A"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_CG_m2p_A")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0

#### GET  APOBEC VARS C>G with n2pos T
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("T"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_CG_m2p_T")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0

#### GET  APOBEC VARS C>G with n2pos G
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("G"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_CG_m2p_G")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0


#### GET  APOBEC VARS C>G with n2pos C
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("G"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT_TCA_TCT, C_REF_ALT_GT_TCA_TCT$minus2pos %in% c("C"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_CG_m2p_C")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0

#### GET CLASSIC APOBEC VARS C>T
C_REF <- subset(var_data, var_data$REF %in% c("C"))
C_REF_ALT_GT <- subset(C_REF, C_REF$ALT %in% c("T"))
C_REF_ALT_GT_TCA_TCT <- subset(C_REF_ALT_GT, C_REF_ALT_GT$tri_nuc_context %in% c("TCA"))
n_apobec_var <- data.frame(table(C_REF_ALT_GT_TCA_TCT$study_id))
names(n_apobec_var) <- c("study_id", "n_classic_apobec_CT")

apobec_expr_data$study_id <- row.names(apobec_expr_data)
expr_with_var_freq <- merge(expr_with_var_freq, n_apobec_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0

### Total Vars
n_var <- data.frame(table(var_data$study_id))
names(n_var) <- c("study_id", "n_vars")
expr_with_var_freq <- merge(expr_with_var_freq, n_var, by = "study_id", all = TRUE)
expr_with_var_freq[is.na(expr_with_var_freq)] <- 0





plot_corr_data_A3A <- expr_with_var_freq[, names(expr_with_var_freq) %in% c("n_classic_apobec_CG", "n_classic_apobec_CG_m2p_A", "n_classic_apobec_CG_m2p_T", "n_classic_apobec_CG_m2p_G", "n_classic_apobec_CG_m2p_C", "n_classic_apobec_CT")] 
plot_corr_data_A3A <- log(plot_corr_data_A3A+1, base = 2)
plot_corr_data_A3A$all_A3A <- expr_with_var_freq$all_A3A
plot_corr_data_A3A$all_A3B <- expr_with_var_freq$all_A3B

pdf(file = "performance_paris_TCA_CG_A3A_minus2.pdf", height = 12, width = 12)
print(chart.Correlation(plot_corr_data_A3A, histogram = TRUE, method = "spearman" ))
dev.off()




plot_corr_data_A3A <- expr_with_var_freq[, names(expr_with_var_freq) %in% c("n_classic_apobec", "n_classic_apobec_m2p_A", "n_classic_apobec_m2p_T", "n_classic_apobec_m2p_G", "n_classic_apobec_m2p_C")]
plot_corr_data_A3A <- log(plot_corr_data_A3A+1, base = 2)
plot_corr_data_A3A$all_A3A <- expr_with_var_freq$all_A3A
plot_corr_data_A3A$all_A3B <- expr_with_var_freq$all_A3B

pdf(file = "performance_paris_classic_TCA_A3A_minus2.pdf", height = 12, width = 12)
print(chart.Correlation(plot_corr_data_A3A, histogram = TRUE, method = "spearman" ))
dev.off()

names(plot_corr_data_A3A) <- c("tCa_TG", "atCa_TG","ttCa_TG","gtCa_TG","ctCa_TG", "Exp.A3A", "Exp.A3B")
plot_corr_data_A3A <- plot_corr_data_A3A[, c("Exp.A3A", "Exp.A3B", "tCa_TG", "atCa_TG","ttCa_TG","gtCa_TG","ctCa_TG")]
library(psych)
pdf(file = "psych_paris_TCA_CG_A3A_minus2.pdf", height = 12, width = 12)
pairs.panels(plot_corr_data_A3A,
             smooth = FALSE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = TRUE,    # If TRUE, draws ellipses
             method = "spearman", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = TRUE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4,       # Histograms color
             stars = FALSE,       # If TRUE, adds significance level with stars
             ci = TRUE)   
dev.off()



library(corrplot)
library(wesanderson)

load("clin_fact_merged_HPVpos_inclusive_TCGA.RData")
clin_data <- clinical_factors_merged
load("integrated_annotaion_df.RData")
integ_data <- integrated_annotaion_df

clin_data$study_id_dot <- gsub(clin_data$sample_id, pattern = "-", replacement = ".")
clin_data$hpv16_pos_opscc <- clin_data$study_id_dot %in% integ_data$study_id_dot

small_clin <- clin_data[, names(clin_data) %in% c("study_id_dot", "hpv16_pos_opscc")]

expr_with_var_freq$study_id_dot <- str_sub(expr_with_var_freq$study_id, start = 1L, end = 12L)

merged <- expr_with_var_freq
merged$hpv16_pos_opscc <- merged$study_id_dot %in% integ_data$study_id_dot



merged$perc_A <- merged$n_classic_apobec_m2p_A/merged$n_classic_apobec
merged$perc_T <- merged$n_classic_apobec_m2p_T/merged$n_classic_apobec
merged$perc_C <- merged$n_classic_apobec_m2p_C/merged$n_classic_apobec
merged$perc_G <- merged$n_classic_apobec_m2p_G/merged$n_classic_apobec

merged$perc_TC <- (merged$n_classic_apobec_m2p_T + merged$n_classic_apobec_m2p_C)/merged$n_classic_apobec

merged$n_TC <- (merged$n_classic_apobec_m2p_T + merged$n_classic_apobec_m2p_C)
merged$n_TCA <- (merged$n_classic_apobec_m2p_T + merged$n_classic_apobec_m2p_C + merged$n_classic_apobec_m2p_A)
merged$n_G <- merged$n_classic_apobec_m2p_G
merged$n_T <- merged$n_classic_apobec_m2p_T
merged$n_C <- merged$n_classic_apobec_m2p_C
merged$n_A <- merged$n_classic_apobec_m2p_A
#merged$n_GA <- (merged$n_classic_apobec_m2p_G + merged$n_classic_apobec)
merged$n_GA <- (merged$n_classic_apobec_m2p_G + merged$n_classic_apobec_m2p_A)


merged$A3B <- merged$n_TC <= merged$n_GA
merged$A3A <- merged$n_TC > merged$n_GA
merged$unclear <- merged$n_TC == merged$n_GA
merged$A3B_perc <- merged$n_classic_apobec_m2p_G > 0


merged$apobec_enhanced <- merged$n_classic_apobec > 4



#$sudo umount /mnt/ramdisk


load("TCGA_NFkB_PCA_data.RData")
NFkB_data <- TCGA_NFkB_PCA_data
merged$sample <- truncate_new(merged$study_id)
merged$hpv16_pos_opscc_61 <- merged$sample  %in% NFkB_data$sample
hpv_neg_data <- subset(merged, merged$hpv16_pos_opscc_61 %in% c(FALSE))


apobec_enhanced <- subset(merged, merged$apobec_enhanced %in% c("TRUE"))

hpv_opscc_data <- merge(merged, NFkB_data, by = "sample")
hpv_opscc_data$NFkB_POS <- hpv_opscc_data$NFkB_mod_PCA >= 0


hpv_neg_data <- subset(merged, merged$hpv16_pos_opscc_61 %in% c("FALSE"))
hpv_neg_data <- subset(hpv_neg_data, hpv_neg_data$hpv16_pos_opscc %in% c("FALSE"))
APOBEC_hpv_neg_data <- hpv_neg_data[, names(hpv_neg_data) %in% c("sample", "n_classic_apobec", "A3A", "A3B", "NFkB_mod_PCA", "NFkB_POS")]




APOBEC_hpv_opscc_data <- hpv_opscc_data[, names(hpv_opscc_data) %in% c("sample", "n_classic_apobec", "A3A", "A3B", "NFkB_mod_PCA", "NFkB_POS", "A3B_perc")]


load("apobec_expr_data.RData")
apobec_expr_data$sample <- truncate_new(row.names(apobec_expr_data))
apobec_expr_data$sample_type <- tail_new(row.names(apobec_expr_data))
apobec_expr_data <- subset(apobec_expr_data, !(apobec_expr_data$sample_type %in% c("11", "06")))


APOBEC_hpv_opscc_data <- merge(APOBEC_hpv_opscc_data, apobec_expr_data, by = "sample")
save(APOBEC_hpv_opscc_data, file = "APOBEC_hpv_opscc_data.RData")
APOBEC_hpv_opscc_data$hpv_status <- rep("pos", times = length(APOBEC_hpv_opscc_data$sample))

APOBEC_hpv_neg_data <- merge(APOBEC_hpv_neg_data, apobec_expr_data, by = "sample")
save(APOBEC_hpv_neg_data, file = "APOBEC_hpv_neg_data.RData")
APOBEC_hpv_neg_data$hpv_status <- rep("neg", times = length(APOBEC_hpv_neg_data$sample))


all_APOBEC_data_HNSCC <- rbind(APOBEC_hpv_neg_data, APOBEC_hpv_opscc_data[, !(names(APOBEC_hpv_opscc_data) %in% c("A3B_perc", "NFkB_mod_PCA", "NFkB_POS"))])

library(ggplot2)
library(wesanderson)
my_plot_A3A <- ggplot(all_APOBEC_data_HNSCC, aes(x = hpv_status, y = all_A3A, fill = A3A, color = A3A))+
  geom_boxplot(show.legend = TRUE, outlier.alpha = 0.0)+
  geom_jitter(position=position_jitterdodge(0.1), show.legend = FALSE, size = 0.5, alpha = 0.75)+
  scale_fill_manual(values = c("white", "white"))+
  scale_color_manual(values = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]))+
  theme_classic()

pdf("A3A_box.pdf", height = 3, width = 3)
my_plot_A3A
dev.off()




my_plot_A3B <- ggplot(all_APOBEC_data_HNSCC, aes(x = hpv_status, y = all_A3B, fill = A3A, color = A3A))+
  geom_boxplot(show.legend = TRUE, outlier.alpha = 0.0)+
  geom_jitter(position=position_jitterdodge(0.1), show.legend = FALSE, size = 0.5, alpha = 0.75)+
  scale_fill_manual(values = c("white", "white"))+
  scale_color_manual(values = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]))+
  ylim(-1.5, 8.5)+
  theme_classic()


pdf("A3B_box.pdf", height = 3, width = 3)
my_plot_A3B
dev.off()



all_APOBEC_data_HNSCC$log_n_apobec <- log(all_APOBEC_data_HNSCC$n_classic_apobec + 1, base = 10)
my_plot_nAPOBEC <- ggplot(all_APOBEC_data_HNSCC, aes(x = hpv_status, y = log_n_apobec, fill = A3A, color = A3A))+
  geom_boxplot(show.legend = TRUE, outlier.alpha = 0.0)+
  geom_jitter(position=position_jitterdodge(0.1), show.legend = FALSE, size = 0.5, alpha = 0.75)+
  scale_fill_manual(values = c("white", "white"))+
  scale_color_manual(values = c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2]))+
  theme_classic()


pdf("nAPOBEC_box.pdf", height = 3, width = 3)
my_plot_nAPOBEC
dev.off()




### VAF in apobec vars ###

apobec_vars_individ <- C_REF_ALT_GT_TCA_TCT
apobec_vars_individ$A3A <- apobec_vars_individ$minus2pos %in% c("T", "C")


