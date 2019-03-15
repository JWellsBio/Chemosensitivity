install.packages('survminer')
source("https://bioconductor.org/biocLite.R")
biocLite("RTCGA.clinical")

library(survminer)
library(RTCGA.clinical)
survivalTCGA(BRCA.clinical, OV.clinical,
             extract.cols = "admin.disease_code") -> BRCAOV.survInfo
library(survival)
fit <- survfit(Surv(times, patient.vital_status) ~ admin.disease_code,
               data = BRCAOV.survInfo)
# Visualize with survminer
ggsurvplot(fit, data = BRCAOV.survInfo, risk.table = TRUE)

ov_surv_times <- ov_clinical_cisplatin_short$PFS
ov_status <- ifelse(ov_clinical_cisplatin_short$PFS == ov_clinical_cisplatin_short$OS, 1, 0)

ov_code <- rep(NA, nrow(ov_clinical_cisplatin_short))
for (i in 1:nrow(ov_clinical_cisplatin_short)) {
  if (ov_clinical_cisplatin_short$most_sensitive[i] == 1) {
    ov_code[i] <- 'least_sensitive'
  }
  else if (ov_clinical_cisplatin_short$least_sensitive[i] == 1) {
    ov_code[i] <- 'most_sensitive'
  }
  else if (ov_clinical_cisplatin_short$most_sensitive[i] == 0 & ov_clinical_cisplatin_short$least_sensitive[i] == 0) {
    ov_code[i] <- 'neither'
  }
}

ov_surv_df <- data.frame(ov_surv_times, ov_status, ov_code)
ov_surv_df <- ov_surv_df[ov_surv_df$ov_code != 'neither', ]
fit <- survfit(Surv(ov_surv_times, ov_status) ~ ov_code,
               data = ov_surv_df)
ggsurvplot(fit, data = ov_surv_df, risk.table = FALSE, pval = TRUE, test.for.trend = TRUE, censor.shape = 124, risk.table.y.text.col = T, risk.table.y.text = F, legend.title = 'Predicted:', legend.labs = c('Resistant', 'Sensitive', 'Neither'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free')
plot(fit, col = c('limegreen', 'darkviolet', 'orange'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 1500, y = 0.8, legend = 'p = 0.0082', bty = 'n', cex = 0.8)
legend(x = 500, y = 0.4, legend = c('predicted sensitive', 'predicted resistant', 'neither'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen', 'orange'), bty = 'n', cex = 0.8)

fit <- pairwise_survdiff(Surv(ov_surv_times, ov_status) ~ ov_code,
               data = ov_surv_df)
fit

lusc_surv_times <- lusc_clinical_gemcitabine_short$PFS
lusc_status <- ifelse(lusc_clinical_gemcitabine_short$PFS == lusc_clinical_gemcitabine_short$OS, 1, 0)

lusc_code <- rep(NA, nrow(lusc_clinical_gemcitabine_short))
for (i in 1:nrow(lusc_clinical_gemcitabine_short)) {
  if (lusc_clinical_gemcitabine_short$most_sensitive[i] == 1) {
    lusc_code[i] <- 'least_sensitive'
  }
  else if (lusc_clinical_gemcitabine_short$least_sensitive[i] == 1) {
    lusc_code[i] <- 'most_sensitive'
  }
  else if (lusc_clinical_gemcitabine_short$most_sensitive[i] == 0 & lusc_clinical_gemcitabine_short$least_sensitive[i] == 0) {
    lusc_code[i] <- 'neither'
  }
}

lusc_surv_df <- data.frame(lusc_surv_times, lusc_status, lusc_code)
fit <- survfit(Surv(lusc_surv_times, lusc_status) ~ lusc_code,
               data = lusc_surv_df)
ggsurvplot(fit, data = lusc_surv_df, risk.table = FALSE, pval = TRUE, test.for.trend = TRUE, censor.shape = 124, risk.table.y.text.col = T, risk.table.y.text = F, legend.title = 'Predicted:', legend.labs = c('Resistant', 'Sensitive', 'Neither'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free')
plot(fit, col = c('limegreen', 'darkviolet', 'orange'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
