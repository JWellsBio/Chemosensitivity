


blca_surv_times <- blca_clinical_cisplatin_short$PFS
blca_status <- ifelse(blca_clinical_cisplatin_short$PFS == blca_clinical_cisplatin_short$OS, 0, 1)

blca_surv_df <- data.frame(blca_surv_times, blca_status, new_blca_tcga_cisplatin)
fit <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_cisplatin,
               data = blca_surv_df)
fit2 <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_cisplatin,
               data = blca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted sensitive', 'predicted resistant'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)



ucs_surv_times <- ucs_clinical_carboplatin_short$PFS
ucs_status <- ifelse(ucs_clinical_carboplatin_short$PFS == ucs_clinical_carboplatin_short$OS, 1, 0)

ucs_code <- rep(NA, nrow(ucs_clinical_carboplatin_short))
for (i in 1:nrow(ucs_clinical_carboplatin_short)) {
  if (new_ucs_tcga_carboplatin_most_sensitive_1se[i] < 0.1970) {
    ucs_code[i] <- 'most_sensitive'
  }
  else if (new_ucs_tcga_carboplatin_least_sensitive_1se[i] > 0.2057) {
    ucs_code[i] <- 'least_sensitive'
  }
  else {
    ucs_code[i] <- 'neither'
  }
}

ucs_surv_df <- data.frame(ucs_surv_times, ucs_status, ucs_code)
fit <- survfit(Surv(ucs_surv_times, ucs_status) ~ ucs_code,
               data = ucs_surv_df)
ggsurvplot(fit, data = ucs_surv_df, risk.table = FALSE, pval = TRUE, test.for.trend = TRUE, censor.shape = 124, risk.table.y.text.col = T, risk.table.y.text = F, legend.title = 'Predicted:', legend.labs = c('Resistant', 'Sensitive', 'Neither'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free')
plot(fit, col = c('limegreen', 'darkviolet', 'orange'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 1000, y = 0.8, legend = 'log-rank\np = 0.0063', bty = 'n', cex = 0.8)
legend(x = 1000, y = 0.3, legend = c('predicted resistant', 'predicted sensitive', 'neither'), lty = c(1,1), lwd = 2, col = c('limegreen', 'darkviolet', 'orange'), bty = 'n', cex = 0.8)

fit <- pairwise_survdiff(Surv(ucs_surv_times, ucs_status) ~ ucs_code,
                         data = ucs_surv_df)
fit
