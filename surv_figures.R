#script for survival curves with log-rank p values aired with distribution curves with ks p values

### load libraries ----
if (!require('survival')) install.packages('survival')
library(survival) #for generating survival curves

if (!require('survminer')) install.packages('survminer')
library(survminer) #for generating customized surviva plots

if (!require('patchwork')) install.packages('patchwork')
library(patchwork) #for putting plots together

#function to swap a few data sets when needed
swap <- function(vec, from, to) {
  tmp <- to[ match(vec, from) ]
  tmp[is.na(tmp)] <- vec[is.na(tmp)]
  return(tmp)
}


## blca w cisplatin 47 ----
blca_cisplatin_surv_df <- read.table('Survival_Data/blca_cisplatin_surv.df.txt', header = TRUE)

#fit survival curves
blca_cisplatin_fit <- survfit(Surv(blca_surv_times, blca_status) ~ X1,
               data = blca_cisplatin_surv_df)

#plot them
surv_plot <- ggsurvplot(blca_cisplatin_fit, data = blca_cisplatin_surv_df, 
           size = 1, #change line size
           palette = c('limegreen', 'purple'), #custom color scheme
           conf.int = TRUE, # add CIs
           pval = TRUE, 
           pval.method = TRUE,
           pval.size = 6, 
           risk.table = TRUE, # add p-value
           legend = 'bottom', 
           legend.labs = c('Sensitive', 'Resistant'), 
           xlab = 'Time (d)', 
           ylab = 'Proportion\nRecurrence-Free', 
           conf.int.style = 'ribbon',
           ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 2500, y = 0.25, 
            label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                    axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                    axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                    axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                    legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                    legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                    panel.grid.major = element_blank(),
                                                    panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
                                          theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                          axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                          axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
                                          axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                          plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank())


#ks test for dist curves
blca_cisplatin_ks_test <- ks.test(blca_cisplatin_surv_df$blca_surv_times[blca_cisplatin_surv_df$X1 == 1], 
                                  blca_cisplatin_surv_df$blca_surv_times[blca_cisplatin_surv_df$X1 == 0])
blca_cisplatin_ks_test_p_value <- round(blca_cisplatin_ks_test$p.value, digits = 6)

#density plot
dens_plot <- ggplot(blca_cisplatin_surv_df, aes(x = blca_cisplatin_surv_df$blca_surv_times, fill = factor(blca_cisplatin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 3500, y = 0.0013, 
           label = paste0('KS\np = ', blca_cisplatin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/blca_cisplatin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
  AAACCC
  AAACCC
  AAACCC
  AAACCC
  AAACCC
  BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()



## cesc w cisplatin 100 ----
cesc_cisplatin_surv_df <- read.table('Survival_Data/cesc_cisplatin_surv.df.txt', header = TRUE)

#fit survival curves
cesc_cisplatin_fit <- survfit(Surv(cesc_surv_times, cesc_status) ~ X1,
                              data = cesc_cisplatin_surv_df)

#plot them
surv_plot <- ggsurvplot(cesc_cisplatin_fit, data = cesc_cisplatin_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 2500, y = 0.25, 
           label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#ks test for dist curves
cesc_cisplatin_ks_test <- ks.test(cesc_cisplatin_surv_df$cesc_surv_times[cesc_cisplatin_surv_df$X1 == 1], 
                                  cesc_cisplatin_surv_df$cesc_surv_times[cesc_cisplatin_surv_df$X1 == 0])
cesc_cisplatin_ks_test_p_value <- round(cesc_cisplatin_ks_test$p.value, digits = 3)

#density plot
dens_plot <- ggplot(cesc_cisplatin_surv_df, aes(x = cesc_cisplatin_surv_df$cesc_surv_times, fill = factor(cesc_cisplatin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 3500, y = 7.5e-4, 
           label = paste0('KS\np = ', cesc_cisplatin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/cesc_cisplatin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## esca w cisplatin 10 ----
esca_cisplatin_surv_df <- read.table('Survival_Data/esca_cisplatin_surv_df.txt', header = TRUE)

esca_cisplatin_surv_df$X1 <- swap(esca_cisplatin_surv_df$X1, c(0, 1), c(1, 0))

#fit survival curves
esca_cisplatin_fit <- survfit(Surv(esca_surv_times, esca_status) ~ X1,
                              data = esca_cisplatin_surv_df)


#plot them
surv_plot <- ggsurvplot(esca_cisplatin_fit, data = esca_cisplatin_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 1*365, linetype = 'dashed') + 
  annotate('text', x = 500, y = 0.25, 
           label = '1-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                 axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                 axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                 axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                 legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                 legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                 panel.grid.major = element_blank(),
                                                                 panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
esca_cisplatin_ks_test <- ks.test(esca_cisplatin_surv_df$esca_surv_times[esca_cisplatin_surv_df$X1 == 1], 
                                  esca_cisplatin_surv_df$esca_surv_times[esca_cisplatin_surv_df$X1 == 0])
esca_cisplatin_ks_test_p_value <- round(esca_cisplatin_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(esca_cisplatin_surv_df, aes(x = esca_cisplatin_surv_df$esca_surv_times, fill = factor(esca_cisplatin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 750, y = 0.015, 
           label = paste0('KS\np = ', esca_cisplatin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/esca_cisplatin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## hnsc w cisplatin 65 ----
hnsc_cisplatin_surv_df <- read.table('Survival_Data/hnsc_cisplatin_surv_df.txt', header = TRUE)

#fit survival curves
hnsc_cisplatin_fit <- survfit(Surv(hnsc_surv_times, hnsc_status) ~ X1,
                              data = hnsc_cisplatin_surv_df)

#plot them
surv_plot <- ggsurvplot(hnsc_cisplatin_fit, data = hnsc_cisplatin_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 2300, y = 0.25, 
           label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
hnsc_cisplatin_ks_test <- ks.test(hnsc_cisplatin_surv_df$hnsc_surv_times[hnsc_cisplatin_surv_df$X1 == 1], 
                                  hnsc_cisplatin_surv_df$hnsc_surv_times[hnsc_cisplatin_surv_df$X1 == 0])
hnsc_cisplatin_ks_test_p_value <- round(hnsc_cisplatin_ks_test$p.value, digits = 3)

#density plot
dens_plot <- ggplot(hnsc_cisplatin_surv_df, aes(x = hnsc_cisplatin_surv_df$hnsc_surv_times, fill = factor(hnsc_cisplatin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 2000, y = 0.0008, 
           label = paste0('KS\np = ', hnsc_cisplatin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/hnsc_cisplatin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## luad w cisplatin 23 ----
luad_cisplatin_surv_df <- read.table('Survival_Data/luad_cisplatin_surv_df.txt', header = TRUE)

#fit survival curves
luad_cisplatin_fit <- survfit(Surv(luad_surv_times, luad_status) ~ X1,
                              data = luad_cisplatin_surv_df)

#plot them
surv_plot <- ggsurvplot(luad_cisplatin_fit, data = luad_cisplatin_surv_df, 
                        size = 1, #change line size
                        palette = c('purple', 'limegreen'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Resistant', 'Sensitive'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 2300, y = 0.25, 
           label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
luad_cisplatin_ks_test <- ks.test(luad_cisplatin_surv_df$luad_surv_times[luad_cisplatin_surv_df$X1 == 1], 
                                  luad_cisplatin_surv_df$luad_surv_times[luad_cisplatin_surv_df$X1 == 0])
luad_cisplatin_ks_test_p_value <- round(luad_cisplatin_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(luad_cisplatin_surv_df, aes(x = luad_cisplatin_surv_df$luad_surv_times, fill = factor(luad_cisplatin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('purple', 'limegreen'), labels = c('Resistant', 'Sensitive')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 2000, y = 0.0012, 
           label = paste0('KS\np = ', luad_cisplatin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/luad_cisplatin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## ov w cisplatin 17 ----
ov_cisplatin_surv_df <- read.table('Survival_Data/ov_cisplatin_surv_df.txt', header = TRUE)

#fit survival curves
ov_cisplatin_fit <- survfit(Surv(ov_surv_times, ov_status) ~ X1,
                              data = ov_cisplatin_surv_df)

#plot them
surv_plot <- ggsurvplot(ov_cisplatin_fit, data = ov_cisplatin_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 1500, y = 0.85, 
           label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#ks test for dist curves
ov_cisplatin_ks_test <- ks.test(ov_cisplatin_surv_df$ov_surv_times[ov_cisplatin_surv_df$X1 == 1], 
                                  ov_cisplatin_surv_df$ov_surv_times[ov_cisplatin_surv_df$X1 == 2])
ov_cisplatin_ks_test_p_value <- round(ov_cisplatin_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(ov_cisplatin_surv_df, aes(x = ov_cisplatin_surv_df$ov_surv_times, fill = factor(ov_cisplatin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 1300, y = 0.0008, 
           label = paste0('KS\np = ', ov_cisplatin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/ov_cisplatin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## tgct w cisplatin 48 ----
tgct_cisplatin_surv_df <- read.table('Survival_Data/tgct_cisplatin_surv_df.txt', header = TRUE)

#fit survival curves
tgct_cisplatin_fit <- survfit(Surv(tgct_surv_times, tgct_status) ~ X1,
                            data = tgct_cisplatin_surv_df)

#plot them
surv_plot <- ggsurvplot(tgct_cisplatin_fit, data = tgct_cisplatin_surv_df, 
                        size = 1, #change line size
                        palette = c('purple', 'limegreen'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Resistant', 'Sensitive'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 3100, y = 0.25, 
           label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#ks test for dist curves
tgct_cisplatin_ks_test <- ks.test(tgct_cisplatin_surv_df$tgct_surv_times[tgct_cisplatin_surv_df$X1 == 1], 
                                tgct_cisplatin_surv_df$tgct_surv_times[tgct_cisplatin_surv_df$X1 == 0])
tgct_cisplatin_ks_test_p_value <- round(tgct_cisplatin_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(tgct_cisplatin_surv_df, aes(x = tgct_cisplatin_surv_df$tgct_surv_times, fill = factor(tgct_cisplatin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('purple', 'limegreen'), labels = c('Resistant', 'Sensitive')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 4000, y = 1.9e-4, 
           label = paste0('KS\np = ', tgct_cisplatin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/tgct_cisplatin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## brca w doxorubicin 45 ----
brca_doxorubicin_surv_df <- read.table('Survival_Data/brca_doxorubicin_surv_df.txt', header = TRUE)
brca_doxorubicin_surv_df$X1 <- swap(brca_doxorubicin_surv_df$X1, c('sensitive', 'resistant'), c('resistant', 'sensitive'))
#fit survival curves
brca_doxorubicin_fit <- survfit(Surv(brca_surv_times, brca_status) ~ X1,
                              data = brca_doxorubicin_surv_df)

#plot them
surv_plot <- ggsurvplot(brca_doxorubicin_fit, data = brca_doxorubicin_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 2300, y = 0.25, 
           label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
brca_doxorubicin_ks_test <- ks.test(brca_doxorubicin_surv_df$brca_surv_times[brca_doxorubicin_surv_df$X1 == 'sensitive'], 
                                  brca_doxorubicin_surv_df$brca_surv_times[brca_doxorubicin_surv_df$X1 == 'resistant'])
brca_doxorubicin_ks_test_p_value <- round(brca_doxorubicin_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(brca_doxorubicin_surv_df, aes(x = brca_doxorubicin_surv_df$brca_surv_times, fill = factor(brca_doxorubicin_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 3000, y = 4e-4, 
           label = paste0('KS\np = ', brca_doxorubicin_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/brca_doxorubicin_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## luad w etoposide 6 ----
luad_etoposide_surv_df <- read.table('Survival_Data/luad_etoposide_surv_df.txt', header = TRUE)
luad_etoposide_surv_df$X1 <- swap(luad_etoposide_surv_df$X1, c('resistant', 'sensitive'), c('sensitive', 'resistant'))
#fit survival curves
luad_etoposide_fit <- survfit(Surv(luad_surv_times, luad_status) ~ X1,
                                data = luad_etoposide_surv_df)

#plot them
surv_plot <- ggsurvplot(luad_etoposide_fit, data = luad_etoposide_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 3*365, linetype = 'dashed') + 
  annotate('text', x = 1300, y = 0.25, 
           label = '3-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
luad_etoposide_ks_test <- ks.test(luad_etoposide_surv_df$luad_surv_times[luad_etoposide_surv_df$X1 == 'sensitive'], 
                                    luad_etoposide_surv_df$luad_surv_times[luad_etoposide_surv_df$X1 == 'resistant'])
luad_etoposide_ks_test_p_value <- round(luad_etoposide_ks_test$p.value, digits = 6)

#density plot
dens_plot <- ggplot(luad_etoposide_surv_df, aes(x = luad_etoposide_surv_df$luad_surv_times, fill = factor(luad_etoposide_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 300, y = 1e-3, 
           label = paste0('KS\np = ', luad_etoposide_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/luad_etoposide_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## blca w gemcitabine 29 ----
blca_gemcitabine_surv_df <- read.table('Survival_Data/blca_gemcitabine_surv_df.txt', header = TRUE)

#fit survival curves
blca_gemcitabine_fit <- survfit(Surv(blca_surv_times, blca_status) ~ X1,
                              data = blca_gemcitabine_surv_df)

#plot them
surv_plot <- ggsurvplot(blca_gemcitabine_fit, data = blca_gemcitabine_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 3*365, linetype = 'dashed') + 
  annotate('text', x = 1700, y = 0.25, 
           label = '3-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
blca_gemcitabine_ks_test <- ks.test(blca_gemcitabine_surv_df$blca_surv_times[blca_gemcitabine_surv_df$X1 == 'sensitive'], 
                                  blca_gemcitabine_surv_df$blca_surv_times[blca_gemcitabine_surv_df$X1 == 'resistant'])
blca_gemcitabine_ks_test_p_value <- round(blca_gemcitabine_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(blca_gemcitabine_surv_df, aes(x = blca_gemcitabine_surv_df$blca_surv_times, fill = factor(blca_gemcitabine_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 2000, y = 0.0008, 
           label = paste0('KS\np = ', blca_gemcitabine_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/blca_gemcitabine_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()



## chol w gemcitabine 5 ----
chol_gemcitabine_surv_df <- read.table('Survival_Data/chol_gemcitabine_surv_df.txt', header = TRUE)

#fit survival curves
chol_gemcitabine_fit <- survfit(Surv(chol_surv_times, chol_status) ~ X1,
                                data = chol_gemcitabine_surv_df)

#plot them
surv_plot <- ggsurvplot(chol_gemcitabine_fit, data = chol_gemcitabine_surv_df, 
                        size = 1, #change line size
                        palette = c('purple', 'limegreen'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Resistant', 'Sensitive'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 3*365, linetype = 'dashed') + 
  annotate('text', x = 1300, y = 0.25, 
           label = '3-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#ks test for dist curves
chol_gemcitabine_ks_test <- ks.test(chol_gemcitabine_surv_df$chol_surv_times[chol_gemcitabine_surv_df$X1 == 'sensitive'], 
                                    chol_gemcitabine_surv_df$chol_surv_times[chol_gemcitabine_surv_df$X1 == 'resistant'])
chol_gemcitabine_ks_test_p_value <- round(chol_gemcitabine_ks_test$p.value, digits = 6)

#density plot
dens_plot <- ggplot(chol_gemcitabine_surv_df, aes(x = chol_gemcitabine_surv_df$chol_surv_times, fill = factor(chol_gemcitabine_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('purple', 'limegreen'), labels = c('Resistant', 'Sensitive')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 400, y = 0.0008, 
           label = paste0('KS\np = ', chol_gemcitabine_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/chol_gemcitabine_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()

## lihc w gemcitabine 8 ----
lihc_gemcitabine_surv_df <- read.table('Survival_Data/lihc_gemcitabine_surv_df.txt', header = TRUE)

#fit survival curves
lihc_gemcitabine_fit <- survfit(Surv(lihc_surv_times, lihc_status) ~ X1,
                                data = lihc_gemcitabine_surv_df)

#plot them
surv_plot <- ggsurvplot(lihc_gemcitabine_fit, data = lihc_gemcitabine_surv_df, 
                        size = 1, #change line size
                        palette = c('purple', 'limegreen'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Resistant', 'Sensitive'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 3*365, linetype = 'dashed') + 
  annotate('text', x = 1300, y = 0.25, 
           label = '3-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
lihc_gemcitabine_ks_test <- ks.test(lihc_gemcitabine_surv_df$lihc_surv_times[lihc_gemcitabine_surv_df$X1 == 'sensitive'], 
                                    lihc_gemcitabine_surv_df$lihc_surv_times[lihc_gemcitabine_surv_df$X1 == 'resistant'])
lihc_gemcitabine_ks_test_p_value <- round(lihc_gemcitabine_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(lihc_gemcitabine_surv_df, aes(x = lihc_gemcitabine_surv_df$lihc_surv_times, fill = factor(lihc_gemcitabine_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('purple', 'limegreen'), labels = c('Resistant', 'Sensitive')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 275, y = 0.0045, 
           label = paste0('KS\np = ', lihc_gemcitabine_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/lihc_gemcitabine_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()


## luad w gemcitabine 13 ----
luad_gemcitabine_surv_df <- read.table('Survival_Data/luad_gemcitabine_surv_df.txt', header = TRUE)
luad_gemcitabine_surv_df$X1 <- swap(luad_gemcitabine_surv_df$X1, c('sensitive', 'resistant'), c('resistant', 'sensitive'))
#fit survival curves
luad_gemcitabine_fit <- survfit(Surv(luad_surv_times, luad_status) ~ X1,
                                data = luad_gemcitabine_surv_df)

#plot them
surv_plot <- ggsurvplot(luad_gemcitabine_fit, data = luad_gemcitabine_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 3*365, linetype = 'dashed') + 
  annotate('text', x = 900, y = 0.25, 
           label = '3-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
luad_gemcitabine_ks_test <- ks.test(luad_gemcitabine_surv_df$luad_surv_times[luad_gemcitabine_surv_df$X1 == 'sensitive'], 
                                    luad_gemcitabine_surv_df$luad_surv_times[luad_gemcitabine_surv_df$X1 == 'resistant'])
luad_gemcitabine_ks_test_p_value <- round(luad_gemcitabine_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(luad_gemcitabine_surv_df, aes(x = luad_gemcitabine_surv_df$luad_surv_times, fill = factor(luad_gemcitabine_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 800, y = 0.0015, 
           label = paste0('KS\np = ', luad_gemcitabine_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/luad_gemcitabine_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()

## lusc w gemcitabine 18 ----
lusc_gemcitabine_surv_df <- read.table('Survival_Data/lusc_gemcitabine_surv_df.txt', header = TRUE)

#fit survival curves
lusc_gemcitabine_fit <- survfit(Surv(lusc_surv_times, lusc_status) ~ X1,
                                data = lusc_gemcitabine_surv_df)

#plot them
surv_plot <- ggsurvplot(lusc_gemcitabine_fit, data = lusc_gemcitabine_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 3*365, linetype = 'dashed') + 
  annotate('text', x = 1400, y = 0.25, 
           label = '3-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
lusc_gemcitabine_ks_test <- ks.test(lusc_gemcitabine_surv_df$lusc_surv_times[lusc_gemcitabine_surv_df$X1 == 'sensitive'], 
                                    lusc_gemcitabine_surv_df$lusc_surv_times[lusc_gemcitabine_surv_df$X1 == 'resistant'])
lusc_gemcitabine_ks_test_p_value <- round(lusc_gemcitabine_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(lusc_gemcitabine_surv_df, aes(x = lusc_gemcitabine_surv_df$lusc_surv_times, fill = factor(lusc_gemcitabine_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 1200, y = 0.0025, 
           label = paste0('KS\np = ', lusc_gemcitabine_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/lusc_gemcitabine_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()

## ov w gemcitabine 17 ----
ov_gemcitabine_surv_df <- read.table('Survival_Data/ov_gemcitabine_surv_df.txt', header = TRUE)

#fit survival curves
ov_gemcitabine_fit <- survfit(Surv(ov_surv_times, ov_status) ~ X1,
                                data = ov_gemcitabine_surv_df)

#plot them
surv_plot <- ggsurvplot(ov_gemcitabine_fit, data = ov_gemcitabine_surv_df, 
                        size = 1, #change line size
                        palette = c('purple', 'limegreen'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Resistant', 'Sensitive'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 3*365, linetype = 'dashed') + 
  annotate('text', x = 1400, y = 0.25, 
           label = '3-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
ov_gemcitabine_ks_test <- ks.test(ov_gemcitabine_surv_df$ov_surv_times[ov_gemcitabine_surv_df$X1 == 'sensitive'], 
                                    ov_gemcitabine_surv_df$ov_surv_times[ov_gemcitabine_surv_df$X1 == 'resistant'])
ov_gemcitabine_ks_test_p_value <- round(ov_gemcitabine_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(ov_gemcitabine_surv_df, aes(x = ov_gemcitabine_surv_df$ov_surv_times, fill = factor(ov_gemcitabine_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('purple', 'limegreen'), labels = c('Resistant', 'Sensitive')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 1200, y = 0.0018, 
           label = paste0('KS\np = ', ov_gemcitabine_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/ov_gemcitabine_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()

## paad w gemcitabine 72 ----
paad_gemcitabine_surv_df <- read.table('Survival_Data/paad_gemcitabine_surv_df.txt', header = TRUE)
paad_gemcitabine_surv_df$X1 <- swap(paad_gemcitabine_surv_df$X1, c('sensitive', 'resistant'), c('resistant', 'sensitive'))
#fit survival curves
paad_gemcitabine_fit <- survfit(Surv(paad_surv_times, paad_status) ~ X1,
                              data = paad_gemcitabine_surv_df)

#plot them
surv_plot <- ggsurvplot(paad_gemcitabine_fit, data = paad_gemcitabine_surv_df, 
                        size = 1, #change line size
                        palette = c('limegreen', 'purple'), #custom color scheme
                        conf.int = TRUE, # add CIs
                        pval = TRUE, 
                        pval.method = TRUE,
                        pval.size = 6, 
                        risk.table = TRUE, # add p-value
                        legend = 'bottom', 
                        legend.labs = c('Sensitive', 'Resistant'), 
                        xlab = 'Time (d)', 
                        ylab = 'Proportion\nRecurrence-Free', 
                        conf.int.style = 'ribbon',
                        ggtheme = theme_bw())

#customize and ready for patchwork
surv_plot_new <- surv_plot$plot + geom_vline(xintercept = 5*365, linetype = 'dashed') + 
  annotate('text', x = 2200, y = 0.25, 
           label = '5-year\nRFS', size = 6, fontface = 2) + theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                                                                  axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                                                                  axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                                                                  legend.title = element_text(size = 1, face = 'bold', colour = 'white'), 
                                                                  panel.grid.major = element_blank(),
                                                                  panel.grid.minor = element_blank())

surv_table_new <- surv_plot$table + labs(title = 'Number at Risk of Recurrence') + 
  theme(axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
        axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
        axis.title.y = element_text(size = 16, face = 'bold', colour = 'white'),
        axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
        plot.title = element_text(size = 16, face = 'bold', colour = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ks test for dist curves
paad_gemcitabine_ks_test <- ks.test(paad_gemcitabine_surv_df$paad_surv_times[paad_gemcitabine_surv_df$X1 == 'sensitive'], 
                                  paad_gemcitabine_surv_df$paad_surv_times[paad_gemcitabine_surv_df$X1 == 'resistant'])
paad_gemcitabine_ks_test_p_value <- round(paad_gemcitabine_ks_test$p.value, digits = 2)

#density plot
dens_plot <- ggplot(paad_gemcitabine_surv_df, aes(x = paad_gemcitabine_surv_df$paad_surv_times, fill = factor(paad_gemcitabine_surv_df$X1))) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = c('limegreen', 'purple'), labels = c('Sensitive', 'Resistant')) + theme_bw() + 
  labs(x = 'Recurrence-Free Survival (d)', y = 'Density', fill = '') + 
  annotate('text', x = 1500, y = 0.0008, 
           label = paste0('KS\np = ', paad_gemcitabine_ks_test_p_value), 
           size = 5) + theme(legend.position = 'bottom', 
                             axis.title.x = element_text(size = 16, face = 'bold', colour = 'black'), 
                             axis.text.x = element_text(size = 14, face = 'bold', colour = 'black'), 
                             axis.title.y = element_text(size = 16, face = 'bold', colour = 'black'),
                             axis.text.y = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.text = element_text(size = 14, face = 'bold', colour = 'black'), 
                             legend.title = element_text(size = 16, face = 'bold', colour = 'black'), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

#plot with patchwork
tiff(filename = 'Images/paad_gemcitabine_surv.tiff', units = 'in', width = 10, height = 6, res = 300)
layout <- '
AAACCC
AAACCC
AAACCC
AAACCC
AAACCC
BBBCCC'
surv_plot_new + surv_table_new + dens_plot + plot_layout(design = layout)
dev.off()

