## plotting CIs around over all accuracy

library(plotrix)

drugs <- 1:11
Accuracy <- c(0.92, 0.94, 0.88, 0.90, 0.89, 0.94, 0.92, 0.92, 0.94, 0.92, 0.93)
lower_limit <- c(0.90, 0.92, 0.85, 0.87, 0.86, 0.92, 0.90, 0.90, 0.92, 0.89, 0.91)
upper_limit <- c(0.94, 0.96, 0.90, 0.92, 0.91, 0.96, 0.94, 0.94, 0.95, 0.93, 0.95)
ci_colors <- c('red', 'blue', 'green', 'purple', 'black', 'orange', 'hotpink', 'magenta', 'dodgerblue', 'black', 'limegreen')

tiff(filename = 'ci_plot.tiff', units = 'in', width = 10, height = 6, res = 300)
plotCI(drugs, Accuracy, ui=upper_limit, li=lower_limit, 
       ylab = 'Accuracy w 95% CI', xlab = '', ylim = c(0.85, 1.0), bty = 'n', xaxt = 'n', lwd = 3, font.axis = 2, cex.axis = 1.4, font.lab = 2, cex.lab = 1.4)
dev.off()
