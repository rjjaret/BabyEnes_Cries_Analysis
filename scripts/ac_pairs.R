# Pair plots - correlations between acoustic characteristics per baby across ages (eg some babies consistently have high pitch relative to their age peers)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Model the effects of age on acoustics

library(reshape2)

df = read.csv('../data/dataset_44605_short.csv', stringsAsFactors = TRUE)
df = droplevels(df[df$cause_stop_engl %in% c('discomfort', 'hunger', 'loneliness'), ])

plot_var = function(v, main = '') {
  b = aggregate(df[, v] ~ baby + age_month, df, mean)
  b = dcast(b, baby ~ age_month)
  colnames(b) = c('baby', '0.5 months', '1.5 months', '2.5 months', '3.5 months')
  # pairs(b[, 2:5])
  # p1 = ggpairs(b, columns = 2:5, lower = list(continuous = "smooth"))
  
  xlim = range(b[, 2:5], na.rm = TRUE)
  ylim = range(b[, 2:5], na.rm = TRUE)
  panel.cor = function(x, y){
    # usr <- par("usr"); on.exit(par(usr))
    par(cex.lab = 1.5, xaxs = 'r')
    points(x, y, pch = 19, cex = .75)
    r = round(cor(x, y, use = 'pairwise.complete.obs'), 2)
    # cex.cor = 0.8/strwidth(txt)
    xl = min(x, na.rm = TRUE) # mean(range(x, na.rm = TRUE))
    yl = max(y, na.rm = TRUE)
    dl = diff(range(y, na.rm = TRUE))
    text(xl, yl - .05 * dl, paste0("r = ", r), cex = 1.25, adj = 0) #  cex = cex.cor * r
    text(xl, yl - .2 * dl, paste0('n = ', nrow(na.omit(data.frame(x=x, y=y)))), cex = 1.25, adj = 0)
    mod = lm(y ~ x)
    coef = summary(mod)$coef
    abline(a = coef[1], b = coef[2])
  }
  png(paste0('../pix/pairs_', v, '.png'), width = 12, height = 12, units = 'cm', res = 600)
  pairs(b[, 2:5], upper.panel = panel.cor, lower.panel = NULL, main = main)
  dev.off()
  # write.csv(b, paste0('../data/fig2c_', v, '.csv'))
}

plot_var('pitch_median_S', 'Pitch (Hz)')
plot_var('specCentroidVoiced_median_S', 'Spectral centroid (Hz)')
plot_var('entropyVoiced_median_S', 'Entropy (0 to 1)')
plot_var('duration_noSilence_S', 'Duration (s)')

