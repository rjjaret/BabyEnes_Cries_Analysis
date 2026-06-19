## Acoustic differences between crying based on what caused it

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(brms)
library(ggplot2)
library(patchwork)
source('zz_formatting.R')

df = read.csv('../data/dataset_44605_short.csv', stringsAsFactors = TRUE)
df = droplevels(df[df$cause_stop_engl %in% c('discomfort', 'hunger', 'loneliness'), ])
audioFolder = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'

# for the flow chart
aggregate(baby ~ age_month, df, function(x) length(unique(x)))
aggregate(file_seq_S ~ age_month, df, function(x) length(unique(x)))
aggregate(file_seq_S ~ cause_stop_engl + age_month, df, function(x) length(unique(x)))
table(df$baby)
table(df$cause_stop_engl)
table(df$baby, df$cause_stop_engl)

vars_to_model = sort(c('duration_segment_P', 'voiced_S', 'pitch_median_S', 'pitch_iqr_S', 'specCentroidVoiced_median_S', 'entropyVoiced_median_S', 'HNRVoiced_median_S', 'roughnessVoiced_median_S', 'jitter_P', 'shimmer_P'))
# vars_to_model[which(!vars_to_model %in% colnames(df))]
# summary(df[, vars_to_model])
vars_for_log = c('duration_segment_P', 'domVoiced_median_S', 'domVoiced_iqr_S', 'pitch_median_S', 'pitch_iqr_S', 'quartile25Voiced_median_S', 'quartile50Voiced_median_S', 'quartile75Voiced_median_S', 'specCentroidVoiced_median_S', 'specCentroidVoiced_iqr_S')

# scale per baby
for (v in vars_for_log) df[, v] = log(df[, v] + 1e-3)
for (b in unique(df$baby)) {
  idx_b = which(df$baby == b)
  for (v in vars_to_model) df[idx_b, v] = as.numeric(scale(df[idx_b, v]))
}

for (v in vars_to_model) boxplot(df[, v] ~ df$cause_stop_engl, main = v)


# observed
plot_v = function(v) {
  ggplot(df, aes(x = cause_stop_engl, y = df[, v])) +
    geom_boxplot() + facet_wrap(~sex) + 
    ylab(v) + ggtitle(v) +
    theme_bw() + theme(axis.text = element_text(angle = 45, hjust = 1))
}
# plot_v('pitch_median_S')
paste0("plot_v('", vars_to_model, "')", collapse = ', ')
source('ggplot_multiple.R')
ggsave(multiplot(plot_v('duration_segment_P'), plot_v('entropyVoiced_median_S'), plot_v('HNRVoiced_median_S'), plot_v('jitter_P'), plot_v('pitch_iqr_S'), plot_v('pitch_median_S'), plot_v('roughnessVoiced_median_S'), plot_v('shimmer_P'), plot_v('specCentroidVoiced_median_S'), plot_v('voiced_S'), cols = ceiling(sqrt(length(vars_to_model)))), filename = '../pix/ac~cause*sex_obs.png', width = 40, height = 25, units = 'cm', dpi = 300)

if (FALSE) {
  # Do we need interaction with age? Check with lme4 for pitch_median
  library(lme4)
  mod1 = lmer(pitch_median_S ~ cause_stop_engl + age_month + (age_month + cause_stop_engl|baby), data = df)
  summary(mod1)
  mod2 = lmer(pitch_median_S ~ cause_stop_engl * age_month + (age_month * cause_stop_engl|baby), data = df)
  summary(mod2)
  anova(mod1, mod2)
  #        npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod1   15 110300 110428 -55135   110270                         
  # mod2   28 110035 110275 -54989   109979 290.92 13  < 2.2e-16 ***
  
  
  ## Do we need to account for bout (multiple calls from the same sequence / recording)?
  mod1 = lmer(pitch_median_S ~ cause_stop_engl * age_month + (age_month * cause_stop_engl|baby), data = df)
  mod2 = lmer(pitch_median_S ~ cause_stop_engl * age_month + (age_month * cause_stop_engl|baby) + (1|file_seq_S), data = df)
  anova(mod1, mod2)
  #        npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # mod1   28 110035 110275 -54989   109979                         
  # mod2   29 108887 109135 -54414   108829 1150.2  1  < 2.2e-16 ***
}

# models
myformula = formula(paste(
  'mvbind(', paste(vars_to_model, collapse = ', '), 
  ') ~ cause_stop_engl + age_month + (age_month + cause_stop_engl|baby) + (1|file_seq_S)'
))
myformula

# run model (84 h / 2000 iter)
# mod_cause = brm(myformula, family = 'gaussian', data = df, cores = 3, chains = 3, warmup = 500, iter = 2000, init_r = 0.1)  # , control = list(adapt_delta = 0.95, max_treedepth = 20)
# saveRDS(mod_cause, '../mod/ac~cause*age.RDS')
mod_cause = readRDS('../mod/ac~cause*age.RDS')

# plot(mod_cause)
# summary(mod_cause)
# conditional_effects(mod_cause)
# waic(mod_cause)

# plot effects
rope = c(-0.1, 0.1)
df_rope = data.frame(xmin = rope[1], xmax = rope[2], ymin = -Inf, ymax = Inf)

newdata = expand.grid(cause_stop_engl = levels(df$cause_stop_engl), 
                      age_month = sort(unique((df$age_month))))
fit = fitted(mod_cause, newdata = newdata, re_formula = NA, summary = TRUE)
# dim(fit)
df_plot = expand.grid(cause_stop_engl = levels(df$cause_stop_engl), 
                      age_month = sort(unique((df$age_month))),
                      outcome = vars_to_model)
for (i in 1:nrow(df_plot)) {
  idx_newdata = which(newdata$cause_stop_engl == df_plot$cause_stop_engl[i] &
                        newdata$age_month == df_plot$age_month[i])
  idx_ac = match(df_plot$outcome[i], vars_to_model)
  df_plot[i, c('fit', 'lwr', 'upr')] = fit[idx_newdata, c(1, 3, 4), idx_ac]  # convert back to natural units?
}
# rope = log2(1.1)
# df_rope = data.frame(xmin = -Inf, xmax = Inf, ymin = -rope, ymax = rope)
# cntr$outside_rope = factor((cntr$lwr > rope) | (cntr$upr < -rope), levels = c(T, F))
df_plot$alpha = ifelse((df_plot$lwr > rope[2] | df_plot$upr < rope[1]), 'visible', 'hidden')

df_plot$outcomeFmt = unlist(acVarsFmt[df_plot$outcome])
df_plot$outcomeFmt = factor(df_plot$outcomeFmt, levels = rev(sort(unique(df_plot$outcomeFmt))))
df_plot$ageFmt = paste(df_plot$age_month, 'months')
df_plot$causeFmt = factor(unlist(contextNames[as.character(df_plot$cause_stop_engl)])) 
df_plot[, c('cause_stop_engl', 'age_month', 'outcomeFmt', 'fit', 'lwr', 'upr')]

p_perAge = ggplot(df_plot, aes(y = outcomeFmt, x = fit, xmin = lwr, xmax = upr, color = causeFmt, shape = causeFmt, alpha = alpha)) +  # , color = outside_rope
  geom_rect(data = df_rope, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = NA, fill = 'gray50', alpha = .2, inherit.aes = FALSE) +
  geom_point(position = position_dodge(.5), size = 2.5) +
  geom_errorbar(width = 0, position = position_dodge(.5)) +
  geom_vline(xintercept = 0, linetype = 3) +
  facet_wrap(~ageFmt) +
  scale_alpha_discrete(range = c(.2, 1), guide = 'none') +
  xlab('Effect size, SD') + ylab('') +
  ggtitle('B. Cries per age group') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom', 
        legend.title = element_blank())
# ggsave('../pix/ac~cause*age.png', width = 20, height = 15, units = 'cm', dpi = 600)



## Acoustic differences between crying averaging across age groups
fit1 = fitted(mod_cause, newdata = newdata, re_formula = NA, summary = FALSE)
# dim(fit1)
df_plot1 = expand.grid(cause_stop_engl = levels(df$cause_stop_engl), 
                       age_month = 'All ages',
                       outcome = vars_to_model)
for (i in 1:nrow(df_plot1)) {
  idx_cause = which(newdata$cause_stop_engl == df_plot1$cause_stop_engl[i])
  idx_ac = which(vars_to_model == df_plot1$outcome[i])
  d = rowMeans(fit1[, idx_cause, idx_ac])
  df_plot1[i, c('fit', 'lwr', 'upr')] = quantile(d, probs = c(.5, .025, .975))  # convert back to natural units?
}
df_plot1$alpha = ifelse((df_plot1$lwr > rope[2] | df_plot1$upr < rope[1]), 'visible', 'hidden')
df_plot1$outcomeFmt = unlist(acVarsFmt[df_plot1$outcome])
df_plot1$outcomeFmt = factor(df_plot1$outcomeFmt, levels = rev(sort(unique(df_plot1$outcomeFmt))))
df_plot1$ageFmt = df_plot1$age_month
df_plot1$causeFmt = factor(unlist(contextNames[as.character(df_plot1$cause_stop_engl)])) 

p_overall = ggplot(df_plot1, aes(y = outcomeFmt, x = fit, xmin = lwr, xmax = upr, color = causeFmt, shape = causeFmt, alpha = alpha)) +  # , color = outside_rope
  geom_rect(data = df_rope, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = NA, fill = 'gray50', alpha = .2, inherit.aes = FALSE) +
  geom_point(position = position_dodge(.5), size = 2.5) +
  geom_errorbar(width = 0, position = position_dodge(.5)) +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_alpha_discrete(range = c(.2, 1), guide = 'none') +
  xlab('Effect size, SD') + ylab('Cry features') +
  ggtitle('A. All cries together') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none', 
        legend.title = element_blank())
# ggsave('../pix/ac~cause.png', width = 10, height = 9, units = 'cm', dpi = 600)

p_overall + p_perAge

pl_comb = rbind(df_plot, df_plot1)
pl_comb$ageFmt = factor(pl_comb$ageFmt, levels = c('All ages', '0.5 months', '1.5 months', '2.5 months', '3.5 months'))
ggplot(pl_comb, aes(y = outcomeFmt, x = fit, xmin = lwr, xmax = upr, color = causeFmt, shape = causeFmt, alpha = alpha)) +  # , color = outside_rope
  geom_rect(data = df_rope, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = NA, fill = 'gray50', alpha = .2, inherit.aes = FALSE) +
  geom_point(position = position_dodge(.5), size = 2.5) +
  geom_errorbar(width = 0, position = position_dodge(.5)) +
  geom_vline(xintercept = 0, linetype = 3) +
  facet_wrap(~ageFmt, nrow = 1) +
  scale_alpha_discrete(range = c(.2, 1), guide = 'none') +
  xlab('Effect size, SD') + ylab('Cry features') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'top', 
        legend.title = element_blank())
# ggsave('../pix/ac~cause_perAge.png', width = 20, height = 8, units = 'cm', dpi = 600)

## plot for each individual baby
length(levels(df$baby))
table(df$baby)
newdata_baby = expand.grid(age_month = sort(unique(df$age_month)), 
                           cause_stop_engl = levels(df$cause_stop_engl), 
                           baby = levels(df$baby))  # or might exclude missing levels
fit_baby = fitted(mod_cause, newdata = newdata_baby, re_formula = '~(age_month * cause_stop_engl|baby)', summary = FALSE)
# dim(fit_baby)
df_plot_baby = expand.grid(cause_stop_engl = levels(df$cause_stop_engl), 
                           baby = levels(df$baby),
                           outcome = vars_to_model)
df_plot_baby[, c('fit', 'lwr', 'upr')] = NA
for (i in 1:nrow(df_plot_baby)) {
  idx = which(newdata_baby$baby == df_plot_baby$baby[i] &
                newdata_baby$cause_stop_engl == df_plot_baby$cause_stop_engl[i])
  idx_ac = which(vars_to_model == df_plot_baby$outcome[i])
  d = rowMeans(fit_baby[idx, c(1, 3, 4), idx_ac])
  # d = fit_baby[idx, c(1, 3, 4), idx_ac]
  df_plot_baby[i, c('fit', 'lwr', 'upr')] = quantile(d, probs = c(.5, .025, .975))  # convert back to natural units?
}
df_plot_baby$outcome = factor(df_plot_baby$outcome, levels = rev(sort(vars_to_model)))

df_plot_baby$alpha = ifelse((df_plot_baby$lwr > rope[2] | df_plot_baby$upr < rope[1]), 'visible', 'hidden')
df_plot_baby$outcomeFmt = unlist(acVarsFmt[df_plot_baby$outcome])
df_plot_baby$outcomeFmt = factor(df_plot_baby$outcomeFmt, levels = rev(sort(unique(df_plot_baby$outcomeFmt))))
df_plot_baby$causeFmt = factor(unlist(contextNames[as.character(df_plot_baby$cause_stop_engl)]))
head(df_plot_baby)
# saveRDS(df_plot_baby, '../data/df_plot_baby.RDS')

p_baby = ggplot(df_plot_baby, aes(y = outcomeFmt, x = fit, xmin = lwr, xmax = upr, color = causeFmt, shape = causeFmt, alpha = alpha)) +  # , color = outside_rope
  geom_rect(data = df_rope, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = NA, fill = 'gray50', alpha = .25, inherit.aes = FALSE) +
  geom_point() +
  geom_point(position = position_dodge(.5), size = 2) +
  geom_errorbar(width = 0, position = position_dodge(.5)) +
  geom_vline(xintercept = 0, linetype = 3) +
  facet_wrap(~baby, ncol = 6) +
  scale_alpha_discrete(range = c(.2, 1), guide = 'none') +
  xlab('') + ylab('') +
  ggtitle('C. Per baby') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        legend.title = element_blank())
# ggsave('../pix/ac~cause_per-baby.png', width = 20, height = 20, units = 'cm', dpi = 600)


## Composite plot
design = "
  1222
  3333
  3333
"
ylim = max(abs(df_plot[, c('lwr', 'upr')]))
p_overall + p_perAge + p_baby +
  plot_layout(nrow = 1, design = design) & 
  scale_x_continuous(limits = c(-ylim, ylim), n.breaks = 4) &
  theme(plot.margin = unit(c(.1, .1, 0, .1), 'cm'), 
        plot.title.position = 'plot', 
        legend.margin = margin(-.75, .1, 0, .1, 'cm'))
# ggsave('../pix/ac~cause_composite.png', width = 18, height = 24, units = 'cm', dpi = 600)



## Do babies cluster by "acoustic strategy" of encoding context?
colnames(df_plot_baby)
bc = expand.grid(cause_stop_engl = unique(df_plot_baby$cause_stop_engl),
                 baby = unique(df_plot_baby$baby))
for (v in unique(df_plot_baby$outcome)) {
  temp = df_plot_baby[df_plot_baby$outcome == v, ]
  bc[, v] = temp$fit
}

# if groups of babies use similar strategies, these babies should form clusters for each context, eg babies 1 & 2 should be closer than babies 2 & 3 in each context.
for (context in unique(bc$cause_stop_engl)) {
  assign(paste0('dist_', context), 
         dist(bc[bc$cause_stop_engl == context, vars_to_model]))
}
dists = cbind(as.numeric(dist_discomfort), as.numeric(dist_hunger), as.numeric(dist_loneliness))
colnames(dists) = c('discomfort', 'hunger', 'loneliness')
cd = cor(dists)
png('../pix/cor_bw_babies.png', width = 12, height = 12, units = 'cm', res = 600)
corrplot::corrplot(cd, method = "color", addCoef.col = "black", tl.srt = 0, type = 'upper') 
dev.off()
#            discomfort      hunger  loneliness
# discomfort 1.00000000  0.03963242  0.05686247
# hunger     0.03963242  1.00000000 -0.02474608
# loneliness 0.05686247 -0.02474608  1.00000000
# write.csv(cd, '../data/fig3d_corMat.csv', row.names = FALSE)

source('umap_helperFun.R')
mydist = dist(bc[, vars_to_model])
# u = as.data.frame(umap(mydist, init = "spectral"))
# saveRDS(u, '../data/umap_baby-cluster-from-dist.RDS')
u = readRDS('../data/umap_baby-cluster-from-dist.RDS')
addVars = c('cause_stop_engl', 'baby')
u[, addVars] = bc[, addVars]

ggplot(u, aes(x = V1, y = V2, color = cause_stop_engl)) +
  geom_point() +
  xlab('') + ylab('') +
  theme_bw()

ggplot(u, aes(x = V1, y = V2, color = cause_stop_engl, label = as.numeric(as.factor(baby)))) +
  geom_text(size = 3) +
  facet_wrap(~cause_stop_engl) +
  xlab('') + ylab('') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none') # no visible clusters
# ggsave('../pix/umap_dist_bw_babies.png', width = 12, height = 5, units = 'cm', dpi = 600)
# write.csv(u, '../data/fig3d_umap.csv', row.names = FALSE)
