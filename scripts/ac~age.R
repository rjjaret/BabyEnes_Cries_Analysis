## Model the effects of age on acoustics

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(brms)
library(ggplot2)
library(randomForest)
source('zz_formatting.R')

df = read.csv('../data/dataset_44605_short.csv', stringsAsFactors = TRUE)
df = droplevels(df[df$cause_stop_engl %in% c('discomfort', 'hunger', 'loneliness'), ])
audioFolder = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'

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


# observed
ggplot(df, aes(x = age_month, y = pitch_iqr_S, linetype = sex)) +
  geom_smooth(method = 'lm') +
  theme_bw()

# observed smooth
plot_v = function(v) {
  ggplot(df, aes(x = age_month, y = df[, v], linetype = sex)) +
    geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'cs', k = 3)) +
    ylab(v) + ggtitle(v) +
    theme_bw()
}
# plot_v('pitch_median_S')
paste0("plot_v('", vars_to_model, "')", collapse = ', ')
source('ggplot_multiple.R')
ggsave(multiplot(plot_v('duration_segment_P'), plot_v('entropyVoiced_median_S'), plot_v('HNRVoiced_median_S'), plot_v('jitter_P'), plot_v('pitch_iqr_S'), plot_v('pitch_median_S'), plot_v('roughnessVoiced_median_S'), plot_v('shimmer_P'), plot_v('specCentroidVoiced_median_S'), plot_v('voiced_S'), cols = ceiling(sqrt(length(vars_to_model)))), filename = '../pix/ac~age*sex_obs.png', width = 40, height = 25, units = 'cm', dpi = 300)

boxplot(pitch_median_S ~ sex + age, df)


# observed boxplots per age group
df$age = as.factor(df$age_month)
plot_v2 = function(v) {
  ggplot(df, aes(x = age, y = df[, v], fill = sex)) +
    geom_boxplot(position = position_dodge(.8), notch = TRUE) +
    ylab(v) + ggtitle(v) +
    theme_bw()
}
# plot_v2('pitch_median_S')
paste0("plot_v2('", vars_to_model, "')", collapse = ', ')
source('ggplot_multiple.R')

library(ggpubr)
ggpubr::ggarrange(plot_v2('duration_segment_P'), plot_v2('entropyVoiced_median_S'), plot_v2('HNRVoiced_median_S'), plot_v2('jitter_P'), plot_v2('pitch_iqr_S'), plot_v2('pitch_median_S'), plot_v2('roughnessVoiced_median_S'), plot_v2('shimmer_P'), plot_v2('specCentroidVoiced_median_S'), plot_v2('voiced_S'), common.legend = TRUE)

ggpubr::ggarrange(plot_v2('duration_segment_P'), plot_v2('entropyVoiced_median_S'))


ggsave(ggpubr::ggarrange(plot_v2('duration_segment_P'), plot_v2('entropyVoiced_median_S'), plot_v2('HNRVoiced_median_S'), plot_v2('jitter_P'), plot_v2('pitch_iqr_S'), plot_v2('pitch_median_S'), plot_v2('roughnessVoiced_median_S'), plot_v2('shimmer_P'), plot_v2('specCentroidVoiced_median_S'), plot_v2('voiced_S'), common.legend = TRUE), filename = '../pix/ac~age*sex_obs_boxplots.png', width = 40, height = 25, units = 'cm', dpi = 300)



# model controlling for context (re-use the model from ac~cause.R)
mod_age = readRDS('../mod/ac~cause.RDS')

summary(mod_age)
conditional_effects(mod_age)

# plot contrasts
newdata = data.frame(age_month = c(1.5, 2.5),  # range(df$age_month),
                     cause_stop_engl = 'discomfort')  # any cause_stop_engl will do - identical results
fit_age = fitted(mod_age, newdata = newdata, re_formula = NA, summary = FALSE)
cntr = data.frame(predictor = vars_to_model)
for (i in 1:nrow(cntr)) {
  d = fit_age[, 2, i] - fit_age[, 1, i]
  cntr[i, c('fit', 'lwr', 'upr')] = quantile(d, probs = c(.5, .025, .975))
}
cntr$predictor = factor(cntr$predictor, levels = rev(sort(unique(cntr$predictor))))

rope = .1 / 3
df_rope = data.frame(xmin = -Inf, xmax = Inf, ymin = -rope, ymax = rope)
cntr$outside_rope = as.factor((cntr$lwr > rope) | (cntr$upr < -rope))
# levels(cntr$outside_rope) = c(T, F)

# contrasts per baby
newdata_baby = expand.grid(age_month = c(1.5, 2.5),  # range(df$age_month),
                           baby = levels(df$baby),
                           cause_stop_engl = 'discomfort')
fit_age_baby = fitted(mod_age, newdata = newdata_baby, summary = FALSE)
dim(fit_age_baby)  # 1 = MCMC, 2 = baby * age, 3 = vars_to_model
cntr_baby = expand.grid(predictor = vars_to_model,
                        baby = levels(df$baby))
for (i in 1:nrow(cntr_baby)) {
  idx_pred = which(sort(levels(cntr$predictor)) == cntr_baby$predictor[i])
  idx_baby = which(newdata_baby$baby == cntr_baby$baby[i])
  d = fit_age_baby[, idx_baby[2], idx_pred] - fit_age_baby[, idx_baby[1], idx_pred]
  cntr_baby[i, c('fit', 'lwr', 'upr')] = quantile(d, probs = c(.5, .025, .975))
}
cntr_baby$predictor = factor(cntr_baby$predictor, levels = rev(sort(unique(cntr_baby$predictor))))

ggplot(cntr, aes(x = predictor, y = fit, ymin = lwr, ymax = upr, color = outside_rope)) +  # , color = outside_rope
  geom_violin(data = cntr_baby, aes(x = predictor, y = fit), inherit.aes = FALSE, color = NA, fill = 'gray50', alpha = .25) +
  geom_rect(data = df_rope, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = NA, fill = 'gray50', alpha = .25, inherit.aes = FALSE) +
  geom_pointrange(position = position_dodge(.3)) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_color_manual(name = NULL, values = c('gray70', 'black'), guide = FALSE) +
  xlab('') +
  ylab('Change with age, SD/month') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave('../pix/ac~age.png', width = 20, height = 10, units = 'cm', dpi = 600)



## Random Forest for classifying into boys/girls
df$age_month = as.factor(df$age_month)

nIter = 100
out = data.frame(iter = 1:nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  idx_train = numeric(0)
  for (b in unique(df$age_month)) {
    idx_age_month = which(df$age_month == b)
    bouts_age_month = as.character(unique(df$file_seq_S[idx_age_month]))
    bouts_train = sample(bouts_age_month, size = length(bouts_age_month) * 2/3)
    idx_train = c(idx_train, which(df$file_seq_S %in% bouts_train))
  }
  train = df[idx_train, ]
  test = df[-idx_train, ]
  out$nAges[iter] = length(levels(train$age_month))
  min_per_age_month = min(table(train$age_month))
  mod_rf_age_month = randomForest(x = train[, vars_to_model], y = train$age_month, strata = train$age_month, sampsize = rep(min_per_age_month, length(levels(train$age_month))))
  # varImpPlot(mod_rf_age_month)
  pr = as.character(predict(mod_rf_age_month, newdata = test[, vars_to_model], type = 'response'))
  out$acc[iter] = mean(pr == test$age_month) 
  soundgen:::reportTime(iter, time_start, nIter)
}
out$OR = odds(out$acc) / odds(1 / out$nAges)
# saveRDS(out, '../data/ageRec.RDS')

out = readRDS('../data/ageRec.RDS')
report(out$acc * 100, 0)  # 40 [35, 43]
report(out$OR)  # 2 [1.6, 2.2]
