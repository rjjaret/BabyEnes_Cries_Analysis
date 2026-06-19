## Model the effects of age on acoustics

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(brms)
library(ggplot2)
library(randomForest)
library(patchwork)
source('zz_formatting.R')

df = read.csv('../data/dataset_44605_short.csv', stringsAsFactors = TRUE)
df = droplevels(df[df$cause_stop_engl %in% c('discomfort', 'hunger', 'loneliness'), ])
orig = df  # copy before scaling
audioFolder = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'

vars_to_model = sort(c('duration_segment_P', 'voiced_S', 'pitch_median_S', 'pitch_iqr_S', 'specCentroidVoiced_median_S', 'entropyVoiced_median_S', 'HNRVoiced_median_S', 'roughnessVoiced_median_S', 'jitter_P', 'shimmer_P'))
# vars_to_model[which(!vars_to_model %in% colnames(df))]
# summary(df[, vars_to_model])
vars_for_log = c('duration_segment_P', 'domVoiced_median_S', 'domVoiced_iqr_S', 'pitch_median_S', 'pitch_iqr_S', 'quartile25Voiced_median_S', 'quartile50Voiced_median_S', 'quartile75Voiced_median_S', 'specCentroidVoiced_median_S', 'specCentroidVoiced_iqr_S')

# scale across all babies
for (v in vars_for_log) df[, v] = log(df[, v] + 1e-3)
for (v in vars_to_model) df[, v] = as.numeric(scale(df[, v]))

# observed
ggplot(df, aes(x = age_month, y = pitch_median_S, linetype = sex)) +
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

library(ggpubr)
ggsave(ggpubr::ggarrange(plot_v('duration_segment_P'), plot_v('entropyVoiced_median_S'), plot_v('HNRVoiced_median_S'), plot_v('jitter_P'), plot_v('pitch_iqr_S'), plot_v('pitch_median_S'), plot_v('roughnessVoiced_median_S'), plot_v('shimmer_P'), plot_v('specCentroidVoiced_median_S'), plot_v('voiced_S'), common.legend = TRUE), filename = '../pix/ac~age*sex_obs.png', width = 40, height = 25, units = 'cm', dpi = 300)

## Model
# Model selection
if (FALSE) {
  # Do we need to account for age, maybe with sex * age interaction? Too slow to test with brms, so here's a hack - take a couple of variables with possible age interaction based on the plots above, test with lme4
  library(lme4)
  mod1 = lmer(specCentroidVoiced_median_S ~ sex + (1|baby), data = df)
  mod2 = lmer(specCentroidVoiced_median_S ~ sex + age_month + (age_month|baby), data = df)
  mod3 = lmer(specCentroidVoiced_median_S ~ sex * age_month + (age_month|baby), data = df)
  anova(mod1, mod2, mod3)  # no interaction with age
  
  mod4 = lmer(specCentroidVoiced_median_S ~ sex + age_month + cause_stop_engl + (1|baby), data = df)
  mod5 = lmer(specCentroidVoiced_median_S ~ sex + age_month + cause_stop_engl + (age_month + cause_stop_engl|baby), data = df)
  anova(mod4, mod5)  # yes, definitely need random slopes... Ah well
  
  # same for pitch_median_S because it's so important
  mod1 = lmer(pitch_median_S ~ sex + (1|baby), data = df)
  mod2 = lmer(pitch_median_S ~ sex + age_month + (age_month|baby), data = df)
  mod3 = lmer(pitch_median_S ~ sex * age_month + (age_month|baby), data = df)
  anova(mod1, mod2, mod3)  # no interaction with age
  
  mod4 = lmer(pitch_median_S ~ sex + age_month + cause_stop_engl + (1|baby), data = df)
  mod5 = lmer(pitch_median_S ~ sex + age_month + cause_stop_engl + (age_month + cause_stop_engl|baby), data = df)
  anova(mod4, mod5)   # same as for specCentroid
}


myformula = formula(paste(
  'mvbind(', paste(vars_to_model, collapse = ', '), 
  ') ~ sex * age_month + (age_month|baby)'
))
myformula

# run model (>48 h)
# mod_sex = brm(myformula, family = 'gaussian', data = df, cores = 2, chains = 2, warmup = 500, iter = 2000)  # , control = list(adapt_delta = 0.95, max_treedepth = 20)  
# saveRDS(mod_sex, '../mod/ac~sex_2022.04.RDS')
mod_sex = readRDS('../mod/ac~sex_2022.04.RDS')

# plot(mod_sex)
# summary(mod_sex)
# conditional_effects(mod_sex, effects = 'sex:age_month')


# plot acoustic changes with age separately for F/M
newdata = expand.grid(age_month = sort(unique(df$age_month)),  # mean(df$age_month),
                      sex = c('M', 'W'))
fit_sex = fitted(mod_sex, newdata = newdata, re_formula = NA, robust = TRUE)
colnames(fit_sex) = c('fit', 'se', 'lwr', 'upr')
df_pl =  expand.grid(age_month = sort(unique(df$age_month)),  # mean(df$age_month),
                     sex = c('M', 'W'),
                     predictor = vars_to_model)
for (v in vars_to_model) {
  idx_df = which(df_pl$predictor == v)
  df_pl[idx_df, c('fit', 'se', 'lwr', 'upr')] = fit_sex[, , which(vars_to_model == v)]
}

ggplot(df_pl, aes(age_month, fit, ymin = lwr, ymax = upr, group = sex, color = sex)) +
  geom_ribbon(color = NA, alpha = .25) +
  geom_line() +
  facet_wrap(~predictor, nrow = 2) +
  scale_color_manual(name = NULL, values = c('gray70', 'black'), guide = FALSE) +
  xlab('Age (months)') +
  ylab('Z-score') +
  theme_bw() +
  theme(panel.grid = element_blank())



# contrasts: acoustic changes over 1 month, averaging M/F
newdata_cntr = expand.grid(age_month = c(0.5, 3.5),  # change over 1 month
                           sex = c('M', 'W'))
fit_sex = fitted(mod_sex, newdata = newdata_cntr, re_formula = NA, summary = FALSE)
cntr = data.frame(predictor = vars_to_model)
for (i in 1:nrow(cntr)) {
  idx_older = which(newdata_cntr$age_month == 3.5)  # both male and female, so we average
  idx_younger = which(newdata_cntr$age_month == 0.5)
  d = rowMeans(fit_sex[, idx_older, i]) - 
    rowMeans(fit_sex[, idx_younger, i])  
  cntr[i, c('fit', 'lwr', 'upr')] = quantile(d, probs = c(.5, .025, .975))
  cntr$pp[i] = mean(d > 0)
}
cntr
cntr$panel = 'Change from 0.5 to 3.5 months'
# write.csv(cntr, '../data/fig1a_ac~age.csv', row.names = FALSE)

# contrasts: differences between F/M, averaging across age groups
newdata_cntr_fm = expand.grid(age_month = unique(df$age_month),  # change over 1 month
                              sex = c('M', 'W'))
fit_fm = fitted(mod_sex, newdata = newdata_cntr_fm, re_formula = NA, summary = FALSE)
cntr_fm = data.frame(predictor = vars_to_model)
for (i in 1:nrow(cntr_fm)) {
  idx_fem = which(newdata_cntr_fm$sex == 'W')
  idx_mal = which(newdata_cntr_fm$sex == 'M')
  d = rowMeans(fit_fm[, idx_fem, i]) - 
    rowMeans(fit_fm[, idx_mal, i])  
  cntr_fm[i, c('fit', 'lwr', 'upr')] = quantile(d, probs = c(.5, .025, .975))
  cntr_fm$pp[i] = mean(d > 0)
}
cntr_fm
cntr_fm$panel = 'Difference "girls - boys"'
# write.csv(cntr_fm, '../data/fig1a_ac~sex.csv', row.names = FALSE)

cntr_pooled = rbind(cntr, cntr_fm)
cntr_pooled$predictor = factor(cntr_pooled$predictor, levels = rev(sort(unique(cntr_pooled$predictor))))
cntr_pooled$cleared_rope = cntr_pooled$upr < -0.1 | cntr_pooled$lwr > 0.1
rope = data.frame(x = 1:2, y = c(-0.1, 0.1))
cntr_pooled$predictor_fmt = as.character(acVarsFmt[as.character(cntr_pooled$predictor)])
cntr_pooled$predictor_fmt = factor(cntr_pooled$predictor_fmt, levels = rev(sort(unique(cntr_pooled$predictor_fmt))))

p1 = ggplot(cntr_pooled, aes(fit, predictor_fmt, xmin = lwr, xmax = upr)) +  # , alpha = cleared_rope
  geom_pointrange() +
  facet_wrap(~panel) +
  # scale_alpha_manual(values = c(.25, 1)) +
  geom_vline(xintercept = 0, linetype = 3) +
  xlab('z-score') +
  ylab('') +
  theme_bw() +
  theme(panel.grid = element_blank())
# ggsave('../pix/ac~sex&age.png', width = 20, height = 10, units = 'cm', dpi = 600)

cntr_pooled$ci = soundgen:::reportCI(cntr_pooled[, c('fit', 'lwr', 'upr')], 2)
cntr_pooled[c('predictor', 'panel', 'ci')]


# interaction plot
cntr_fm_perAge = expand.grid(predictor = vars_to_model,
                             age_month = unique(df$age_month))
for (i in 1:nrow(cntr_fm_perAge)) {
  idx_pred = which(vars_to_model == cntr_fm_perAge$predictor[i])
  idx_fem = which(newdata_cntr_fm$sex == 'W' & 
                    newdata_cntr_fm$age_month == cntr_fm_perAge$age_month[i])
  idx_mal = which(newdata_cntr_fm$sex == 'M' & 
                    newdata_cntr_fm$age_month == cntr_fm_perAge$age_month[i])
  d = fit_fm[, idx_fem, idx_pred] - fit_fm[, idx_mal, idx_pred]
  cntr_fm_perAge[i, c('fit', 'lwr', 'upr')] = quantile(d, probs = c(.5, .025, .975))
  cntr_fm_perAge$pp[i] = mean(d > 0)
}
cntr_fm_perAge$predictor_fmt = as.character(acVarsFmt[as.character(cntr_fm_perAge$predictor)])
cntr_fm_perAge$predictor_fmt = factor(cntr_fm_perAge$predictor_fmt, levels = rev(sort(unique(cntr_fm_perAge$predictor_fmt))))

p2 = ggplot(cntr_fm_perAge, aes(fit, predictor_fmt, xmin = lwr, xmax = upr)) +  # , alpha = cleared_rope
  geom_pointrange() +
  facet_wrap(~age_month, nrow = 1) +
  geom_vline(xintercept = 0, linetype = 3) +
  xlab('Difference "girls - boys", z-score') +
  ylab('') +
  theme_bw() +
  theme(panel.grid = element_blank())
# ggsave('../pix/ac~sex*age.png', width = 20, height = 7, units = 'cm', dpi = 600)
# write.csv(cntr_fm_perAge, '../data/fig1a_ac~sex*age.csv', row.names = FALSE)

p1 / p2 + plot_annotation(tag_levels = 'A')
ggsave('../pix/ac~sex*age.png', width = 20, height = 15, units = 'cm', dpi = 600)


###
# another way to group these plots (July 21, 2022)
cntr_fm$predictor_fmt = as.character(acVarsFmt[as.character(cntr_fm$predictor)])
cntr_fm$predictor_fmt = factor(cntr_fm$predictor_fmt, levels = rev(sort(unique(cntr_fm$predictor_fmt))))
pl1 = ggplot(cntr_fm, aes(fit, predictor_fmt, xmin = lwr, xmax = upr)) +  # , alpha = cleared_rope
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_x_continuous('', n.breaks = 3) +
  ylab('') +
  theme_bw() +
  theme(panel.grid = element_blank())

cntr_fm_perAge$age_month_fmt = paste(cntr_fm_perAge$age_month, 'months')
pl2 = ggplot(cntr_fm_perAge, aes(fit, predictor_fmt, xmin = lwr, xmax = upr)) +  # , alpha = cleared_rope
  geom_pointrange() +
  facet_wrap(~age_month_fmt, nrow = 1) +
  geom_vline(xintercept = 0, linetype = 3) +
  xlab('Difference "girls - boys", z-score') +
  ylab('') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())

pl1 + pl2 + plot_layout(widths = c(1, 3))
ggsave('../pix/ac~sex*age_sex-dif.png', width = 20, height = 6, units = 'cm', dpi = 600)

cntr$predictor_fmt = as.character(acVarsFmt[as.character(cntr$predictor)])
cntr$predictor_fmt = factor(cntr$predictor_fmt, levels = rev(sort(unique(cntr$predictor_fmt))))
ggplot(cntr, aes(fit, predictor_fmt, xmin = lwr, xmax = upr)) +  # , alpha = cleared_rope
  geom_pointrange() +
  geom_vline(xintercept = 0, linetype = 3) +
  xlab('z-score') +
  ylab('') +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave('../pix/ac~sex*age_change-per-age.png', width = 6, height = 8, units = 'cm', dpi = 600)



## Random Forest for classifying into boys/girls
df$sexMF = as.factor(ifelse(df$sex == 'M', 'M', 'F'))
table(df$sexMF)

nIter = 100
out = data.frame(iter = 1:nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  idx_train = numeric(0)
  for (b in unique(df$sexMF)) {
    idx_sexMF = which(df$sexMF == b)
    bouts_sexMF = as.character(unique(df$file_seq_S[idx_sexMF]))
    bouts_train = sample(bouts_sexMF, size = length(bouts_sexMF) * 2/3)
    idx_train = c(idx_train, which(df$file_seq_S %in% bouts_train))
  }
  train = df[idx_train, ]
  test = df[-idx_train, ]
  out$nSexes[iter] = length(levels(train$sexMF))
  min_per_sexMF = min(table(train$sexMF))
  mod_rf_sexMF = randomForest(x = train[, vars_to_model], y = train$sexMF, strata = train$sexMF, sampsize = rep(min_per_sexMF, length(levels(train$sexMF))))
  # varImpPlot(mod_rf_sexMF)
  pr = as.character(predict(mod_rf_sexMF, newdata = test[, vars_to_model], type = 'response'))
  # table(pr, test$sexMF)
  out$acc[iter] = mean(pr == test$sexMF) 
  soundgen:::reportTime(iter, time_start, nIter)
}
out$OR = odds(out$acc) / odds(1 / out$nSexes)
# saveRDS(out, '../data/sexRec_global.RDS')

out = readRDS('../data/sexRec_global.RDS')
report(out$acc * 100, 0)  # 70 [66, 73]
report(out$OR)  # 2.4 [2, 2.8]


## within age groups
table(df$sexMF)
nIter = 100
outAge = vector('list', nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  acc_age = data.frame(age = sort(unique(df$age_month)), acc = NA, nLevels = NA)
  for (i in 1:nrow(acc_age)) {
    temp = droplevels(df[df$age_month == acc_age$age[i], ])
    idx_train = numeric(0)
    for (b in unique(temp$sexMF)) {
      idx_sexMF = which(temp$sexMF == b)
      bouts_sexMF = as.character(unique(temp$file_seq_S[idx_sexMF]))
      bouts_train = sample(bouts_sexMF, size = length(bouts_sexMF) * 2/3)
      idx_train = c(idx_train, which(temp$file_seq_S %in% bouts_train))
    }  
    train = droplevels(temp[idx_train, ])
    test = temp[-idx_train, ]
    test = test[test$sexMF %in% unique(train$sexMF), ]
    min_per_sexMF = min(table(train$sexMF))
    mod_temp = randomForest(x = train[, vars_to_model], y = train$sexMF, strata = train$sexMF, sampsize = rep(min_per_sexMF, length(levels(train$sexMF))))
    # varImpPlot(mod_temp)
    pr = as.character(predict(mod_temp, newdata = test[, vars_to_model], type = 'response'))
    acc_age$acc[i] = mean(pr == test$sexMF) 
    acc_age$nLevels[i] = length(unique(temp$sexMF))
  }
  acc_age$OR = odds(acc_age$acc) / odds(1 / acc_age$nLevels)
  outAge[[iter]] = acc_age
  soundgen:::reportTime(iter, time_start, nIter)
}
# saveRDS(outAge, '../data/sexRec_perAge.RDS')

outAge = readRDS('../data/sexRec_perAge.RDS')
# Reduce('+', outAge) / length(outAge)

sumtab = outAge[[1]]
for (c in 2:4) {
  for (r in 1:4) {
    if (c == 2) {
      # accuracy 
      sumtab[r, c] = report(unlist(lapply(outAge, function(x) x[r, c])) * 100, 0)
    } else {
      # OR & nLevels
      sumtab[r, c] = report(unlist(lapply(outAge, function(x) x[r, c])), 1)
    }
  }
}
sumtab
#   age         acc  nLevels             OR
# 1 0.5 68 [55, 75] 2 [2, 2] 2.1 [1.2, 2.9]
# 2 1.5 70 [64, 74] 2 [2, 2] 2.3 [1.8, 2.9]
# 3 2.5 64 [56, 72] 2 [2, 2] 1.8 [1.3, 2.5]
# 4 3.5 75 [63, 85] 2 [2, 2] 3.1 [1.7, 5.5]


