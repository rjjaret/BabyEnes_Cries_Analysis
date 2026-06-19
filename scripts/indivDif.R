## Individual recognition with Random Forest models

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(randomForest)
library(brms)
library(ggplot2)
library(corrplot)
library(soundgen)
source('zz_formatting.R')

df = read.csv('../data/dataset_44605_short.csv', stringsAsFactors = TRUE)
audioFolder = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'
odds = function(x) x / (1 - x)

vars_to_model = sort(c('duration_segment_P', 'voiced_S', 'pitch_median_S', 'pitch_iqr_S', 'specCentroidVoiced_median_S', 'entropyVoiced_median_S', 'HNRVoiced_median_S', 'roughnessVoiced_median_S', 'jitter_P', 'shimmer_P'))
# vars_to_model[which(!vars_to_model %in% colnames(df))]
# summary(df[, vars_to_model])
vars_for_log = c('duration_segment_P', 'domVoiced_median_S', 'domVoiced_iqr_S', 'pitch_median_S', 'pitch_iqr_S', 'quartile25Voiced_median_S', 'quartile50Voiced_median_S', 'quartile75Voiced_median_S', 'specCentroidVoiced_median_S', 'specCentroidVoiced_iqr_S')

# scale
summary(df[, vars_to_model])
for (v in vars_to_model) {
  if (v %in% vars_for_log) df[, v] = log(df[, v] + 1e-3)
  # if (!v %in% c('duration', 'duration_noSilence', 'voiced'))
  df[, v] = as.numeric(scale(df[, v]))
} 
summary(df[, vars_to_model])


## Hs 
library(IDmeasurer)
pc = prcomp(df[, vars_to_model])
df_pc = data.frame(baby = df$baby)
df_pc = cbind(df_pc, pc$x)
head(df_pc)
IDmeasurer::calcHS(df_pc)  # 0.63; 2^0.63 = 1.55
# same as: calcHS(calcPCA(df[, c('baby', vars_to_model)]), sumHS=T)
# doesn't really depend on the number of obs: calcHS(calcPCA(df[sample(1:nrow(df), size = 1500), c('baby', vars_to_model)]), sumHS=T)

# as a function of dur
qn = unique(quantile(df$duration_noSilence_S, probs = seq(0, 1, length.out = 26)))  # every 4%
length(unique(qn))
hs = data.frame(duration = qn, n = NA, Hs = NA)
for (i in 1:nrow(hs)) {
  df_temp = df[df$duration_noSilence_S > hs$duration[i] & df$duration_noSilence_S < hs$duration[i + 1], c('baby', vars_to_model)]
  hs$n[i] = nrow(df_temp)
  if (hs$n[i] > (24 * 20)) {
    temp = try(calcHS(calcPCA(df_temp), sumHS=T)[2])
    if (class(temp) != 'try-error') hs$Hs[i] = temp
  }
}
# hs = na.omit(hs)
plot(hs, log = 'x')
library(cowplot)  # for adding customized marginal distributions: https://www.lreding.com/nonstandard_deviations/2017/08/19/cowmarg/
pl_hs = ggplot(hs, aes(duration, 2 ^ Hs)) +
  geom_point() +
  geom_smooth() + # method = 'lm'
  scale_x_continuous('Duration (s)', trans = 'log2', n.breaks = 7, limits = c(.2, 3)) + # limits = range(df$duration_noSilence_S)
  theme_bw() +
  theme(panel.grid = element_blank())
xhist = axis_canvas(pl_hs, axis = "x") +
  geom_density(data = df, aes(duration_noSilence_S)) +
  scale_x_continuous('Duration (s)', trans = 'log2', limits = c(.2, 3))
pl_hs2 = insert_xaxis_grob(pl_hs, xhist, position = "top")
ggdraw(pl_hs2)
# ggsave('../pix/hs_baby_dur.png', width = 10, height = 8, units = 'cm', dpi = 600)
summary(lm(I(2 ^ Hs) ~ duration, hs))  # R2 = .36

# for more variables
calcHS(calcPCA(df[, c('baby', colnames(df)[13:35])]), sumHS=T)  # 1.1

# per age
age = data.frame(age = sort(unique(df$age_month)), Hs = NA)
for (i in 1:nrow(age)) {
  age$Hs[i] = calcHS(calcPCA(df[df$age_month == age$age[i], c('baby', vars_to_model)]), sumHS=T)[2]
}
#   age   Hs
# 1 0.5 1.14
# 2 1.5 0.83
# 3 2.5 0.81
# 4 3.5 1.68



## Individual recognition by Random Forest
min_per_baby = min(table(df$baby))
mod_rf_baby = randomForest(x = df[, vars_to_model], y = df$baby, strata = df$baby, sampsize = rep(min_per_baby, length(levels(df$baby))))
mod_rf_baby  # 36% correct
varImpPlot(mod_rf_baby)

idx_train = sample(1:nrow(df), size = nrow(df) * .66)
train = df[idx_train, ]
test = df[-idx_train, ]

min_per_baby = min(table(train$baby))
mod_rf_baby = randomForest(x = train[, vars_to_model], y = train$baby, strata = train$baby, sampsize = rep(min_per_baby, length(levels(train$baby))))
pr = as.character(predict(mod_rf_baby, newdata = test[, vars_to_model], type = 'response'))
mean_acc = mean(pr == test$baby)  # 34%
odds(mean_acc) / odds(1 / 24)  # ~12

## same, but using different bouts of crying for training and testing
nIter = 100
out = data.frame(iter = 1:nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  idx_train = numeric(0)
  for (b in unique(df$baby)) {
    idx_baby = which(df$baby == b)
    bouts_baby = as.character(unique(df$file_seq_S[idx_baby]))
    bouts_train = sample(bouts_baby, size = length(bouts_baby) * 2/3)
    idx_train = c(idx_train, which(df$file_seq_S %in% bouts_train))
  }
  train = df[idx_train, ]
  test = df[-idx_train, ]
  out$nBabies[iter] = length(levels(train$baby))
  min_per_baby = min(table(train$baby))
  mod_rf_baby = randomForest(x = train[, vars_to_model], y = train$baby, strata = train$baby, sampsize = rep(min_per_baby, length(levels(train$baby))))
  # varImpPlot(mod_rf_baby)
  pr = as.character(predict(mod_rf_baby, newdata = test[, vars_to_model], type = 'response'))
  out$acc[iter] = mean(pr == test$baby) 
  soundgen:::reportTime(iter, time_start, nIter)
}
out$OR = odds(out$acc) / odds(1 / out$nBabies)
# saveRDS(out, '../data/indRec_global.RDS')

out = readRDS('../data/indRec_global.RDS')
report(out$acc * 100, 0)  # 28% [23, 31]
report(out$OR)  # 8.7 [6.8, 10.2]


## within age groups
table(df$age_month)
nIter = 100
outAge = vector('list', nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  acc_age = data.frame(age = sort(unique(df$age_month)), acc = NA, nLevels = NA)
  for (i in 1:nrow(acc_age)) {
    temp = droplevels(df[df$age_month == acc_age$age[i], ])
    idx_train = numeric(0)
    for (b in unique(temp$baby)) {
      idx_baby = which(temp$baby == b)
      bouts_baby = as.character(unique(temp$file_seq_S[idx_baby]))
      bouts_train = sample(bouts_baby, size = length(bouts_baby) * 2/3)
      idx_train = c(idx_train, which(temp$file_seq_S %in% bouts_train))
    }  
    train = droplevels(temp[idx_train, ])
    test = temp[-idx_train, ]
    test = test[test$baby %in% unique(train$baby), ]
    min_per_baby = min(table(train$baby))
    mod_temp = randomForest(x = train[, vars_to_model], y = train$baby, strata = train$baby, sampsize = rep(min_per_baby, length(levels(train$baby))))
    # varImpPlot(mod_temp)
    pr = as.character(predict(mod_temp, newdata = test[, vars_to_model], type = 'response'))
    acc_age$acc[i] = mean(pr == test$baby) 
    acc_age$nLevels[i] = length(unique(temp$baby))
  }
  acc_age$OR = odds(acc_age$acc) / odds(1 / acc_age$nLevels)
  outAge[[iter]] = acc_age
  soundgen:::reportTime(iter, time_start, nIter)
}
saveRDS(outAge, '../data/indRec_perAge.RDS')

outAge = readRDS('../data/indRec_perAge.RDS')
Reduce('+', outAge) / length(outAge)

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
#   age         acc     nLevels              OR
# 1 0.5 38 [28, 45] 17 [17, 17] 9.7 [6.2, 13.1]
# 2 1.5 26 [21, 30] 23 [23, 23]  7.8 [5.8, 9.6]
# 3 2.5 26 [19, 33] 23 [23, 23] 7.7 [5.1, 10.7]
# 4 3.5 32 [21, 44] 12 [12, 12]  5.2 [2.9, 8.6]


# Can a model trained on a baby 0.5 months old recognize the same baby when they are 3.5 months old? Etc, with all possible age group combinations
ag = aggregate(age_month ~ baby, df, function(x) length(unique(x)))
table(ag$age_month)  # 6 babies recorded at all 4 ages
babies_allAges = ag$baby[ag$age_month == 4]

df_allAges = droplevels(df[df$baby %in% babies_allAges, ])
nIter = 100
outAge = vector('list', nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  out = expand.grid(train = c(0.5, 1.5, 2.5, 3.5),
                    test = c(0.5, 1.5, 2.5, 3.5))
  out = out[order(out$train), ]
  out[, c('acc', 'OR')] = NA
  for (i in 1:nrow(out)) {
    if (out$train[i] == out$test[i]) {
      # same age group - need to train and test on different recording sessions
      temp = df_allAges[df_allAges$age_month == out$train[i], ]
      sessions = data.frame(file_seq_S = unique(temp$file_seq_S))
      sessions$baby = temp$baby[match(sessions$file_seq_S, temp$file_seq_S)]
      # stratified sampling to make sure all causes are represented in the training and testing sample, while training and testing on different recording sessions
      idx = unlist(tapply(1:nrow(sessions), sessions$baby, function(x) sample(x, size = .66 * length(x))))
      # table(sessions$baby[idx])
      # table(sessions$baby[-idx])
      train = df[df$file_seq_S %in% sessions$file_seq_S[idx], ]
      test = df[!df$file_seq_S %in% sessions$file_seq_S[idx], ]
    } else {
      # different groups - non-overlapping recording sessions
      train = df_allAges[df_allAges$age_month == out$train[i], ]
      test = df_allAges[df_allAges$age_month == out$test[i], ]
    }
    mod_rf_ab = randomForest(x = train[, vars_to_model], y = train$baby, strata = train$baby, sampsize = rep(sum(train$baby == levels(train$baby)[which.min(table(train$baby))]), length(levels(train$baby))))
    # mod_rf_ab
    pr = as.character(predict(mod_rf_ab, newdata = test[, vars_to_model], type = 'response'))
    out$acc[i] = mean(pr == test$baby)  # prop correct
    out$OR[i] = odds(out$acc[i]) / odds(1 / length(babies_allAges))  # OR
  }
  outAge[[iter]] = out
  reportTime(iter, time_start, nIter)
}
saveRDS(outAge, '../data/RF_indRec_per_ageGroup.RDS')

outAge = readRDS('../data/RF_indRec_per_ageGroup.RDS')

out = Reduce("+", outAge) / length(outAge)
out$age = outAge[[1]]$age_month

out$acc = round(out$acc, 2)
out$OR = round(out$OR, 1)
out
m = reshape2::acast(out, test~train, value.var = 'OR')
m  # columns = train, rows = test
png('../pix/corrplot_baby~ageGr.png', width = 10, height = 10, units = 'cm', res = 600)
corrplot(m, method = "color", addCoef.col = "black", tl.srt = 0, is.corr = FALSE) 
dev.off()
# write.csv(m, '../data/fig2b_RF-matrix.csv', row.names = FALSE)

ac = reshape2::acast(out, test~train, value.var = 'acc')
ac = ac * 100  # to %; columns = train, rows = test
png('../pix/corrplot_baby~ageGr_acc.png', width = 10, height = 10, units = 'cm', res = 600)
corrplot(ac, method = "color", addCoef.col = "black", tl.srt = 0, is.corr = FALSE) 
dev.off()


## Replicate with rejected syllables as control (should be just noise)
df0 = read.csv('../data/old/babySeparate_soundgen_sumTab.csv', stringsAsFactors = TRUE)

vars_to_model0 = sort(c('duration_noSilence', 'voiced', 'pitch_median', 'pitch_iqr', 'specCentroidVoiced_median', 'entropyVoiced_median', 'HNRVoiced_median', 'roughnessVoiced_median'))  # praat vars jitter & shimmer NA for the full sample of 78K sounds
idx_na = which(apply(df0[, vars_to_model0], 1, function(x) any(is.na(x))))
df0 = df0[-idx_na, ]
# vars_to_model0[which(!vars_to_model0 %in% colnames(df0))]
# summary(df0[, vars_to_model0])
vars_for_log0 = c('duration_noSilence', 'pitch_median', 'pitch_iqr', 'specCentroidVoiced_median')

# scale
summary(df0[, vars_to_model0])
for (v in vars_to_model0) {
  if (v %in% vars_for_log0) df0[, v] = log(df0[, v] + 1e-3)
  # if (!v %in% c('duration', 'duration_noSilence', 'voiced'))
  df0[, v] = as.numeric(scale(df0[, v]))
} 
summary(df0[, vars_to_model0])

min_per_baby0 = min(table(df0$baby))
mod_rf_baby0 = randomForest(x = df0[, vars_to_model0], y = df0$baby, strata = df0$baby, sampsize = rep(min_per_baby0, length(levels(df0$baby))))
mod_rf_baby0  # 31% correct for all syllables

df2 = df0[!df0$file %in% unique(df$file), ]
min_per_baby2 = min(table(df2$baby))
mod_rf_baby2 = randomForest(x = df2[, vars_to_model0], y = df2$baby, strata = df2$baby, sampsize = rep(min_per_baby2, length(levels(df2$baby))))
mod_rf_baby2  # 25% correct for rejected syllables only (n = 25K)
odds(.25) / odds(1 / 24)  # 7.7 times better than chance

## same, but using different bouts of crying for training and testing
nIter = 2
out = data.frame(iter = 1:nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  idx_train = numeric(0)
  for (b in unique(df2$baby)) {
    idx_baby = which(df2$baby == b)
    bouts_baby = as.character(unique(df2$file_seq[idx_baby]))
    bouts_train = sample(bouts_baby, size = length(bouts_baby) * 2/3)
    idx_train = c(idx_train, which(df2$file_seq %in% bouts_train))
  }
  train = df2[idx_train, ]
  test = df2[-idx_train, ]
  out$nBabies[iter] = length(levels(train$baby))
  min_per_baby = min(table(train$baby))
  mod_rf_baby = randomForest(x = train[, vars_to_model0], y = train$baby, strata = train$baby, sampsize = rep(min_per_baby, length(levels(train$baby))))
  pr = as.character(predict(mod_rf_baby, newdata = test[, vars_to_model0], type = 'response'))
  out$acc[iter] = mean(pr == test$baby) 
  soundgen:::reportTime(iter, time_start, nIter)
}
out$OR = out$acc / odds(1 / out$nBabies)  # ~4 times better than chance, so presumably noise is not just noise!
