# Random Forest models for predicting the cause of crying from acoustics

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(randomForest)
library(soundgen)
library(corrplot)
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


## Global
reportRF = function(predictors, nIter = 100, propTrain = .66) {
  perIter = data.frame(iter = 1:nIter)
  time_start = proc.time()
  for (i in 1:nIter) {
    sessions = data.frame(file_seq_S = unique(df$file_seq_S))
    sessions$cause_stop_engl = df$cause_stop_engl[match(sessions$file_seq_S, df$file_seq_S)]
    # stratified sampling to make sure all causes are represented in the training and testing sample, while training and testing on different recording sessions
    idx = unlist(tapply(1:nrow(sessions), sessions$cause_stop_engl, function(x) sample(x, size = propTrain * length(x))))
    # table(sessions$cause_stop_engl[idx])
    # table(sessions$cause_stop_engl[-idx])
    train = df[df$file_seq_S %in% sessions$file_seq_S[idx], ]
    test = df[!df$file_seq_S %in% sessions$file_seq_S[idx], ]
    mod_rf_bouts = randomForest(x = train[, predictors], y = train$cause_stop_engl, strata = train$cause_stop_engl, sampsize = rep(sum(train$cause_stop_engl == levels(train$cause_stop_engl)[which.min(table(train$cause_stop_engl))]), length(levels(train$cause_stop_engl))))
    perIter$acc_train[i] = 1 - mean(mod_rf_bouts$confusion[, ncol(mod_rf_bouts$confusion)]) 
    test$pr = predict(mod_rf_bouts, newdata = test)
    perIter$acc_test[i] = mean(test$pr == test$cause_stop_engl) 
    perIter[i, c('acc_test_discomfort', 'acc_test_hunger', 'acc_test_loneliness')] = 
      aggregate((pr == cause_stop_engl) ~ cause_stop_engl, test, mean)[, 2]
    soundgen:::reportTime(i, time_start, nIter)
  }
  perIter$OR_train = odds(perIter$acc_train) / odds(1 / length(unique(df$cause_stop_engl))) 
  perIter$OR_test = odds(perIter$acc_test) / odds(1 / length(unique(df$cause_stop_engl))) 
  out = list(
    perIter = perIter,
    summary = list(
      acc_train = report(perIter$acc_train * 100, 0),  
      OR_train = report(perIter$OR_train), 
      acc_test = report(perIter$acc_test * 100, 0), 
      OR_test = report(perIter$OR_test),
      discomfort = report(perIter$acc_test_discomfort * 100),
      hunger = report(perIter$acc_test_hunger * 100),
      loneliness = report(perIter$acc_test_loneliness * 100)
    )
  )
  return(out)
}

global = reportRF(vars_to_model, nIter = 100)
# saveRDS(global, '../data/RF_global.RDS')
global = readRDS('../data/RF_global.RDS')

global$summary
# $acc_train
# [1] "45 [44, 46]"
# 
# $OR_train
# [1] "1.6 [1.5, 1.7]"
# 
# $acc_test
# [1] "36 [33, 38]"
# 
# $OR_test
# [1] "1.1 [1, 1.2]"
# 
# $discomfort
# [1] "29.8 [25.2, 34.6]"
# 
# $hunger
# [1] "33.5 [26.2, 40]"
# 
# $loneliness
# [1] "42.3 [37, 48.1]"


# +age group as a predictor
# withAge = reportRF(c('age_month', vars_to_model), nIter = 100)
# saveRDS(withAge, '../data/RF_withAge.RDS')
withAge = readRDS('../data/RF_withAge.RDS')
withAge$summary
# $acc_train
# [1] "52 [50, 55]"
# 
# $OR_train
# [1] "2.2 [2, 2.4]"
# 
# $acc_test
# [1] "38 [33, 41]"
# 
# $OR_test
# [1] "1.2 [1, 1.4]"

# +baby as a predictor
# withBaby = reportRF(c('baby', vars_to_model), nIter = 100)
# saveRDS(withBaby, '../data/RF_withBaby.RDS')
withBaby = readRDS('../data/RF_withBaby.RDS')
withBaby$summary
# $acc_train
# [1] "66 [64, 69]"
# 
# $OR_train
# [1] "3.9 [3.5, 4.4]"
# 
# $acc_test
# [1] "43 [36, 50]"
# 
# $OR_test
# [1] "1.5 [1.1, 2]"

# +baby and age
# withAgeBaby = reportRF(c('baby', 'age_month', vars_to_model), nIter = 100)
# saveRDS(withAgeBaby, '../data/RF_withAgeBaby.RDS')
withAgeBaby = readRDS('../data/RF_withAgeBaby.RDS')
withAgeBaby$summary
# $acc_train
# [1] "75 [73, 77]"
# 
# $OR_train
# [1] "5.9 [5.4, 6.7]"
# 
# $acc_test
# [1] "43 [35, 48]"
# 
# $OR_test
# [1] "1.5 [1.1, 1.9]"


## One baby at a time
babies = data.frame(baby = sort(unique(df$baby)), acc = NA, acc_otherBabies = NA, nLevels = NA)
for (i in 1:nrow(babies)) {
  baby_i = droplevels(df[which(df$baby == babies$baby[i]), ])
  babies$nLevels[i] = length(unique(baby_i$cause_stop_engl))
  if (babies$nLevels[i] == 3) {
    # train on either just "vars_to_model" or on "c(vars_to_model, 'age_month')"
    mod_rf_baby_i = randomForest(x = baby_i[, c(vars_to_model, 'age_month')], y = baby_i$cause_stop_engl, strata = baby_i$cause_stop_engl, sampsize = rep(sum(baby_i$cause_stop_engl == levels(baby_i$cause_stop_engl)[which.min(table(baby_i$cause_stop_engl))]), length(levels(baby_i$cause_stop_engl))))
    babies$acc[i] = 1 - mean(mod_rf_baby_i$confusion[, ncol(mod_rf_baby_i$confusion)])
    temp_not_i = droplevels(df[which(df$baby != babies$baby[i]), ])
    babies$acc_otherBabies[i] = mean(predict(mod_rf_baby_i, newdata = temp_not_i) == temp_not_i$cause_stop_engl)
  }
}
babies$OR = odds(babies$acc) / odds(1 / babies$nLevels)
babies$OR_otherBabies = odds(babies$acc_otherBabies) / odds(1 / babies$nLevels)
mean(babies$acc, na.rm = TRUE)  # 62%
mean(babies$OR, na.rm = TRUE)   # OR = 3.5
mean(babies$acc_otherBabies, na.rm = TRUE)  # 34%
mean(babies$OR_otherBabies, na.rm = TRUE)   # OR = 1.0
babies$nObs = aggregate(get(vars_to_model[1]) ~ baby, df, function(x) length(unique(x)))[, 2]
plot(acc ~ nObs, babies)  # no obvious patterns
cor(babies$acc, babies$nObs)
babies$n_obs = aggregate(cause_stop_engl ~ baby, df, length)[, 2]
bp = babies[, c('baby', 'sex', 'n_obs')]
bp$same = paste0(round(babies$acc * 100, 0), ' (', round(babies$OR, 1), ')')
bp$dif = paste0(round(babies$acc_otherBabies * 100, 0), ' (', round(babies$OR_otherBabies, 1), ')')
head(bp)
# write.csv(bp, '../data/RF_per_baby.csv', row.names = FALSE)
bp_pl = data.frame(
  Train = rep(babies$baby, 2),
  Test = rep(c('Same', 'Different'), each = nrow(babies)),
  Accuracy = 100 * c(babies$acc, babies$acc_otherBabies),
  OR = c(babies$OR, babies$OR_otherBabies)
)
ggplot(bp_pl, aes(Test, Accuracy)) +
  geom_violin()


babies$sex = df$sex[match(babies$baby, df$baby)]
aggregate(acc ~ sex, babies, mean)  # no difference between boys and girls in terms of model accuracy 
#   sex       acc
# 1   M 0.6698051
# 2   W 0.6543060


## One baby at a time, with different cry bouts in train/test
nIter = 100
babies_list = vector('list', nIter)
time_start = proc.time()
for (iter in 1:nIter) {
  babies = data.frame(baby = sort(unique(df$baby)), acc_train = NA, acc_test = NA, 
                      acc_otherBabies = NA, nLevels = NA)
  for (i in 1:nrow(babies)) {
    baby_i = droplevels(df[which(df$baby == babies$baby[i]), ])
    bouts_baby_i = as.character(unique(baby_i$file_seq_S))
    bouts_train = sample(bouts_baby_i, size = length(bouts_baby_i) * 2/3)
    train_baby_i = droplevels(baby_i[baby_i$file_seq_S %in% bouts_train, ])
    test_baby_i = try(baby_i[!baby_i$file_seq_S %in% bouts_train & 
                               baby_i$cause_stop_engl %in% unique(train_baby_i$cause_stop_engl), ])
    if (nrow(test_baby_i) < 1) next
    # NB: we need the same factor levels in train, test, and not_i (2 or 3 causes of crying)
    babies$nLevels[i] = length(unique(train_baby_i$cause_stop_engl))
    if (babies$nLevels[i] > 1) {
      # train on either just "vars_to_model" or on "c(vars_to_model, 'age_month')"
      mod_rf_baby_i = randomForest(x = train_baby_i[, vars_to_model], y = train_baby_i$cause_stop_engl, strata = train_baby_i$cause_stop_engl, sampsize = rep(sum(train_baby_i$cause_stop_engl == levels(train_baby_i$cause_stop_engl)[which.min(table(train_baby_i$cause_stop_engl))]), length(levels(train_baby_i$cause_stop_engl))))
      babies$acc_train[i] = 1 - mean(mod_rf_baby_i$confusion[, ncol(mod_rf_baby_i$confusion)])
      babies$acc_test[i] = mean(predict(mod_rf_baby_i, newdata = test_baby_i) == as.character(test_baby_i$cause_stop_engl))
      temp_not_i = droplevels(df[which(
        df$baby != babies$baby[i] &
          df$cause_stop_engl %in% unique(train_baby_i$cause_stop_engl)), ])
      babies$acc_otherBabies[i] = mean(predict(mod_rf_baby_i, newdata = temp_not_i) == temp_not_i$cause_stop_engl)
    }
  }
  babies$OR_train = odds(babies$acc_train) / odds(1 / babies$nLevels)
  babies$OR_test = odds(babies$acc_test) / odds(1 / babies$nLevels)
  babies$OR_otherBabies = odds(babies$acc_otherBabies) / odds(1 / babies$nLevels)
  babies_list[[iter]] = babies
  soundgen:::reportTime(iter, time_start, nIter)
}
# saveRDS(babies_list, '../data/RF_cause_per_baby.RDS')

babies_list = readRDS('../data/RF_cause_per_baby.RDS')
bm = Reduce("+", babies_list) / length(babies_list)
bm$baby = babies_list[[1]]$baby

sumtab = babies_list[[1]]
for (c in 2:ncol(sumtab)) {
  for (r in 1:nrow(sumtab)) {
    if (c %in% c(2:4)) {
      # accuracy 
      sumtab[r, c] = report(unlist(lapply(babies_list, function(x) x[r, c])) * 100, 0)
    } else {
      # OR & nLevels
      sumtab[r, c] = report(unlist(lapply(babies_list, function(x) x[r, c])), 1)
    }
  }
}

pl = rbind(
  data.frame(test = 'Same', baby = bm$baby, acc = bm$acc_test, OR = bm$OR_test),
  data.frame(test = 'Other', baby = bm$baby, acc = bm$acc_otherBabies, OR = bm$OR_otherBabies)
)
pl$test = factor(pl$test, levels = c('Same', 'Other'))
ag = aggregate(file ~ baby, df, function(x) length(unique(x)))
pl$n_cries = ag$file[match(pl$baby, ag$baby)]
ags = aggregate(file_seq_S ~ baby, df, function(x) length(unique(x)))
pl$n_sessions = ags$file[match(pl$baby, ag$baby)]
ggplot(pl, aes(test, OR, size = n_sessions)) +  # or size = n_cries
  geom_violin(fill = 'gray75', alpha = .5) +
  geom_point() +
  scale_size_continuous(breaks = (1:4)*10, range = c(.5, 3), trans = scales::trans_new('square', function(x)x^2, function(x)sqrt(x))) +
  xlab('Test babies') +
  ylab('Odds Ratio to chance') +
  geom_hline(yintercept = 1, linetype = 3) +
  theme_bw() +
  theme(panel.grid = element_blank())
# ggsave('../pix/RF_cause_by_baby_violin.png', width = 15, height = 10, units = 'cm', dpi = 600)
# write.csv(pl, '../data/fig3c_violins.csv', row.names = FALSE)

plot(OR ~ n_sessions, pl)
cor(pl$OR, pl$n_sessions, use = 'complete.obs') # r = .1
plot(OR ~ n_cries, pl)
cor(pl$OR, pl$n_cries, use = 'complete.obs') # r = .01

all_acc = unlist(lapply(babies_list, function(x) x$acc_test))
report(all_acc * 100, 0)  # 38% [17, 65]
all_OR = unlist(lapply(babies_list, function(x) x$OR_test))
report(all_OR, 1)  # 1.2 [0.4, 3.2]
all_acc_other = unlist(lapply(babies_list, function(x) x$acc_other))
report(all_acc_other * 100, 0)  # 34% [30, 54]
all_OR_other = unlist(lapply(babies_list, function(x) x$OR_other))
report(all_OR_other)  # 1.0 [0.9, 1.2]

babies$nObs = aggregate(get(vars_to_model[1]) ~ baby, df, function(x) length(unique(x)))[, 2]
plot(acc ~ nObs, babies)  # no obvious patterns
cor(babies$acc, babies$nObs)
babies$n_obs = aggregate(cause_stop_engl ~ baby, df, length)[, 2]
bp = babies[, c('baby', 'sex', 'n_obs')]
bp$same = paste0(round(babies$acc * 100, 0), ' (', round(babies$OR, 1), ')')
bp$dif = paste0(round(babies$acc_otherBabies * 100, 0), ' (', round(babies$OR_otherBabies, 1), ')')
head(bp)
# write.csv(bp, '../data/RF_per_baby.csv', row.names = FALSE)
bp_pl = data.frame(
  Train = rep(babies$baby, 2),
  Test = rep(c('Same', 'Different'), each = nrow(babies)),
  Accuracy = 100 * c(babies$acc, babies$acc_otherBabies),
  OR = c(babies$OR, babies$OR_otherBabies)
)
ggplot(bp_pl, aes(Test, Accuracy)) +
  geom_violin()



## try for each age group separately + train on one age group, test on all other ages,  with all possible age group combinations
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
      temp = df[df$age_month == out$train[i], ]
      sessions = data.frame(file_seq_S = unique(temp$file_seq_S))
      sessions$cause_stop_engl = temp$cause_stop_engl[match(sessions$file_seq_S, temp$file_seq_S)]
      # stratified sampling to make sure all causes are represented in the training and testing sample, while training and testing on different recording sessions
      idx = unlist(tapply(1:nrow(sessions), sessions$cause_stop_engl, function(x) sample(x, size = .66 * length(x))))  
      # table(sessions$cause_stop_engl[idx])
      # table(sessions$cause_stop_engl[-idx])
      train = temp[temp$file_seq_S %in% sessions$file_seq_S[idx], ]
      test = temp[!temp$file_seq_S %in% sessions$file_seq_S[idx], ]
    } else {
      # different groups - non-overlapping recording sessions
      train = df[df$age_month == out$train[i], ]
      test = df[df$age_month == out$test[i], ]
    }
    mod_rf_ab = randomForest(x = train[, vars_to_model], y = train$cause_stop_engl, strata = train$cause_stop_engl, sampsize = rep(sum(train$cause_stop_engl == levels(train$cause_stop_engl)[which.min(table(train$cause_stop_engl))]), length(levels(train$cause_stop_engl))))
    pr = as.character(predict(mod_rf_ab, newdata = test[, vars_to_model], type = 'response'))
    out$acc[i] = mean(pr == test$cause_stop_engl)  # prop correct
    out$OR[i] = odds(out$acc[i]) / odds(1 / length(unique(temp$cause_stop_engl)))  # OR
  }
  outAge[[iter]] = out
  reportTime(iter, time_start, nIter)
}
saveRDS(outAge, '../data/RF_cause_per_ageGroup.RDS')

outAge = readRDS('../data/RF_cause_per_ageGroup.RDS')

out = Reduce("+", outAge) / length(outAge)
out$age = outAge[[1]]$age_month

out$acc = round(out$acc, 2)
out$OR = round(out$OR, 1)
out
m = reshape2::acast(out, test~train, value.var = 'OR')
m  # columns = train, rows = test
png('../pix/corrplot_cause~ageGr.png', width = 10, height = 10, units = 'cm', res = 600)
corrplot(m, method = "color", addCoef.col = "black", tl.srt = 0, is.corr = FALSE) 
dev.off()

write.csv(m, '../data/fig3b_RF-matrix.csv', row.names = FALSE)
