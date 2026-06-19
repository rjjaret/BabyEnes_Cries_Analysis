## Helper functions for UMAP, especially for computing a distance matrix based on DTW 
# of frame-by-frame features, running in parallel on multiple cores

library(soundgen)
library(ggplot2)
# library(patchwork)  # for arranging multiple ggplots
library(uwot)  # for umap
library(dbscan)  # for hierarchical clustering
library(randomForest)
library(foreach)
library(doParallel)
# This version uses multiple cores. See:
# https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r  
# https://stackoverflow.com/questions/5423760/how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r 

audioFolder = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'
files = list.files(audioFolder, pattern = '.wav')
df = read.csv('../data/dataset_44605_full.csv', stringsAsFactors = TRUE)
df = droplevels(df[df$cause_stop_engl %in% c('discomfort', 'hunger', 'loneliness'), ])

acAn_full = readRDS('../data/old/babySeparate_soundgen_list.RDS')
length(acAn_full)
# df = df[1:200, ]
acAn_full = acAn_full[which(names(acAn_full) %in% df$file)]

### helper functions
# identify points and play the corresponding files by clicking the plot
# (adapted from ex. for identify())
identifyPch = function(x, y = NULL, n = length(x), plot = FALSE, pch = 19, data = df, ...) {
  xy = xy.coords(x, y); x <- xy$x; y <- xy$y
  sel = rep(FALSE, length(x))
  answers = numeric(0)
  while(sum(sel) < n) {
    ans = identify(x[!sel], y[!sel], labels = which(!sel), n = 1, plot = plot, ...)
    if(!length(ans)) break
    ans = which(!sel)[ans]
    answers = c(answers, ans)
    points(x[ans], y[ans], pch = pch)
    ## play the selected point
    pl(ans, data = data)
  }
  ## return indices of selected points
  return(answers)
}

pl = function(idx, data = df) {
  for (i in idx)
    playme(paste0(audioFolder, '/', data[i, 'file']))
}
# pl(1:3)

# Can't deal with all 44K? Reduce the dataset
if (TRUE) {
  # take a random sample
  idx = sample(1:length(acAn_full), size = 10000)
  acAn = acAn_full[idx]
} else {
  # or analyze one baby at a time
  table(df$baby)
  idx_baby = which(df$baby == 'BM29')
  idx = which(names(acAn_full) %in% unique(df$file[idx_baby]))
  acAn = acAn_full[idx]
}


# distance matrix based on frame-by-frame features
add_vars = c('duration_noSilence', 'voiced', 'ampl', 'dom', 'entropy',  'HNR', 'pitch', 'harmEnergy', 'harmHeight', 'loudness', 'peakFreq', 'quartile25', 'quartile50', 'quartile75', 'roughness', 'specCentroid', 'specSlope')
vars_for_log = c('duration_noSilence', 'dom', 'pitch', 'harmHeight', 'peakFreq', 'quartile25', 'quartile50', 'quartile75', 'specCentroid', 'specCentroid')

getDistMat = function(acAn, add_vars, vars_for_log) {
  # normalize
  for (v in add_vars) {
    for (i in 1:length(acAn)) {
      if (v %in% vars_for_log) acAn[[i]][, v] = log(acAn[[i]][, v] + 1e-3)
    }
    all_v = as.numeric(unlist(lapply(acAn, function(x) x[, v])))
    mean_v = mean(all_v, na.rm = TRUE); sd_v = sd(all_v, na.rm = TRUE)
    for (i in 1:length(acAn)) {
      acAn[[i]][, v] = (acAn[[i]][, v] - mean_v) / sd_v
    }
  } 
  # distance matrix
  n = length(acAn) 
  mydist = matrix(NA, nrow = n, ncol = n)
  time_start = proc.time()
  n_done = 0; n_total = n * (n - 1) / 2
  for (i in 1:(n - 1)) {
    for (j in ((i + 1):n)) {
      d = dtw::dtw(
        x = acAn[[i]][, add_vars], 
        y = acAn[[j]][, add_vars], 
        distance.only = TRUE)
      mydist[j, i] = d$normalizedDistance 
      if (is.na(d$normalizedDistance)) print(paste('NA in ij', i, j))
      n_done = n_done + 1
      reportTime(n_done, time_start, n_total)
    }
    # reportTime(i, time_start, n, 1:n)
  }  # ~1 h / 5K, 12 h / 10K
  mydist = as.dist(mydist)
  # any(is.na(mydist) | !is.finite(mydist))
  return(mydist)
  # saveRDS(mydist, '../data/distMatrix_features_10K.RDS')
  # saveRDS(mydist, '../data/distMatrix_features_BM29.RDS')
}
# getDistMat(acAn = acAn_full[1:15], add_vars, vars_for_log)


## same with parallel multi-core processing
getDistMat_multiCore = function(acAn, add_vars, vars_for_log) {
  # normalize
  for (v in add_vars) {
    for (i in 1:length(acAn)) {
      if (v %in% vars_for_log) acAn[[i]][, v] = log(acAn[[i]][, v] + 1e-3)
    }
    all_v = as.numeric(unlist(lapply(acAn, function(x) x[, v])))
    mean_v = mean(all_v, na.rm = TRUE); sd_v = sd(all_v, na.rm = TRUE)
    for (i in 1:length(acAn)) {
      acAn[[i]][, v] = (acAn[[i]][, v] - mean_v) / sd_v
    }
  } 
  # distance matrix
  n = length(acAn) 
  cl = makeCluster(4, outfile = '') # set the number of cores to use
  registerDoParallel(cl)
  time_start = proc.time()
  n_total = n * (n - 1) / 2
  mydist = foreach(i = 1:(n - 1), .combine = rbind) %dopar% {  # .combine = rbind_reportTime()
    out_col = rep(NA, n)
    for (j in ((i + 1):n)) {
      d = dtw::dtw(
        x = acAn[[i]][, add_vars], 
        y = acAn[[j]][, add_vars], 
        distance.only = TRUE)
      out_col[j] = d$normalizedDistance 
      # n_done = n_done + 1
      # soundgen::reportTime(n_done, time_start, n_total)
    }
    # print(soundgen::reportTime(i, time_start, n, 1:n))
    out_col
  }  # ~3 h / 5K, 12 h / 10K
  stopCluster(cl)
  mydist = t(rbind(mydist, rep(NA, n)))  # mydist[1:5, 1:5]
  # dim(mydist)
  mydist = as.dist(mydist)
  # any(is.na(mydist) | !is.finite(mydist))
  return(mydist)
  # saveRDS(mydist, '../data/distMatrix_features_10K.RDS')
  # saveRDS(mydist, '../data/distMatrix_features_BM29.RDS')
}

rbind_reportTime = function() {
  # works, but for some reason not if there is another loop within foreach :(
  time_start = proc.time()
  count = 0
  pb = txtProgressBar(min = 1, max = n - 1, style = 3)
  function(...) {  # m1 = m1, m2 = m2, i = i, n_total = n_total, time_start = time_start
    # browser()
    count <<- count + 1   # use "<<-" instead of "=" to modify a value in parent environment
    soundgen::reportTime(i = count, time_start = time_start, nIter = n - 1, jobs = n:1)
    # count <<- count + length(list(...)) - 1
    # setTxtProgressBar(pb, count)
    flush.console()
    rbind(...)
  }
}
