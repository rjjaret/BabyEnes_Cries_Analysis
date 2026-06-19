# Script for extracting individual cries (syllables) from long recordings
# AND for analyzing them acoustically

# NB: use cause_stop for analyses instead of cause_parents

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the necessary libraries
library(soundgen)

# Prepare meta-info
info = read.csv('../data/info_676_modifyV3.csv', stringsAsFactors = FALSE)
table(info$cause_parents)
table(info$cause_stop)
table(info$cause_parents, info$cause_stop)
mean(info$cause_parents == info$cause_stop)
# double categories - just use the first letter
info$cause_parents_first = substr(info$cause_parents, 1, 1)
info$cause_stop_first = substr(info$cause_stop, 1, 1)
for (i in 1:nrow(info)) {
  info$cause_parents_engl[i] = switch(
    info$cause_parents_first[i],
    'D' = 'pain',
    'F' = 'hunger',
    'I' = 'discomfort',
    'S' = 'loneliness',
    'J' = 'DK',
    'N' = 'DK',
    'other'  # residual category - everything else goes in here
  )
  info$cause_stop_engl[i] = switch(
    info$cause_stop_first[i],
    'D' = 'pain',
    'F' = 'hunger',
    'I' = 'discomfort',
    'S' = 'loneliness',
    'J' = 'DK',
    'N' = 'DK',
    'other'  # residual category - everything else goes in here
  )
}
info$cause_stop_engl[info$cause == 'SOIR'] = 'other'
table(info$cause_parents_engl)
table(info$cause_stop_engl)
write.csv(info, '../data/info_676_edited.csv', row.names = FALSE)

## Prepare the audio
# Relative path to the folder with your files (check syntax, eg \ on Windows PCs)
master_path = '../audio/Longcry_Database_quality_concatenation/'
all_folders = list.files(master_path)
newpath = paste0(master_path, '/00_pooled_preprocessed/')
time_start = proc.time()
for (folder_n in 2:length(all_folders)) {
  mypath = all_folders[folder_n]
  # Preprocessing
  files = list.files(paste0(master_path, mypath), pattern = '.wav')
  for (f in files) {
    # Open the file
    file_f = try(tuneR::readWave(paste0(master_path, mypath, '/', f)))
    if (class(file_f) == 'try-error') next
    
    # Downsample to 22050 Hz to speed up processing and remove high-frequency noise
    file_f_downs = seewave::resamp(file_f, g = 22050, output = 'Wave')
    
    # High-pass above 100 Hz to remove low-frequency noise
    file_f_filt = seewave::ffilter(file_f_downs, from = 100, output = 'Wave')
    
    # Check the spectrum and save
    # seewave::meanspec(file_f_filt, f = 22050)
    # spectrogram(file_f_filt[,1], samplingRate = 22050, ylim = c(0, 2))
    seewave::savewav(file_f_filt, filename = paste0(newpath, f))
  }
  reportTime(i = folder_n, nIter = length(all_folders), time_start = time_start)
}


# Normalize for peak ampl, overwriting the files (100 at a time to avoid running out of RAM)
p = paste0(newpath, '/07')
normalizeFolder(p, savepath = p, type = 'peak')
# if needed, normalize for max or mean RMS:
# normalizeFolder(newpath, savepath = newpath, type = 'rms', summaryFun = 'max')
# getRMSFolder(newpath, summaryFun = c('mean', 'max'))

# try to analyze a few in pitch_app to find reasonable settings
# pitch_app()


## Segment into syllables and save each as a separate audio file
# preview - use to fine-tune the settings
seg = segment(
  '../audio/Longcry_Database_quality_concatenation/00_pooled_preprocessed/BM29_M_2_F_F_10102019_2058.wav',
  from = 10, to = 30,  # 10-30 s, just as an example
  shortestSyl = 50,  # min accepted duration of syllable
  shortestPause = 100,  # min accepted pause (if smaller, merge)
  method = 'mel',  # segment based on mel-spectrogram
  propNoise = NULL,  # expected proportion of noise - detect automatically
  SNR = 3,           # expeted signal-to-noise ratio
  windowLength = 40,
  step = 5,
  interburst = 200,
  summaryFun = NULL,  # or run w/o it and save summaries as well
  plot = TRUE
)

# now segment all the files and save each syllable as a separate file
seg = segment(
  '../audio/Longcry_Database_quality_concatenation/00_pooled_preprocessed',
  shortestSyl = 50,  # min accepted duration of syllable
  shortestPause = 100,  # min accepted pause (if smaller, merge)
  method = 'mel',  # segment based on mel-spectrogram
  propNoise = NULL,  # expected proportion of noise - detect automatically
  SNR = 3,           # expected signal-to-noise ratio
  windowLength = 40,
  step = 5,
  interburst = 200,
  summaryFun = NULL,  # or run w/o it and save summaries as well
  plot = FALSE,
  saveAudio = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'
)
saveRDS(seg, '../data/syllables&bursts.RDS')
seg_syl = lapply(seg, function(x)
  soundgen:::summarizeAnalyze(x$syllables[, c('sylLen', 'pauseLen')],
                              summaryFun = c('median', 'sd'),
                              var_noSummary = NULL))
seg_syl_sum = do.call('rbind', seg_syl)
write.csv(seg_syl_sum, '../data/syllables_perSeq.csv')

seg_bursts = lapply(seg, function(x)
  soundgen:::summarizeAnalyze(x$bursts[, c('interburst'), drop = FALSE],
                              summaryFun = c('median', 'sd'),
                              var_noSummary = NULL))
seg_bursts_sum = do.call('rbind', seg_syl)
write.csv(seg_bursts_sum, '../data/bursts_perSeq.csv')


# Analyze all preprocessed files
newpath = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'
files = list.files(newpath, pattern = '.wav')
out = vector('list', length(files))
time_start = proc.time()
for (f in 1:length(files)) {
  temp_anal = try(analyze(
    paste0(newpath, '/', files[f]),
    pitchFloor = 100, pitchCeiling = 1750,  # absolute pitch limits
    priorMean = 400, priorSD = 12,          # "soft" pitch priors: 500 Hz ± 1 octave
    nFormants = 0,     # don't analyze formants
    roughness = list(amRes = 30),
    summaryFun = NULL,
    plot = FALSE
  ))
  if (class(temp_anal) == 'try-error') {
    warning(paste('Failed to analyze file', files[f]))
    next
  } else {
    out[[f]] = temp_anal
  }
  reportTime(i = f, nIter = length(files), reportEvery = 100, time_start = time_start)
}  # takes about 11 h with default win and step
names(out) = files

iqr = function(x) IQR(x, na.rm = TRUE)
out = analyze(
  '../audio/Longcry_Database_quality_concatenation/00_pooled_separate', 
  pitchFloor = 100, pitchCeiling = 1750,  # absolute pitch limits
  priorMean = 400, priorSD = 12,          # "soft" pitch priors: 500 Hz ± 1 octave
  nFormants = 0,     # don't analyze formants
  roughness = list(amRes = 30),
  summaryFun = c('median', 'iqr'),
  plot = FALSE, cores = 4
)
saveRDS(out, '../data/babySeparate_soundgen_list.RDS')


## summary table
timeMax = function(x) which.max(x) / length(x)  # without omitting NAs
timeMin = function(x) which.min(x) / length(x)
meanAbsDif = function(x) mean(abs(diff(x)), na.rm = TRUE)
iqr = function(x) IQR(x, na.rm = TRUE)
# summaryFun = c('mean', 'sd', 'median', 'iqr', 'max', 'timeMax', 'min', 'timeMin', 'meanAbsDif')
summaryFun = c('median', 'iqr')
nFiles = length(out)
out_sum_list = vector('list', nFiles)

time_start = proc.time()
for (i in 1:nFiles) {
  temp = soundgen:::summarizeAnalyze(
    out[[i]],
    summaryFun = summaryFun
  )
  cln = colnames(temp)
  temp$file = names(out)[i]
  out_sum_list[[i]] = temp[, c('file', cln)]
  reportTime(i = i, nIter = nFiles, reportEvery = 1000, time_start = time_start)
}

df = do.call('rbind', out_sum_list)
for (i in 1:nrow(df)) {
  temp = df$file[i]
  temp1 = unlist(strsplit(temp, '_'))
  temp2 = paste0(temp1[1:(length(temp1) - 1)], collapse = '_')
  df$file_seq[i] = paste0(temp2, '.wav')
}


# Add meta-info
df = read.csv('../data/babySeparate_soungen_praat_sumtab.csv')
info = read.csv('../data/info_676_edited.csv', stringsAsFactors = FALSE)
table(info$cause_parents_engl)
table(info$cause_stop_engl)
info$file[info$cause_parents_engl == '']
table(info$cause_stop)

addVars = c('baby', 'sex', 'age_month', 'cause_parents', 'cause_stop', 'cause_parents_engl', 'cause_stop_engl')
idx = match(df$file_seq_S, info$file)
any(is.na(idx))
# unique(df$file_seq_S[which(is.na(idx))])

df[, addVars] = info[idx, addVars]
write.csv(df, '../data/dataset_44605_full.csv', row.names = FALSE)


# add segment start/end + dur of original sequence
for (i in 1:nrow(df)) {
  temp = as.character(df$file[i])
  temp1 = unlist(strsplit(temp, '_'))
  temp2 = temp1[length(temp1)]
  temp3 = substr(temp2, 1, nchar(temp2) - 4)
  temp4 = unlist(strsplit(temp3, '-'))
  df$start[i] = as.numeric(temp4[1]) / 1000
  df$end[i] = as.numeric(temp4[2]) / 1000
}
df$seq_dur = info$duree[match(df$file_seq, info$file)]
df$duration_segment = df$end - df$start
df[1:10, c('file', 'start', 'end', 'seq_dur', 'duration_segment', 'duration', 'duration_noSilence')]


# short version of the dataset (only the most relevant columns)
df = read.csv('../data/dataset_44605_full.csv')
colnames(df)
plot(df$pitch_median_S, df$f0_mean_P)
cor(df$pitch_median_S, df$f0_mean_P)  # r = .67

plot(df$pitch_iqr_S, df$f0_stdev_P)
cor(df$pitch_iqr_S, df$f0_stdev_P)  # r = .41

plot(df$HNR_median_S, df$harm_P)
cor(df$HNR_median_S, df$harm_P)  # r = .95

# check a few sounds manually

# short-list the most relevant columns
df = read.csv('../data/dataset_44605_full.csv', stringsAsFactors = TRUE)
colnames(df)
df_short = df[, c('file', 'file_seq_S', 'baby', 'sex', 'age_month', 'cause_parents', 'cause_parents_engl', 'cause_stop', 'cause_stop_engl', 'start_P', 'end_P', 'duration_segment_P', 'duration_noSilence_S', 'voiced_S', 'domVoiced_median_S', 'domVoiced_iqr_S', 'entropyVoiced_median_S', 'entropyVoiced_iqr_S', 'harmEnergy_median_S', 'harmHeight_median_S', 'HNRVoiced_median_S', 'HNRVoiced_iqr_S', 'pitch_median_S', 'pitch_iqr_S', 'quartile25Voiced_median_S', 'quartile50Voiced_median_S', 'quartile75Voiced_median_S', 'roughnessVoiced_median_S', 'roughnessVoiced_iqr_S', 'specCentroidVoiced_median_S', 'specCentroidVoiced_iqr_S', 'specSlopeVoiced_median_S', 'specSlopeVoiced_iqr_S', 'jitter_P', 'shimmer_P')]
write.csv(df_short, '../data/dataset_44605_short.csv', row.names = FALSE)


# add novelty
if (FALSE) {
  nov = data.frame(file = files)
  time_start = proc.time()
  for (i in 1:length(files)) {
    # settings optimized for salience detection (see the salience paper)
    s = ssm(paste0(newpath, '/', files[35]),
            windowLength = 20, step = 5,
            ssmWin = 40, input = 'audiogram',
            kernelLen = 200, plot = T)
    nov$novelty_mean[i] = mean(s$novelty, na.rm = TRUE)
    reportTime(i = i, nIter = nrow(out), time_start = time_start)
  }
  df$novelty_mean = nov$novelty_mean
  # probably doesn't make sense b/c these sounds are too short
}

# write.csv(df, '../data/babySeparate_soundgen_sumTab.csv', row.names = FALSE)


## save modulation spectra
dim_ms = c(32, 32)
newpath = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'
files = list.files(newpath, pattern = '.wav')
out_ms = vector('list', length(files))
names(out_ms) = files
time_start = proc.time()
for (f in 1:length(files)) {
  temp_ms = try(modulationSpectrum(
    paste0(newpath, '/', files[f]),
    amRes = NULL,
    maxDur = 15,
    plot = FALSE
  ))
  if (class(temp_ms) == 'try-error') {
    warning(paste('Failed to analyze file', files[f]))
    next
  } else {
    out_ms[[f]] =  as.numeric(soundgen:::interpolMatrix(
      temp_ms$original,
      nr = dim_ms[1], nc = dim_ms[2],
      interpol = 'approx'))
  }
  reportTime(i = f, nIter = length(files), time_start = time_start)
}  # takes about 1 h
# saveRDS(out_ms, '../data/babySeparate_modSpec.RDS')
# out_ms = readRDS('../data/babySeparate_modSpec.RDS')
ms = as.data.frame(do.call('rbind', out_ms))
ms = round(ms, 6)
ms$file = files
ms = ms[, c('file', colnames(ms)[2:(ncol(ms) - 1)])]
write.csv(ms, '../data/babySeparate_modSpec.csv', row.names = FALSE)


## save mel-spectrograms
dim_melSpec = c(32, 32)
newpath = '../audio/Longcry_Database_quality_concatenation/00_pooled_separate'
files = list.files(newpath, pattern = '.wav')
df = read.csv('../data/babySeparate_soundgen_sumTab.csv', stringsAsFactors = TRUE)
durs = df$duration
log_durs = log(df$duration * 1000)
nSamples = round(log_durs / max(log_durs) * dim_melSpec[2])  # range(nSamples)
out_melSpec = vector('list', length(files))
names(out_melSpec) = files
time_start = proc.time()
for (f in 1:length(files)) {  # 1:length(files)
  sound_wave = tuneR::readWave(paste0(newpath, '/', files[f]))
  temp_mel = try(tuneR::melfcc(
    sound_wave, wintime = .04, hoptime = .01, preemph = 0,
    lifterexp = 0, nbands = dim_melSpec[1], spec_out = TRUE
  ))
  if (class(temp_mel) == 'try-error') {
    warning(paste('Failed to analyze file', files[f]))
    next
  } else {
    temp_as = t(temp_mel$aspectrum)
    out_melSpec[[f]] = as.numeric(cbind(
      temp_as[, seq(1, ncol(temp_as), length.out = nSamples[f])],
      matrix(0, nrow = dim_melSpec[1], ncol = dim_melSpec[2] - nSamples[f])
    ))
  }
  reportTime(i = f, nIter = length(files), time_start = time_start)
}  # takes about 10 min
# saveRDS(out_melSpec, '../data/babySeparate_melSpec.RDS')

range(lapply(out_melSpec[1:5], length))

as = as.data.frame(do.call('rbind', out_melSpec))
as$file = files
as = as[, c('file', colnames(as)[2:(ncol(as) - 1)])]
# write.csv(as, '../data/babySeparate_melSpec.csv', row.names = FALSE)

