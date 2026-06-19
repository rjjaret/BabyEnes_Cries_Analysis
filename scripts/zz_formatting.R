
# a table of formatted names of acoustic variables
acVarsFmt = list(
  'duration_segment_P' = 'Duration',
  'entropyVoiced_median_S' = 'Entropy',
  'HNRVoiced_median_S' = 'Harmonicity',
  'jitter_P' = 'Jitter',
  'pitch_iqr_S' = 'Pitch (IQR)',
  'pitch_median_S' = 'Pitch (median)',
  'roughnessVoiced_median_S' = 'Roughness',
  'shimmer_P' = 'Shimmer',
  'specCentroidVoiced_median_S' = 'Spectral centroid',
  'voiced_S' = 'Voicing'
)

# formatted names of contexts
contextNames = list(
  'discomfort' = 'Discomfort', 
  'hunger' = 'Hunger', 
  'loneliness' = 'Isolation'
)

odds = function(x) x / (1 - x)

report = function(x, digits = 1) soundgen:::reportCI(
  quantile(x, probs = c(.5, .025, .975), na.rm = TRUE), 
  digits = digits)
