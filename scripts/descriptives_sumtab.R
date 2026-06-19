# Summary of the dataset

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('zz_formatting.R')

df = read.csv('../data/dataset_44605_short.csv', stringsAsFactors = TRUE)
df = droplevels(df[df$cause_stop_engl %in% c('discomfort', 'hunger', 'loneliness'), ])

sumtab = as.data.frame(aggregate(file_seq_S ~ baby, df, function(x) length(unique(x))))
sumtab$baby = as.character(sumtab$baby)
sumtab$sex = as.character(df$sex[match(sumtab$baby, df$baby)])
sumtab$sex[sumtab$sex == 'W'] = 'F'
sumtab$nCries = aggregate(file ~ baby, df, function(x) length(unique(x)))[, 2]
sumtab$nByCause = NA
sessions = data.frame(session = unique(df$file_seq_S))
sessions[, c('cause_stop_engl', 'age_month', 'baby')] = df[match(sessions$session, df$file_seq_S), c('cause_stop_engl', 'age_month', 'baby')]
sessions$age_month = as.factor(session$age_month)  # otherwise drops empty levels
df$age_month = as.factor(df$age_month) 
for (i in 1:nrow(sumtab)) {
  sumtab$nSesByCause[i] = paste(table(sessions$cause_stop_engl[sessions$baby == sumtab$baby[i]]), collapse = ' / ')
  sumtab$nSesByAge[i] = paste(table(sessions$age_month[sessions$baby == sumtab$baby[i]]), collapse = ' / ')
  
  sumtab$nByCause[i] = paste(table(df$cause_stop_engl[df$baby == sumtab$baby[i]]), collapse = ' / ')
  sumtab$nByAge[i] = paste(table(df$age_month[df$baby == sumtab$baby[i]]), collapse = ' / ')
}
sumtab = sumtab[, c('baby', 'sex', 'file_seq_S', 'nSesByAge', 'nSesByCause', 'nCries', 'nByAge', 'nByCause')]
sumtab = rbind(
  sumtab,
  c(
    'Total (24 babies)',
     '',
    length(unique(df$file_seq_S)),
    paste(table(sessions$age_month), collapse = ' / '),
    paste(table(sessions$cause_stop_engl), collapse = ' / '),
    length(unique(df$file)),
    paste(table(df$age_month), collapse = ' / '),
    paste(table(df$cause_stop_engl), collapse = ' / ')
  )
)
colnames(sumtab) = c('Baby', 'Sex', 'Sessions', 'Sessions by age (0.5/1.5/2.5/3.5)', 'Sessions by cause (discomfort / hunger / loneliness)', 'Cries', 'Cries by age (0.5/1.5/2.5/3.5)', 'Cries by cause (discomfort / hunger / loneliness)')
write.csv(sumtab, '../data/summary_nCries.csv', row.names = FALSE)
