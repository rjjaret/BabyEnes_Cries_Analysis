# A separate distance matrix for each cause of crying. Input to UMAP is distance matrix based on DTW of frame-by-frame features

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('umap_helperFun.R')
source('zz_formatting.R')

### All babies, but separate UMAP for each cause or age group
nMax = Inf  # max sample size (5K ~= 3.5 h of machine time to get the full distance matrix)

dm_filename = paste0('../data/distMatrix_features_byCause_', nMax / 1000, 'K.RDS')

if (file.exists(dm_filename)) {
  dm = readRDS(dm_filename)
} else {
  # ~24 h on 4 cores
  table(df$cause_stop_engl)
  causes = levels(df$cause_stop_engl)
  dm_by_cause = vector('list', length(causes))
  names(dm_by_cause) = causes
  for (i in 1:length(causes)) {
    idx_i = which(df$cause_stop_engl == causes[i])
    if (length(idx_i) > nMax) idx_i = sample(idx_i, size = nMax)
    dm_by_cause[[i]]$sample = df$file[idx_i]
    dm_by_cause[[i]]$dm = getDistMat_multiCore(acAn = acAn_full[idx_i], add_vars, vars_for_log)
  }
  saveRDS(dm_by_cause, dm_filename)
  dm = dm_by_cause
}

## UMAP
if (file.exists('../data/umap_plot_per_cause.RDS')) {
  df_plot = readRDS('../data/umap_plot_per_cause.RDS')
} else {
  n = length(dm)
  out = vector('list', n)
  names(out) = names(dm)
  # a separate UMAP for each level of our grouping factor (age / cause)
  for (i in 1:length(dm)) {  
    mydist = dm[[i]]$dm
    idx = match(dm[[i]]$sample, df$file)
    
    out[[i]] = as.data.frame(umap(mydist, init = "spectral", n_neighbors = 15))
    addVars = c('file_seq_S', 'baby', 'age_month', 'cause_stop_engl')
    out[[i]][, addVars] = df[idx, addVars]
    
    if (FALSE) {
      plot(out[[i]]$V1, out[[i]]$V2, cex = .5)
      ip = identifyPch(out[[i]]$V1, out[[i]]$V2, data = df)  # if it works, clicking plays the corresponding sound
      u[ip, ]
      
      ## cluster the umap representation - find the optimal number of clusters
      # hierarchical dbscan
      if (TRUE) {
        u_cl = out[[i]][, c('V1', 'V2')]
        idx_cl = 1:nrow(u_cl)
      } else {
        # if >30 K, take a sample
        idx_cl = sample(1:nrow(df_umap), size = 30000)  # take only 30K points to avoid running out of RAM
        u_cl = u[idx_cl, c('V1', 'V2')]
      }
      cl = hdbscan(u_cl[, c('V1', 'V2')], minPts = 25)  # NB: very RAM-intensive! May crash with 40k points
      # cl
      plot(cl, show_flat = TRUE)
      # cl$cluster
      
      u_cl$cl = as.factor(cl$cluster)
      u_cl[, addVars] = u[, addVars]
      table(u_cl$cl, u_cl$baby)
      
      ggplot(u_cl, aes(x = V1, y = V2, color = cl, label = as.character(cl))) +  # 
        geom_text(size = 1.5) +
        # geom_point(shape = 16, size = .5, alpha = .25) +
        theme_bw()
      # ggsave('../pix/hdbscan.png', width = 15, height = 15, unit = 'cm', dpi = 300)
    }
  }
  
  df_plot = NULL
  for (i in 1:length(dm)) {  
    temp = out[[i]]
    temp$Context = names(out)[i]
    if (is.null(df_plot)) {
      df_plot = temp
    } else {
      df_plot = rbind(df_plot, temp)
    }
  }
  df_plot$Baby = df_plot$baby
  df_plot$Age = df_plot$age
  df_plot$Context = factor(unlist(contextNames[as.character(df_plot$Context)]))
  
  saveRDS(df_plot, '../data/umap_plot_per_cause.RDS')
}

ggplot(df_plot, aes(x = V1, y = V2, color = Baby)) +
  geom_point(size = .1, alpha = .5) +
  facet_wrap(~Context) +
  xlab('') + ylab('') +
  ggtitle('Separate UMAP for each context, labeled by baby') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(.01, 'npc'))
# ggsave(filename = '../pix/umap_cause~baby.png', width = 18, height = 8, units = 'cm', dpi = 600)

ggplot(df_plot, aes(x = V1, y = V2, color = Age)) +
  geom_point(size = .1, alpha = .5) +
  facet_wrap(~Context) +
  xlab('') + ylab('') +
  ggtitle('Separate UMAP for each context, labeled by age') +
  theme_bw() +
  theme(panel.grid = element_blank())
# ggsave(filename = '../pix/umap_cause~age.png', width = 20, height = 10, units = 'cm', dpi = 600)
