# A separate distance matrix for each age group. Input to UMAP is distance matrix based on DTW of frame-by-frame features

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('umap_helperFun.R')
source('zz_formatting.R')
library(ggplot2)
library(patchwork)

### All babies, but separate UMAP for each cause or age group
nMax = Inf  # max sample size (5K ~= 3.5 h of machine time to get the full distance matrix)

dm_filename = paste0('../data/distMatrix_features_byAge_', nMax / 1000, 'K.RDS')

if (file.exists(dm_filename)) {
  dm = readRDS(dm_filename)
} else {
  table(df$age_month)
  ages = sort(unique(df$age_month))
  dm_by_age = vector('list', length(ages))
  names(dm_by_age) = ages
  for (i in 1:length(ages)) {
    idx_i = which(df$age_month == ages[i])
    if (length(idx_i) > nMax) idx_i = sample(idx_i, size = nMax)
    dm_by_age[[i]]$sample = df$file[idx_i]
    dm_by_age[[i]]$dm = getDistMat_multiCore(acAn = acAn_full[idx_i], add_vars, vars_for_log)
  }
  saveRDS(dm_by_age, dm_filename)
  dm = dm_by_age
}

## UMAP
if (file.exists('../data/umap_plot_per_age.RDS')) {
  df_plot = readRDS('../data/umap_plot_per_age.RDS')
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
    temp$age = paste(names(out)[i], ' months')
    if (is.null(df_plot)) {
      df_plot = temp
    } else {
      df_plot = rbind(df_plot, temp)
    }
  }
  df_plot$Baby = df_plot$baby
  df_plot$Context = factor(unlist(contextNames[as.character(df_plot$cause_stop_engl)]))
  saveRDS(df_plot, '../data/umap_plot_per_age.RDS')
}

ggplot(df_plot, aes(x = V1, y = V2, color = Baby)) +
  geom_point(size = .03, alpha = .5) +
  facet_wrap(~age) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  xlab('') + ylab('') +
  # ggtitle('A. UMAP per age group, labeled by baby') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'right')
# ggsave(filename = '../pix/umap_age~baby.png', width = 14, height = 10, units = 'cm', dpi = 600)

# n per age and per baby
ag_age = aggregate(V1 ~ age, df_plot, length)
range(ag_age$V1) # 4580 12399
ag_baby = aggregate(V1 ~ baby, df_plot, length)
range(ag_baby$V1)  # 327 4140

p_all = ggplot(df_plot, aes(x = V1, y = V2, color = Context)) +
  geom_point(size = .1, alpha = .5) +
  facet_wrap(~age) +
  xlab('') + ylab('') +
  # ggtitle('B. UMAP per age group, labeled by context') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'right')
p_all
# ggsave(filename = '../pix/umap_age~cause.png', width = 14, height = 10, units = 'cm', dpi = 600)

df_plot$sex = df$sex[match(df_plot$baby, df$baby)]
df_plot$sex_fmt = ifelse(df$sex == 'M', 'Boy', 'Girl')
ggplot(df_plot, aes(x = V1, y = V2, color = sex_fmt)) +
  geom_point(size = .03, alpha = .25) +
  facet_wrap(~age, scales = 'free', nrow = 1) +
  xlab('') + ylab('') +
  # ggtitle('UMAP per age group, labeled by sex') +
  scale_color_manual('', values = c('blue', 'yellow')) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'top')
ggsave(filename = '../pix/umap_age~sex.png', width = 20, height = 7, units = 'cm', dpi = 600)
# write.csv(df_plot, '../data/fig1a_umap_age~sex.csv', row.names = FALSE)

## clustering quality
library(IDmeasurer)
hs_baby = data.frame(age_month = unique(df_plot$age_month))
for (i in 1:nrow(hs_baby))
  hs_baby$hs[i] = IDmeasurer::calcHS(df_plot[df_plot$age_month == hs_baby$age_month[i], c('baby', 'V1', 'V2')])[2]

hs_context = data.frame(age_month = unique(df_plot$age_month))
for (i in 1:nrow(hs_context))
  hs_context$hs[i] = IDmeasurer::calcHS(df_plot[df_plot$age_month == hs_context$age_month[i], c('Context', 'V1', 'V2')])[2]


## A separate plot just for 1.5-m.o.
d15 = df_plot[df_plot$age_month == 1.5, ]
p15 = ggplot(d15, aes(x = V1, y = V2, color = Baby)) +
  geom_point(size = .25, alpha = .5) +
  facet_wrap(~Context) +
  xlab('') + ylab('') +
  ggtitle('(B) UMAP for 1.5-m.o., labeled by baby') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'right')
# ggsave(filename = '../pix/umap_1.5~cause.png', width = 20, height = 8, units = 'cm', dpi = 600)


## Plots per baby from this 1.5-m.o. group, coloring by context
p15_perBaby = ggplot(d15, aes(x = V1, y = V2, color = Context)) +
  geom_point(size = .5, alpha = .5) +
  facet_wrap(~Baby, ncol = 6) +
  xlab('') + ylab('') +
  scale_x_continuous(n.breaks = 4) +
  ggtitle('(C) UMAP for 1.5-m.o. per baby') +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'right')
# ggsave(filename = '../pix/umap_1.5~cause_per_baby.png', width = 20, height = 20, units = 'cm', dpi = 600)


## composite plot
p_all + p15 + p15_perBaby +
  plot_layout(nrow = 3, heights = c(.8, 1.2, 2)) &
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        plot.margin = unit(c(.1, 0, 0, 0), 'cm'))
# ggsave(filename = '../pix/umap_composite_per-age.png', width = 18, height = 24, units = 'cm', dpi = 600)

