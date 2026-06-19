# Analysis of the playbacks results

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(ggplot2)
library(brms)

### Analyses for implicit design ###

df_implicit = read.table(file="../data/df_implicit.csv", header=T, sep=";",dec=".")

## select only the trials 72 to 102 (=test session) for the analysis
df_implicit <- df_implicit %>% filter (num_Row > 72)
df_implicit <- df_implicit %>% mutate_if(is.character,as.factor)
summary(df_implicit)

## Bayesian analysis ##
model4_file = 'path_to_file_model4.RDS'
if (file.exists(model4_file)) {
  model4 = readRDS(model4_file) 
} else {
  model4 <- brm(
    success ~ nSameSession + real_cause + sex + parentality +
      real_cause:parentality + sex:parentality + 
      (1 | subjectID + babyID),  
    data = df_implicit, family = 'bernoulli',  
    chains = 2, cores = 2, warmup = 1000, iter = 3000,
    init_r = 0.1,  
    control = list(adapt_delta = 0.99),
    prior = c(
      set_prior('normal(0, 1)'),
      set_prior('normal(-.84, 1)', class = 'Intercept')  # log(.3/.7) = -.84 (chance at 1/3)
    )
  )
  saveRDS(model4, 'path_to_file_model4.RDS')
}

summary(model4)
pp_check(model4)
conditional_effects(model4)

# extract fitted values
newdata = expand.grid(sex = levels(df_implicit$sex),
                      parentality = levels(df_implicit$parentality),
                      real_cause = levels(df_implicit$real_cause),
                      nSameSession = 0)  # set nSameSession to 0 assuming subjects were not exposed to the same session during the implicit training
fit = fitted(model4, newdata = newdata, robust = TRUE, re_formula = NA) * 100
colnames(fit) = c('fit', 'se', 'lwr', 'upr')
pl = cbind(newdata, fit)


sex_label <- c(
  "female" = "Female adult listeners",
  "male" = "Male adult listeners"
)

ggplot(pl, aes(real_cause, fit, ymin = lwr, ymax = upr, color=parentality, shape=parentality)) +
  geom_pointrange(position = position_dodge(.3)) +
  geom_hline(yintercept = 100/3, linetype = 3) +
  xlab('') +
  ylab('Recognition performance (Accuracy, %)') +
  facet_grid(cols=vars(sex), labeller = labeller(sex=sex_label)) +
  #facet_wrap(~nSameSession + sex, nrow = 2) +
  ylim(20,50) +
  theme_bw()+
  ggtitle("a")+
  labs(color="Parentality", shape="Parentality") +
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=15, color="black"), 
        axis.text.x = element_text(size=12, color="black"),
        axis.title.y = element_text(size=13, color="black"),
        plot.title = element_text(size=15, color="black"), 
        strip.text.x = element_text(size=15, color="black"))



### Analyses for explicit design ###
df_explicit = read.table(file="../data/df_explicit.csv", header=T, sep=";",dec=".")
df_explicit <- df_explicit %>% mutate_if(is.character,as.factor)
df_explicit$subjectID <- as.factor(df_explicit$subjectID)
summary(df_explicit)

# observed values for the 30 last trials of each participant (designed to have 10 syllables of each condition)
df_test <- df_explicit %>% filter (trial_nr >= 72)

## Bayesian analysis ##
model8_file = 'path_to_file_model8.RDS'
if (file.exists(model8_file)) {
  model8 = readRDS(model8_file)  # load this one instead of re-running the model every time
} else {
  model8 <- brm(
    success ~ nSameSession + real_cause + sex + parentality +
      real_cause:parentality + sex:parentality + 
      (1 | subjectID) + (1|babyID),  # AA: no random slopes? Eg real_cause|baby_id
    data = df_test, family = 'bernoulli',  
    chains = 2, cores = 2, warmup = 1000, iter = 3000,
    init_r = 0.1,  # a way to ensure convergence with a lot of random slopes - see https://discourse.mc-stan.org/t/model-with-many-correlated-varying-slopes-sampling-not-done/6189)
    # sample_prior = 'only',  # first just preview the prior
    control = list(adapt_delta = 0.99),
    prior = c(
      set_prior('normal(0, 1)'),
      set_prior('normal(-.84, 1)', class = 'Intercept')  # log(.3/.7) = -.84 (chance at 1/3)
    )
  )
  saveRDS(model8, 'path_to_file_model8.RDS')
}

summary(model8)
pp_check(model8)
conditional_effects(model8)

# extract fitted values
newdata = expand.grid(sex = levels(df_explicit$sex),
                      parentality = levels(df_explicit$parentality),
                      real_cause = levels(df_explicit$real_cause),
                      nSameSession = 1) # set nSameSession to 1 (1st presentation of this given recording session)
fit = fitted(model8, newdata = newdata, robust = TRUE, re_formula = NA) * 100  # *100 = convert to %
colnames(fit) = c('fit', 'se', 'lwr', 'upr')
pl = cbind(newdata, fit)

sex_label <- c(
  "female" = "Female adult listeners",
  "male" = "Male adult listeners"
)

ggplot(pl, aes(real_cause, fit, ymin = lwr, ymax = upr, color=parentality, shape=parentality)) +
  geom_pointrange(position = position_dodge(.3)) +
  geom_hline(yintercept = 100/3, linetype = 3) +
  xlab('') +
  ylab('') +
  facet_grid(cols=vars(sex), labeller = labeller(sex=sex_label)) +
  #facet_wrap(~nSameSession + sex, nrow = 2) +
  # ylim(20,55) +
  theme_bw()+
  ggtitle("b")+
  labs(color="Parentality", shape="Parentality") +
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=15, color="black"), 
        axis.text.x = element_text(size=12, color="black"),
        axis.title.y = element_text(size=13, color="black"),
        plot.title = element_text(size=15, color="black"), 
        strip.text.x = element_text(size=15, color="black"))
