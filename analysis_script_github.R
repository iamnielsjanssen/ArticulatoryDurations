# Data analysis of Articulation Data
# 
# (c) Niels Janssen and Anna Siyanova-Chanturia 2017
# 
# If you use this code or data for analysis that is published in an indexed journal or repository, please cite the following article:
# Siyanova-Chanturia & Janssen "Production of familiar phrases: Frequency effects in native speakers and second language learners"
#
# Run with R (3.4.1), lme4 (1.1-13), afex (0.17-8)
# For questions email: njanssen@ull.es or anna.siyanova@vuw.ac.nz

#-------------------------------------------------------------------
# LOAD DATA AND INITIAL CLEAN UP
#-------------------------------------------------------------------
# read data matrix
rm(list=ls())
setwd('/path/to/') # CHANGE TO FOLDER WITH DATA MATRIX
dt = read.table('data_matrix_github.txt', header=TRUE)

# quick look at data, note this is the cleaned data set
head(dt) 
str(dt)

#-------------------------------------------------------------------
# PRE-PROCESSING OF DATA
#-------------------------------------------------------------------

# subject is a factor
dt$subject = as.factor(dt$subject)

# log transform relevant vars
dt$logphrase_freq = log(dt$phrasal_freq + 1) # +1 to avoid log(0)

dt$logw1_freq = log(dt$w1_bnc_freq + 1)
dt$logw2_freq = log(dt$w2_bnc_freq + 1)

dt$logw1_and_freq = log(dt$w1_and_bnc_freq + 1)
dt$logand_w2_freq = log(dt$and_w2_bnc_freq + 1)

# also log transform the dependent variable Articulation Time
dt$logAT = log(dt$AT)

# center variables
dt$phrase_lenC = scale(dt$phrasal_length, scale=FALSE)
dt$trialC = scale(dt$trial, scale=FALSE)
dt$nativeC = scale(as.integer(dt$native), scale=FALSE)
dt$stim_typeC = scale(as.integer(dt$stim_type), scale=FALSE)

dt$logphrase_freqC = scale(dt$logphrase_freq, scale=FALSE)
dt$logw1_freqC = scale(dt$logw1_freq, scale=FALSE)
dt$logw2_freqC = scale(dt$logw2_freq, scale=FALSE)
dt$logw1_and_freqC = scale(dt$logw1_and_freq, scale = FALSE)
dt$logand_w2_freqC = scale(dt$logand_w2_freq, scale = FALSE)

#-------------------------------------------------------------------
# DATA MODELING
# Following modeling strategy outlined in Bates, Kliegl, Vasishth & Baayen (2015)
#-------------------------------------------------------------------

# load lme4
require(lme4)

# 1. *maximal model fails to converge*
# note only native var for item random effect structure (other vars are between items)
# (this takes a while to end...)
model1 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC)|subject) + (1+ nativeC|stim), data=dt, REML=FALSE)
summary(model1)

# 2. attempting to find model that converges by deleting terms from random effect structure i.e., those with low variance in non-converged model
# Specifically we take out the interaction with nativeC in the by-subject random slopes
model1 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ nativeC + phrase_lenC +logphrase_freqC + logw1_freqC + logw2_freqC|subject) + (1+ nativeC|stim), data=dt, REML=FALSE)
summary(model1)

# Per Bates et al. we next run the rePCA script to see if the RE structure is too complex
E = rePCA(model1)
summary(E)
# PCA suggests that there may be four random effect components in subject structure that do not contribute variance

# next, model without correlation parameters (|| for subject and item RE)
model1 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ nativeC + phrase_lenC +logphrase_freqC + logw1_freqC + logw2_freqC||subject) + (1+ nativeC||stim), data=dt, REML=FALSE)
summary(model1)

# suggests zero variance for RE by-item nativeC, and by-subject logw2_freqC, logw1_freqC, and intercept.
model2 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (0+ nativeC + phrase_lenC +logphrase_freqC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model1, model2) # not significant, prefer simpler model2
summary(model2)

# remove by-subject slope for logphrase_freqC
model3 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (0+ nativeC + phrase_lenC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model2, model3) # not significant, prefer simpler model3
summary(model3)

# remove by-subject slope for phrase_lenC
model4 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (0+ nativeC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model3, model4) # not significant, prefer simpler model3
summary(model4)

# remove by-subject slope for nativeC
model5 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1|stim), data=dt, REML=FALSE)
anova(model4, model5) # significant! prefer model4
summary(model5)

# remove stim
model6 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (0+ nativeC ||subject), data=dt, REML=FALSE)
anova(model4, model6) # significant! prefer model4
summary(model6)

# We end with model 4 as the best model.
anova(model1,model4) # not significant, keep model4 (no accumulation)
# note cannot put correlation back in because subject has random slope but not random intercept and so correlation impossible

# SAME FOR FIXED EFFECTS

# remove main effect logw2_freq
model7 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC) + nativeC:logw2_freqC+ (0+ nativeC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model4, model7) # not significant, prefer simpler model7
summary(model7)

# remove stim_typeC
model8 = lmer(logAT ~ trialC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC) + nativeC:logw2_freqC+ (0+ nativeC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model7, model8) # not significant, prefer simpler model8
summary(model8)

# remove nativeC:logw1_freqC
model9 = lmer(logAT ~ trialC + phrase_lenC + nativeC*(logphrase_freqC)+ logw1_freqC + nativeC:logw2_freqC+ (0+ nativeC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model8, model9) # not significant, prefer simpler model9
summary(model9)

# remove nativeC:logw2_freqC
model10 = lmer(logAT ~ trialC + phrase_lenC + nativeC*(logphrase_freqC)+ logw1_freqC +  (0+ nativeC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model9, model10) # not significant, prefer simpler model10
summary(model10)

# remove logw1_freqC
model11 = lmer(logAT ~ trialC + phrase_lenC + nativeC*(logphrase_freqC)+  (0+ nativeC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model10, model11) # not significant, prefer simpler model11
summary(model11)

# MODEL 11 IS FINAL MODEL
anova(model1,model11) # no accumulation

# check condition index
source("mer-utils.R") # from Jaeger website
kappa.mer(model11)

# obtain p-values
require(afex)
mixed(logAT ~ trialC + phrase_lenC + nativeC*(logphrase_freqC)+  (0+ nativeC ||subject) + (1|stim), data=dt, REML=FALSE, method="S")

# Further examination of the interaction between nativeC:log_phrase_freqC

dtNS = dt[dt$native=="Native speaker",]
dtNNS = dt[dt$native=="Non-native speaker",]

modelNS = lmer(logAT ~ trialC + phrase_lenC + logphrase_freqC + (1|subject) + (1|stim), data=dtNS, REML=FALSE)
summary(modelNS)

modelNNS = lmer(logAT ~ trialC + phrase_lenC + logphrase_freqC  + (1|subject) + (1|stim), data=dtNNS, REML=FALSE)
summary(modelNNS)

mixed(logAT ~ trialC + phrase_lenC + logphrase_freqC + (1|subject) + (1|stim), data=dtNS, REML=FALSE,method="S")
mixed(logAT ~ trialC + phrase_lenC + logphrase_freqC + (1|subject) + (1|stim), data=dtNNS, REML=FALSE,method="S")
