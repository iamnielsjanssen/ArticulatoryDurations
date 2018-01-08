# Data analysis of Articulation Data
# 
# Niels Janssen and Anna Siyanova-Chanturia 2017
# 
# If you use this code for analysis that is published in an indexed journal or repository, please cite the following article:
# Siyanova-Chanturia & Janssen "Production of familiar phrases: Frequency effects in native speakers and second language learners"
#
# For questions email: njanssen@ull.es or anna.siyanova@vuw.ac.nz
#
# Update: Dec 2017, added analysis Length of stay, age Exposure, and speech rate


#-------------------------------------------------------------------
# LOAD DATA AND INITIAL CLEAN UP
#-------------------------------------------------------------------
# read data matrix
rm(list=ls())
setwd('/Data/Dropbox/Documents/R_analyses/La_Laguna/Binomial reading/') # CHANGE TO FOLDER WITH DATA MATRIX
#dt = read.table('data_matrix_github.txt', header=TRUE)
dt = read.table("data_matrix_DEC2017.txt", header=TRUE, sep="\t")

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

# compute by-participant speech rate (UPDATE DEC-2017)
# speech rate defined as the mean (item phrase duration / item phrase length in phonemes) for each participant

dt$itemrate = dt$AT / dt$phrasal_length
dt$speechrate = ave(dt$itemrate, dt$subject)
dt$speechrateC = scale(log(dt$speechrate), scale=FALSE)

#-------------------------------------------------------------------
# DATA MODELING
# Following modeling strategy outlined in Bates, Kliegl, Vasishth & Baayen (2015)
#-------------------------------------------------------------------

# load lme4
require(lme4)

# 1. *maximal model fails to converge*
# note only native var for item random effect structure (other vars are between items)
# (this takes a while to end...)
model1 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ trialC + stim_typeC + phrase_lenC + logphrase_freqC + logw1_freqC + logw2_freqC|subject) + (1+ nativeC|stim), data=dt, REML=FALSE)
summary(model1)

# attempt to find model that converges by deleting random-effect terms that produce near zero variance in non-converged model
# delete by-subject trialC
model1 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ stim_typeC + phrase_lenC + logphrase_freqC + logw1_freqC + logw2_freqC|subject) + (1+ nativeC|stim), data=dt, REML=FALSE)
summary(model1)

# Per Bates et al. we next run the rePCA script to see if the RE structure is too complex
# devtools::install_github("dmbates/RePsychLing")
require(RePsychLing)
E = rePCA(model1)
summary(E)
# PCA suggests that there may be three random effect components in subject structure that do not contribute variance

# next, model without correlation parameters (|| for subject and item RE)
model1 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ stim_typeC + phrase_lenC + logphrase_freqC + logw1_freqC + logw2_freqC||subject) + (1+ nativeC||stim), data=dt, REML=FALSE)
summary(model1)

# suggests zero variance for RE by-item nativeC, and by-subject logw2_freqC, logw1_freqC, logphrase_freqC.
model2 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ stim_typeC + phrase_lenC||subject) + (1|stim), data=dt, REML=FALSE)
anova(model1, model2) # not significant, prefer simpler model2
summary(model2)

# remove by-subject slope for phrase_lenC
model3 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1+ stim_typeC ||subject) + (1|stim), data=dt, REML=FALSE)
anova(model2, model3) # not significant, prefer simpler model3
summary(model3)

# remove by-subject slope for stim_typeC
model4 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC*(logphrase_freqC + logw1_freqC + logw2_freqC) + (1|subject) + (1|stim), data=dt, REML=FALSE)
anova(model3, model4) # not significant, prefer simpler model4
summary(model4)

# We end with model 4 as the best model.
anova(model1,model4) # not significant, keep model4 (no accumulation)
# note cannot put correlation back in because subject has random slope but not random intercept and so correlation impossible

# SAME FOR FIXED EFFECTS

# remove main effect logw2_freq
model5 = lmer(logAT ~ trialC + stim_typeC + phrase_lenC + nativeC:logw2_freqC + nativeC*(logphrase_freqC + logw1_freqC) + (1|subject) + (1|stim), data=dt, REML=FALSE)
anova(model4, model5) # not significant, prefer simpler model7
summary(model5)

# remove stim_typeC
model6 = lmer(logAT ~ trialC + phrase_lenC + nativeC:logw2_freqC + nativeC*(logphrase_freqC + logw1_freqC) + (1|subject) + (1|stim), data=dt, REML=FALSE)
anova(model5, model6) # not significant, prefer simpler model8
summary(model6)

# remove nativeC:logw1_freqC
model7 = lmer(logAT ~ trialC + phrase_lenC + nativeC:logw2_freqC + nativeC*logphrase_freqC + logw1_freqC + (1|subject) + (1|stim), data=dt, REML=FALSE)
anova(model6, model7) # not significant, prefer simpler model9
summary(model7)

# remove nativeC:logw2_freqC
model8 = lmer(logAT ~ trialC + phrase_lenC + nativeC*logphrase_freqC + logw1_freqC + (1|subject) + (1|stim), data=dt, REML=FALSE)
anova(model7, model8) # not significant, prefer simpler model10
summary(model8)

# remove logw1_freqC
model9 = lmer(logAT ~ trialC + phrase_lenC + nativeC*logphrase_freqC + (1|subject) + (1|stim), data=dt, REML=FALSE)
anova(model8, model9) # not significant, prefer simpler model11
summary(model9)

# remove nativeC:logphrase_freqC
model10 = lmer(logAT ~ trialC + phrase_lenC + nativeC+logphrase_freqC + (1|subject) + (1|stim), data=dt, REML=FALSE)
anova(model9, model10) # significant!, prefer model9
summary(model10)

# MODEL 9 IS FINAL MODEL
anova(model1,model9) # no accumulation

# check condition index
source("mer-utils.R") # from Jaeger website
kappa.mer(model9) # = 1.12

# obtain p-values
require(afex)
mixed(logAT ~ trialC + phrase_lenC + nativeC*logphrase_freqC + (1|subject) + (1|stim), data=dt, REML=FALSE, method="S")


# Further examination of the interaction between nativeC:log_phrase_freqC

dtNS = dt[dt$native=="Native speaker",]
dtNNS = dt[dt$native=="Non-native speaker",]

dtNNS$length_stayC = scale(dtNNS$length_stay, scale=FALSE)
dtNNS$age_exposeC =  scale(dtNNS$age_expose, scale=FALSE)

modelNS = lmer(logAT ~  trialC + phrase_lenC + logphrase_freqC + (1|subject) + (1|stim), data=dtNS, REML=FALSE)
summary(modelNS)

modelNNS = lmer(logAT ~ trialC + phrase_lenC + logphrase_freqC  + (1|subject) + (1|stim), data=dtNNS, REML=FALSE)
summary(modelNNS)

mixed(logAT ~ trialC + phrase_lenC + logphrase_freqC + (1|subject) + (1|stim), data=dtNS, REML=FALSE,method="S")
mixed(logAT ~ trialC + phrase_lenC + logphrase_freqC + (1|subject) + (1|stim), data=dtNNS, REML=FALSE,method="S")

# Looking at Length of Stay in ESC and Age of Eposure

# Length of stay in NNS and interaction with phrase_freq
modelNNS = lmer(logAT ~ trialC + phrase_lenC + length_stayC*logphrase_freqC  + (1|subject) + (1|stim), data=dtNNS, REML=FALSE)
summary(modelNNS)

# Age of exposure in NNS and interaction with phrase_freq
modelNNS = lmer(logAT ~ trialC + phrase_lenC + age_exposeC*logphrase_freqC  + (1|subject) + (1|stim), data=dtNNS, REML=FALSE)
summary(modelNNS)
