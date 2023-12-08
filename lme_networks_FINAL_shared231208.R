require(MuMIn)
library(tidyverse)
library(caret)
library(leaps)
library(MASS)
library(psych)
library(modelsummary)
library(dplyr)
library(kableExtra)
library(effectsize)
library(lmtest)

df <- read.csv('/Users/yutakatsumi/Dropbox (Partners HealthCare)/Dickerson lab/Papers_chapters/Papers_submitted/Katsumi_FTP-DMN_CDRSB/Manuscript/NatComms/df_network_analysis_N48.csv', sep=',')

options(na.action = "na.fail")

# Standardization
df$CDRSB_annualized_change_Z <- (df$CDRSB_annualized_change-mean(df$CDRSB_annualized_change))/sd(df$CDRSB_annualized_change)
df$age_Z <- (df$age-mean(df$age))/sd(df$age)
df$CDRSB_t1_Z <- (df$CDRSB_t1-mean(df$CDRSB_t1))/sd(df$CDRSB_t1)
df$ftp_DMN_Z <- (df$ftp_DMN-mean(df$ftp_DMN))/sd(df$ftp_DMN)
df$ftp_FPN_Z <- (df$ftp_FPN-mean(df$ftp_FPN))/sd(df$ftp_FPN)
df$ftp_LIM_Z <- (df$ftp_LIM-mean(df$ftp_LIM))/sd(df$ftp_LIM)
df$ftp_VIS_Z <- (df$ftp_VIS-mean(df$ftp_VIS))/sd(df$ftp_VIS)
df$ftp_SAL_Z <- (df$ftp_SAL-mean(df$ftp_SAL))/sd(df$ftp_SAL)
df$ftp_DAN_Z <- (df$ftp_DAN-mean(df$ftp_DAN))/sd(df$ftp_DAN)
df$ftp_SOM_Z <- (df$ftp_SOM-mean(df$ftp_SOM))/sd(df$ftp_SOM)

## Basic and simple network models
# Basic model with age, sex, and CDR-SB at T1
mBasic <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z, data=df)
summary(mBasic)
extractAIC(mBasic)

# Full model
mFull <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DMN_Z + ftp_FPN_Z + ftp_LIM_Z + ftp_VIS_Z + ftp_DAN_Z + ftp_SAL_Z + ftp_SOM_Z, data=df)
summary(mFull)
extractAIC(mFull)

# Basic models
mDMN <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DMN_Z, data=df)
summary(mDMN)
extractAIC(mDMN)
anova(mBasic, mDMN)

mFPN <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_FPN_Z, data=df)
summary(mFPN)
extractAIC(mFPN)
anova(mBasic, mFPN)

mLIM <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_LIM_Z, data=df)
summary(mLIM)
extractAIC(mLIM)
anova(mBasic, mLIM)

mVIS <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_VIS_Z, data=df)
summary(mVIS)
extractAIC(mVIS)
anova(mBasic, mVIS)

mDAN <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DAN_Z, data=df)
summary(mDAN)
extractAIC(mDAN)
anova(mBasic, mDAN)

mSAL <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_SAL_Z, data=df)
summary(mSAL)
extractAIC(mSAL)
anova(mBasic, mSAL)

mSOM <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_SOM_Z, data=df)
summary(mSOM)
extractAIC(mSOM)
anova(mBasic, mSOM)


## Multi-model inference via an automated data-driven model selection  
m0 <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DMN_Z + ftp_LIM_Z + ftp_SAL_Z + ftp_VIS_Z + ftp_FPN_Z, data=df)
summary(m0)
extractAIC(m0)
ms1 <- dredge(m0)


## Multimodal DMN models
df$thickness_DMN_Z <- (df$thickness_DMN-mean(df$thickness_DMN))/sd(df$thickness_DMN)
df$amyloid_DMN_Z <- (df$amyloid_DMN-mean(df$amyloid_DMN))/sd(df$amyloid_DMN)

mDMN_ftp_thickness_amyloid <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DMN_Z + thickness_DMN_Z + amyloid_DMN_Z, data=df)
summary(mDMN_ftp_thickness_amyloid)
extractAIC(mDMN_ftp_thickness_amyloid)
ms2 <- dredge(mDMN_ftp_thickness_amyloid)

mDMN_ftp_amyloid <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DMN_Z + amyloid_DMN_Z, data=df)
summary(mDMN_ftp_amyloid)
extractAIC(mDMN_ftp_amyloid)

mDMN_ftp_thickness <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DMN_Z + thickness_DMN_Z, data=df)
summary(mDMN_ftp_thickness)
extractAIC(mDMN_ftp_thickness)

mDMN_amyloid <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + amyloid_DMN_Z, data=df)
summary(mDMN_amyloid)
extractAIC(mDMN_amyloid)
anova(mBasic,mDMN_amyloid)

mDMN_thickness <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + thickness_DMN_Z, data=df)
summary(mDMN_thickness)
extractAIC(mDMN_thickness)
anova(mBasic,mDMN_thickness)

mDMN_ftp <- lm(CDRSB_annualized_change_Z ~ age_Z + factor(sex) + CDRSB_t1_Z + ftp_DMN_Z, data=df)
summary(mDMN_ftp)
extractAIC(mDMN_ftp)