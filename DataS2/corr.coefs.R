#####################################################################
#### This script is for creating Figures in Chapter 1 Manuscript ####
#####################################################################

rm(list=ls())
library(here)
library(rstan)
library(coda)

#### Load data for all figures ####
# ARCA 

#load in data, subset out C treatment and columns of groups then remove NAs
ARCA<-read.csv("Data/groups.ARCA.csv")
ARCA<-subset(ARCA, treatment=="C")
ARCA<-subset(ARCA, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
ARCA<-as.matrix(ARCA)
# load in the posteriors 
load(here("Bayes data", "ARCA_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

ARCA_nat <- posteriors$alpha_NatForb
ARCA_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(ARCA)
ARCAF_intra <- matrix(NA, length(posteriors$alpha_NatForb), 44)
ARCAF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 44)
ARCAF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 44)
ARCAF <-matrix(NA, length(posteriors$alpha_NatForb), 44)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:44){
  if (ARCA[i,3]>0){
    ARCAF_nat[,i]<-posteriors$alpha_NatForb*ARCA[i,3]
  }
  if (ARCA[i,2]+ARCA[i,4]>0){
    ARCAF_ex[,i]<-posteriors$alpha_InvGraminoid*ARCA[i,2]+posteriors$alpha_InvForb*ARCA[i,4]
  }
  if (sum(ARCA[i,1:4]>0)){
    ARCAF[,i]<-posteriors$alpha_intra*ARCA[i,1]+posteriors$alpha_InvGraminoid*ARCA[i,2]+posteriors$alpha_NatForb*ARCA[i,3]+posteriors$alpha_InvForb*ARCA[i,4]
  }
}

# calculate credible intervals 
HPDARCA<-HPDinterval((as.mcmc(as.numeric(ARCAF_nat))))
lowerARCAnat<-HPDARCA[,1]
upperARCAnat<-HPDARCA[,2] 

HPDARCA<-HPDinterval((as.mcmc(as.numeric(ARCAF_ex))))
lowerARCAex<-HPDARCA[,1]
upperARCAex<-HPDARCA[,2] 


remove(HPDARCA)

ARCAmeanF<-mean(ARCAF, na.rm=T)
ArcaL <- posteriors$lambda
ARCAmeanL<-mean(posteriors$lambda)
HPDARCA<-HPDinterval((as.mcmc(as.numeric(ARCAF))))
HPDARCAlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerARCA<-HPDARCA[,1]
upperARCA<-HPDARCA[,2] 
lowerARCAlam<-HPDARCAlam[,1]
upperARCAlam<-HPDARCAlam[,2]

remove(PrelimFit)

# DAGL

#load in data, subset out C treatment and columns of groups then remove NAs
DAGL<-read.csv("Data/groups.DAGL.csv")
DAGL<-subset(DAGL, treatment=="C")
DAGL<-subset(DAGL, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
DAGL<-as.matrix(DAGL)
# load in the posteriors 
load(here("Bayes data", "DAGL_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

DAGL_nat <- posteriors$alpha_NatForb
DAGL_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(DAGL)
DAGLF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 43)
DAGLF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 43)
DAGLF <-matrix(NA, length(posteriors$alpha_NatForb), 43)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:43){
  
  if (DAGL[i,3]>0){
    DAGLF_nat[,i]<-posteriors$alpha_NatForb*DAGL[i,3]
  }
  if (DAGL[i,2]+DAGL[i,4]>0){
    DAGLF_ex[,i]<-posteriors$alpha_InvGraminoid*DAGL[i,2]+posteriors$alpha_InvForb*DAGL[i,4]
  }
  if (sum(DAGL[i,1:4]>0)){
    DAGLF[,i]<-posteriors$alpha_intra*DAGL[i,1]+posteriors$alpha_InvGraminoid*DAGL[i,2]+posteriors$alpha_NatForb*DAGL[i,3]+posteriors$alpha_InvForb*DAGL[i,4]
  }
}

# calculate credible intervals 
HPDDAGL<-HPDinterval((as.mcmc(as.numeric(DAGLF_nat))))
lowerDAGLnat<-HPDDAGL[,1]
upperDAGLnat<-HPDDAGL[,2] 

HPDDAGL<-HPDinterval((as.mcmc(as.numeric(DAGLF_ex))))
lowerDAGLex<-HPDDAGL[,1]
upperDAGLex<-HPDDAGL[,2] 

remove(HPDDAGL)

DAGLmeanF<-mean(DAGLF, na.rm=T)
DaglL <- posteriors$lambda
DAGLmeanL<-mean(posteriors$lambda)
HPDDAGL<-HPDinterval((as.mcmc(as.numeric(DAGLF))))
HPDDAGLlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerDAGL<-HPDDAGL[,1]
upperDAGL<-HPDDAGL[,2] 
lowerDAGLlam<-HPDDAGLlam[,1]
upperDAGLlam<-HPDDAGLlam[,2]

remove(PrelimFit)

# GITE

#load in data, subset out C treatment and columns of groups then remove NAs
GITE<-read.csv("Data/groups.GITE.csv")
GITE<-subset(GITE, treatment=="C")
GITE<-subset(GITE, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
GITE<-as.matrix(GITE)
# load in the posteriors 
load(here("Bayes data", "GITE_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

GITE_nat <- posteriors$alpha_NatForb
GITE_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(GITE)
GITEF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 48)
GITEF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 48)
GITEF <-matrix(NA, length(posteriors$alpha_NatForb), 48)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:48){
  if (GITE[i,3]>0){
    GITEF_nat[,i]<-posteriors$alpha_NatForb*GITE[i,3]
  }
  if (GITE[i,2]+GITE[i,4]>0){
    GITEF_ex[,i]<-posteriors$alpha_InvGraminoid*GITE[i,2]+posteriors$alpha_InvForb*GITE[i,4]
  }
  if (sum(GITE[i,1:4]>0)){
    GITEF[,i]<-posteriors$alpha_intra*GITE[i,1]+posteriors$alpha_InvGraminoid*GITE[i,2]+posteriors$alpha_NatForb*GITE[i,3]+posteriors$alpha_InvForb*GITE[i,4]
  }
}

# calculate credible intervals 
HPDGITE<-HPDinterval((as.mcmc(as.numeric(GITEF_nat))))
lowerGITEnat<-HPDGITE[,1]
upperGITEnat<-HPDGITE[,2] 

HPDGITE<-HPDinterval((as.mcmc(as.numeric(GITEF_ex))))
lowerGITEex<-HPDGITE[,1]
upperGITEex<-HPDGITE[,2] 

remove(HPDGITE)

GITEmeanF<-mean(GITEF, na.rm=T)
GiteL <- posteriors$lambda
GITEmeanL<-mean(posteriors$lambda)
HPDGITE<-HPDinterval((as.mcmc(as.numeric(GITEF))))
HPDGITElam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerGITE<-HPDGITE[,1]
upperGITE<-HPDGITE[,2] 
lowerGITElam<-HPDGITElam[,1]
upperGITElam<-HPDGITElam[,2]

remove(PrelimFit)

# HYGL
#load in data, subset out C treatment and columns of groups then remove NAs
HYGL<-read.csv("Data/groups.HYGL.csv")
HYGL<-subset(HYGL, treatment=="C")
HYGL<-subset(HYGL, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
HYGL<-as.matrix(HYGL)
# load in the posteriors 
load(here("Bayes data", "HYGL_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

HYGL_nat <- posteriors$alpha_NatForb
HYGL_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(HYGL)
HYGLF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 46)
HYGLF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 46)
HYGLF <-matrix(NA, length(posteriors$alpha_NatForb), 46)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:46){
  if (HYGL[i,3]>0){
    HYGLF_nat[,i]<-posteriors$alpha_NatForb*HYGL[i,3]
  }
  if (HYGL[i,2]+HYGL[i,4]>0){
    HYGLF_ex[,i]<-posteriors$alpha_InvGraminoid*HYGL[i,2]+posteriors$alpha_InvForb*HYGL[i,4]
  }
  if (sum(HYGL[i,1:4]>0)){
    HYGLF[,i]<-posteriors$alpha_intra*HYGL[i,1]+posteriors$alpha_InvGraminoid*HYGL[i,2]+posteriors$alpha_NatForb*HYGL[i,3]+posteriors$alpha_InvForb*HYGL[i,4]
  }
}

# calculate credible intervals 
HPDHYGL<-HPDinterval((as.mcmc(as.numeric(HYGLF_nat))))
lowerHYGLnat<-HPDHYGL[,1]
upperHYGLnat<-HPDHYGL[,2] 

HPDHYGL<-HPDinterval((as.mcmc(as.numeric(HYGLF_ex))))
lowerHYGLex<-HPDHYGL[,1]
upperHYGLex<-HPDHYGL[,2] 

remove(HPDHYGL)

HYGLmeanF<-mean(HYGLF, na.rm=T)
HyglL <- posteriors$lambda
HYGLmeanL<-mean(posteriors$lambda)
HPDHYGL<-HPDinterval((as.mcmc(as.numeric(HYGLF))))
HPDHYGLlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerHYGL<-HPDHYGL[,1]
upperHYGL<-HPDHYGL[,2] 
lowerHYGLlam<-HPDHYGLlam[,1]
upperHYGLlam<-HPDHYGLlam[,2]

remove(PrelimFit)

# MEDI
#load in data, subset out C treatment and columns of groups then remove NAs
MEDI<-read.csv("Data/groups.MEDI.csv")
MEDI<-subset(MEDI, treatment=="C")
MEDI<-subset(MEDI, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
MEDI<-as.matrix(MEDI)
# load in the posteriors 
load(here("Bayes data", "MEDI_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

MEDI_nat <- posteriors$alpha_NatForb
MEDI_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(MEDI)
MEDIF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 42)
MEDIF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 42)
MEDIF <-matrix(NA, length(posteriors$alpha_NatForb), 42)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:42){
  if (MEDI[i,3]>0){
    MEDIF_nat[,i]<-posteriors$alpha_NatForb*MEDI[i,3]
  }
  if (MEDI[i,2]+MEDI[i,4]>0){
    MEDIF_ex[,i]<-posteriors$alpha_InvGraminoid*MEDI[i,2]+posteriors$alpha_InvForb*MEDI[i,4]
  }
  if (sum(MEDI[i,1:4]>0)){
    MEDIF[,i]<-posteriors$alpha_intra*MEDI[i,1]+posteriors$alpha_InvGraminoid*MEDI[i,2]+posteriors$alpha_NatForb*MEDI[i,3]+posteriors$alpha_InvForb*MEDI[i,4]
  }
}

# calculate credible intervals 
HPDMEDI<-HPDinterval((as.mcmc(as.numeric(MEDIF_nat))))
lowerMEDInat<-HPDMEDI[,1]
upperMEDInat<-HPDMEDI[,2] 

HPDMEDI<-HPDinterval((as.mcmc(as.numeric(MEDIF_ex))))
lowerMEDIex<-HPDMEDI[,1]
upperMEDIex<-HPDMEDI[,2] 

remove(HPDMEDI)

MEDImeanF<-mean(MEDIF, na.rm=T)
MediL <- posteriors$lambda
MEDImeanL<-mean(posteriors$lambda)
HPDMEDI<-HPDinterval((as.mcmc(as.numeric(MEDIF))))
HPDMEDIlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerMEDI<-HPDMEDI[,1]
upperMEDI<-HPDMEDI[,2] 
lowerMEDIlam<-HPDMEDIlam[,1]
upperMEDIlam<-HPDMEDIlam[,2]

remove(PrelimFit)

# MOMO
#load in data, subset out C treatment and columns of groups then remove NAs
MOMO<-read.csv("Data/groups.MOMO.csv")
MOMO<-subset(MOMO, treatment=="C")
MOMO<-subset(MOMO, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
MOMO<-as.matrix(MOMO)
# load in the posteriors 
load(here("Bayes data", "MOMO_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

MOMO_nat <- posteriors$alpha_NatForb
MOMO_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(MOMO)
MOMOF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 45)
MOMOF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 45)
MOMOF <-matrix(NA, length(posteriors$alpha_NatForb), 45)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:45){
  if (MOMO[i,3]>0){
    MOMOF_nat[,i]<-posteriors$alpha_NatForb*MOMO[i,3]
  }
  if (MOMO[i,2]+MOMO[i,4]>0){
    MOMOF_ex[,i]<-posteriors$alpha_InvGraminoid*MOMO[i,2]+posteriors$alpha_InvForb*MOMO[i,4]
  }
  if (sum(MOMO[i,1:4]>0)){
    MOMOF[,i]<-posteriors$alpha_intra*MOMO[i,1]+posteriors$alpha_InvGraminoid*MOMO[i,2]+posteriors$alpha_NatForb*MOMO[i,3]+posteriors$alpha_InvForb*MOMO[i,4]
  }
}

# calculate credible intervals 
HPDMOMO<-HPDinterval((as.mcmc(as.numeric(MOMOF_nat))))
lowerMOMOnat<-HPDMOMO[,1]
upperMOMOnat<-HPDMOMO[,2] 

HPDMOMO<-HPDinterval((as.mcmc(as.numeric(MOMOF_ex))))
lowerMOMOex<-HPDMOMO[,1]
upperMOMOex<-HPDMOMO[,2] 

remove(HPDMOMO)

MOMOmeanF<-mean(MOMOF, na.rm=T)
MomoL <- posteriors$lambda
MOMOmeanL<-mean(posteriors$lambda)
HPDMOMO<-HPDinterval((as.mcmc(as.numeric(MOMOF))))
HPDMOMOlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerMOMO<-HPDMOMO[,1]
upperMOMO<-HPDMOMO[,2] 
lowerMOMOlam<-HPDMOMOlam[,1]
upperMOMOlam<-HPDMOMOlam[,2]

remove(PrelimFit)

# PEAI
#load in data, subset out C treatment and columns of groups then remove NAs
PEAI<-read.csv("Data/groups.PEAI.csv")
PEAI<-subset(PEAI, treatment=="C")
PEAI<-subset(PEAI, select=c(intra, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
PEAI<-as.matrix(PEAI)
# load in the posteriors 
load(here("Bayes data", "PEAI_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

dim(PEAI)
PEAIF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 41)
PEAIF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 41)
PEAIF <-matrix(NA, length(posteriors$alpha_NatForb), 41)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:41){
  if (PEAI[i,2]>0){
    PEAIF_nat[,i]<-posteriors$alpha_NatForb*PEAI[i,2]
  }
  if (PEAI[i,3]>0){
    PEAIF_ex[,i]<-posteriors$alpha_InvForb*PEAI[i,3]
  }
  if (sum(PEAI[i,1:3]>0)){
    PEAIF[,i]<-posteriors$alpha_intra*PEAI[i,1]+posteriors$alpha_NatForb*PEAI[i,2]+posteriors$alpha_InvForb*PEAI[i,3]
  }
}

# calculate credible intervals 
HPDPEAI<-HPDinterval((as.mcmc(as.numeric(PEAIF_nat))))
lowerPEAInat<-HPDPEAI[,1]
upperPEAInat<-HPDPEAI[,2] 

HPDPEAI<-HPDinterval((as.mcmc(as.numeric(PEAIF_ex))))
lowerPEAIex<-HPDPEAI[,1]
upperPEAIex<-HPDPEAI[,2] 

remove(HPDPEAI)

PEAImeanF<-mean(PEAIF, na.rm=T)
PeaiL <- posteriors$lambda
PEAImeanL<-mean(posteriors$lambda)
HPDPEAI<-HPDinterval((as.mcmc(as.numeric(PEAIF))))
HPDPEAIlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerPEAI<-HPDPEAI[,1]
upperPEAI<-HPDPEAI[,2] 
lowerPEAIlam<-HPDPEAIlam[,1]
upperPEAIlam<-HPDPEAIlam[,2]

remove(PrelimFit)

# PLDE
#load in data, subset out C treatment and columns of groups then remove NAs
PLDE<-read.csv("Data/groups.PLDE.csv")
PLDE<-subset(PLDE, treatment=="C")
PLDE<-subset(PLDE, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
PLDE<-as.matrix(PLDE)
# load in the posteriors 
load(here("Bayes data", "PLDE_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

PLDE_nat <- posteriors$alpha_NatForb
PLDE_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(PLDE)
PLDEF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 45)
PLDEF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 45)
PLDEF <-matrix(NA, length(posteriors$alpha_NatForb), 45)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:45){
  if (PLDE[i,3]>0){
    PLDEF_nat[,i]<-posteriors$alpha_NatForb*PLDE[i,3]
  }
  if (PLDE[i,2]+PLDE[i,4]>0){
    PLDEF_ex[,i]<-posteriors$alpha_InvGraminoid*PLDE[i,2]+posteriors$alpha_InvForb*PLDE[i,4]
  }
  if (sum(PLDE[i,1:4]>0)){
    PLDEF[,i]<-posteriors$alpha_intra*PLDE[i,1]+posteriors$alpha_InvGraminoid*PLDE[i,2]+posteriors$alpha_NatForb*PLDE[i,3]+posteriors$alpha_InvForb*PLDE[i,4]
  }
}

# calculate credible intervals 
HPDPLDE<-HPDinterval((as.mcmc(as.numeric(PLDEF_nat))))
lowerPLDEnat<-HPDPLDE[,1]
upperPLDEnat<-HPDPLDE[,2] 

HPDPLDE<-HPDinterval((as.mcmc(as.numeric(PLDEF_ex))))
lowerPLDEex<-HPDPLDE[,1]
upperPLDEex<-HPDPLDE[,2] 

remove(HPDPLDE)

PLDEmeanF<-mean(PLDEF, na.rm=T)
PldeL <- posteriors$lambda
PLDEmeanL<-mean(posteriors$lambda)
HPDPLDE<-HPDinterval((as.mcmc(as.numeric(PLDEF))))
HPDPLDElam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerPLDE<-HPDPLDE[,1]
upperPLDE<-HPDPLDE[,2] 
lowerPLDElam<-HPDPLDElam[,1]
upperPLDElam<-HPDPLDElam[,2]

remove(PrelimFit)

# POCA
#load in data, subset out C treatment and columns of groups then remove NAs
POCA<-read.csv("Data/groups.POCA.csv")
POCA<-subset(POCA, treatment=="C")
POCA<-subset(POCA, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
POCA<-as.matrix(POCA)
# load in the posteriors 
load(here("Bayes data", "POCA_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

dim(POCA)
POCAF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 44)
POCAF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 44)
POCAF <-matrix(NA, length(posteriors$alpha_NatForb), 44)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:44){
  if (POCA[i,3]>0){
    POCAF_nat[,i]<-posteriors$alpha_NatForb*POCA[i,3]
  }
  if (POCA[i,2]+POCA[i,4]>0){
    POCAF_ex[,i]<-posteriors$alpha_InvGraminoid*POCA[i,2]+posteriors$alpha_InvForb*POCA[i,4]
  }
  if (sum(POCA[i,1:4]>0)){
    POCAF[,i]<-posteriors$alpha_intra*POCA[i,1]+posteriors$alpha_InvGraminoid*POCA[i,2]+posteriors$alpha_NatForb*POCA[i,3]+posteriors$alpha_InvForb*POCA[i,4]
  }
}

# calculate credible intervals 
HPDPOCA<-HPDinterval((as.mcmc(as.numeric(POCAF_nat))))
lowerPOCAnat<-HPDPOCA[,1]
upperPOCAnat<-HPDPOCA[,2] 

HPDPOCA<-HPDinterval((as.mcmc(as.numeric(POCAF_ex))))
lowerPOCAex<-HPDPOCA[,1]
upperPOCAex<-HPDPOCA[,2] 

remove(HPDPOCA)

POCAmeanF<-mean(POCAF, na.rm=T)
PocaL <- posteriors$lambda
POCAmeanL<-mean(posteriors$lambda)
HPDPOCA<-HPDinterval((as.mcmc(as.numeric(POCAF))))
HPDPOCAlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerPOCA<-HPDPOCA[,1]
upperPOCA<-HPDPOCA[,2] 
lowerPOCAlam<-HPDPOCAlam[,1]
upperPOCAlam<-HPDPOCAlam[,2]

remove(PrelimFit)

# TRCY
#load in data, subset out C treatment and columns of groups then remove NAs
TRCY<-read.csv("Data/groups.TRCY.csv")
TRCY<-subset(TRCY, treatment=="C")
TRCY<-subset(TRCY, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
TRCY<-as.matrix(TRCY)
# load in the posteriors 
load(here("Bayes data", "TRCY_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

TRCY_nat <- posteriors$alpha_NatForb
TRCY_ex <- c(posteriors$alpha_InvGraminoid, posteriors$InvForb)

dim(TRCY)
TRCYF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 51)
TRCYF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 51)
TRCYF <-matrix(NA, length(posteriors$alpha_NatForb), 51)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:51){
  if (TRCY[i,3]>0){
    TRCYF_nat[,i]<-posteriors$alpha_NatForb*TRCY[i,3]
  }
  if (TRCY[i,2]+TRCY[i,4]>0){
    TRCYF_ex[,i]<-posteriors$alpha_InvGraminoid*TRCY[i,2]+posteriors$alpha_InvForb*TRCY[i,4]
  }
  if (sum(TRCY[i,1:4]>0)){
    TRCYF[,i]<-posteriors$alpha_intra*TRCY[i,1]+posteriors$alpha_InvGraminoid*TRCY[i,2]+posteriors$alpha_NatForb*TRCY[i,3]+posteriors$alpha_InvForb*TRCY[i,4]
  }
}

# calculate credible intervals 
HPDTRCY<-HPDinterval((as.mcmc(as.numeric(TRCYF_nat))))
lowerTRCYnat<-HPDTRCY[,1]
upperTRCYnat<-HPDTRCY[,2] 

HPDTRCY<-HPDinterval((as.mcmc(as.numeric(TRCYF_ex))))
lowerTRCYex<-HPDTRCY[,1]
upperTRCYex<-HPDTRCY[,2] 

remove(HPDTRCY)

TRCYmeanF<-mean(TRCYF, na.rm=T)
TrcyL <- posteriors$lambda
TRCYmeanL<-mean(posteriors$lambda)
HPDTRCY<-HPDinterval((as.mcmc(as.numeric(TRCYF))))
HPDTRCYlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerTRCY<-HPDTRCY[,1]
upperTRCY<-HPDTRCY[,2] 
lowerTRCYlam<-HPDTRCYlam[,1]
upperTRCYlam<-HPDTRCYlam[,2]

remove(PrelimFit)

# VERO
#load in data, subset out C treatment and columns of groups then remove NAs
VERO<-read.csv("Data/groups.VERO.csv")
VERO<-subset(VERO, treatment=="C")
VERO<-subset(VERO, select=c(intra, exotic.grass, native.forb, exotic.forb, unknown))

# make a matrix of the abundances 
VERO<-as.matrix(VERO)
# load in the posteriors 
load(here("Bayes data", "VERO_posteriors.rdata"))
posteriors<-rstan::extract(PrelimFit)

dim(VERO)
VEROF_nat <-matrix(NA, length(posteriors$alpha_NatForb), 44)
VEROF_ex <-matrix(NA, length(posteriors$alpha_NatForb), 44)
VEROF <-matrix(NA, length(posteriors$alpha_NatForb), 44)
# calculate fhat/lambda for each posterior of intra and then not intra (all of ther alphas and other groups)
for (i in 1:44){
  #if (VERO[i,3]>0){
    VEROF_nat[,i]<-posteriors$alpha_NatForb*VERO[i,3] # present as ln(VeroF_nat)
  #}
  if (VERO[i,2]+VERO[i,4]>0){
    VEROF_ex[,i]<-posteriors$alpha_InvGraminoid*VERO[i,2]+posteriors$alpha_InvForb*VERO[i,4]
  }
  if (sum(VERO[i,1:4]>0)){
    VEROF[,i]<-posteriors$alpha_intra*VERO[i,1]+posteriors$alpha_InvGraminoid*VERO[i,2]+posteriors$alpha_NatForb*VERO[i,3]+posteriors$alpha_InvForb*VERO[i,4]
  }
}

log(VEROF_nat+1)
# need to log the net neighbourhood effect values 
#VEROF_nat <- log(VEROF_nat+1) #??


# calculate credible intervals 
HPDVERO<-HPDinterval((as.mcmc(as.numeric(VEROF_nat))))
lowerVEROnat<-HPDVERO[,1]
upperVEROnat<-HPDVERO[,2] 

HPDVERO<-HPDinterval((as.mcmc(as.numeric(VEROF_ex))))
lowerVEROex<-HPDVERO[,1]
upperVEROex<-HPDVERO[,2] 

remove(HPDVERO)

VEROmeanF<-mean(VEROF, na.rm=T)
VeroL <- posteriors$lambda
VEROmeanL<-mean(posteriors$lambda)
HPDVERO<-HPDinterval((as.mcmc(as.numeric(VEROF))))
HPDVEROlam<-HPDinterval((as.mcmc(as.numeric(posteriors$lambda))))

lowerVERO<-HPDVERO[,1]
upperVERO<-HPDVERO[,2] 
lowerVEROlam<-HPDVEROlam[,1]
upperVEROlam<-HPDVEROlam[,2]

remove(PrelimFit)


# try correlation coeffs 
lambdas <- c(DaglL, HyglL, GiteL, PldeL, PocaL, TrcyL, VeroL, ArcaL, MediL, MomoL, PeaiL)
nne.1 <- c(DAGLF_nat[,1], HYGLF_nat[,1], GITEF_nat[,1], PLDEF_nat[,1], POCAF_nat[,1], TRCYF_nat[,1], VEROF_nat[,1], ARCAF_nat[,1], MEDIF_nat[,1], MOMOF_nat[,1], PEAIF_nat[,1]) # per plot (would have to do up to min number of plots (41)) 
nne.2 <- c(DAGLF_nat[,3], HYGLF_nat[,3], GITEF_nat[,3], PLDEF_nat[,3], POCAF_nat[,3], TRCYF_nat[,3], VEROF_nat[,3], ARCAF_nat[,3], MEDIF_nat[,3], MOMOF_nat[,3], PEAIF_nat[,3])
nne.list <- list(nne.1, nne.2)
cor.list <- list()

for (p in 1:2){ # just do min number of plots? 41 - and only the ones without NA's
cor.list[p] <- cor(nne.list[[p]], lambdas, use="complete.obs", method = "spearman")
}
cor.list


