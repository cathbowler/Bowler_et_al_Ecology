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
  if (VERO[i,3]>0){
    VEROF_nat[,i]<-posteriors$alpha_NatForb*VERO[i,3]
  }
  if (VERO[i,2]+VERO[i,4]>0){
    VEROF_ex[,i]<-posteriors$alpha_InvGraminoid*VERO[i,2]+posteriors$alpha_InvForb*VERO[i,4]
  }
  if (sum(VERO[i,1:4]>0)){
    VEROF[,i]<-posteriors$alpha_intra*VERO[i,1]+posteriors$alpha_InvGraminoid*VERO[i,2]+posteriors$alpha_NatForb*VERO[i,3]+posteriors$alpha_InvForb*VERO[i,4]
  }
}

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



#### Calculate all alphas figure parameters ####
library(here)
library(coda)

sp <- 11

load(here("Bayes data/Arca_posteriors.rdata"))
arca <- rstan::extract(PrelimFit)

arca_intra <- HPDinterval(as.mcmc(as.numeric(arca$alpha_intra)))
arca_intra_mean <- mean(arca$alpha_intra)

arca_ig <- HPDinterval(as.mcmc(as.numeric(arca$alpha_InvGraminoid)))
arca_ig_mean <- mean(arca$alpha_InvGraminoid)

arca_nf <- HPDinterval(as.mcmc(as.numeric(arca$alpha_NatForb)))
arca_nf_mean <- mean(arca$alpha_NatForb)

arca_if <- HPDinterval(as.mcmc(as.numeric(arca$alpha_InvForb)))
arca_if_mean <- mean(arca$alpha_InvForb)
# arca_if <- rep(NA, 2)
# arca_if_mean <- NA

arca_mean <- c(arca_intra_mean, arca_nf_mean, arca_if_mean, arca_ig_mean)
arca_low <- c(arca_intra[1], arca_nf[1], arca_if[1], arca_ig[1])
arca_high <- c(arca_intra[2], arca_nf[2], arca_if[2], arca_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Dagl_posteriors.rdata"))
dagl <- rstan::extract(PrelimFit)

dagl_intra <- HPDinterval(as.mcmc(as.numeric(dagl$alpha_intra)))
dagl_intra_mean <- mean(dagl$alpha_intra)

dagl_ig <- HPDinterval(as.mcmc(as.numeric(dagl$alpha_InvGraminoid)))
dagl_ig_mean <- mean(dagl$alpha_InvGraminoid)

dagl_nf <- HPDinterval(as.mcmc(as.numeric(dagl$alpha_NatForb)))
dagl_nf_mean <- mean(dagl$alpha_NatForb)

dagl_if <- HPDinterval(as.mcmc(as.numeric(dagl$alpha_InvForb)))
dagl_if_mean <- mean(dagl$alpha_InvForb)

dagl_mean <- c(dagl_intra_mean, dagl_nf_mean, dagl_if_mean, dagl_ig_mean)
dagl_low <- c(dagl_intra[1], dagl_nf[1], dagl_if[1], dagl_ig[1])
dagl_high <- c(dagl_intra[2], dagl_nf[2], dagl_if[2], dagl_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Gite_posteriors.rdata"))
gite <- rstan::extract(PrelimFit)

gite_intra <- HPDinterval(as.mcmc(as.numeric(gite$alpha_intra)))
gite_intra_mean <- mean(gite$alpha_intra)

gite_ig <- HPDinterval(as.mcmc(as.numeric(gite$alpha_InvGraminoid)))
gite_ig_mean <- mean(gite$alpha_InvGraminoid)

gite_nf <- HPDinterval(as.mcmc(as.numeric(gite$alpha_NatForb)))
gite_nf_mean <- mean(gite$alpha_NatForb)

gite_if <- HPDinterval(as.mcmc(as.numeric(gite$alpha_InvForb)))
gite_if_mean <- mean(gite$alpha_InvForb)


gite_mean <- c(gite_intra_mean, gite_nf_mean, gite_if_mean, gite_ig_mean)
gite_low <- c(gite_intra[1], gite_nf[1], gite_if[1], gite_ig[1])
gite_high <- c(gite_intra[2], gite_nf[2], gite_if[2], gite_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Hygl_posteriors.rdata"))
hygl <- rstan::extract(PrelimFit)

hygl_intra <- HPDinterval(as.mcmc(as.numeric(hygl$alpha_intra)))
hygl_intra_mean <- mean(hygl$alpha_intra)

hygl_ig <- HPDinterval(as.mcmc(as.numeric(hygl$alpha_InvGraminoid)))
hygl_ig_mean <- mean(hygl$alpha_InvGraminoid)

hygl_nf <- HPDinterval(as.mcmc(as.numeric(hygl$alpha_NatForb)))
hygl_nf_mean <- mean(hygl$alpha_NatForb)

hygl_if <- HPDinterval(as.mcmc(as.numeric(hygl$alpha_InvForb)))
hygl_if_mean <- mean(hygl$alpha_InvForb)

#hygl_o <- rep(NA, 2)
#hygl_o_mean <- NA

hygl_mean <- c(hygl_intra_mean, hygl_nf_mean, hygl_if_mean, hygl_ig_mean)
hygl_low <- c(hygl_intra[1], hygl_nf[1], hygl_if[1], hygl_ig[1])
hygl_high <- c(hygl_intra[2], hygl_nf[2], hygl_if[2], hygl_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Medi_posteriors.rdata"))
medi <- rstan::extract(PrelimFit)

medi_intra <- HPDinterval(as.mcmc(as.numeric(medi$alpha_intra)))
medi_intra_mean <- mean(medi$alpha_intra)

medi_ig <- HPDinterval(as.mcmc(as.numeric(medi$alpha_InvGraminoid)))
medi_ig_mean <- mean(medi$alpha_InvGraminoid)

medi_nf <- HPDinterval(as.mcmc(as.numeric(medi$alpha_NatForb)))
medi_nf_mean <- mean(medi$alpha_NatForb)

medi_o <- rep(NA, 2)
medi_o_mean <- NA

medi_if <- HPDinterval(as.mcmc(as.numeric(medi$alpha_InvForb)))
medi_if_mean <- mean(medi$alpha_InvForb)

medi_mean <- c(medi_intra_mean, medi_nf_mean, medi_if_mean, medi_ig_mean)
medi_low <- c(medi_intra[1], medi_nf[1], medi_if[1], medi_ig[1])
medi_high <- c(medi_intra[2], medi_nf[2], medi_if[2], medi_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Momo_posteriors.rdata"))
momo <- rstan::extract(PrelimFit)

momo_intra <- HPDinterval(as.mcmc(as.numeric(momo$alpha_intra)))
momo_intra_mean <- mean(momo$alpha_intra)

momo_ig <- HPDinterval(as.mcmc(as.numeric(momo$alpha_InvGraminoid)))
momo_ig_mean <- mean(momo$alpha_InvGraminoid)

momo_nf <- HPDinterval(as.mcmc(as.numeric(momo$alpha_NatForb)))
momo_nf_mean <- mean(momo$alpha_NatForb)

momo_if <- HPDinterval(as.mcmc(as.numeric(momo$alpha_InvForb)))
momo_if_mean <- mean(momo$alpha_InvForb)

momo_mean <- c(momo_intra_mean, momo_nf_mean, momo_if_mean, momo_ig_mean)
momo_low <- c(momo_intra[1], momo_nf[1], momo_if[1], momo_ig[1])
momo_high <- c(momo_intra[2], momo_nf[2], momo_if[2], momo_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Peai_posteriors.rdata"))
peai <- rstan::extract(PrelimFit)

peai_intra <- HPDinterval(as.mcmc(as.numeric(peai$alpha_intra)))
peai_intra_mean <- mean(peai$alpha_intra)

peai_ig <- rep(NA, 2)
peai_ig_mean <- NA

peai_nf <- HPDinterval(as.mcmc(as.numeric(peai$alpha_NatForb)))
peai_nf_mean <- mean(peai$alpha_NatForb)

peai_if <- HPDinterval(as.mcmc(as.numeric(peai$alpha_InvForb)))
peai_if_mean <- mean(peai$alpha_InvForb)

peai_mean <- c(peai_intra_mean, peai_nf_mean, peai_if_mean, peai_ig_mean)
peai_low <- c(peai_intra[1], peai_nf[1], peai_if[1], peai_ig[1])
peai_high <- c(peai_intra[2], peai_nf[2], peai_if[2], peai_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Plde_posteriors.rdata"))
plde <- rstan::extract(PrelimFit)

plde_intra <- HPDinterval(as.mcmc(as.numeric(plde$alpha_intra)))
plde_intra_mean <- mean(plde$alpha_intra)

plde_ig <- HPDinterval(as.mcmc(as.numeric(plde$alpha_InvGraminoid)))
plde_ig_mean <- mean(plde$alpha_InvGraminoid)

plde_nf <- HPDinterval(as.mcmc(as.numeric(plde$alpha_NatForb)))
plde_nf_mean <- mean(plde$alpha_NatForb)

plde_if <- HPDinterval(as.mcmc(as.numeric(plde$alpha_InvForb)))
plde_if_mean <- mean(plde$alpha_InvForb)

#plde_o <- rep(NA,2)
#plde_o_mean <- NA

plde_mean <- c(plde_intra_mean, plde_nf_mean, plde_if_mean, plde_ig_mean)
plde_low <- c(plde_intra[1], plde_nf[1], plde_if[1], plde_ig[1])
plde_high <- c(plde_intra[2], plde_nf[2], plde_if[2], plde_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Poca_posteriors.rdata"))
poca <- rstan::extract(PrelimFit)

poca_intra <- HPDinterval(as.mcmc(as.numeric(poca$alpha_intra)))
poca_intra_mean <- mean(poca$alpha_intra)

poca_ig <- HPDinterval(as.mcmc(as.numeric(poca$alpha_InvGraminoid)))
poca_ig_mean <- mean(poca$alpha_InvGraminoid)

poca_nf <- HPDinterval(as.mcmc(as.numeric(poca$alpha_NatForb)))
poca_nf_mean <- mean(poca$alpha_NatForb)

poca_if <- HPDinterval(as.mcmc(as.numeric(poca$alpha_InvForb)))
poca_if_mean <- mean(poca$alpha_InvForb)

#poca_o <- rep(NA, 2)
#poca_o_mean <- NA

poca_mean <- c(poca_intra_mean, poca_nf_mean, poca_if_mean, poca_ig_mean)
poca_low <- c(poca_intra[1], poca_nf[1], poca_if[1], poca_ig[1])
poca_high <- c(poca_intra[2], poca_nf[2], poca_if[2], poca_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Trcy_posteriors.rdata"))
trcy <- rstan::extract(PrelimFit)

trcy_intra <- HPDinterval(as.mcmc(as.numeric(trcy$alpha_intra)))
trcy_intra_mean <- mean(trcy$alpha_intra)

trcy_ig <- HPDinterval(as.mcmc(as.numeric(trcy$alpha_InvGraminoid)))
trcy_ig_mean <- mean(trcy$alpha_InvGraminoid)

trcy_nf <- HPDinterval(as.mcmc(as.numeric(trcy$alpha_NatForb)))
trcy_nf_mean <- mean(trcy$alpha_NatForb)

trcy_if <- HPDinterval((as.mcmc(as.numeric(trcy$alpha_InvForb))))
trcy_if_mean <- mean(trcy$alpha_InvForb)

trcy_mean <- c(trcy_intra_mean, trcy_nf_mean, trcy_if_mean, trcy_ig_mean)
trcy_low <- c(trcy_intra[1], trcy_nf[1], trcy_if[1], trcy_ig[1])
trcy_high <- c(trcy_intra[2], trcy_nf[2], trcy_if[2], trcy_ig[2])

remove(PrelimFit)

load(here("Bayes data", "Vero_posteriors.rdata"))
vero <- rstan::extract(PrelimFit)

vero_intra <- HPDinterval(as.mcmc(as.numeric(vero$alpha_intra)))
vero_intra_mean <- mean(vero$alpha_intra)

vero_ig <- HPDinterval(as.mcmc(as.numeric(vero$alpha_InvGraminoid)))
vero_ig_mean <- mean(vero$alpha_InvGraminoid)

vero_nf <- HPDinterval(as.mcmc(as.numeric(vero$alpha_NatForb)))
vero_nf_mean <- mean(vero$alpha_NatForb)

vero_if <- HPDinterval(as.mcmc(as.numeric(vero$alpha_InvForb)))
vero_if_mean <- mean(vero$alpha_InvForb)

vero_mean <- c(vero_intra_mean, vero_nf_mean, vero_if_mean, vero_ig_mean)
vero_low <- c(vero_intra[1], vero_nf[1], vero_if[1], vero_ig[1])
vero_high <- c(vero_intra[2], vero_nf[2], vero_if[2], vero_ig[2])

remove(PrelimFit)

# for supp occurrance vs posterior width figure ####
ARCA.O <- ifelse(ARCA>0, 1,0)
DAGL.O <- ifelse(DAGL>0, 1,0)
GITE.O <- ifelse(GITE>0, 1,0)
HYGL.O <- ifelse(HYGL>0, 1,0)
MEDI.O <- ifelse(MEDI>0, 1,0)
MOMO.O <- ifelse(MOMO>0, 1,0)
PLDE.O <- ifelse(PLDE>0, 1,0)
POCA.O <- ifelse(POCA>0, 1,0)
PEAI.O <- ifelse(PEAI>0, 1,0)
TRCY.O <- ifelse(TRCY>0, 1,0)
VERO.O <- ifelse(VERO>0, 1,0)


pdf("Figures/posterior.width.supp.pdf")
par(mfrow=c(2,2), mar=c(4,4,2,2))
ex.forb.pw <- c(abs(arca_if[2]-arca_if[1]), abs(medi_if[2]-medi_if[1]), abs(momo_if[2]-momo_if[1]), abs(peai_if[2]-peai_if[1]),
                abs(dagl_if[2]-dagl_if[1]), abs(gite_if[2]-gite_if[1]), abs(hygl_if[2]-hygl_if[1]), abs(plde_if[2]-plde_if[1]),
                abs(poca_if[2]-poca_if[1]), abs(trcy_if[2]-trcy_if[1]), abs(vero_if[2]-vero_if[1]))
ex.forb.ab <- c(sum(ARCA.O[,4]), sum(MEDI.O[,4]), sum(MOMO.O[,4]), sum(PEAI.O[,3]), sum(DAGL.O[,4]), sum(GITE.O[,4]), sum(HYGL.O[,4]),
                sum(PLDE.O[,4]), sum(POCA.O[,4]), sum(TRCY.O[,4]), sum(VERO.O[,4]))
plot(ex.forb.ab, ex.forb.pw, ylab="Posterior width", xlab="")
mtext(side=3, "a) Neighbor: exotic forb", cex=0.8, line=1)

ex.grass.pw <- c(abs(arca_ig[2]-arca_ig[1]), abs(medi_ig[2]-medi_ig[1]), abs(momo_ig[2]-momo_ig[1]),
                abs(dagl_ig[2]-dagl_ig[1]), abs(gite_ig[2]-gite_ig[1]), abs(hygl_ig[2]-hygl_ig[1]), abs(plde_ig[2]-plde_ig[1]),
                abs(poca_ig[2]-poca_ig[1]), abs(trcy_ig[2]-trcy_ig[1]), abs(vero_ig[2]-vero_ig[1]))
ex.grass.ab <- c(sum(ARCA.O[,2]), sum(MEDI.O[,2]), sum(MOMO.O[,2]), sum(DAGL.O[,2]), sum(GITE.O[,2]), sum(HYGL.O[,2]),
                sum(PLDE.O[,2]), sum(POCA.O[,2]), sum(TRCY.O[,2]), sum(VERO.O[,2]))
plot(ex.grass.ab, ex.grass.pw, ylab="", xlab="")
mtext(side=3, "b) Neighbor: exotic grass", cex=0.8, line=1)

nat.forb.pw <- c(abs(arca_nf[2]-arca_nf[1]), abs(medi_nf[2]-medi_nf[1]), abs(momo_nf[2]-momo_nf[1]), abs(peai_nf[2]-peai_nf[1]),
                 abs(dagl_nf[2]-dagl_nf[1]), abs(gite_nf[2]-gite_nf[1]), abs(hygl_nf[2]-hygl_nf[1]), abs(plde_nf[2]-plde_nf[1]),
                 abs(poca_nf[2]-poca_nf[1]), abs(trcy_nf[2]-trcy_nf[1]), abs(vero_nf[2]-vero_nf[1]))
nat.forb.ab <- c(sum(ARCA.O[,3]), sum(MEDI.O[,3]), sum(MOMO.O[,3]), sum(PEAI.O[,2]), sum(DAGL.O[,3]), sum(GITE.O[,3]), sum(HYGL.O[,3]),
                 sum(PLDE.O[,3]), sum(POCA.O[,3]), sum(TRCY.O[,3]), sum(VERO.O[,3]))
plot(nat.forb.ab, nat.forb.pw, ylab="Posterior width", xlab="Data points")
mtext(side=3, "c) Neighbor: native forb", cex=0.8, line=1)

intra.pw <- c(abs(arca_intra[2]-arca_intra[1]), abs(medi_intra[2]-medi_intra[1]), abs(momo_intra[2]-momo_intra[1]), abs(peai_intra[2]-peai_intra[1]),
                 abs(dagl_intra[2]-dagl_intra[1]), abs(gite_intra[2]-gite_intra[1]), abs(hygl_intra[2]-hygl_intra[1]), abs(plde_intra[2]-plde_intra[1]),
                 abs(poca_intra[2]-poca_intra[1]), abs(trcy_intra[2]-trcy_intra[1]), abs(vero_intra[2]-vero_intra[1]))
intra.ab <- c(sum(ARCA.O[,1]), sum(MEDI.O[,1]), sum(MOMO.O[,1]), sum(PEAI.O[,1]), sum(DAGL.O[,1]), sum(GITE.O[,1]), sum(HYGL.O[,1]),
                 sum(PLDE.O[,1]), sum(POCA.O[,1]), sum(TRCY.O[,1]), sum(VERO.O[,1]))
plot(intra.ab, intra.pw, ylab="", xlab="Data points")
mtext(side=3, "d) Neighbor: conspecific", cex=0.8, line=1)
dev.off()








#### Fig 3 and 4 stuff ####
intra_alphas <- c(dagl_intra_mean, gite_intra_mean, hygl_intra_mean, plde_intra_mean, poca_intra_mean, 
                  trcy_intra_mean, vero_intra_mean, arca_intra_mean, medi_intra_mean, momo_intra_mean, 
                  peai_intra_mean)

nf_alphas<-c(dagl_nf_mean, gite_nf_mean, hygl_nf_mean, plde_nf_mean, poca_nf_mean, trcy_nf_mean, vero_nf_mean, arca_nf_mean, medi_nf_mean, momo_nf_mean, peai_nf_mean)

if_alphas<-c(dagl_if_mean, gite_if_mean, hygl_if_mean, plde_if_mean, poca_if_mean, trcy_if_mean, vero_if_mean, arca_if_mean, medi_if_mean, momo_if_mean, peai_if_mean)

ig_alphas<-c(dagl_ig_mean, gite_ig_mean, hygl_ig_mean, plde_ig_mean, poca_ig_mean, trcy_ig_mean, vero_if_mean, arca_ig_mean, medi_ig_mean, momo_ig_mean, peai_ig_mean)
focal.sp<-c("dagl", "gite", "hygl", "plde", "poca", "trcy", "vero", "arca", "medi", "momo", "poca")
alphas_together<-rbind(focal.sp, intra_alphas, nf_alphas, if_alphas, ig_alphas)
alphas_together<-as.data.frame(alphas_together)

Fmean<-c(DAGLmeanF, GITEmeanF, HYGLmeanF, PLDEmeanF, POCAmeanF, TRCYmeanF, VEROmeanF, ARCAmeanF, MEDImeanF, MOMOmeanF, PEAImeanF)

meanLambdas<-c(DAGLmeanL, GITEmeanL, HYGLmeanL,  PLDEmeanL, POCAmeanL, TRCYmeanL, VEROmeanL, ARCAmeanL, MEDImeanL, MOMOmeanL, PEAImeanL)

meanFnats <- c(mean(DAGLF_nat, na.rm = T), mean(GITEF_nat, na.rm = T), mean(HYGLF_nat, na.rm = T),
               mean(PLDEF_nat, na.rm = T), mean(POCAF_nat, na.rm = T), mean(TRCYF_nat, na.rm = T), 
               mean(VEROF_nat, na.rm = T), mean(ARCAF_nat, na.rm = T), mean(MEDIF_nat, na.rm = T), 
               mean(MOMOF_nat, na.rm = T), mean(PEAIF_nat, na.rm = T))

meanFexs <- c(mean(DAGLF_ex, na.rm = T), mean(GITEF_ex, na.rm = T), mean(HYGLF_ex, na.rm = T),
              mean(PLDEF_ex, na.rm = T), mean(POCAF_ex, na.rm = T), mean(TRCYF_ex, na.rm = T), 
              mean(VEROF_ex, na.rm = T), mean(ARCAF_ex, na.rm = T), mean(MEDIF_ex, na.rm = T), 
              mean(MOMOF_ex, na.rm = T), mean(PEAIF_ex, na.rm = T))


# Fnat lower HPD intervals
lower_nat<-c(lowerDAGLnat, lowerGITEnat, lowerHYGLnat, lowerPLDEnat, lowerPOCAnat,
             lowerTRCYnat, lowerVEROnat, lowerARCAnat, lowerMEDInat, lowerMOMOnat,
             lowerPEAInat)
# Fnat upper HPD intervals
upper_nat<-c(upperDAGLnat, upperGITEnat, upperHYGLnat, upperPLDEnat, upperPOCAnat,
             upperTRCYnat, upperVEROnat, upperARCAnat, upperMEDInat, upperMOMOnat,
             upperPEAInat)
interval_nat<-(upper_nat-lower_nat)/2

# same for Finters
lower_ex<-c(lowerDAGLex, lowerGITEex, lowerHYGLex, lowerPLDEex, lowerPOCAex,
            lowerTRCYex, lowerVEROex, lowerARCAex, lowerMEDIex, lowerMOMOex,
            lowerPEAIex)
upper_ex<-c(upperDAGLex, upperGITEex, upperHYGLex, upperPLDEex, upperPOCAex,
            upperTRCYex, upperVEROex, upperARCAex, upperMEDIex, upperMOMOex,
            upperPEAIex)
interval_ex<-(upper_ex-lower_ex)/2

# same but overall F
lowerF<-c(lowerDAGL, lowerGITE, lowerHYGL, lowerPLDE, lowerPOCA,
          lowerTRCY, lowerVERO, lowerARCA, lowerMEDI, lowerMOMO,
          lowerPEAI)
upperF<-c(upperDAGL, upperGITE, upperHYGL, upperPLDE, upperPOCA,
          upperTRCY, upperVERO, upperARCA, upperMEDI, upperMOMO,
          upperPEAI)
lowerlam<-c(lowerDAGLlam, lowerGITElam, lowerHYGLlam, lowerPLDElam, lowerPOCAlam,
            lowerTRCYlam, lowerVEROlam, lowerARCAlam, lowerMEDIlam, lowerMOMOlam,
            lowerPEAIlam)
upperlam<-c(upperDAGLlam, upperGITElam, upperHYGLlam, upperPLDElam, upperPOCAlam,
            upperTRCYlam, upperVEROlam, upperARCAlam, upperMEDIlam, upperMOMOlam,
            upperPEAIlam)
intervalF<-(upperF-lowerF)/2
intervalL<-(upperlam-lowerlam)/2

# for intra the inter neighbours split by native forbs and exotic (forbs and grass) 
intra_low <- c(dagl_low[1], gite_low[1], hygl_low[1], plde_low[1], poca_low[1], trcy_low[1], 
               vero_low[1], arca_low[1],  medi_low[1], momo_low[1], peai_low[1])
intra_high <- c(dagl_high[1], gite_high[1], hygl_high[1], plde_high[1], poca_high[1], trcy_high[1], 
                vero_high[1], arca_high[1],  medi_high[1], momo_high[1], peai_high[1])
intra_comp <- length(which(intra_low < 0 & intra_high < 0))/length(intra_low)
intra_fac <- length(which(intra_low > 0 & intra_high > 0))/length(intra_low)
intra_no <- 1-(intra_comp+intra_fac)

native_low <- c(dagl_low[2], gite_low[2], hygl_low[2], plde_low[2], poca_low[2], trcy_low[2], 
                vero_low[2], arca_low[2],  medi_low[2], momo_low[2], peai_low[2])
native_high <- c(dagl_high[2], gite_high[2], hygl_high[2], plde_high[2], poca_high[2], trcy_high[2], 
                 vero_high[2], arca_high[2],  medi_high[2], momo_high[2], peai_high[2])
native_comp <- length(which(native_low < 0 & native_high < 0))/length(native_low)
native_fac <- length(which(native_low > 0 & native_high > 0))/length(native_low)
native_no <- 1-(native_comp+native_fac)

exotic_low <- c(dagl_low[3:4], gite_low[3:4], hygl_low[3:4], plde_low[3:4], poca_low[3:4], 
                trcy_low[3:4], vero_low[3:4], arca_low[3:4],  medi_low[3:4], momo_low[3:4], 
                peai_low[3:4])
exotic_low <- exotic_low[is.na(exotic_low)==FALSE]
exotic_high <- c(dagl_high[3:4], gite_high[3:4], hygl_high[3:4], plde_high[3:4], poca_high[3:4], 
                 trcy_high[3:4], vero_high[3:4], arca_high[3:4],  medi_high[3:4], momo_high[3:4], 
                 peai_high[3:4])
exotic_high <- exotic_high[is.na(exotic_high)==FALSE]
exotic_comp <- length(which(exotic_low < 0 & exotic_high < 0))/length(exotic_low)
exotic_fac <- length(which(exotic_low > 0 & exotic_high > 0))/length(exotic_low)
exotic_no <- 1-(exotic_comp+exotic_fac)

intra <- c(intra_comp, intra_fac, intra_no)
native <- c(native_comp, native_fac, native_no)
exotic <- c(exotic_comp, exotic_fac, exotic_no)


# FOR PANEL B OF REALISED FECUNDITY PLOT
# percent realised fecundity greater than 1 (focal facilitated by neighbourhood)
x1 <- 0


 dagl.f <- ecdf(DAGLF) 
 dagl.perc <- 1-dagl.f(x1)
# 
 gite.f <- ecdf(GITEF)
 gite.perc <- 1-gite.f(x1)
# 
 hygl.f <- ecdf(HYGLF)
 hygl.perc <- 1-hygl.f(x1)
# 
 plde.f <- ecdf(PLDEF)
 plde.perc <- 1-plde.f(x1)
# 
 poca.f <- ecdf(POCAF)
 poca.perc <- 1-poca.f(x1)
# 
 trcy.f <- ecdf(TRCYF)
 trcy.perc <- 1-trcy.f(x1)
# 
 vero.f <- ecdf(VEROF)
 vero.perc <- 1-vero.f(x1)
# 
 arca.f <- ecdf(ARCAF)
 arca.perc <- 1-arca.f(x1)
# 
 medi.f <- ecdf(MEDIF)
 medi.perc <- 1-medi.f(x1)
# 
 momo.f <- ecdf(MOMOF)
 momo.perc <- 1-momo.f(x1)
# 
 peai.f <- ecdf(PEAIF)
 peai.perc <- 1-peai.f(x1)
# 

