# load in the environmental data and subset to the most prominant environmental parameters 
env <- read.csv("raw data/environmental.data.complete.csv")
head(env)
env.sub <- subset(env, select=c("plot", "focal.species", "P.mg.kg", "canopy.cover.percentage", "litter.cover.percentage"))

# DAGL
env.sub.dagl <- subset(env.sub, focal.species=="Daucus glochidiatus")

# load in plot level effects from Bayes models 
load("Bayes data/DAGL_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
dagl.plot.eff <- as.data.frame(posteriors$epsilon)# do i need to +1?
dagl.plot.eff <- log(dagl.plot.eff)

int.low.dagl <- rep(NA,16)
int.up.dagl <- rep(NA, 16)
for (i in 1:length(dagl.plot.eff)){
  int.low.dagl[i] <- HPDinterval((as.mcmc(as.numeric(dagl.plot.eff[,i]))))[,1]
  int.up.dagl[i] <- HPDinterval((as.mcmc(as.numeric(dagl.plot.eff[,i]))))[,2]
}
int.dagl <- (int.up.dagl-int.low.dagl)/2

dagl.plot.eff <- colMeans(dagl.plot.eff)
dagl.env <- cbind(env.sub.dagl,dagl.plot.eff)

# HYGL
env.sub.hygl <- subset(env.sub, focal.species=="Hyalosperma glutinosum")

# load in plot level effects from Bayes models 
load("Bayes data/HYGL_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
hygl.plot.eff <- as.data.frame(posteriors$epsilon)
hygl.plot.eff <- log(hygl.plot.eff)

int.low.hygl <- rep(NA,16)
int.up.hygl <- rep(NA, 16)
for (i in 1:length(hygl.plot.eff)){
  int.low.hygl[i] <- HPDinterval((as.mcmc(as.numeric(hygl.plot.eff[,i]))))[,1]
  int.up.hygl[i] <- HPDinterval((as.mcmc(as.numeric(hygl.plot.eff[,i]))))[,2]
}
int.hygl <- (int.up.hygl-int.low.hygl)/2

hygl.plot.eff <- colMeans(hygl.plot.eff)
hygl.env <- cbind(env.sub.hygl,hygl.plot.eff)

# GITE
env.sub.gite <- subset(env.sub, focal.species=="Gilberta tenuifolia")

# load in plot level effects from Bayes models 
load("Bayes data/GITE_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
gite.plot.eff <- as.data.frame(posteriors$epsilon)
gite.plot.eff <- log(gite.plot.eff)

int.low.gite <- rep(NA,16)
int.up.gite <- rep(NA, 16)
for (i in 1:length(gite.plot.eff)){
  int.low.gite[i] <- HPDinterval((as.mcmc(as.numeric(gite.plot.eff[,i]))))[,1]
  int.up.gite[i] <- HPDinterval((as.mcmc(as.numeric(gite.plot.eff[,i]))))[,2]
}
int.gite <- (int.up.gite-int.low.gite)/2

gite.plot.eff <- colMeans(gite.plot.eff)
gite.env <- cbind(env.sub.gite,gite.plot.eff)

# PLDE
env.sub.plde <- subset(env.sub, focal.species=="Plantago debilis")

# load in plot level effects from Bayes models 
load("Bayes data/PLDE_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
plde.plot.eff <- as.data.frame(posteriors$epsilon)
plde.plot.eff <- log(plde.plot.eff)

int.low.plde <- rep(NA,16)
int.up.plde <- rep(NA, 16)
for (i in 1:length(plde.plot.eff)){
  int.low.plde[i] <- HPDinterval((as.mcmc(as.numeric(plde.plot.eff[,i]))))[,1]
  int.up.plde[i] <- HPDinterval((as.mcmc(as.numeric(plde.plot.eff[,i]))))[,2]
}
int.plde <- (int.up.plde-int.low.plde)/2

plde.plot.eff <- colMeans(plde.plot.eff)
plde.env <- cbind(env.sub.plde,plde.plot.eff)

# POCA
env.sub.poca <- subset(env.sub, focal.species=="Podolepis canescens")

# load in plot level effects from Bayes models 
load("Bayes data/POCA_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
poca.plot.eff <- as.data.frame(posteriors$epsilon)
poca.plot.eff <- log(poca.plot.eff)

int.low.poca <- rep(NA,16)
int.up.poca <- rep(NA, 16)
for (i in 1:length(poca.plot.eff)){
  int.low.poca[i] <- HPDinterval((as.mcmc(as.numeric(poca.plot.eff[,i]))))[,1]
  int.up.poca[i] <- HPDinterval((as.mcmc(as.numeric(poca.plot.eff[,i]))))[,2]
}
int.poca <- (int.up.poca-int.low.poca)/2

poca.plot.eff <- colMeans(poca.plot.eff)
poca.env <- cbind(env.sub.poca,poca.plot.eff)

# TRCY
env.sub.trcy <- subset(env.sub, focal.species=="Trachymene cyanopetala")

# load in plot level effects from Bayes models 
load("Bayes data/TRCY_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
trcy.plot.eff <- as.data.frame(posteriors$epsilon)
trcy.plot.eff <- log(trcy.plot.eff)

int.low.trcy <- rep(NA,24)
int.up.trcy <- rep(NA, 24)
for (i in 1:length(trcy.plot.eff)){
  int.low.trcy[i] <- HPDinterval((as.mcmc(as.numeric(trcy.plot.eff[,i]))))[,1]
  int.up.trcy[i] <- HPDinterval((as.mcmc(as.numeric(trcy.plot.eff[,i]))))[,2]
}
int.trcy <- (int.up.trcy-int.low.trcy)/2

trcy.plot.eff <- colMeans(trcy.plot.eff)
trcy.env <- cbind(env.sub.trcy,trcy.plot.eff)


# VERO
env.sub.vero <- subset(env.sub, focal.species=="Velleia rosea")

# load in plot level effects from Bayes models 
load("Bayes data/VERO_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
vero.plot.eff <- as.data.frame(posteriors$epsilon)
vero.plot.eff <- log(vero.plot.eff)

int.low.vero <- rep(NA,16)
int.up.vero <- rep(NA, 16)
for (i in 1:length(vero.plot.eff)){
  int.low.vero[i] <- HPDinterval((as.mcmc(as.numeric(vero.plot.eff[,i]))))[,1]
  int.up.vero[i] <- HPDinterval((as.mcmc(as.numeric(vero.plot.eff[,i]))))[,2]
}
int.vero <- (int.up.vero-int.low.vero)/2

vero.plot.eff <- colMeans(vero.plot.eff)
vero.env <- cbind(env.sub.vero,vero.plot.eff)


# ARCA
env.sub.arca <- subset(env.sub, focal.species=="Arctotheca calendula")

# load in plot level effects from Bayes models 
load("Bayes data/ARCA_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
arca.plot.eff <- as.data.frame(posteriors$epsilon)
arca.plot.eff <- log(arca.plot.eff)

int.low.arca <- rep(NA,16)
int.up.arca <- rep(NA, 16)
for (i in 1:length(arca.plot.eff)){
  int.low.arca[i] <- HPDinterval((as.mcmc(as.numeric(arca.plot.eff[,i]))))[,1]
  int.up.arca[i] <- HPDinterval((as.mcmc(as.numeric(arca.plot.eff[,i]))))[,2]
}
int.arca <- (int.up.arca-int.low.arca)/2

arca.plot.eff <- colMeans(arca.plot.eff)
arca.env <- cbind(env.sub.arca,arca.plot.eff)


# MEDI
env.sub.medi <- subset(env.sub, focal.species=="Medicago minima")

# load in plot level effects from Bayes models 
load("Bayes data/MEDI_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
medi.plot.eff <- as.data.frame(posteriors$epsilon)
medi.plot.eff <- log(medi.plot.eff)

int.low.medi <- rep(NA,16)
int.up.medi <- rep(NA, 16)
for (i in 1:length(medi.plot.eff)){
  int.low.medi[i] <- HPDinterval((as.mcmc(as.numeric(medi.plot.eff[,i]))))[,1]
  int.up.medi[i] <- HPDinterval((as.mcmc(as.numeric(medi.plot.eff[,i]))))[,2]
}
int.medi <- (int.up.medi-int.low.medi)/2

medi.plot.eff <- colMeans(medi.plot.eff)
medi.env <- cbind(env.sub.medi,medi.plot.eff)

# MOMO
env.sub.momo <- subset(env.sub, focal.species=="Monoculus monstrosus")

# load in plot level effects from Bayes models 
load("Bayes data/MOMO_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
momo.plot.eff <- as.data.frame(posteriors$epsilon)
momo.plot.eff <- log(momo.plot.eff)

int.low.momo <- rep(NA,16)
int.up.momo <- rep(NA, 16)
for (i in 1:length(momo.plot.eff)){
  int.low.momo[i] <- HPDinterval((as.mcmc(as.numeric(momo.plot.eff[,i]))))[,1]
  int.up.momo[i] <- HPDinterval((as.mcmc(as.numeric(momo.plot.eff[,i]))))[,2]
}
int.momo <- (int.up.momo-int.low.momo)/2

momo.plot.eff <- colMeans(momo.plot.eff)
momo.env <- cbind(env.sub.momo,momo.plot.eff)

# PEAI
env.sub.peai <- subset(env.sub, focal.species=="Pentameris airoides")

# load in plot level effects from Bayes models 
load("Bayes data/PEAI_posteriors.rdata")
posteriors<-rstan::extract(PrelimFit)
library(coda)
library(tidyverse)
peai.plot.eff <- as.data.frame(posteriors$epsilon)
peai.plot.eff <- log(peai.plot.eff)

int.low.peai <- rep(NA,16)
int.up.peai <- rep(NA, 16)
for (i in 1:length(peai.plot.eff)){
  int.low.peai[i] <- HPDinterval((as.mcmc(as.numeric(peai.plot.eff[,i]))))[,1]
  int.up.peai[i] <- HPDinterval((as.mcmc(as.numeric(peai.plot.eff[,i]))))[,2]
}
int.peai <- (int.up.peai-int.low.peai)/2

peai.plot.eff <- colMeans(peai.plot.eff)
peai.env <- cbind(env.sub.peai,peai.plot.eff)



# MAKE THE FIGURES
pdf("Figures/Supp.env.nat.pdf") #####
par(mfrow=c(7,3))
par(mar=c(1.5,1.5,1.5,1.5), oma=c(4,4,4,3))
# dagl
plot(dagl.env$P.mg.kg, dagl.env$dagl.plot.eff, ylim=c(min(int.low.dagl), max(int.up.dagl)), ylab="", xlab="", pch=16,
     arrows(dagl.env$P.mg.kg, int.up.dagl, dagl.env$P.mg.kg, int.low.dagl, length=0.05, angle=90, code=3))
mtext(expression(italic("Daucus glochidiatus")), side=3, line=0, cex = 0.7, adj = 0)

plot(dagl.env$canopy.cover.percentage, dagl.env$dagl.plot.eff, ylim=c(min(int.low.dagl), max(int.up.dagl)), ylab="", xlab="", pch=16,
     arrows(dagl.env$canopy.cover.percentage, int.up.dagl, dagl.env$canopy.cover.percentage, int.low.dagl, length=0.05, angle=90, code=3))

plot(dagl.env$litter.cover.percentage, dagl.env$dagl.plot.eff,  ylim=c(min(int.low.dagl), max(int.up.dagl)), ylab="", xlab="", pch=16,
     arrows(dagl.env$litter.cover.percentage, int.up.dagl, dagl.env$litter.cover.percentage, int.low.dagl, length=0.05, angle=90, code=3))

# gite
plot(gite.env$P.mg.kg, gite.env$gite.plot.eff,  ylim=c(min(int.low.gite), max(int.up.gite)), ylab="", xlab="", pch=16,
     arrows(gite.env$P.mg.kg, int.up.gite, gite.env$P.mg.kg, int.low.gite, length=0.05, angle=90, code=3))
mtext(expression(italic("Gilberta tenuifolia")), side=3, line=0, cex = 0.7, adj = 0)

plot(gite.env$canopy.cover.percentage, gite.env$gite.plot.eff,  ylim=c(min(int.low.gite), max(int.up.gite)), ylab="", xlab="", pch=16,
     arrows(gite.env$canopy.cover.percentage, int.up.gite, gite.env$canopy.cover.percentage, int.low.gite, length=0.05, angle=90, code=3))

plot(gite.env$litter.cover.percentage, gite.env$gite.plot.eff, ylim=c(min(int.low.gite), max(int.up.gite)), ylab="", xlab="", pch=16,
     arrows(gite.env$litter.cover.percentage, int.up.gite, gite.env$litter.cover.percentage, int.low.gite, length=0.05, angle=90, code=3))

# hygl
plot(hygl.env$P.mg.kg, hygl.env$hygl.plot.eff, ylim=c(min(int.low.hygl), max(int.up.hygl)), ylab="", xlab="", pch=16,
     arrows(hygl.env$P.mg.kg, int.up.hygl, hygl.env$P.mg.kg, int.low.hygl, length=0.05, angle=90, code=3))
mtext(expression(italic("Hyalosperma glutinosum")), side=3, line=0, cex = 0.7, adj = 0)

plot(hygl.env$canopy.cover.percentage, hygl.env$hygl.plot.eff, ylim=c(min(int.low.hygl), max(int.up.hygl)), ylab="", xlab="", pch=16,
     arrows(hygl.env$canopy.cover.percentage, int.up.hygl, hygl.env$canopy.cover.percentage, int.low.hygl, length=0.05, angle=90, code=3))

plot(hygl.env$litter.cover.percentage, hygl.env$hygl.plot.eff, ylim=c(min(int.low.hygl), max(int.up.hygl)), ylab="", xlab="", pch=16,
     arrows(hygl.env$litter.cover.percentage, int.up.hygl, hygl.env$litter.cover.percentage, int.low.hygl, length=0.05, angle=90, code=3))

# plde
plot(plde.env$P.mg.kg, plde.env$plde.plot.eff, ylim=c(min(int.low.plde), max(int.up.plde)), ylab="", xlab="", pch=16,
     arrows(plde.env$P.mg.kg, int.up.plde, plde.env$P.mg.kg, int.low.plde, length=0.05, angle=90, code=3))
mtext(expression(italic("Plantago debilis")), side=3, line=0, cex = 0.7, adj = 0)

plot(plde.env$canopy.cover.percentage, plde.env$plde.plot.eff, ylim=c(min(int.low.plde), max(int.up.plde)), ylab="", xlab="", pch=16,
     arrows(plde.env$canopy.cover.percentage, int.up.plde, plde.env$canopy.cover.percentage, int.low.plde, length=0.05, angle=90, code=3))

plot(plde.env$litter.cover.percentage, plde.env$plde.plot.eff, ylim=c(min(int.low.plde), max(int.up.plde)), ylab="", xlab="", pch=16,
     arrows(plde.env$litter.cover.percentage, int.up.plde, plde.env$litter.cover.percentage, int.low.plde, length=0.05, angle=90, code=3))

# poca
plot(poca.env$P.mg.kg, poca.env$poca.plot.eff, ylim=c(min(int.low.poca), max(int.up.poca)), ylab="", xlab="", pch=16,
     arrows(poca.env$P.mg.kg, int.up.poca, poca.env$P.mg.kg, int.low.poca, length=0.05, angle=90, code=3))
mtext(expression(italic("Podalepis canescens")), side=3, line=0, cex = 0.7, adj = 0)

plot(poca.env$canopy.cover.percentage, poca.env$poca.plot.eff, ylim=c(min(int.low.poca), max(int.up.poca)), ylab="", xlab="", pch=16,
     arrows(poca.env$canopy.cover.percentage, int.up.poca, poca.env$canopy.cover.percentage, int.low.poca, length=0.05, angle=90, code=3))

plot(poca.env$litter.cover.percentage, poca.env$poca.plot.eff,  ylim=c(min(int.low.poca), max(int.up.poca)), ylab="", xlab="", pch=16,
     arrows(poca.env$litter.cover.percentage, int.up.poca, poca.env$litter.cover.percentage, int.low.poca, length=0.05, angle=90, code=3))

# trcy
plot(trcy.env$P.mg.kg, trcy.env$trcy.plot.eff,  ylim=c(min(int.low.trcy), max(int.up.trcy)), ylab="", xlab="", pch=16,
     arrows(trcy.env$P.mg.kg, int.up.trcy, trcy.env$P.mg.kg, int.low.trcy, length=0.05, angle=90, code=3))
mtext(expression(italic("Trachymene cyanopetala")), side=3, line=0, cex = 0.7, adj = 0)

plot(trcy.env$canopy.cover.percentage, trcy.env$trcy.plot.eff,  ylim=c(min(int.low.trcy), max(int.up.trcy)), ylab="", xlab="", pch=16,
     arrows(trcy.env$canopy.cover.percentage, int.up.trcy, trcy.env$canopy.cover.percentage, int.low.trcy, length=0.05, angle=90, code=3))

plot(trcy.env$litter.cover.percentage, trcy.env$trcy.plot.eff,  ylim=c(min(int.low.trcy), max(int.up.trcy)), ylab="", xlab="", pch=16,
     arrows(trcy.env$litter.cover.percentage, int.up.trcy, trcy.env$litter.cover.percentage, int.low.trcy, length=0.05, angle=90, code=3))

# vero 
plot(vero.env$P.mg.kg, vero.env$vero.plot.eff,  ylim=c(min(int.low.vero), max(int.up.vero)), ylab="", xlab="", pch=16,
     arrows(vero.env$P.mg.kg, int.up.vero, vero.env$P.mg.kg, int.low.vero, length=0.05, angle=90, code=3))
mtext(side=1, line=3, "Phosphorus mg/kg%")
mtext(side=2, line=3, adj=-4, "ln(Plot Effect)") #??
mtext(expression(italic("Velleia rosea")), side=3, line=0, cex = 0.7, adj = 0)
plot(vero.env$canopy.cover.percentage, vero.env$vero.plot.eff, ylim=c(min(int.low.vero), max(int.up.vero)), ylab="", xlab="", pch=16,
     arrows(vero.env$canopy.cover.percentage, int.up.vero, vero.env$canopy.cover.percentage, int.low.vero, length=0.05, angle=90, code=3))
mtext(side=1, line=3, "Canopy cover %")
plot(vero.env$litter.cover.percentage, vero.env$vero.plot.eff, ylim=c(min(int.low.vero), max(int.up.vero)), ylab="", xlab="", pch=16,
     arrows(vero.env$litter.cover.percentage, int.up.vero, vero.env$litter.cover.percentage, int.low.vero, length=0.05, angle=90, code=3))
mtext(side=1, line=3, "Litter cover %")
dev.off()


pdf("Figures/Supp.env.ex.pdf") ##### 
par(mfrow=c(4,3))
par(mar=c(1.5,1.5,1.5,1.5), oma=c(3,3,3,3))
# arca 
plot(arca.env$P.mg.kg, arca.env$arca.plot.eff, ylim=c(min(int.low.arca), max(int.up.arca)), ylab="", xlab="", pch=16,
     arrows(arca.env$P.mg.kg, int.up.arca, arca.env$P.mg.kg, int.low.arca, length=0.05, angle=90, code=3))
mtext(expression(italic("Arctotheca calendula")), side=3, line=0, cex = 0.7, adj = 0)

plot(arca.env$canopy.cover.percentage, arca.env$arca.plot.eff, ylim=c(min(int.low.arca), max(int.up.arca)), ylab="", xlab="", pch=16,
     arrows(arca.env$canopy.cover.percentage, int.up.arca, arca.env$canopy.cover.percentage, int.low.arca, length=0.05, angle=90, code=3))

plot(arca.env$litter.cover.percentage, arca.env$arca.plot.eff, ylim=c(min(int.low.arca), max(int.up.arca)), ylab="", xlab="", pch=16,
     arrows(arca.env$litter.cover.percentage, int.up.arca, arca.env$litter.cover.percentage, int.low.arca, length=0.05, angle=90, code=3))

# medi
plot(medi.env$P.mg.kg, medi.env$medi.plot.eff, ylim=c(min(int.low.medi), max(int.up.medi)), ylab="", xlab="", pch=16,
     arrows(medi.env$P.mg.kg, int.up.medi, medi.env$P.mg.kg, int.low.medi, length=0.05, angle=90, code=3))
mtext(expression(italic("Medicago minima")), side=3, line=0, cex = 0.7, adj = 0)
mtext(side=2, line=3,"ln(Plot effect)", adj = -5) #??

plot(medi.env$canopy.cover.percentage, medi.env$medi.plot.eff, ylim=c(min(int.low.medi), max(int.up.medi)), ylab="", xlab="", pch=16,
     arrows(medi.env$canopy.cover.percentage, int.up.medi, medi.env$canopy.cover.percentage, int.low.medi, length=0.05, angle=90, code=3))

plot(medi.env$litter.cover.percentage, medi.env$medi.plot.eff,  ylim=c(min(int.low.medi), max(int.up.medi)), ylab="", xlab="", pch=16,
     arrows(medi.env$litter.cover.percentage, int.up.medi, medi.env$litter.cover.percentage, int.low.medi, length=0.05, angle=90, code=3))

# momo
plot(momo.env$P.mg.kg, momo.env$momo.plot.eff, ylim=c(min(int.low.momo), max(int.up.momo)), ylab="", xlab="", pch=16,
     arrows(momo.env$P.mg.kg, int.up.momo, momo.env$P.mg.kg, int.low.momo, length=0.05, angle=90, code=3))
mtext(expression(italic("Monoculus monstrosus")), side=3, line=0, cex = 0.7, adj = 0)

plot(momo.env$canopy.cover.percentage, momo.env$momo.plot.eff, ylim=c(min(int.low.momo), max(int.up.momo)), ylab="", xlab="", pch=16,
     arrows(momo.env$canopy.cover.percentage, int.up.momo, momo.env$canopy.cover.percentage, int.low.momo, length=0.05, angle=90, code=3))

plot(momo.env$litter.cover.percentage, momo.env$momo.plot.eff, ylim=c(min(int.low.momo), max(int.up.momo)), ylab="", xlab="", pch=16,
     arrows(momo.env$litter.cover.percentage, int.up.momo, momo.env$litter.cover.percentage, int.low.momo, length=0.05, angle=90, code=3))

# peai
plot(peai.env$P.mg.kg, peai.env$peai.plot.eff, ylim=c(min(int.low.peai), max(int.up.peai)), ylab="", xlab="", pch=16,
     arrows(peai.env$P.mg.kg, int.up.peai, peai.env$P.mg.kg, int.low.peai, length=0.05, angle=90, code=3))
mtext(side=1, line=3, "Phosphorus mg/kg")
mtext(expression(italic("Pentameris airoides")), side=3, line=0, cex = 0.7, adj = 0)
plot(peai.env$canopy.cover.percentage, peai.env$peai.plot.eff, ylim=c(min(int.low.peai), max(int.up.peai)), ylab="", xlab="", pch=16,
     arrows(peai.env$canopy.cover.percentage, int.up.peai, peai.env$canopy.cover.percentage, int.low.peai, length=0.05, angle=90, code=3))
mtext(side=1, line=3, "Canopy cover %")
plot(peai.env$litter.cover.percentage, peai.env$peai.plot.eff, ylim=c(min(int.low.peai), max(int.up.peai)), ylab="", xlab="", pch=16,
     arrows(peai.env$litter.cover.percentage, int.up.peai, peai.env$litter.cover.percentage, int.low.peai, length=0.05, angle=90, code=3))
mtext(side=1, line=3, "Litter cover %")

dev.off()




