###########################################################################################
# This code runs Bayes models (sourcing from stan models in "Bayes models" folder)
# for all focal species 
###########################################################################################
library("rstan")
# This line of code detects the number of cores available to your machine so that
#    chains can be run in parallel. If you will be doing other things on your 
#    computer while the model is running, I suggest not running this line of code
#    as it can slow down other applications on your computer
options(mc.cores = parallel::detectCores())

# This line allows a compiled model to be saved directly, so that stan does
#    not need to recompile it if it hasn't changed which can save time
rstan_options(auto_write = TRUE)


# ARCA ####

# Load in the data
SpData <- read.csv("Data/groups.ARCA.csv")
SpData <- subset(SpData, select = c(total.fecundity, plot, intra, exotic.grass, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

# Run the model 
initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit

save(PrelimFit, file = "Bayes data/Arca_posteriors.rdata")


# DAGL --------------------------------------------------
SpData <- read.csv("Data/groups.DAGL.csv")
SpData <- subset(SpData, select = c(plot, total.fecundity, intra, exotic.grass, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}


initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
#control = list(adapt_delta = 0.9
PrelimFit 

save(PrelimFit, file = "Bayes data/Dagl_posteriors.rdata")

# GITE -------------------------------------
SpData <- read.csv("Data/groups.Gite.csv")

SpData <- subset(SpData, select = c(plot, total.fecundity, intra, exotic.grass, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c( intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}


initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
#control = list(adapt_delta = 0.9
PrelimFit

save(PrelimFit, file = "Bayes data/Gite_posteriors.rdata")

# HYGL ------------------------------------
SpData <- read.csv("Data/groups.HYGL.csv")
SpData <- subset(SpData, select = c(plot, total.fecundity, intra, exotic.grass, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit

save(PrelimFit, file = "Bayes data/Hygl_posteriors.rdata")

# MEDI ---------------------------------------------
SpData <- read.csv("Data/groups.MEDI.csv")

SpData <- subset(SpData, select = c(plot, total.fecundity, intra, exotic.grass, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}


initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)

PrelimFit

save(PrelimFit, file = "Bayes data/Medi_posteriors.rdata")

# MOMO -----------------------------------

SpData <- read.csv("Data/groups.MOMO.csv")

SpData <- subset(SpData, select = c(plot, total.fecundity, intra, exotic.grass, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}


initials <- list(epsilon=rep(1,P), sigma = .01) 
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)

PrelimFit

save(PrelimFit, file = "Bayes data/Momo_posteriors.rdata")

# PEAI ------------------------------------

SpData <- read.csv("Data/groups.PEAI.csv")
SpData <- subset(SpData, select = c(plot, total.fecundity, intra, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(16)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}


initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = ("Stan models/FecundityModel_random2.stan"), data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)

PrelimFit

save(PrelimFit, file = "Bayes data/Peai_posteriors.rdata")


# PLDE -------------------------------------
PldeData <- read.csv("Data/groups.PLDE.csv")
PldeData <- subset(PldeData, select = c(total.fecundity, plot, intra, exotic.grass, exotic.forb, native.forb, unknown))
PldeData <- na.omit(PldeData)
N <- as.integer(dim(PldeData)[1])
P <- as.integer(16)
Fecundity <- as.integer(PldeData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(PldeData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(PldeData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}


initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)

PrelimFit
save(PrelimFit, file = "Bayes data/Plde_posteriors.rdata")


# POCA ---------------------------------------
PocaData <- read.csv("Data/groups.POCA.csv")
PocaData <- subset(PocaData, select = c(total.fecundity, plot, intra, exotic.grass, exotic.forb, native.forb, unknown))
PocaData <- na.omit(PocaData)
N <- as.integer(dim(PocaData)[1])
P <- as.integer(16)
Fecundity <- as.integer(PocaData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(PocaData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(PocaData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit
save(PrelimFit, file = "Bayes data/Poca_posteriors.rdata")

# TRCY ----------------------------------
SpData <- read.csv("Data/groups.TRCY.csv")

SpData <- subset(SpData, select = c(total.fecundity, plot, intra, exotic.grass, exotic.forb, native.forb, unknown))
SpData <- na.omit(SpData)
N <- as.integer(dim(SpData)[1])
P <- as.integer(24)
Fecundity <- as.integer(SpData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(SpData, select = c(intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(SpData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}


initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)

PrelimFit
save(PrelimFit, file = "Bayes data/Trcy_posteriors.rdata")


# VERO ------------------

VeroData <- read.csv("Data/groups.VERO.csv")
VeroData <- subset(VeroData, select = c(total.fecundity, plot, intra, exotic.grass, exotic.forb, native.forb, unknown))
VeroData <- na.omit(VeroData)
N <- as.integer(dim(VeroData)[1])
P <- as.integer(16)
Fecundity <- as.integer(VeroData$total.fecundity)

# Create the model matrix
ModelMatrix <- subset(VeroData, select = c( intra, exotic.grass, exotic.forb, native.forb, unknown))
ModelMatrix <- as.matrix(ModelMatrix)
PlotData <- subset(VeroData, select = plot)
Plot <- rep(NA, nrow(ModelMatrix))
for(i in 1:nrow(ModelMatrix)){
  Plot[i] <- PlotData[i,1]
}

initials <- list(epsilon=rep(1,P), sigma = .01)
initials1<- list(initials, initials, initials)

PrelimFit <- stan(file = 'Stan models/FecundityModel_random1.stan', data = c("Fecundity", "N", "ModelMatrix", "P", "Plot"),
                  iter = 6000, chains = 3, thin = 2, init = initials1)
PrelimFit
save(PrelimFit, file = "Bayes data/Vero_posteriors.rdata")









