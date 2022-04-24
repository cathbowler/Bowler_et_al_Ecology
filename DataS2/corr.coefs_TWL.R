#####################################################################
###### This script is for calculating correlation coefficients ######
#####################################################################

rm(list=ls())
#setwd("~/Desktop/Mentoring/Cath/Stochasticity_chapter-/Manuscript1 code/")
library(rstan)
library(coda)

# Make a function to load in and process the data for each species. This function
#  two versions of the species' acronym for its primary arguments: the version
#  used for the .csv files and the version used to save the posteriors. Then it
#  will also a take Boolean (TRUE/FALSE) argument specifying whether or not
#  neighborhoods with 0 individuals should be considered.
DataProcessing <- function(DataPrefix, PostPrefix, Status){
  # Load in and subset the neighborhood data
  GroupData <- read.csv(paste("Data/groups.", DataPrefix, ".csv", sep = ""))
  if(Status == "native"){
    GroupData_nat <- subset(GroupData, treatment == "C", select = c(native.forb, intra))
    GroupData_ex <- subset(GroupData, treatment == "C", select = c(exotic.grass, exotic.forb))
  }else{
    GroupData_nat <- subset(GroupData, treatment == "C", select = c(native.forb))
    GroupData_ex <- subset(GroupData, treatment == "C", select = c(exotic.grass, exotic.forb, intra))
  }
  GroupData_all <- subset(GroupData, treatment == "C", select = c(native.forb, exotic.grass, exotic.forb, intra))
  
  # Load and extract the model posteriors
  load(paste("Bayes data/", PostPrefix, "_posteriors.rdata", sep = ""))
  #load(here("Manuscript1 code/Bayes data", paste(PostPrefix, "posteriors.rdata", sep = "_")))
  posteriors <- rstan::extract(PrelimFit)
  # Save the lambda posterior and the posterior length for making the matrices below
  Lambda <- posteriors$lambda
  PostLength <- length(Lambda)
  
  # If we are excluding 0's, we need to first determine the number of non-zero
  # neighborhoods
  NonZero_nat <- which(rowSums(GroupData_nat) > 0)
  NonZero_ex <- which(rowSums(GroupData_ex) > 0)
  NonZero_all <- which(rowSums(GroupData_all) > 0)
  N_nat <- length(NonZero_nat)
  N_ex <- length(NonZero_ex)
  N_all <- length(NonZero_all)
    
  # Now calculate and save ln(F_hat/lambda) and for each neighborhood
  NetNeighbor_nat <- matrix(NA, nrow = PostLength, ncol = N_nat)
  NetNeighbor_ex <- matrix(NA, nrow = PostLength, ncol = N_ex)
  NetNeighbor_all <- matrix(NA, nrow = PostLength, ncol = N_all)
  
  if(Status == "native"){
    for(i in 1:max(N_nat, N_ex, N_all)){ # Loop through whichever is larger here
      if(i <= N_nat){
        NetNeighbor_nat[,i] <- posteriors$alpha_NatForb * GroupData_nat$native.forb[NonZero_nat[i]] + 
          posteriors$alpha_intra * GroupData_nat$intra[NonZero_nat[i]]
      }
      if(i <= N_ex){
        if(!is.null(posteriors$alpha_InvGraminoid)){
          NetNeighbor_ex[,i] <- posteriors$alpha_InvGraminoid * GroupData_ex$exotic.grass[NonZero_ex[i]] +
            posteriors$alpha_InvForb * GroupData_ex$exotic.forb[NonZero_ex[i]]
        }else{
          NetNeighbor_ex[,i] <- posteriors$alpha_InvForb * GroupData_ex$exotic.forb[NonZero_ex[i]]
        }
      }
      if(i <= N_all){
        if(!is.null(posteriors$alpha_InvGraminoid)){
          NetNeighbor_all[,i] <- posteriors$alpha_InvGraminoid * GroupData_all$exotic.grass[NonZero_all[i]] +
            posteriors$alpha_InvForb * GroupData_all$exotic.forb[NonZero_all[i]] + 
            posteriors$alpha_NatForb * GroupData_all$native.forb[NonZero_all[i]] +
            posteriors$alpha_intra * GroupData_all$intra[NonZero_all[i]]
        }else{
          NetNeighbor_all[,i] <- posteriors$alpha_InvForb * GroupData_all$exotic.forb[NonZero_all[i]] + 
            posteriors$alpha_NatForb * GroupData_all$native.forb[NonZero_all[i]] +
            posteriors$alpha_intra * GroupData_all$intra[NonZero_all[i]]
        }
      }
    }
  }else{
    for(i in 1:max(N_nat, N_ex, N_all)){ # Loop through whichever is larger here
      if(i <= N_nat){
        NetNeighbor_nat[,i] <- posteriors$alpha_NatForb * GroupData_nat$native.forb[NonZero_nat[i]]
      }
      if(i <= N_ex){
        if(!is.null(posteriors$alpha_InvGraminoid)){
          NetNeighbor_ex[,i] <- posteriors$alpha_InvGraminoid * GroupData_ex$exotic.grass[NonZero_ex[i]] +
            posteriors$alpha_InvForb * GroupData_ex$exotic.forb[NonZero_ex[i]] + 
            posteriors$alpha_intra * GroupData_ex$intra[NonZero_ex[i]]
        }else{
          NetNeighbor_ex[,i] <- posteriors$alpha_InvForb * GroupData_ex$exotic.forb[NonZero_ex[i]] + 
            posteriors$alpha_intra * GroupData_ex$intra[NonZero_ex[i]]
        }
      }
      if(i <= N_all){
        if(!is.null(posteriors$alpha_InvGraminoid)){
          NetNeighbor_all[,i] <- posteriors$alpha_InvGraminoid * GroupData_all$exotic.grass[NonZero_all[i]] +
            posteriors$alpha_InvForb * GroupData_all$exotic.forb[NonZero_all[i]] + 
            posteriors$alpha_NatForb * GroupData_all$native.forb[NonZero_all[i]] +
            posteriors$alpha_intra * GroupData_all$intra[NonZero_all[i]]
        }else{
          NetNeighbor_all[,i] <- posteriors$alpha_InvForb * GroupData_all$exotic.forb[NonZero_all[i]] + 
            posteriors$alpha_NatForb * GroupData_all$native.forb[NonZero_all[i]] +
            posteriors$alpha_intra * GroupData_all$intra[NonZero_all[i]]
        }
      }
    }
  }

  # Finally, return all the necessary values for the function
  Results <- list(N_ex = N_ex, N_nat = N_nat, N_all = N_all, NetNeighbor_ex = NetNeighbor_ex, 
                  NetNeighbor_nat = NetNeighbor_nat, NetNeighbor_all = NetNeighbor_all,
                  Lambda = Lambda, PostLength = PostLength)
  return(Results)
}

# Now use the above function to load in and process the data for all 11 species
S <- 11
SpeciesStatus <- c("exotic", "native", "native", "native", "exotic", "exotic", "exotic", "native",
                   "native", "native", "native")
AllDataPrefixes <- c("ARCA", "Dagl", "GITE", "HYGL", "MEDI", "MOMO", "PEAI", "PLDE",
                     "POCA", "TRCY", "VERO")
AllPostPrefixes <- c("Arca", "Dagl", "Gite", "Hygl", "Medi", "Momo", "Peai", "Plde",
                     "Poca", "Trcy", "Vero")
AllSpeciesData <- vector(mode = "list", length = S)
PostLengths <- rep(NA, S)
for(s in 1:S){
  AllSpeciesData[[s]] <- DataProcessing(DataPrefix = AllDataPrefixes[s], PostPrefix = AllPostPrefixes[s],
                                        Status = SpeciesStatus[s])
  PostLengths[s] <- AllSpeciesData[[s]]$PostLength
}

# Now make a function to calculate the correlation coefficient for varying
# combinations of focal and neighbor species
CalcCorCoef <- function(focals, neighbors, PostLength, AllSpeciesData){
  SpeciesStatus <- c("exotic", "native", "native", "native", "exotic", "exotic", "exotic", "native",
                     "native", "native", "native")
  S <- length(SpeciesStatus)
  LambdaVals <- NULL
  NeighborEffects <- NULL
  
  if(focals != "all"){
    CurSpecies <- which(SpeciesStatus == focals)
  }else{
    CurSpecies <- 1:S
  }
  CurPost <- sample(1:PostLength, size = S, replace = TRUE)
  for(s in CurSpecies){
    if(neighbors == "exotic"){
      LambdaVals <- c(LambdaVals, rep(AllSpeciesData[[s]]$Lambda[CurPost[s]], AllSpeciesData[[s]]$N_ex))
      NeighborEffects <- c(NeighborEffects, AllSpeciesData[[s]]$NetNeighbor_ex[CurPost[s],])
    } else if(neighbors == "native"){
      LambdaVals <- c(LambdaVals, rep(AllSpeciesData[[s]]$Lambda[CurPost[s]], AllSpeciesData[[s]]$N_nat))
      NeighborEffects <- c(NeighborEffects, AllSpeciesData[[s]]$NetNeighbor_nat[CurPost[s],])
    }else{
      LambdaVals <- c(LambdaVals, rep(AllSpeciesData[[s]]$Lambda[CurPost[s]], AllSpeciesData[[s]]$N_all))
      NeighborEffects <- c(NeighborEffects, AllSpeciesData[[s]]$NetNeighbor_all[CurPost[s],])
    }
  }
  CorCoef <- cor(x = LambdaVals, y = NeighborEffects, method = "spearman")
  return(CorCoef)
}


# Now loop through each entry of the posterior to calculate the corresponding
# Spearman's Rank Correlation Coefficient
PostLength <- min(PostLengths)
CorCoef <- array(NA, dim = c(3, 3, PostLength))
for(i in 1:PostLength){
  CorCoef[1,1,i] <- CalcCorCoef(focals = "native", neighbors = "native", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[1,2,i] <- CalcCorCoef(focals = "native", neighbors = "exotic", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[1,3,i] <- CalcCorCoef(focals = "native", neighbors = "all", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[2,1,i] <- CalcCorCoef(focals = "exotic", neighbors = "native", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[2,2,i] <- CalcCorCoef(focals = "exotic", neighbors = "exotic", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[2,3,i] <- CalcCorCoef(focals = "exotic", neighbors = "all", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[3,1,i] <- CalcCorCoef(focals = "all", neighbors = "native", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[3,2,i] <- CalcCorCoef(focals = "all", neighbors = "exotic", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
  CorCoef[3,3,i] <- CalcCorCoef(focals = "all", neighbors = "all", PostLength = PostLength, AllSpeciesData = AllSpeciesData)
}

CorMeans <- matrix(NA, nrow = 3, ncol = 3)
CorInts <- array(NA, dim = c(3,3,2))
for(i in 1:3){
  for(j in 1:3){
    CorMeans[i,j] <- mean(CorCoef[i,j,])
    CorInts[i,j,] <- HDInterval::hdi(CorCoef[i,j,])
  }
}

xRange <- range(CorCoef)
RowSeq <- c("Native focals", "Exotic focals", "All focals")
ColSeq <- c("Native neighbors", "Exotic neighbors", "All neighbors")
pdf(file = "Figures/CorCoefsDistributions.pdf", width = 8, height = 8, onefile = FALSE, paper = "special")
  par(mfrow = c(3,3), mar = c(1.25,1.25,1.25,1.25), oma = c(3,4,4,1))
  for(i in 1:3){
    for(j in 1:3){
      plot(density(CorCoef[i,j,]), xlab = "", ylab = "", main = "", xlim = xRange)
      abline(v = 0, lty = 2)
      Message <- paste(round(CorMeans[i,j], digits = 2), "(", round(CorInts[i,j,1], digits = 2),
                       " - ", round(CorInts[i,j,2], digits = 2), ")", sep = "")
      mtext(Message, side = 3, adj = 0.1, line = -1.5, cex = 0.75)
      if(j == 1){
        mtext(RowSeq[i], side = 2, line = 2.5)
      }
      if(i == 1){
        mtext(ColSeq[j], side = 3, line = 1)
      }
    }
  }
  mtext("Spearman's rank correlation coefficient", side = 1, outer = TRUE, line = 1.5)
dev.off()

# Make a single graph with all of them and 95% credible intervals
library(RColorBrewer)
FocalCols <- brewer.pal(n = 3, name = "Dark2")
xSeqNat <- seq(0.85, 2.85, by = 1)
xSeqEx <- 1:3
xSeqAll <- seq(1.15, 3.15, by = 1)
yRange <- c(-1, 0.5)
pdf(file = "Figures/CorCoefsAll.pdf", width = 8, height = 5, onefile = FALSE, paper = "special")
  # x axis shows neighbor groups, different points/colors for focals
  plot(NA, NA, xlim = c(0.5,3.5), ylim = yRange, las = 1, xaxt = "n",
       xlab = "", ylab = expression(paste("Spearman's ", rho, sep = "")))
  points(x = xSeqNat, y = CorMeans[1,], pch = 19, col = "orange")
  #segments(x0 = xSeqNat, y0 = CorInts[1,,1], x1 = xSeqNat, y1 = CorInts[1,,2], col = FocalCols[1])
  arrows(xSeqNat, y0 = CorInts[1,,1], xSeqNat, y1 = CorInts[1,,2], col = "orange", code=3, length=0.05, angle=90)
  
  points(x = xSeqEx, y = CorMeans[2,], pch = 15, col = "skyblue")
  #segments(x0 = xSeqEx, y0 = CorInts[2,,1], x1 = xSeqEx, y1 = CorInts[2,,2], col = FocalCols[2])
  arrows(xSeqEx, y0 = CorInts[2,,1], xSeqEx, y1 = CorInts[2,,2], col = "skyblue", code=3, length=0.05, angle=90)
  
  points(x = xSeqAll, y = CorMeans[3,], pch = 17, col = "black")
  #segments(x0 = xSeqAll, y0 = CorInts[3,,1], x1 = xSeqAll, y1 = CorInts[3,,2], col = FocalCols[3])
  arrows(xSeqAll, y0 = CorInts[3,,1], xSeqAll, y1 = CorInts[3,,2], col = "black", code=3, length=0.05, angle=90)
  
  abline(h = 0, lty = 2)
  #abline(v = c(1.5, 2.5), lty = 1, col = "grey")
  
  mtext("Types of Neighbours", side = 1, line = 3)
  axis(1, at = 1:3, labels = c("Native", "Exotic", "All"))
  legend("topleft", horiz = FALSE, legend = c("Native focals", "Exotic focals", "All focals"), col = c("orange", "skyblue", "black"), pch = c(19, 15, 17),
         bty = "n")
dev.off()


