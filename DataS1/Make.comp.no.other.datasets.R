# Code to convert the community composition data from long form to wide form 
#        and then join it with the fecundity data for each species 
#        and then group neighbours into categories (without 'other' category)
library(tidyverse)

rm(list=ls())
# the latest data on RDM
Fecundity <- read.csv("raw data/Fecundity_2018.csv") 
Community <- read.csv("raw data/Community_composition.csv") 
Community <- na.omit(Community)
env <- read.csv("raw data/environmental.data.complete.csv")
env.sub <- subset(env, select=c("plot", "focal.species", "P.mg.kg", "canopy.cover.percentage", "litter.cover.percentage"))

# remove notes columns
f <- subset(Fecundity, select = c(focal.species, plot, treatment, subplot, mature, immature))
comm <- subset(Community, select = c(focal.species, plot, treatment, subplot, neighbour.species,
                                          neighbour.abundance))

# remove rows where mature and immature both have NA's
library(data.table)
f=data.table(f)
setkeyv(f, c("mature", "immature"))
dontwant <- f[J(NA, NA)] 
f <- f[!J(NA, NA)]  # Pull everything mature=na and immature=na
f=as.data.frame(f)

# then change the NA when its only in one of these to a zero (hence not losing these rows completely by initially removing NAs)
f$mature[is.na(f$mature)] <- 0
f$immature[is.na(f$immature)] <- 0


f <- tbl_df(f)
comm <- tbl_df(comm)
aa <- comm %>% group_by(subplot) %>% summarize(m = mean(neighbour.abundance))


# change community comp data to wide format 
spread.comm <- comm
spread.comm <- spread.comm %>% 
  group_by_at(vars(-neighbour.abundance)) %>% # group by everything other than the number of neighbours
  mutate(row_id=1:n()) %>% ungroup() %>%  # build group index column
  spread(neighbour.species, neighbour.abundance, fill = 0) %>% 
  select(-row_id)  # drop the index


mean(rowSums(spread.comm[,5:87]))
sd(rowSums(spread.comm[,5:87]))

# combine with fecundity data (need to keep the E subplots that aren't in the comm data)
spread.data<- right_join(spread.comm, f) 

# idensity the mismatch in row number between spread.data and f (they're duplicates so need to remove them with 'distinct')
spread.data <- spread.data %>% distinct(focal.species, plot, treatment, subplot, mature, immature, .keep_all = TRUE)

# make a total fecundity column because we are treating immature same as mature - presumably would also be mature if plant harvested later
spread.data$total.fecundity <- spread.data$mature + spread.data$immature

# HERE MAKE A COLUMN OF GROUP (NATIVE/EXOTIC/FORB/GRASS/UNKNOWN)
#unique(comm$neighbour.species)
spread.data <- as.data.frame(spread.data)
spread.data$unknown <- rowSums(spread.data[,53:85])

spread.data$native.forb <- rowSums(subset(spread.data, select=c(`Brachyscome iberidifolia`, `Brachyscome Iberidifolia`, `Brachyscome perpusilla`, `Bulbine semibarbata`,
                                                                `Calandrinia sp`, `Calotis hispidula`, `Ceratogyne obionoides`, `Cheilanthes austrotenuifolia`,
                                                                `Chthonocephalus pseudevax`, `Daucus glochidiatus`, `Dichopogon sp.`, `Drosera macrantha`,
                                                                `Drosera sp`, `Erodium cygnorum`, `Gilberta tenuifolia`, `Gonocarpus nodulosus`, 
                                                                `Goodenia berardiana`, `Goodenia pusilliflora`, `Haloragis odontocarpa`, `Haloragus sp2`, `Hyalosperma glutinosum`,
                                                                `Hydrocotyle pilifera`, `Lawrencella rosea`, `Nicotiana rotundifolia`, `Omphalolappula concava`,
                                                                `Ophioglossum lusitanicum`, `Plantago debilis`, `Podolepis canescens`, `Podolepis lessonii`, `Pogonolepis muelleriana`,
                                                                `Ptilotus gaudichaudii`, `Rhodanthe laevis`, `Schoenia cassiniana`, `Schoenus nanus`, `Stenopetalum filifolium`, 
                                                                `Trachymene cyanopetala`,`Trachymene ornata`, `Velleia rosea`, `Waitzia acuminata`)))

spread.data$exotic.forb <- rowSums(subset(spread.data, select=c(`Arctotheca calendula`, `Boragia sp`, `Cuscuta campestris`, `Hypochaeris glabra`, `Medicago minima`, 
                                                                `Medicago polymorpha`, `Monoculus monstrosus`, `Sonchus oleraceus`, `Spergula arvensis`)))
spread.data$exotic.grass <- spread.data$`Pentameris airoides`

# change the NA's in competitor columns to 0's so they won't be included as a neighbour in downstream analysis 
#spread.data[rowSums(is.na(spread.data)) > 0,]
spread.data[ , 5:94][is.na(spread.data[ , 5:94])] = 0 
                           
dat <- merge(spread.data,env.sub)

#################################################################################
arca <- subset(spread.data, spread.data$focal.species=="Arctotheca calendula")
arca <- subset(arca, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Arctotheca calendula", "native.forb", "exotic.forb", "exotic.grass", "unknown",  "P.mg.kg", "canopy.cover.percentage"))
arca$intra <- arca$`Arctotheca calendula`
arca$exotic.forb <- arca$exotic.forb-arca$`Arctotheca calendula`
write.csv(arca, "Data/groups.ARCA.csv")

# medi
medi <- subset(spread.data, spread.data$focal.species=="Medicago minima")
medi <- subset(medi, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Medicago minima", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
medi$intra <- medi$`Medicago minima`
medi$exotic.forb <- medi$exotic.forb-medi$`Medicago minima`
write.csv(medi, "Data/groups.MEDI.csv")

# momo
momo <- subset(spread.data, spread.data$focal.species=="Monoculus monstrosus")
momo <- subset(momo, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Monoculus monstrosus", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
momo$intra <- momo$`Monoculus monstrosus`
momo$exotic.forb <- momo$exotic.forb-momo$`Monoculus monstrosus`
write.csv(momo, "Data/groups.MOMO.csv")

# peai
peai <- subset(spread.data, spread.data$focal.species=="Pentameris airoides")
peai <- subset(peai, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Pentameris airoides", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
peai$intra <- peai$`Pentameris airoides`
write.csv(peai, "Data/groups.PEAI.csv")

# dagl
dagl <- subset(spread.data, spread.data$focal.species=="Daucus glochidiatus")
dagl <- subset(dagl, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Daucus glochidiatus", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
dagl$intra <- dagl$`Daucus glochidiatus`
dagl$native.forb <- dagl$native.forb-dagl$`Daucus glochidiatus`
write.csv(dagl, "Data/groups.DAGL.csv")

# gite
gite <- subset(spread.data, spread.data$focal.species=="Gilberta tenuifolia")
gite <- subset(gite, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Gilberta tenuifolia", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
gite$intra <- gite$`Gilberta tenuifolia`
gite$native.forb <- gite$native.forb-gite$`Gilberta tenuifolia`
write.csv(gite, "Data/groups.GITE.csv")

# hygl
hygl <- subset(spread.data, spread.data$focal.species=="Hyalosperma glutinosum")
hygl <- subset(hygl, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Hyalosperma glutinosum", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
hygl$intra <- hygl$`Hyalosperma glutinosum`
hygl$native.forb <- hygl$native.forb-hygl$`Hyalosperma glutinosum`
write.csv(hygl, "Data/groups.HYGL.csv")

# plde
plde <- subset(spread.data, spread.data$focal.species=="Plantago debilis")
plde <- subset(plde, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Plantago debilis", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
plde$intra <- plde$`Plantago debilis`
plde$native.forb <- plde$native.forb-plde$`Plantago debilis`
write.csv(plde, "Data/groups.PLDE.csv")

# poca
poca <- subset(spread.data, spread.data$focal.species=="Podolepis canescens")
poca <- subset(poca, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Podolepis canescens", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
poca$intra <- poca$`Podolepis canescens`
poca$native.forb <- poca$native.forb-poca$`Podolepis canescens`
write.csv(poca, "Data/groups.POCA.csv")

# vero
vero <- subset(spread.data, spread.data$focal.species=="Velleia rosea")
vero <- subset(vero, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Velleia rosea", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
vero$intra <- vero$`Velleia rosea`
vero$native.forb <- vero$native.forb-vero$`Velleia rosea`
write.csv(vero, "Data/groups.VERO.csv")

# trcy
trcy <- subset(spread.data, spread.data$focal.species=="Trachymene cyanopetala")
trcy <- subset(trcy, select = c("focal.species", "plot", "treatment", "subplot", "total.fecundity", "Trachymene cyanopetala", "native.forb", "exotic.forb", "exotic.grass", "unknown"))
trcy$intra <- trcy$`Trachymene cyanopetala`
trcy$native.forb <- trcy$native.forb-trcy$`Trachymene cyanopetala`
write.csv(trcy, "Data/groups.TRCY.csv")

#-------------------------------------------------------------
# supp figure for linear response to neighbour abundance ####
pdf("Figures/supp.linear.abundance.pdf")
par(mfrow=c(4,3), mar=c(2,2,2,2), oma=c(4,4,1,1))
plot(arca$intra, arca$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(arca$native.forb, arca$total.fecundity, col="skyblue", pch=2)
points(arca$exotic.forb, arca$total.fecundity, col="red", pch=3)
points(arca$exotic.grass, arca$total.fecundity, col="orange", pch=4)
legend("topright", legend=c("intra", "native forb", "exotic forb", "exotic grass"), 
       col = c("black", "skyblue", "red", "orange"), pch = c(1,2,3,4), cex = 0.7, bty = "n")
mtext(side=3, expression(italic("A. calendula")), adj = 0, cex=0.8)
mtext(side=2, "Fecundity", line=3)


plot(medi$intra, medi$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(medi$native.forb, medi$total.fecundity, col="skyblue", pch=2)
points(medi$exotic.forb, medi$total.fecundity, col="red", pch=3)
points(medi$exotic.grass, medi$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("M. minima")), adj = 0, cex=0.8)

plot(momo$intra, momo$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(momo$native.forb, momo$total.fecundity, col="skyblue", pch=2)
points(momo$exotic.forb, momo$total.fecundity, col="red", pch=3)
points(momo$exotic.grass, momo$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("M. monstrosus")), adj = 0, cex=0.8)


plot(peai$intra, peai$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(peai$native.forb, peai$total.fecundity, col="skyblue", pch=2)
points(peai$exotic.forb, peai$total.fecundity, col="red", pch=3)
mtext(side=3, expression(italic("P. airoides")), adj = 0, cex=0.8)
mtext(side=2, "Fecundity", line=3)

plot(dagl$intra, dagl$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(dagl$native.forb, dagl$total.fecundity, col="skyblue", pch=2)
points(dagl$exotic.forb, dagl$total.fecundity, col="red", pch=3)
points(dagl$exotic.grass, dagl$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("D. glochidiatus")), adj = 0, cex=0.8)

plot(gite$intra, gite$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(gite$native.forb, gite$total.fecundity, col="skyblue", pch=2)
points(gite$exotic.forb, gite$total.fecundity, col="red", pch=3)
points(gite$exotic.grass, gite$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("G. tenuifolia")), adj = 0, cex=0.8)

plot(hygl$intra, hygl$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(hygl$native.forb, hygl$total.fecundity, col="skyblue", pch=2)
points(hygl$exotic.forb, hygl$total.fecundity, col="red", pch=3)
points(hygl$exotic.grass, hygl$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("H. glutinosum")), adj = 0, cex=0.8)
mtext(side=2, "Fecundity", line=3)

plot(plde$intra, plde$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(plde$native.forb, plde$total.fecundity, col="skyblue", pch=2)
points(plde$exotic.forb, plde$total.fecundity, col="red", pch=3)
points(plde$exotic.grass, plde$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("P. debilis")), adj = 0, cex=0.8)

plot(poca$intra, poca$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(poca$native.forb, poca$total.fecundity, col="skyblue", pch=2)
points(poca$exotic.forb, poca$total.fecundity, col="red", pch=3)
points(poca$exotic.grass, poca$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("P. canescens")), adj = 0, cex=0.8)

plot(trcy$intra, trcy$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(trcy$native.forb, trcy$total.fecundity, col="skyblue", pch=2)
points(trcy$exotic.forb, trcy$total.fecundity, col="red", pch=3)
points(trcy$exotic.grass, trcy$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("T. cyanopetala")), adj = 0, cex=0.8)
mtext(side=2, "Fecundity", line=3)
mtext(side=1, "Neighbour abundance", line=3)

plot(vero$intra, vero$total.fecundity, ylab = "fecundity", xlab = "neighbour abundance")
points(vero$native.forb, vero$total.fecundity, col="skyblue", pch=2)
points(vero$exotic.forb, vero$total.fecundity, col="red", pch=3)
points(vero$exotic.grass, vero$total.fecundity, col="orange", pch=4)
mtext(side=3, expression(italic("V. rosea")), adj = 0, cex=0.8)
mtext(side=1, "Neighbour abundance", line=3)

dev.off()
 




####





