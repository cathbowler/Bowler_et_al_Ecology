source("Figures/Make_data_for_figures.R")

# Set colour pallete 
cols <- c("#E69F00", "tan1", "tan4", "gold", "yellow2", "lightgoldenrod", "chocolate",
          "#56B4E9", "steelblue", "blue", "#009E73")

# Figure 2 ####
sp <- 11
x <- seq(1:((sp-1)*7+4))
blanks <- rep(NA, 3)
vline <- c(7*seq(2:(sp))-1)
cbbPalette <- c("#000000", "brown", "darkblue", "#009E73", "#F0E442", "red", "red")

pdf("Figures/AllAlphas.pdf", width = 11, height=6)
par(mfrow=c(1,2))
par(fig=c(0,.65,0.1,1))
par(mar=c(6,5,5,1))
plot(x, y=c(dagl_mean, blanks, gite_mean, blanks, hygl_mean, blanks, plde_mean, blanks, 
            poca_mean, blanks, trcy_mean, blanks, vero_mean, blanks, arca_mean, blanks, 
            medi_mean, blanks, momo_mean, blanks, peai_mean), 
     xlab="", ylab="", 
     pch=16, col=cbbPalette, xaxt="n", ylim=c(-0.45, 0.25), cex=1)
mtext("Interaction coefficient", side=2, line=2)
abline(h=0, col="black", lwd=1)
abline(v=vline, col="grey")
arrows(x0=x, y0=c(dagl_low, blanks, gite_low, blanks, hygl_low, blanks, plde_low, blanks, 
                  poca_low, blanks, trcy_low, blanks, vero_low, blanks, arca_low, blanks, 
                  medi_low, blanks, momo_low, blanks, peai_low), 
       x1=x, y1=c(dagl_high, blanks, gite_high, blanks, hygl_high, blanks, plde_high, blanks, 
                  poca_high, blanks, trcy_high, blanks, vero_high, blanks, arca_high, blanks, 
                  medi_high, blanks, momo_high, blanks, peai_high), length=.05,
       angle=90, col=cbbPalette, code=3)
axis(side=1, at=(c(vline, 77)-4), tick=TRUE, labels = FALSE)
mtext("(a)", side=3, adj=0)
text(x=(c(vline, 77)-5), y=-.51, srt=90, cex=0.8, adj = c(1,1), labels=c(expression(italic("D. glochidiatus")), expression(italic("G. tenuifolia")), 
                                                            expression(italic("H. glutinosum")), expression(italic("P. debilis")), 
                                                            expression(italic("P. canescens")), expression(italic("T. cyanopetala")), 
                                                            expression(italic("V. rosea")), expression(italic("A. calendula")),
                                                            expression(italic("M. minima")), expression(italic("M. monstrosus")), 
                                                            expression(italic( "P. airoides"))),
     xpd = TRUE, col=c(rep("brown", 7), rep("darkblue",3), "#009E73"))
legend(x=53, y=-.23, c("intraspecific", "native forb", "exotic forb", "exotic grass"), 
       col=c("#000000", "brown", "darkblue", "#009E73", "#F0E442"), lwd=2, bg="white", box.lty=0)
rect(-2,0,100,0.3, col = rgb(0.5,0.5,0.5,1/6), border=NA)

###
par(fig=c(.65,1,0,1), new=TRUE)
par(mar=c(5,3,5,1))
barplot(100*c(intra, NA, native, NA, exotic), ylim=c(0,100), 
        col=c("#000000", "plum", "steelblue1", "black", "#000000", "plum", "steelblue1","black", "#000000", "plum", "steelblue1"), yaxt="n")
abline(h=0, col="black", lwd=2)
abline(v=c(4.3, 9.1), lty=1, col="grey")
axis(side=2, at=c(0,25,50,75,100), tick=T, labels = T)
mtext("Percentage of interactions", side=2, line=2)
legend(6, 100, c("negative effect", "positive effect", "overlap zero"), fill=c("#000000", "plum", "steelblue1"), bg="white", box.lty=0)
text(.5, 100, "(b)", xpd=TRUE)
axis(side=1, at=c(1.9, 6.7, 11.5), tick=TRUE, tck=-0.01, cex.axis=0.8, lab=c("Conspecific \n neighbours", "Native \n neighbours", "Exotic \n neighbours"))

dev.off()

# Figure 3 ####
pdf("Figures/RealisedFecundity.pdf", width=8, height=5)
par(mfrow=c(1,2))
par(fig=c(0,.6,0.1,1))
par(mar=c(5,5,5,1))
plot(Fmean, col=cols, cex.lab=1,
     ylim=c(-2.5,.8), xaxt="n", xlab="", ylab="Net neighbourhood effect", pch=c(19,19,19,19,19,19,19,17,17,17,15))
axis(side=1, at=1:11, las=2, cex.axis=0.8, lab=c(expression(italic("D. glochidiatus")), expression(italic("G. tenuifolia")), 
                                                 expression(italic("H. glutinosum")), expression(italic("P. debilis")), 
                                                 expression(italic("P. canescens")), expression(italic("T. cyanopetala")), 
                                                 expression(italic("V. rosea")), expression(italic("A. calendula")),
                                                 expression(italic("M. minima")), expression(italic("M. monstrosus")), 
                                                 expression(italic( "P. airoides"))))
abline(h=0, lty=2)
arrows(1:11, lowerF, 1:11, upperF, length=0.05, angle=90, code=3, col=adjustcolor(cols))
#legend(1, 1.9, xpd=TRUE, ncol=3, legend=c("Native forb", "Invasive forb", "Invasive grass"), pch=c(19,17,15), col="black", bty="n", cex=0.8)
#text(x=c(6,6), y=c(0.1, 1.5), labels=c("inhibited", "facilitated"), cex=1)
mtext(side=3, adj=0, "(a)")
rect(0.5,0,12,4, col = rgb(0.5,0.5,0.5,1/6), border=NA)


natives.per.fac <- c(dagl.perc, gite.perc, hygl.perc, plde.perc, poca.perc, trcy.perc, vero.perc)
exotics.per.fac <- c(arca.perc, medi.perc, momo.perc, peai.perc)
perc.fac <- c(dagl.perc, gite.perc, hygl.perc, plde.perc, poca.perc, trcy.perc, vero.perc, arca.perc, medi.perc, momo.perc, peai.perc)
perc.fac <- as.data.frame(perc.fac)
perc.fac$X <- c(rep(1,7), rep(2,4))

# calculate quartiles 
nat.upper <- summary(natives.per.fac)[5]-summary(natives.per.fac)[4]
nat.lower <- summary(natives.per.fac)[4]-summary(natives.per.fac)[2]
ex.upper <- summary(exotics.per.fac)[5]-summary(exotics.per.fac)[4]
ex.lower <- summary(exotics.per.fac)[4]-summary(exotics.per.fac)[2]

par(fig=c(.6,1,0.1,1), new=TRUE)
par(mar=c(5,5,5,1))
plot(jitter(perc.fac$X,0.8), perc.fac$perc.fac,
    ylim=c(0,.4), col = cols,  pch=c(19,19,19,19,19,19,19,17,17,17,15), xlab="",ylab="Probability of positive neighbourhood effect", xaxt="n")
points(x=1, y=mean(natives.per.fac), pch=18)
points(x=2, y=mean(exotics.per.fac), pch=18)
arrows(1, mean(natives.per.fac)-nat.lower, 1, mean(natives.per.fac)+nat.upper,  length=0.05, angle=90, code=3)
arrows(2, mean(exotics.per.fac)-ex.lower, 2, mean(exotics.per.fac)+ex.upper,  length=0.05, angle=90, code=3)
axis(side=1, at=c(1,2), lab=c("Native focals", "Exotic focals"))
mtext(side=3, adj=0, "(b)")
dev.off()


# Figure 4 ####

# data for table 
sp <- list("D. glochidiatus", "G. tenuifolia", "H. glutinosum",
           "P. debilis", "P. canescens", "T. cyanopetala", 
           "V. rosea", "A. calendula", "M. minima", "M. monstrosus", 
           "P. airoides")

actualFmean<-Fmean*meanLambdas
difference <- meanLambdas-actualFmean
data <- cbind(sp, round(meanLambdas,2), round(Fmean,2), round(actualFmean,2), round(difference,2))
colnames(data) <- c("Species", "Intrinsic fecundity", "Realised fecundity", "Calculated fecundity", "Difference")
write.csv(data, "Mean_lambdas_realisedfecundity.csv")

#PLOT 

pdf("Figures/FullRatio.pdf", width = 7, height=4)
par(mar=c(4,5,2,2))
plot(meanLambdas, Fmean, col=cols,
     ylim=c(0,1.8), xlim=c(0, 620), ylab="Net neighborhood \n effect", xlab="Intrinsic fecundity", pch=c(19,19,19,19,19,19,19,17,17,17,15))
abline(h=1, lty=2)
arrows(meanLambdas, lowerF, meanLambdas, upperF, length=0.05, angle=90, code=3, col=adjustcolor(cols))
arrows(lowerlam, Fmean, upperlam, Fmean, length=0.05, angle=90, code=3, col=adjustcolor(cols))
legend(20, 1.8, pch=c(19,19,19,19,19,19,19,17,17,17,15), ncol=5,
       legend=c(expression(italic("D. glochidiatus")), expression(italic("G. tenuifolia")), 
                expression(italic("H. glutinosum")), expression(italic("P. debilis")), 
                expression(italic("P. canescens")), expression(italic("T. cyanopetala")), 
                expression(italic("V. rosea")), expression(italic("A. calendula")),
                expression(italic("M. minima")), expression(italic("M. monstrosus")), 
                expression(italic( "P. airoides"))), col=cols, cex=0.8, xpd=T, bty="n")
rect(-100,1,700,2.5, col = rgb(0.5,0.5,0.5,1/6), border=NA)

dev.off()

# NEW FIGURE 4 - NATIVE VERSUS EXOTIC ####

ARCAF_natmeans <- rowMeans(ARCAF_nat, na.rm=T) # or could do correlation per subplot 

cor(ARCAF_natmeans, arca$lambda)
# need a distribution of correlation coeffs between F and lambda for arca
# then a distribtion of corr coeffs across species?


pdf("Figures/Ratio.pdf", width = 7, height=7)
par(mfrow=c(2,1))

par(mar=c(4,5,2,2))
plot(meanLambdas, meanFnats, col=cols, xlim=c(0,620), ylim=c(0,2.8),
     ylab="Net neighborhood \n effect", xlab="Intrinsic fecundity", pch=c(19,19,19,19,19,19,19,17,17,17,15))
arrows(meanLambdas, lower_nat, meanLambdas, upper_nat, length=0.05, angle=90, code=3, col=adjustcolor(cols))
arrows(lowerlam, meanFnats, upperlam, meanFnats, length=0.05, angle=90, code=3, col=adjustcolor(cols))
abline(h=1, lty=2)
mtext(side=3, text=c("(a) Native neighbors"), adj=0)
legend(0, 2.5, pch=c(19,19,19,19,19,19,19,17,17,17,15), ncol=5,
       legend=c(expression(italic("D. glochidiatus")), expression(italic("G. tenuifolia")), 
                expression(italic("H. glutinosum")), expression(italic("P. debilis")), 
                expression(italic("P. canescens")), expression(italic("T. cyanopetala")), 
                expression(italic("V. rosea")), expression(italic("A. calendula")),
                expression(italic("M. minima")), expression(italic("M. monstrosus")), 
                expression(italic( "P. airoides"))), col=cols, cex=0.8, xpd=T, bty="n")
abline(m1, col="grey")
rect(-100,1,1000,4, col = rgb(0.5,0.5,0.5,1/6), border=NA)
#----
par(mar=c(4,5,2,2))
plot(meanLambdas, meanFexs, col=cols, xlim=c(0,620), ylim=c(0,2.8), 
     ylab="Net neighborhood \n effect", xlab="Intrinsic fecundity", pch=c(16,16,16,16,16,16,16,17,17,17,15))
arrows(meanLambdas, lower_ex, meanLambdas, upper_ex, length=0.05, angle=90, code=3, col=adjustcolor(cols))
arrows(lowerlam, meanFexs, upperlam, meanFexs, length=0.05, angle=90, code=3, col=adjustcolor(cols))
abline(h=1, lty=2)
abline(m2, col="grey")
mtext(side=3, text=c("(b) Exotic neighbors"), adj=0)
rect(-100,1,1000,4, col = rgb(0.5,0.5,0.5,1/6), border = NA)


dev.off()


