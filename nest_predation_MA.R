# R script by Malgorzata Lagisz 

# for "Does urbanization affect predation of bird nests? A meta-analysis"
# by Ernő Vincze, Gábor Seress, Malgorzata Lagisz, Shinichi Nakagawa, Niels Dingemanse, Philipp Sprau
# Accepted for journal: "Frontiers in Ecology and Evolution", 2017

# LICERNSE
# The MIT License (MIT)
# 
# Copyright (c) 2017 Malgorzata Lagisz
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Final dataset in the file: nest_predation_MA_data.csv
# #Variables:
# Paper_ID -	Identification number of the paper - P + 3 digits (i.e. P001, P035)
# Group_ID -	Identification number of the research group (NOTE: this will overlap with Species!)
# Journal_ID -	Six-character identification code for the journals (for analysing publication bias - see full list on a separate sheet)
# Pub_year -	Year of publication (for analysing publication bias)
# rho -	Positive value means higher survival in more urbanized habitat; negative value means higher mortality in more urbanized habitat
# total_nests -	The total number of nests used in calculating rho (sample size)
# variance -	Variance calculated from sample size (1 / rho-3)
# source -	The values rho is calculated from. "numbers" = the exact numbers of predated and non-predated nests is given in the text; "percentage" = the numbers are calculated from percentages of predation/survival given in the text or table; "figure" = the numbers are calculated from measuring plotted data points on a figure; "DPR" = the numbers are converted from daily predation rates; other sources?
# death_cause -	0 = predation and other sources of mortality (weather, abandonment) not differentiated; 1 = only nests that died from predation are included
# partial_predation -	How nests that were partially predated are treated. "pred" = nests where at least one egg/offspring died are considered predated, only completely intact nests count as survived; "alive" = nests where at least one egg/offspring survided are considered survived, only nests where no offspring survived count as predated; none = there was only 1 egg/offspring per nest
# Nest_type -	Natural nest population (nat) or artificial nests (art)
# Egg_number -	Number of eggs in a nest (including real and fake eggs)
# Fake eggs - how many fake eggs per nest
# Nest_position -	1 = on the ground; 2 = elevated (shrubs or trees); 1,5 = both positions alternating, without separate effect sizes given for the two
# Nest_openness -	1 = open cup-nest; 2 = closed nest in a box; 1,5 = orb-like nest on tree
# Study_days -	The number of days the eggs were exposed to predation
# Study_year -	The year(s) belonging to the predation rate
# Species -	The species the nest / eggs belong to

sessionInfo()
#R version 3.2.4 (2016-03-10)
#Platform: x86_64-apple-darwin13.4.0 (64-bit)
#Running under: OS X 10.9.5 (Mavericks)
options(scipen=100)
rm(list=ls())

library(metafor)
library(ape)

# ##################  Custom Functions ##################  
# ## Function to convert r (rho) to Zr
# r.to.Zr <- function(r){
#     Zr <- 0.5*(log(1+r)-log(1-r))
#     print(Zr)}
# ## Function to convert Zr to r        
# Zr.to.r <- function(Zr){
#     r <- (exp(2*Zr)-1)/(exp(2*Zr)+1)
#     print(r)}
# ######################################################  
# ### Calculating effect sizes - SKIP, as this was already done and is included in the dataset:
# ## Calculate Zr using rho (original urbanisation levels accoording to the study authors, adjusted to 0-1 scale)
# Data$Zr <- r.to.Zr(Data$rho)
# ## Calculate Zr (VZr will be the same) using rho2 (uniformized urbanisation level scores, 1-5 scale)
# Data$Zr2 <- r.to.Zr(Data$rho2)
# ## Calculate variance for Zr, basing on sample sizes (total number of nests)
# Data$VZr <- 1/(Data$total_nests-3)
# #write.csv(Data, file="nest_predation_MA_data.csv", row.names = FALSE)


###########################################################
# DATA SUMMARY
###########################################################

################### full data set 

Data <- read.csv(file="nest_predation_MA_data.csv", stringsAsFactors = TRUE)
dim(Data) 
str(Data) # 117 ES, 29 columns, 51 papers, 45 groups, 32 species
head(Data)

## Check for missing data ('NA'):
table(is.na(Data)) # NA in total
unlist(lapply(Data, function(x) any(is.na(x)))) #columns that contain NA
unlist(lapply(Data, function(x) sum(is.na(x)))) #how many NA per column: Study_days and Species approx 50% NA, Nest_height 31 missing

## Quick data checks
hist(Data$Pub_year)
range(Data$Pub_year)
table(Data$Study_year) #not usable for analyses due to ranges present, use median
hist(Data$Median_year) #ok
hist(Data$rho)
hist(Data$rho2)
plot(Data$rho, Data$rho2) #just a few different rho vs. rho2 values

par(mfrow=c(1,1))
plot(Data$Median_year, Data$Pub_year) #some delay with publishing relatively to median data collection year, as expected
cor.test(Data$Median_year, Data$Pub_year) # cor=0.79, t=13.822, p<0.0001 - strongly correlated


###########################################################
# ANALYSES on rho2 (uniformized urbanisation scores)
###########################################################

###########################################################
#### Simple meta-analysis - all data

## MA without random effects
MA_0 <- rma(yi=Zr2, vi=VZr, method="REML", data=Data)
summary(MA_0) #I^2 (total heterogeneity / total variability):   91.1%

#funnel plot
plot(jitter(Data$Zr2), jitter(sqrt(1/Data$VZr)), cex=0.75, xlim=c(-1.1,1.1), xlab="Zr", ylab="Precision (1/SE.Zr)", main="all data")
abline(v=0, lwd=0.5)

## MA with random effects only
MA_all <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), method="REML", data=Data)
summary(MA_all) #estimate ns Zr2 -0.0034 CI -0.0806   0.0738 
funnel(MA_all) 
#forest(MA_all)
# calculate I2 values
s2m <- sum(1/Data$VZr) * (MA_all$k-1) / (sum(1/Data$VZr)^2 - sum((1/Data$VZr)^2)) # typical sampling error
s2t <- sum(MA_all$sigma2)+s2m # total sampling error
I2t <- sum(MA_all$sigma2)/s2t
I2t*100 # 92.7% - total heterogeneity
I2s <- MA_all$sigma2[1]/s2t
I2s*100 # 81.7% varaince due to Study
I2e <- MA_all$sigma2[2]/s2t
I2e*100 # 11.0% - residuals against sampling error
transf.ztor(MA_all$b[1,]) # intrecept as r = -0.003
transf.ztor(MA_all$ci.lb[1]) # CI.lb for the intercept as r = -0.080
transf.ztor(MA_all$ci.ub[1]) # CI.ub for the intercept as r = 0.074
res <- data.frame(Model="Meta-analytic mean, all data", M=MA_all$b, CI.lb=MA_all$ci.lb, CI.ub=MA_all$ci.ub, pch=18) #harvest results (Zr)


#### Meta-regression on all data with nest_type as moderator
MR_nest_intc <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~Nest_type-1, method="REML", data=Data)
summary(MR_nest_intc) # in artificial nests signif negative Zr -0.1163 (more predation with incr urbanization); 
# in natural nests ns positive Zr 0.0811 (less predation with incr urbanization)
funnel(MR_nest_intc) # 1 outlier ES < -1 (same point as in the full data set)
forest(MR_nest_intc)
summ <- summary(MR_nest_intc)
forest(x=summ$b, ci.lb=summ$ci.lb, ci.ub=summ$ci.ub, slab=dimnames(summ$b)[[1]],  xlab="Zr2")
# calculate I2 values
s2m <- sum(1/Data$VZr) * (MR_nest_intc$k-1) / (sum(1/Data$VZr)^2 - sum((1/Data$VZr)^2)) # typical sampling error
s2t <- sum(MR_nest_intc$sigma2)+s2m # total sampling error
I2t <- sum(MR_nest_intc$sigma2)/s2t
I2t*100 # 91.9% - total heterogeneity
I2s <- MR_nest_intc$sigma2[1]/s2t
I2s*100 # 79.9% varaince due to Study
I2e <- MR_nest_intc$sigma2[2]/s2t
I2e*100 # 12.0% - residuals against sampling error
transf.ztor(MR_nest_intc$b[1,]) #-0.116
transf.ztor(MR_nest_intc$ci.lb[1]) #-0.224
transf.ztor(MR_nest_intc$ci.ub[1]) #-0.005
transf.ztor(MR_nest_intc$b[2,]) #0.081
transf.ztor(MR_nest_intc$ci.lb[2]) #-0.015
transf.ztor(MR_nest_intc$ci.ub[2]) #0.176
res <- rbind(res, data.frame(Model=c("Artificial vs. Natural nests:","   Artificial nests *","   Natural nests"), M=c(NA,MR_nest_intc$b), CI.lb=c(NA,MR_nest_intc$ci.lb), CI.ub=c(NA,MR_nest_intc$ci.ub), pch=c(20,20,20))) #harvest results (Zr)

MR_nest_diff <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~Nest_type, method="REML", data=Data)
summary(MR_nest_diff)  # diff between artificial and natural nests signif Zr 0.197 CI 0.050  0.345
transf.ztor(MR_nest_diff$b[2,]) #0.195
transf.ztor(MR_nest_diff$ci.lb[2]) #0.050
transf.ztor(MR_nest_diff$ci.ub[2]) #0.332
res <- rbind(res, data.frame(Model="   Difference: Natural - Artificial nests *", M=MR_nest_diff$b[2], CI.lb=MR_nest_diff$ci.lb[2], CI.ub=MR_nest_diff$ci.ub[2], pch=20)) #harvest results (Zr)

#tidy up extracted results table
res$M <- round(transf.ztor(res$M),3)
res$CI.lb <- round(transf.ztor(res$CI.lb),3)
res$CI.ub <- round(transf.ztor(res$CI.ub),3)
write.csv(res,"MA_MR_alldata_main_res_rho2.csv")
revres <- res[rev(rownames(res)),] #reverse for plotting (from bottom to top)


### PLOT - MA and MR models on all data
#opar <- par()      # make a copy of current settings
#par(opar)          # restore original settings

pdf(file="Fig_MA_MR_alldata_rho2.pdf",width=4,height=3,pointsize=8)
par(mfrow=c(1,1))
par(mar=c(5,14,1,1))
plot(revres$M, 1:length(revres$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.5,0.5), ylim=c(0.25, length(revres$M)+.5), xlab="effect size [r]", pch=revres$pch, cex=1, cex.axis=.9)
abline(v=0,lty=3)
mtext(revres$Model, 2, 13, at=1:length(revres$M), las=2, cex=1, font=1, adj=0)
segments(revres$CI.lb, 1:length(revres$M), revres$CI.ub, 1:length(revres$M), lwd=1.25)
dev.off()


#################### ARTIFICIAL NESTS ################

DataA <- Data[Data$Nest_type == "art", ]
dim(DataA) # 59 31
DataA <- droplevels(DataA)

unlist(lapply(DataA, function(x) sum(is.na(x)))) #how many NA per column: Study_days 18, Nest_height 1

#### MA
#funnel plot
plot(jitter(DataA$Zr2), jitter(sqrt(1/DataA$VZr)), cex=0.75, xlim=c(-1.1,1.1), xlab="Zr2", ylab="Precision (1/SE.Zr)", main="all data")
abline(v=0, lwd=0.5)

## MA with random effects only
MA_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), method="REML", data=DataA)
summary(MA_A) #estimate ns Zr -0.1183 CI  -0.2423   0.0056 
funnel(MA_A)
# calculate I2 values
s2m <- sum(1/DataA$VZr) * (MA_A$k-1) / (sum(1/DataA$VZr)^2 - sum((1/DataA$VZr)^2)) # typical sampling error
s2t <- sum(MA_A$sigma2)+s2m # total sampling error
I2t <- sum(MA_A$sigma2)/s2t
I2t*100 # 93.2% - total heterogeneity
I2s <- MA_A$sigma2[1]/s2t
I2s*100 # 83.8% varaince due to Study
I2e <- MA_A$sigma2[2]/s2t
I2e*100 # 9.4% - residuals against sampling error
transf.ztor(MA_A$b[1,]) #-0.118
transf.ztor(MA_A$ci.lb[1]) #-0.2376
transf.ztor(MA_A$ci.ub[1]) #0.0056
resA <- data.frame(Model="Meta-analytic mean", M=MA_A$b, CI.lb=MA_A$ci.lb, CI.ub=MA_A$ci.ub, pch=18) #harvest results (Zr)


#### MR - Meta-Regressions (univariate)

# univariate meta-regression with Nest_openness 
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_openness)-1, method="REML", data=DataA)
summary(MR_A) # cup ns Zr2  -0.1196  -0.2435  -0.0042 
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.2%
plot(DataA$Zr2 ~ DataA$Nest_openness)
resA <- rbind(resA, data.frame(Model=c("Nest openness:","   Cup","   Hole"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20))) #harvest results (Zr)

# univariate meta-regression with Nest_openness Hole-Cup difference
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_openness), method="REML", data=DataA)
summary(MR_A) # diff ns 0.0834  -0.0510   0.2179 

# univariate meta-regression with Nest_position
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_position)-1, method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.0%
plot(DataA$Zr2 ~ DataA$Nest_position)
resA <- rbind(resA, data.frame(Model=c("Nest position:","   Elevated","   Ground","   Mix"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# univariate meta-regression with Egg_number
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(DataA$Egg_number), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.5%
plot(DataA$Zr2 ~ DataA$Egg_number)
resA <- rbind(resA, data.frame(Model=c("Egg number (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with Study_days (NOTE: 18 NA)
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(Study_days), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 95.6%
plot(DataA$Zr2 ~ DataA$Study_days)
resA <- rbind(resA, data.frame(Model=c("Study duration (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with Median_year
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(Median_year), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.3%
plot(DataA$Zr2 ~ DataA$Median_year)
resA <- rbind(resA, data.frame(Model=c("Median study year (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with Pub_year
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(DataA$Pub_year), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.5%
plot(DataA$Zr2 ~ DataA$Pub_year)
resA <- rbind(resA, data.frame(Model=c("Publication year (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with urbmin as a factor
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(urbmin)-1, method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 94.3
plot(DataA$Zr2 ~ DataA$urbmin)
resA <- rbind(resA, data.frame(Model=c("Min urbanisation score:","   1","   2","   3","   4"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20,20,20))) #harvest results (Zr)
# univariate meta-regression with urbmin as continuous - not used becouse the relationship is driven by just 2 data points at 4

# univariate meta-regression with urbmax as a factor
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(urbmax)-1, method="REML", data=DataA)
summary(MR_A) # urbmax=4 signif Zr2 -0.3232 -0.6290  -0.0174, but driven by only 4 data points! (most are at urbanmax=3)
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.2%
plot(DataA$Zr2 ~ DataA$urbmax)
resA <- rbind(resA, data.frame(Model=c("Max urbanisation score:","   3","   4 *","   5"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# meta-regression with urbmin_scaled*urbmax_scaled interaction (as continuous predictors)
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(urbmin) * scale(urbmax), method="REML", data=DataA)
summary(MR_A) # all ns
# meta-regression with Nest_openness*Nest_position interaction
MR_A <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_openness) * as.factor(Nest_position), method="REML", data=DataA)
summary(MR_A) # ns
table(DataA$Nest_openness, DataA$Nest_position) #not enough data for estimating interaction: 0 hole-like nests at ground and mix levels
boxplot(DataA$Zr2 ~  DataA$Nest_position + DataA$Nest_openness)

# meta-regression with Nest_openness + Nest_position without interaction
MR_A1 <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_position) + as.factor(Nest_openness) -1, method="REML", data=DataA)
summary(MR_A1) # all ns
#resA <- rbind(resA, data.frame(Model=c("Multivariate:","   Cup - elevated","   Cup - ground","   Cup - mix","   Hole - elevated"), M=c(NA,MR_A1$b), CI.lb=c(NA,MR_A1$ci.lb), CI.ub=c(NA,MR_A1$ci.ub), pch=c(20,20,20,20,20))) #harvest results (Zr)
resA <- rbind(resA, data.frame(Model=c("Multivariate:","   Cup - elevated","   Cup - ground","   Cup - mix"), M=c(NA,MR_A1$b[1:3]), CI.lb=c(NA,MR_A1$ci.lb[1:3]), CI.ub=c(NA,MR_A1$ci.ub[1:3]), pch=c(20,20,20,20))) #harvest results (Zr)
MR_A2 <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_position) + relevel(as.factor(Nest_openness), ref="hole") -1, method="REML", data=DataA)
summary(MR_A2) # all ns
resA <- rbind(resA, data.frame(Model="   Hole - elevated", M=c(MR_A2$b[1]), CI.lb=c(MR_A2$ci.lb[1]), CI.ub=c(MR_A2$ci.ub[1]), pch=c(20))) #harvest results (Zr)

# tidy up extracted results table
resA$M <- round(transf.ztor(resA$M),3)
resA$CI.lb <- round(transf.ztor(resA$CI.lb),3)
resA$CI.ub <- round(transf.ztor(resA$CI.ub),3)
write.csv(resA,"MA_MR_dataA_resA_rho2.csv")
revresA <- resA[rev(rownames(resA)),] #reverse for plotting (from bottom to top)


### PLOT - MA and MR models on all data

pdf(file="Fig_MA_MR_dataA_rho2.pdf",width=4,height=6,pointsize=10)
par(mfrow=c(1,1))
par(mar=c(4,10,2,0))
plot(revresA$M, 1:length(revresA$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.6,0.6), ylim=c(0.25, length(revresA$M)+.5), xlab="Zr", pch=revresA$pch, cex=1.1, cex.axis=.9)
abline(v=0,lty=3)
mtext(revresA$Model, 2, 9, at=1:length(revresA$M), las=2, cex=.8, font=1, adj=0)
segments(revresA$CI.lb, 1:length(revresA$M), revresA$CI.ub, 1:length(revresA$M), lwd=1.25)
dev.off()



#################### NATURAL NESTS ################

DataN <- Data[Data$Nest_type == "nat", ]
dim(DataN) # 58 31
DataN <- droplevels(DataN)

unlist(lapply(DataN, function(x) sum(is.na(x)))) #how many NA per column: Study_days 31 NA, Nest_height 30 NA, Egg_number 16 NA

#### MA
#funnel plot
plot(jitter(DataN$Zr2), jitter(sqrt(1/DataN$VZr)), cex=0.75, xlim=c(-1.1,1.1), xlab="Zr2", ylab="Precision (1/SE.Zr)", main="all data")
abline(v=0, lwd=0.5)

## MA with 2 random effects
MA_N <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), method="REML", data=DataN)
summary(MA_N) #estimate Zr ns 0.0794 CI  -0.0072   0.1664 
# overall ns positive relationship between urbanisation and survival - less predation with increasing urbanization
funnel(MA_N)
forest(MA_N)
# calculate I2 values
s2m <- sum(1/DataN$VZr) * (MA_N$k-1) / (sum(1/DataN$VZr)^2 - sum((1/DataN$VZr)^2)) # typical sampling error
s2t <- sum(MA_N$sigma2)+s2m # total sampling error
I2t <- sum(MA_N$sigma2)/s2t
I2t*100 # 89.92% - total heterogeneity
I2s <- MA_N$sigma2[1]/s2t
I2s*100 # 73.19% varaince due to Study
I2e <- MA_N$sigma2[2]/s2t
I2e*100 # 16.73% - residuals against sampling error
transf.ztor(MA_N$b[1,]) #0.0793
transf.ztor(MA_N$ci.lb[1]) #-0.0071
transf.ztor(MA_N$ci.ub[1]) #0.1646
resN <- data.frame(Model="Meta-analytic mean", M=MA_N$b, CI.lb=MA_N$ci.lb, CI.ub=MA_N$ci.ub, pch=18) #harvest results (Zr)

## MA with 3 random effects, including Species identity
MA_N <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|Species,~1|ES_ID), method="REML", data=DataN)
summary(MA_N) #estimate Zr ns 0.0657 CI  -0.0139   0.1452 
# overall ns positive relationship between urbanisation and survival - less predation with increasing urbanization
funnel(MA_N)
forest(MA_N)
# calculate I2 values
s2m <- sum(1/DataN$VZr) * (MA_N$k-1) / (sum(1/DataN$VZr)^2 - sum((1/DataN$VZr)^2)) # typical sampling error
s2t <- sum(MA_N$sigma2)+s2m # total sampling error
I2t <- sum(MA_N$sigma2)/s2t
I2t*100 # 88.53% - total heterogeneity
I2s <- MA_N$sigma2[1]/s2t
I2s*100 # 0% varaince due to Study
I2e <- MA_N$sigma2[2]/s2t
I2e*100 # 74.75% varaince due to Species identity - equivalent to using just Study before, usually 1 Species per study, and 1 study per species
#table(DataN$Species,DataN$Paper_ID)
I2e <- MA_N$sigma2[3]/s2t
I2e*100 # 13.77% - residuals against sampling error

## Phylogeny
birds <- unique(DataN$Species_latin) #get list of unique bird species latin names
birds_stree <- read.tree("Ericson.tre") # load birds supertree from a file with Ericson backbone
birds_steer <- collapse.singles(birds_stree)
bird_tree_species <- as.character(birds_stree$tip.label) # extract list of species from the tree
setdiff(birds, bird_tree_species) #Cyanistes_caeruleus not matching 
birds_stree$tip.label <- sub("Parus_caeruleus","Cyanistes_caeruleus",birds_stree$tip.label) #replace with synonym Parus_caeruleus
bird_tree_species <- as.character(birds_stree$tip.label) # extract list of species from the tree
intersect(bird_tree_species, birds) #32 species matching, ok
tree <- drop.tip(birds_stree, birds_stree$tip.label[-match(birds, birds_stree$tip.label)]) # prune the supertree tree to a list of taxa from our list
is.binary.tree(tree) #TRUE
is.ultrametric(tree) #TRUE
plot(tree, cex=0.8) #plot with branch lengths
tree2 <- rotate(tree, 37) #rotate some branches

pdf("tree.pdf",width = 4, height = 5, pointsize=8)
par(mfcol=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(tree2)
dev.off()

write.tree(tree2, file = "birds_32sp_tree.tre", append = FALSE, digits = 10, tree.names = FALSE)

### phylogenetic meta-analysis (does not converge if Species ID added)
tree <- read.tree("birds_32sp_tree.tre")# upload cleaned-up and prepeprocessed phylogenetic tree file
CorMatrix <- vcv(tree, corr=TRUE) # make a phylogenetic correlation matrix

MA_N_phylo <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MA_N_phylo) #intc = 0.0337 ns, CI -0.1645   0.2319 
sum(MA_N_phylo$sigma2)/(sum(MA_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MA_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 92.8%
s2t <- sum(MA_N_phylo$sigma2) + sum(1/MA_N_phylo$vi) * (MA_N_phylo$k-1) / (sum(1/MA_N_phylo$vi)^2 - sum((1/MA_N_phylo$vi)^2))
MA_N_phylo$sigma2[1]/s2t*100 # I^2 between study (id) = 3.67%
MA_N_phylo$sigma2[2]/s2t*100 # I^2 phylogeny = 78.46%
MA_N_phylo$sigma2[3]/s2t*100 # I^2 within study or residual = 10.63%
resN <- rbind(resN, data.frame(Model=c("Phylogenetic meta-analytic mean"), M=c(MA_N_phylo$b), CI.lb=c(MA_N_phylo$ci.lb), CI.ub=c(MA_N_phylo$ci.ub), pch=c(18))) #harvest results (Zr)

# plots
Names <- paste(DataN$Paper, DataN$Species, sep="_")
forest(MA_N_phylo, slab=Names)
funnel(MA_N_phylo, yaxis="seinv")

#### Meta-regression 

# species meta-regression
MA_N_species <- rma.mv(yi=Zr2, V=VZr, mod=~Species_latin-1, random=list(~1|Paper_ID,~1|ES_ID), data=DataN, method="REML") #without phylogeny
summary(MA_N_species) # some species very signif + or -
sum(MA_N_species$sigma2)/(sum(MA_N_species$sigma2)+(sum(1/DataN$VZr)*(MA_N_species$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 57.2%
s2t <- sum(MA_N_species$sigma2) + sum(1/MA_N_species$vi) * (MA_N_species$k-1) / (sum(1/MA_N_species$vi)^2 - sum((1/MA_N_species$vi)^2))
MA_N_species$sigma2[1]/s2t # I^2 between study (id) = 0.00% (studies are almost equivalent to the species)
MA_N_species$sigma2[2]/s2t # I^2 within study or residual = 57.2%
# plots
Names <- paste(DataN$Paper, DataN$Species, sep="_")
forest(MA_N_phylo, slab=Names)
funnel(MA_N_phylo, yaxis="seinv")
#plot(tree)
# extract point estimates and their CI for each species
ssp <- summary(MA_N_species)
sp_df <- data.frame(sp = substr(attr(ssp$b,"dimnames")[[1]], 14, 36), M = ssp$b, CI.lb = ssp$ci.lb, CI.ub = ssp$ci.ub)
sp_df <- sp_df[match(tree$tip.label, sp_df$sp),] #reorder dataframe to match order of the tip labels on the tree
# tidy up extracted results table
sp_df$M <- round(transf.ztor(sp_df$M),3)
sp_df$CI.lb <- round(transf.ztor(sp_df$CI.lb),3)
sp_df$CI.ub <- round(transf.ztor(sp_df$CI.ub),3)
sp_df$Signif <- ifelse(sp_df$CI.lb>0 | sp_df$CI.ub<0, "*", "") #add column with stars for estimates that significantly differ from zero
write.csv(sp_df,"MA_MR_dataN_sp_df_rho2.csv")
revres_sp <- sp_df[rev(rownames(sp_df)),] #reverse for plotting (from bottom to top)

### PLOT - species phylogeny and forest plot
pdf(file="Fig_species_rho2.pdf",width=6,height=4,pointsize=9)
par(mfrow=c(1,2))
par(mar=c(4.5,2,1.35,0))
plot(tree, font=1, cex=0.8, x.lim=90, show.tip.label = FALSE)
par(mar=c(4,6,1,0))
plot(sp_df$M, 1:length(sp_df$M), ylab=NA, yaxt="n", bty="n", xlim=c(-1.5,1.5), ylim=c(0.25, length(sp_df$M)+.5), xlab="effect size [r]", pch=16, cex=0.8, cex.axis=.9)
abline(v=0,lty=3)
mtext(sp_df$sp, 2, -1, at=1:length(sp_df$M), las=2, cex=.8, font=3)
segments(sp_df$CI.lb, 1:length(sp_df$M), sp_df$CI.ub, 1:length(sp_df$M), lwd=1.25)
for (i in 1:length(rownames(sp_df))) mtext(sp_df$Signif[i], 2, -1.5, at=i, las=2, cex=.7) #add stars
dev.off()


# Phylogenetic meta-regression with death_cause
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 88.4
plot(DataN$Zr2 ~ DataN$death_cause) # positive values (less predation in urbanised areas) when "yes" = only nests that died from predation are included
resN <- rbind(resN, data.frame(Model=c("Predation as only source of mortality:","   No","   Yes"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20))) #harvest results (Zr)

MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # ns diff
#MR_N <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause), random=list(~1|Paper_ID,~1|ES_ID), data=DataN, method="REML") #without phylogeny
#summary(MR_N) # signif diff 0.1712 CI 0.0096  0.3328  * - result without phylogeny

# Phylogenetic meta-regression with Nest_openness
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(Nest_openness)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
#MR_N <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(Nest_openness)-1, random=list(~1|Paper_ID,~1|ES_ID), data=DataN, method="REML") #without phylogeny
#summary(MR_N) # hole signif 0.2285 CI 0.0406  0.4163 * - result without phylogeny
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 91.9
plot(DataN$Zr2 ~ DataN$Nest_openness) 
resN <- rbind(resN, data.frame(Model=c("Nest openness:","   Cup","   Hole","   Orb"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# Phylogenetic meta-regression with Nest_position
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(Nest_position)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 92.8
plot(DataN$Zr2 ~ DataN$Nest_position) # more positive values in ground-located nests
resN <- rbind(resN, data.frame(Model=c("Nest position:","   Elevated","   Ground","   Mix"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# Phylogenetic meta-regression with log(Nest_height+1)
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~scale(log(Nest_height+1)), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # similar ns trend as above -0.0909 -0.1999  0.0182
plot(DataN$Zr2 ~ log(DataN$Nest_height+1))
resN <- rbind(resN, data.frame(Model="Nest height above ground (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

hist(DataN$Nest_height[DataN$death_cause=="yes"], col="blue", xlim=c(0,25))
hist(DataN$Nest_height[DataN$death_cause=="no"], col=rgb(0, 0, 0, 0.5), add=TRUE)
table(DataN$death_cause, DataN$Nest_height) # measure predation = "no" in nest height > 10
t.test(DataN$Nest_height ~ DataN$death_cause, var.equal = TRUE) #t = 3.3704, df = 26, p-value = 0.002354 
#high nests usually report overall mortrality, not by predation
table(is.na(DataN$Nest_height), DataN$Nest_position) # 30 missing height values are from elevated nests, i.e. only a subset of heights known

# Phylogenetic meta-regression with Egg_number (NOTE: 16 NA) 
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~scale(DataN$Egg_number), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 94.4
plot(DataN$Zr2 ~ DataN$Egg_number) 
resN <- rbind(resN, data.frame(Model="Egg number (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

# Phylogenetic meta-regression with Study_days (NOTE: 28 NA) 
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~scale(DataN$Study_days), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 98.3
plot(DataN$Zr2 ~ DataN$Study_days) 
resN <- rbind(resN, data.frame(Model="Study duration (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

# Phylogenetic meta-regression with Median_year as moderator
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~scale(DataN$Median_year), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 92.7
plot(DataN$Zr2 ~ DataN$Median_year) # all ns
resN <- rbind(resN, data.frame(Model="Median study year (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

hist(DataN$Median_year[DataN$death_cause=="yes"], col="blue")
hist(DataN$Median_year[DataN$death_cause=="no"], col=rgb(0, 0, 0, 0.5), add=TRUE)
table(DataN$death_cause, DataN$Median_year) # recent studies more likely to only measure predation = "yes"
t.test(DataN$Median_year ~ DataN$death_cause, var.equal = FALSE) #t = -2.3125, df = 26.758, p-value = 0.02868

# Phylogenetic meta-regression with Pub_year as moderator
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~scale(DataN$Pub_year), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 92.8
plot(DataN$Zr2 ~ DataN$Pub_year) 
resN <- rbind(resN, data.frame(Model="Publication year (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

# Phylogenetic meta-regression with urbmin as a factor
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(urbmin)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 93.8
plot(DataN$Zr2 ~ DataN$urbmin) # only 1 data point at urbmin=4, and 2 at 3
resN <- rbind(resN, data.frame(Model=c("Min urbanisation score:","   1","   2","   3","   4"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20,20))) #harvest results (Zr)

# Phylogenetic meta-regression with urbmax as a factor
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(urbmax)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 93.6
plot(DataN$Zr2 ~ DataN$urbmax) 
resN <- rbind(resN, data.frame(Model=c("Min urbanisation score:","   3","   4","   5"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# meta-regression with urbmin_scaled*urbmax_scaled interaction (as continuous predictors)
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(urbmin) * scale(urbmax), method="REML", data=DataN)
summary(MR_N_phylo) # ns interaction
# meta-regression with Death_cause and Nest_openness, without interaction 
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause) + as.factor(Nest_openness) -1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # (death_cause)yes signif 0.1668  0.0528  0.2808 **
# meta-regression with Death_cause*Nest_openness interaction
MR_N_phylo <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause) * as.factor(Nest_openness), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # signif interaction
MR_N_phylo1 <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause) * as.factor(Nest_openness) - 1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo1) # signif interaction
resN <- rbind(resN, data.frame(Model=c("Multivariate meta-regression:","   Cup - No","   Cup - Yes *"), M=c(NA,MR_N_phylo1$b[1:2]), CI.lb=c(NA,MR_N_phylo1$ci.lb[1:2]), CI.ub=c(NA,MR_N_phylo1$ci.ub[1:2]), pch=c(20,20,20))) #harvest results (Zr)
MR_N_phylo2 <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause) * as.factor(relevel(Nest_openness, ref="hole")) - 1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo2) # signif interaction
resN <- rbind(resN, data.frame(Model=c("   Hole - No","   Hole - Yes *"), M=c(MR_N_phylo2$b[1:2]), CI.lb=c(MR_N_phylo2$ci.lb[1:2]), CI.ub=c(MR_N_phylo2$ci.ub[1:2]), pch=c(20,20))) #harvest results (Zr)
MR_N_phylo3 <- rma.mv(yi=Zr2, V=VZr, mod= ~as.factor(death_cause) * as.factor(relevel(Nest_openness, ref="orb")) - 1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo3) # signif interaction
resN <- rbind(resN, data.frame(Model=c("   Orb - Yes"), M=c(MR_N_phylo3$b[2]), CI.lb=c(MR_N_phylo3$ci.lb[2]), CI.ub=c(MR_N_phylo3$ci.ub[2]), pch=c(20))) #harvest results (Zr)
resN <- rbind(resN, data.frame(Model=c("  Nest openness - Death cause interaction *"), M=c(MR_N_phylo3$b[5]), CI.lb=c(MR_N_phylo3$ci.lb[5]), CI.ub=c(MR_N_phylo3$ci.ub[5]), pch=c(20))) #harvest results (Zr)

table(DataN$death_cause, DataN$Nest_openness) # only 5 data points are hole/yes, 2 are are orb/yes, 0 are orb/no
boxplot(DataN$Zr2 ~ DataN$death_cause * DataN$Nest_openness, varwidth=TRUE) #most positive values in hole/yes subest - less predation with increasing urbanisation

# tidy up extracted results table
resN$M <- round(transf.ztor(resN$M),3)
resN$CI.lb <- round(transf.ztor(resN$CI.lb),3)
resN$CI.ub <- round(transf.ztor(resN$CI.ub),3)
write.csv(resN,"MA_MR_dataN_resN_rho2.csv")
revresN <- resN[rev(rownames(resN)),] #reverse for plotting (from bottom to top)


### PLOT - MA and MR models on all data
opar <- par()      # make a copy of current settings
#par(opar)          # restore original settings

pdf(file="Fig_MA_MR_dataN_rho2.pdf",width=4,height=6,pointsize=10)
par(mfrow=c(1,1))
par(mar=c(4,10,2,0))
plot(revresN$M, 1:length(revresN$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.8,0.8), ylim=c(0.25, length(revresN$M)+.5), xlab="effect size [r]", pch=revresN$pch, cex=1.1, cex.axis=.9)
abline(v=0,lty=3)
mtext(revresN$Model, 2, 9, at=1:length(revresN$M), las=2, cex=.8, font=1, adj=0)
segments(revresN$CI.lb, 1:length(revresN$M), revresN$CI.ub, 1:length(revresN$M), lwd=1.25)
dev.off()

## OVERALL for natural nests: tendency for less predation with increasing urbanisation (positive Zr) - small overall effect size, 
# this effect is more pronounced when only nest lost due to predation included
# and no much effect of urbanisation when other causes of mortality (confounding) are potentially present.
# Tendency for more negative valus in nests higher above ground  - more mortality in urbanised areas, 
# in lower nests more positive values more likely - less predation in urbanised areas, 
# but this result is likely to be related on death_cause variable.



### PLOT - bubble plots for natural nests: r-death_cause, r-nest_height
plot(DataN$Zr2 ~ DataN$death_cause) # positive values (less predation in urbanised areas) when "yes" = only nests that died from predation are included
plot(DataN$Zr2 ~ DataN$Nest_height) # K=23

pdf(file="Fig_bubble2_dataN_rho2.pdf",width=6,height=4,pointsize=10)
par(mfrow=c(1,2)) 
par(mar=c(4,4,4,2))
#A
symbols(DataN$death_cause, DataN$rho2, circles=sqrt(1/DataN$VZr),inches=0.4, xlab="Predation as only source of mortality",ylab="effect size [r]",main="A",xlim=c(0.5,2.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2), labels=c("No","Yes"))
abline(h=0,lty=3)
#B
symbols(DataN$Nest_openness, DataN$rho2, circles=sqrt(1/DataN$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="B",xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
dev.off()

pdf(file="Fig_bubble2_dataN_v2_rho2.pdf",width=6,height=4,pointsize=10)
par(mfrow=c(1,2)) 
par(mar=c(4,4,4,2))
DataNy <- subset(DataN, DataN$death_cause=="yes", select=c(Nest_openness, rho2, VZr))
DataNn <- subset(DataN, DataN$death_cause=="no", select=c(Nest_openness, rho2, VZr))
#A
symbols(DataNy$Nest_openness, DataNy$rho2, circles=sqrt(1/DataNy$VZr),inches=0.4, xlab="Nest openness",ylab="effect size [r]",main="A. Only mortality from predation",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
#B
symbols(DataNn$Nest_openness, DataNn$rho2, circles=sqrt(1/DataNn$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="B. Mortality from all sources",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
dev.off()


pdf(file="Fig_bubble4_dataN_rho2.pdf",width=6,height=6,pointsize=10)
par(mfrow=c(2,2)) 
par(mar=c(4,4,4,2))
#A
symbols(DataN$death_cause, DataN$rho2, circles=sqrt(1/DataN$VZr),inches=0.4, xlab="Predation as only source of mortality",ylab="effect size [r]",main="A. Mortality sources",xlim=c(0.5,2.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2), labels=c("No","Yes"))
abline(h=0,lty=3)
#B
symbols(DataN$Nest_openness, DataN$rho2, circles=sqrt(1/DataN$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="B. Nest openness",xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
#C
DataNy <- subset(DataN, DataN$death_cause=="yes", select=c(Nest_openness, rho2, VZr))
symbols(DataNy$Nest_openness, DataNy$rho2, circles=sqrt(1/DataNy$VZr),inches=0.4, xlab="Nest openness",ylab="effect size [r]",main="C. Only mortality from predation",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
#D
DataNn <- subset(DataN, DataN$death_cause=="no", select=c(Nest_openness, rho2, VZr))
symbols(DataNn$Nest_openness, DataNn$rho2, circles=sqrt(1/DataNn$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="D. Mortality from all sources",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
dev.off()


### PLOT - Fig_MA_MR_dataA and Fig_MA_MR_dataN subsets in one figure (skip min and max scores)

pdf(file="Fig_MA_MR_dataA_dataN_rho2.pdf",width=9,height=5,pointsize=10)
par(mfcol=c(1,2)) 
par(mar=c(4,10,2,0))

revresAA <- revresA[c(15:length(revresA$M)), ]
plot(revresAA$M, 1:length(revresAA$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.6,0.6), ylim=c(0.25, length(revresAA$M)+.5), xlab="effect size [r]", pch=revresAA$pch, cex=1.1, cex.axis=.9, main="A.   Artificial nests")
abline(v=0,lty=3)
mtext(revresAA$Model, 2, 9, at=1:length(revresAA$M), las=2, cex=.8, font=1, adj=0)
segments(revresAA$CI.lb, 1:length(revresAA$M), revresAA$CI.ub, 1:length(revresAA$M), lwd=1.25)

revresNN <- revresN[c(17:length(revresN$M)), ]
plot(revresNN$M, 1:length(revresNN$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.6,0.6), ylim=c(0.25, length(revresNN$M)+.5), xlab="effect size [r]", pch=revresNN$pch, cex=1.1, cex.axis=.9, main="B.   Natural nests")
abline(v=0,lty=3)
mtext(revresNN$Model, 2, 9, at=1:length(revresNN$M), las=2, cex=.8, font=1, adj=0)
segments(revresNN$CI.lb, 1:length(revresNN$M), revresNN$CI.ub, 1:length(revresNN$M), lwd=1.25)
dev.off()

#mtext("a)",side=2,line=7,at=14,las=2)

###### PUBLICATION BIAS

MR_nest_intc <- rma.mv(yi=Zr2, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~Nest_type-1, method="REML", data=Data)
summary(MR_nest_intc) # in artificial nests signif negative r (more predation with incr urbanization)
# for natural nests ns positive r (less predation with incr urbanization)
Residuals <- residuals(MR_nest_intc)
Precision <- sqrt(1/MR_nest_intc$vi)
plot(Residuals,Precision, xlim=c(-1,1), xlab="Residuals", ylab="Precision [1/SE]")
abline(v=0,lty=3)

model <- rma(yi=Residuals,sei=1/Precision) 
summary(model) # ns est -0.0087 CI -0.0564 0.0390  
funnel(model,yaxi="seinv")
#Trim and fill
TF <- trimfill(model) 
TF #Estimated number of missing studies on the right side: 0 (SE = 6.2815)
funnel(TF)

#Egger's regression test
regtest(model,model="lm") #test for funnel plot asymmetry: t = 0.7158, df = 115, p = 0.4756
ranktest(model) #Kendall's tau = -0.0540, p = 0.3897 (warning about ties, do not use)

### PLOT - funnel plots
pdf(file="Fig_funnels2_alldata_rho2.pdf",width=6,height=4,pointsize=10)
par(mfcol=c(1,2)) 
par(mar=c(4,4,2,1))
plot(Data$rho2, sqrt(1/Data$variance), xlim=c(-1.2,1.2), xlab = "Effect size [r]", ylab="Precision [1/SE]", main="A")
abline(v=0,lty=3)
plot(Residuals,Precision, xlim=c(-1.2,1.2), xlab = "Residuals", ylab="Precision [1/SE]", main="B")
abline(v=0,lty=3)
dev.off()









###########################################################
# ANALYSES on rho (original urbanisation scores)
###########################################################

###########################################################
#### Simple meta-analysis - all data

## MA without random effects
MA_0 <- rma(yi=Zr, vi=VZr, method="REML", data=Data)
summary(MA_0) #I^2 (total heterogeneity / total variability):   91.0%

#funnel plot
plot(jitter(Data$Zr), jitter(sqrt(1/Data$VZr)), cex=0.75, xlim=c(-1.1,1.1), xlab="Zr", ylab="Precision (1/SE.Zr)", main="all data")
abline(v=0, lwd=0.5)

## MA with random effects only
MA_all <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), method="REML", data=Data)
summary(MA_all) #estimate ns Zr0.0063 CI -0.0723   0.0849 
funnel(MA_all) # 1 outlier with large error
#forest(MA_all)
# calculate I2 values
s2m <- sum(1/Data$VZr) * (MA_all$k-1) / (sum(1/Data$VZr)^2 - sum((1/Data$VZr)^2)) # typical sampling error
s2t <- sum(MA_all$sigma2)+s2m # total sampling error
I2t <- sum(MA_all$sigma2)/s2t
I2t*100 # 92.8% - total heterogeneity
I2s <- MA_all$sigma2[1]/s2t
I2s*100 # 84.5% varaince due to Study
I2e <- MA_all$sigma2[2]/s2t
I2e*100 # 8.3% - residuals against sampling error
transf.ztor(MA_all$b[1,]) # intrecept as r = -0.007
transf.ztor(MA_all$ci.lb[1]) # CI.lb for the intercept as r = -0.085
transf.ztor(MA_all$ci.ub[1]) # CI.ub for the intercept as r = 0.071
res <- data.frame(Model="Meta-analytic mean, all data", M=MA_all$b, CI.lb=MA_all$ci.lb, CI.ub=MA_all$ci.ub, pch=18) #harvest results (Zr)


#### Meta-regression on all data with nest_type as moderator
MR_nest_intc <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~Nest_type-1, method="REML", data=Data)
summary(MR_nest_intc) # in artificial nests signif negative Zr -0.1238 (more predation with incr urbanization); 
# in natural nests ns positive Zr 0.0802 (less predation with incr urbanization)
funnel(MR_nest_intc) # 1 outlier ES < -1 (same point as in the full data set)
forest(MR_nest_intc)
summ <- summary(MR_nest_intc)
forest(x=summ$b, ci.lb=summ$ci.lb, ci.ub=summ$ci.ub, slab=dimnames(summ$b)[[1]],  xlab="Zr")
# calculate I2 values
s2m <- sum(1/Data$VZr) * (MR_nest_intc$k-1) / (sum(1/Data$VZr)^2 - sum((1/Data$VZr)^2)) # typical sampling error
s2t <- sum(MR_nest_intc$sigma2)+s2m # total sampling error
I2t <- sum(MR_nest_intc$sigma2)/s2t
I2t*100 # 91.9% - total heterogeneity
I2s <- MR_nest_intc$sigma2[1]/s2t
I2s*100 # 82.7% varaince due to Study
I2e <- MR_nest_intc$sigma2[2]/s2t
I2e*100 # 9.1% - residuals against sampling error
transf.ztor(MR_nest_intc$b[1,]) #-0.123
transf.ztor(MR_nest_intc$ci.lb[1]) #-0.232
transf.ztor(MR_nest_intc$ci.ub[1]) #-0.012
transf.ztor(MR_nest_intc$b[2,]) #0.080
transf.ztor(MR_nest_intc$ci.lb[2]) #-0.017
transf.ztor(MR_nest_intc$ci.ub[2]) #0.175
res <- rbind(res, data.frame(Model=c("Artificial vs. Natural nests:","   Artificial nests *","   Natural nests"), M=c(NA,MR_nest_intc$b), CI.lb=c(NA,MR_nest_intc$ci.lb), CI.ub=c(NA,MR_nest_intc$ci.ub), pch=c(20,20,20))) #harvest results (Zr)

MR_nest_diff <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~Nest_type, method="REML", data=Data)
summary(MR_nest_diff)  # diff between artificial and natural nests signif Zr 0.203 CI 0.056  0.352
transf.ztor(MR_nest_diff$b[2,]) #0.201
transf.ztor(MR_nest_diff$ci.lb[2]) #0.056
transf.ztor(MR_nest_diff$ci.ub[2]) #0.338
res <- rbind(res, data.frame(Model="   Difference: Natural - Artificial nests *", M=MR_nest_diff$b[2], CI.lb=MR_nest_diff$ci.lb[2], CI.ub=MR_nest_diff$ci.ub[2], pch=20)) #harvest results (Zr)

#tidy up extracted results table
res$M <- round(transf.ztor(res$M),3)
res$CI.lb <- round(transf.ztor(res$CI.lb),3)
res$CI.ub <- round(transf.ztor(res$CI.ub),3)
write.csv(res,"MA_MR_alldata_main_res.csv")
revres <- res[rev(rownames(res)),] #reverse for plotting (from bottom to top)


### PLOT - MA and MR models on all data
#opar <- par()      # make a copy of current settings
#par(opar)          # restore original settings

pdf(file="Fig_MA_MR_alldata.pdf",width=4,height=3,pointsize=8)
par(mfrow=c(1,1))
par(mar=c(5,14,1,1))
plot(revres$M, 1:length(revres$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.5,0.5), ylim=c(0.25, length(revres$M)+.5), xlab="effect size [r]", pch=revres$pch, cex=1, cex.axis=.9)
abline(v=0,lty=3)
mtext(revres$Model, 2, 13, at=1:length(revres$M), las=2, cex=1, font=1, adj=0)
segments(revres$CI.lb, 1:length(revres$M), revres$CI.ub, 1:length(revres$M), lwd=1.25)
dev.off()


#################### ARTIFICIAL NESTS ################

DataA <- Data[Data$Nest_type == "art", ]
dim(DataA) # 59 31
DataA <- droplevels(DataA)

unlist(lapply(DataA, function(x) sum(is.na(x)))) #how many NA per column: Study_days 18, Nest_height 1

#### MA
#funnel plot
plot(jitter(DataA$Zr), jitter(sqrt(1/DataA$VZr)), cex=0.75, xlim=c(-1.1,1.1), xlab="Zr", ylab="Precision (1/SE.Zr)", main="all data")
abline(v=0, lwd=0.5)

## MA with random effects only
MA_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), method="REML", data=DataA)
summary(MA_A) #estimate signif Zr -0.1256 CI  -0.2498   -0.0015 
funnel(MA_A)
# calculate I2 values
s2m <- sum(1/DataA$VZr) * (MA_A$k-1) / (sum(1/DataA$VZr)^2 - sum((1/DataA$VZr)^2)) # typical sampling error
s2t <- sum(MA_A$sigma2)+s2m # total sampling error
I2t <- sum(MA_A$sigma2)/s2t
I2t*100 # 93.1% - total heterogeneity
I2s <- MA_A$sigma2[1]/s2t
I2s*100 # 86.0% varaince due to Study
I2e <- MA_A$sigma2[2]/s2t
I2e*100 # 7.1% - residuals against sampling error
transf.ztor(MA_A$b[1,]) #-0.1250
transf.ztor(MA_A$ci.lb[1]) #-0.2447
transf.ztor(MA_A$ci.ub[1]) #-0.0015
#sum(MA_A$sigma2)/(sum(MA_A$sigma2)+(sum(1/DataA$VZr)*(MA_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity
resA <- data.frame(Model="Meta-analytic mean*", M=MA_A$b, CI.lb=MA_A$ci.lb, CI.ub=MA_A$ci.ub, pch=18) #harvest results (Zr)


#### MR - Meta-Regressions (univariate)

# univariate meta-regression with Nest_openness 
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_openness)-1, method="REML", data=DataA)
summary(MR_A) # cup signif Zr  -0.1269  -0.2509  -0.0029 
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.1%
plot(DataA$Zr ~ DataA$Nest_openness)
resA <- rbind(resA, data.frame(Model=c("Nest openness:","   Cup *","   Hole"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20))) #harvest results (Zr)

# univariate meta-regression with Nest_openness Hole-Cup difference
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_openness), method="REML", data=DataA)
summary(MR_A) # diff Zr ns 0.0805  -0.0459   0.2070 

# univariate meta-regression with Nest_position
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_position)-1, method="REML", data=DataA)
summary(MR_A) # only ground Zr signif -0.1376 -0.2697  -0.0054
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.0%
plot(DataA$Zr ~ DataA$Nest_position)
resA <- rbind(resA, data.frame(Model=c("Nest position:","   Elevated","   Ground *","   Mix"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# univariate meta-regression with Egg_number
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(DataA$Egg_number), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.5%
plot(DataA$Zr ~ DataA$Egg_number)
resA <- rbind(resA, data.frame(Model=c("Egg number (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with Study_days (NOTE: 18 NA)
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(Study_days), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 95.4%
plot(DataA$Zr ~ DataA$Study_days)
resA <- rbind(resA, data.frame(Model=c("Study duration (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with Median_year
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(Median_year), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.2%
plot(DataA$Zr ~ DataA$Median_year)
resA <- rbind(resA, data.frame(Model=c("Median study year (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with Pub_year
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(DataA$Pub_year), method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.5%
plot(DataA$Zr ~ DataA$Pub_year)
resA <- rbind(resA, data.frame(Model=c("Publication year (slope)"), M=c(MR_A$b[2]), CI.lb=c(MR_A$ci.lb[2]), CI.ub=c(MR_A$ci.ub[2]), pch=c(20))) #harvest results (Zr)

# univariate meta-regression with urbmin as a factor
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(urbmin)-1, method="REML", data=DataA)
summary(MR_A) # all ns
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 94.2
plot(DataA$Zr ~ DataA$urbmin)
resA <- rbind(resA, data.frame(Model=c("Min urbanisation score:","   1","   2","   3","   4"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20,20,20))) #harvest results (Zr)
# univariate meta-regression with urbmin as continuous - not used becouse the relationship is driven by just 2 data points at 4

# univariate meta-regression with urbmax as a factor
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(urbmax)-1, method="REML", data=DataA)
summary(MR_A) # urbmax=4 signif Zr -0.3363 -0.6413  -0.0314, but driven by only 4 data points! (most are at urbanmax=3)
sum(MR_A$sigma2)/(sum(MR_A$sigma2)+(sum(1/DataA$VZr)*(MR_A$k-1)/(sum(1/DataA$VZr)^2-sum((1/DataA$VZr)^2))))*100 # total heterogeneity 93.2%
plot(DataA$Zr ~ DataA$urbmax)
resA <- rbind(resA, data.frame(Model=c("Max urbanisation score:","   3","   4 *","   5"), M=c(NA,MR_A$b), CI.lb=c(NA,MR_A$ci.lb), CI.ub=c(NA,MR_A$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# meta-regression with urbmin_scaled*urbmax_scaled interaction (as continuous predictors)
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(urbmin) * scale(urbmax), method="REML", data=DataA)
summary(MR_A) # ns
# meta-regression with Nest_openness*Nest_position interaction
MR_A <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_openness) * as.factor(Nest_position), method="REML", data=DataA)
summary(MR_A) # ns
table(DataA$Nest_openness, DataA$Nest_position) #not enough data for estimating interaction: 0 hole-like nests at ground and mix levels
boxplot(DataA$Zr ~  DataA$Nest_position + DataA$Nest_openness)

# meta-regression with Nest_openness + Nest_position without interaction
MR_A1 <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_position) + as.factor(Nest_openness) -1, method="REML", data=DataA)
summary(MR_A1) # all ns
#resA <- rbind(resA, data.frame(Model=c("Multivariate:","   Cup - elevated","   Cup - ground","   Cup - mix","   Hole - elevated"), M=c(NA,MR_A1$b), CI.lb=c(NA,MR_A1$ci.lb), CI.ub=c(NA,MR_A1$ci.ub), pch=c(20,20,20,20,20))) #harvest results (Zr)
resA <- rbind(resA, data.frame(Model=c("Multivariate:","   Cup - elevated","   Cup - ground","   Cup - mix"), M=c(NA,MR_A1$b[1:3]), CI.lb=c(NA,MR_A1$ci.lb[1:3]), CI.ub=c(NA,MR_A1$ci.ub[1:3]), pch=c(20,20,20,20))) #harvest results (Zr)
MR_A2 <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~as.factor(Nest_position) + relevel(as.factor(Nest_openness), ref="hole") -1, method="REML", data=DataA)
summary(MR_A2) # all ns
resA <- rbind(resA, data.frame(Model="   Hole - elevated", M=c(MR_A2$b[1]), CI.lb=c(MR_A2$ci.lb[1]), CI.ub=c(MR_A2$ci.ub[1]), pch=c(20))) #harvest results (Zr)


# tidy up extracted results table
resA$M <- round(transf.ztor(resA$M),3)
resA$CI.lb <- round(transf.ztor(resA$CI.lb),3)
resA$CI.ub <- round(transf.ztor(resA$CI.ub),3)
write.csv(resA,"MA_MR_dataA_resA.csv")
revresA <- resA[rev(rownames(resA)),] #reverse for plotting (from bottom to top)



### PLOT - MA and MR models on all data
opar <- par()      # make a copy of current settings
#par(opar)          # restore original settings

pdf(file="Fig_MA_MR_dataA.pdf",width=4,height=6,pointsize=10)
par(mfrow=c(1,1))
par(mar=c(4,10,2,0))
plot(revresA$M, 1:length(revresA$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.6,0.6), ylim=c(0.25, length(revresA$M)+.5), xlab="Zr", pch=revresA$pch, cex=1.1, cex.axis=.9)
abline(v=0,lty=3)
mtext(revresA$Model, 2, 9, at=1:length(revresA$M), las=2, cex=.8, font=1, adj=0)
segments(revresA$CI.lb, 1:length(revresA$M), revresA$CI.ub, 1:length(revresA$M), lwd=1.25)
dev.off()



#################### NATURAL NESTS ################

DataN <- Data[Data$Nest_type == "nat", ]
dim(DataN) # 58 31
DataN <- droplevels(DataN)

unlist(lapply(DataN, function(x) sum(is.na(x)))) #how many NA per column: Study_days 31 NA, Nest_height 30 NA, Egg_number 16 NA

#### MA
#funnel plot
plot(jitter(DataN$Zr), jitter(sqrt(1/DataN$VZr)), cex=0.75, xlim=c(-1.1,1.1), xlab="Zr", ylab="Precision (1/SE.Zr)", main="all data")
abline(v=0, lwd=0.5)

## MA with 2 random effects
MA_N <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), method="REML", data=DataN)
summary(MA_N) #estimate Zr ns 0.0790 CI  -0.0083   0.1664 
# overall ns positive relationship between urbanisation and survival - less predation with increasing urbanization
funnel(MA_N)
forest(MA_N)
# calculate I2 values
s2m <- sum(1/DataN$VZr) * (MA_N$k-1) / (sum(1/DataN$VZr)^2 - sum((1/DataN$VZr)^2)) # typical sampling error
s2t <- sum(MA_N$sigma2)+s2m # total sampling error
I2t <- sum(MA_N$sigma2)/s2t
I2t*100 # 89.97% - total heterogeneity
I2s <- MA_N$sigma2[1]/s2t
I2s*100 # 77.43% varaince due to Study
I2e <- MA_N$sigma2[2]/s2t
I2e*100 # 12.54% - residuals against sampling error
transf.ztor(MA_N$b[1,]) #0.0788
transf.ztor(MA_N$ci.lb[1]) #-0.0083
transf.ztor(MA_N$ci.ub[1]) #0.1648
resN <- data.frame(Model="Meta-analytic mean", M=MA_N$b, CI.lb=MA_N$ci.lb, CI.ub=MA_N$ci.ub, pch=18) #harvest results (Zr)

## MA with 3 random effects, including Species identity
MA_N <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|Species,~1|ES_ID), method="REML", data=DataN)
summary(MA_N) #estimate Zr ns 0.0671 CI  -0.0118   0.1461 
# overall ns positive relationship between urbanisation and survival - less predation with increasing urbanization
funnel(MA_N)
forest(MA_N)
# calculate I2 values
s2m <- sum(1/DataN$VZr) * (MA_N$k-1) / (sum(1/DataN$VZr)^2 - sum((1/DataN$VZr)^2)) # typical sampling error
s2t <- sum(MA_N$sigma2)+s2m # total sampling error
I2t <- sum(MA_N$sigma2)/s2t
I2t*100 # 88.36% - total heterogeneity
I2s <- MA_N$sigma2[1]/s2t
I2s*100 # 0% varaince due to Study
I2e <- MA_N$sigma2[2]/s2t
I2e*100 # 74.0% varaince due to Species identity - equivalent to using just Study before, usually 1 Species per study, and 1 study per species
I2e <- MA_N$sigma2[3]/s2t
I2e*100 # 14.36% - residuals against sampling error
table(DataN$Species,DataN$Paper_ID)

## Phylogeny
birds <- unique(DataN$Species_latin) #get list of unique bird species latin names
str(birds)
birds_stree <- read.tree("Ericson.tre") # load birds supertree from a file
birds_stree #9993 tips = species
str(birds_stree) #has edge (branch) lengths
birds_steer <- collapse.singles(birds_stree)
bird_tree_species <- as.character(birds_stree$tip.label) # extract list of species from the tree
intersect(bird_tree_species, birds) #31 matching
setdiff(birds, bird_tree_species) #Cyanistes_caeruleus not matching 
birds_stree$tip.label <- sub("Parus_caeruleus","Cyanistes_caeruleus",birds_stree$tip.label) #replace with synonym Parus_caeruleus
bird_tree_species <- as.character(birds_stree$tip.label) # extract list of species from the tree
intersect(bird_tree_species, birds) #32 matching, ok
tree <- drop.tip(birds_stree, birds_stree$tip.label[-match(birds, birds_stree$tip.label)]) # prune the supertree tree to a list of taxa from our list
is.binary.tree(tree) #TRUE
is.ultrametric(tree) #TRUE
plot(tree, cex=0.8) #plot with branch lengths
nodelabels()
tree2 <- rotate(tree, 37) #rotate some branches

pdf("tree.pdf",width = 4, height = 5, pointsize=8)
par(mfcol=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(tree2)
dev.off()

write.tree(tree2, file = "birds_32sp_tree.tre", append = FALSE, digits = 10, tree.names = FALSE)

### phylogenetic meta-analysis (does not converge if Species ID added)
tree <- read.tree("birds_32sp_tree.tre")# upload cleaned-up and prepeprocessed phylogenetic tree file
CorMatrix <- vcv(tree, corr=TRUE) # make a phylogenetic correlation matrix

MA_N_phylo <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MA_N_phylo) #intc = 0.0446 ns, CI -0.1264   0.2155 
sum(MA_N_phylo$sigma2)/(sum(MA_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MA_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 91.5%
s2t <- sum(MA_N_phylo$sigma2) + sum(1/MA_N_phylo$vi) * (MA_N_phylo$k-1) / (sum(1/MA_N_phylo$vi)^2 - sum((1/MA_N_phylo$vi)^2))
MA_N_phylo$sigma2[1]/s2t*100 # I^2 between study (id) = 23.9%
MA_N_phylo$sigma2[2]/s2t*100 # I^2 phylogeny = 56.9%
MA_N_phylo$sigma2[3]/s2t*100 # I^2 within study or residual = 10.67%
resN <- rbind(resN, data.frame(Model=c("Phylogenetic meta-analytic mean"), M=c(MA_N_phylo$b), CI.lb=c(MA_N_phylo$ci.lb), CI.ub=c(MA_N_phylo$ci.ub), pch=c(18))) #harvest results (Zr)

# plots
Names <- paste(DataN$Paper, DataN$Species, sep="_")
forest(MA_N_phylo, slab=Names)
funnel(MA_N_phylo, yaxis="seinv")

#### Meta-regression 

# species meta-regression
MA_N_species <- rma.mv(yi=Zr, V=VZr, mod=~Species_latin-1, random=list(~1|Paper_ID,~1|ES_ID), data=DataN, method="REML") #without phylogeny
summary(MA_N_species) # some species very signif + or -
sum(MA_N_species$sigma2)/(sum(MA_N_species$sigma2)+(sum(1/DataN$VZr)*(MA_N_species$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 58.1%
s2t <- sum(MA_N_species$sigma2) + sum(1/MA_N_species$vi) * (MA_N_species$k-1) / (sum(1/MA_N_species$vi)^2 - sum((1/MA_N_species$vi)^2))
MA_N_species$sigma2[1]/s2t # I^2 between study (id) = 0.00% (studies are almost equivalent to the species)
MA_N_species$sigma2[2]/s2t # I^2 within study or residual = 58.1%
# plots
Names <- paste(DataN$Paper, DataN$Species, sep="_")
forest(MA_N_phylo, slab=Names)
funnel(MA_N_phylo, yaxis="seinv")
plot(tree)
# extract point estimates and their CI for each species
ssp <- summary(MA_N_species)
sp_df <- data.frame(sp = substr(attr(ssp$b,"dimnames")[[1]], 14, 36), M = ssp$b, CI.lb = ssp$ci.lb, CI.ub = ssp$ci.ub)
sp_df <- sp_df[match(tree$tip.label, sp_df$sp),] #reorder dataframe to match order of the tip labels on the tree
# tidy up extracted results table
sp_df$M <- round(transf.ztor(sp_df$M),3)
sp_df$CI.lb <- round(transf.ztor(sp_df$CI.lb),3)
sp_df$CI.ub <- round(transf.ztor(sp_df$CI.ub),3)
sp_df$Signif <- ifelse(sp_df$CI.lb>0 | sp_df$CI.ub<0, "*", "") #add column with stars for estimates that significantly differ from zero
write.csv(sp_df,"MA_MR_dataN_sp_df.csv")
revres_sp <- sp_df[rev(rownames(sp_df)),] #reverse for plotting (from bottom to top)

opar <- par()      # make a copy of current settings
par(opar)          # restore original settings

### PLOT - species phylogeny and forest plot
pdf(file="Fig_species.pdf",width=6,height=4,pointsize=9)
par(mfrow=c(1,2))
par(mar=c(4.5,2,1.35,0))
plot(tree, font=1, cex=0.8, x.lim=90, show.tip.label = FALSE)
par(mar=c(4,6,1,0))
plot(sp_df$M, 1:length(sp_df$M), ylab=NA, yaxt="n", bty="n", xlim=c(-1.5,1.5), ylim=c(0.25, length(sp_df$M)+.5), xlab="effect size [r]", pch=16, cex=0.8, cex.axis=.9)
abline(v=0,lty=3)
mtext(sp_df$sp, 2, -1, at=1:length(sp_df$M), las=2, cex=.8, font=3)
segments(sp_df$CI.lb, 1:length(sp_df$M), sp_df$CI.ub, 1:length(sp_df$M), lwd=1.25)
for (i in 1:length(rownames(sp_df))) mtext(sp_df$Signif[i], 2, -1.5, at=i, las=2, cex=.7) #add stars
dev.off()


# Phylogenetic meta-regression with death_cause
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # yes signif 0.1615 CI 0.0501  0.2729
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 88.4
plot(DataN$Zr ~ DataN$death_cause) # positive values (less predation in urbanised areas) when "yes" = only nests that died from predation are included
resN <- rbind(resN, data.frame(Model=c("Predation as only source of mortality:","   No","   Yes *"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20))) #harvest results (Zr)

MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # signif diff 0.1780 CI 0.0155  0.3404  *
#MR_N <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause), random=list(~1|Paper_ID,~1|ES_ID), data=DataN, method="REML") #without phylogeny
#summary(MR_N) # signif diff 0.1780 CI 0.0155  0.3404  * - same result without phylogeny

# Phylogenetic meta-regression with Nest_openness
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(Nest_openness)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # hole signif 0.2227 CI 0.0320  0.4133
#MR_N <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(Nest_openness)-1, random=list(~1|Paper_ID,~1|ES_ID), data=DataN, method="REML") #without phylogeny
#summary(MR_N) # hole signif 0.2227 CI 0.0320  0.4133 - same result without phylogeny
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 89.6
plot(DataN$Zr ~ DataN$Nest_openness) 
resN <- rbind(resN, data.frame(Model=c("Nest openness:","   Cup","   Hole *","   Orb"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# Phylogenetic meta-regression with Nest_position
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(Nest_position)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 90.6
plot(DataN$Zr ~ DataN$Nest_position) # more positive values in ground-located nests
resN <- rbind(resN, data.frame(Model=c("Nest position:","   Elevated","   Ground","   Mix"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# Phylogenetic meta-regression with Nest_height (NOTE: 31 NA, 23 present) 
# MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~scale(Nest_height), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
# summary(MR_N_phylo) # scale(Nest_height) signif slope -0.1347 CI -0.2677  -0.0018 - in low-placed nests tendency for positive Zr, in high nests negative Zr
# sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 91.92
# plot(DataN$Zr ~ DataN$Nest_height) # positive values more likely in lower nests, negative values in nests higher than 6m

# Phylogenetic meta-regression with log(Nest_height+1)
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~scale(log(Nest_height+1)), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # similar ns trend as above -0.0991 -0.2099  0.0117
plot(DataN$Zr ~ log(DataN$Nest_height+1))
resN <- rbind(resN, data.frame(Model="Nest height above ground (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

hist(DataN$Nest_height[DataN$death_cause=="yes"], col="blue", xlim=c(0,25))
hist(DataN$Nest_height[DataN$death_cause=="no"], col=rgb(0, 0, 0, 0.5), add=TRUE)
table(DataN$death_cause, DataN$Nest_height) # measure predation = "no" in nest height > 10
t.test(DataN$Nest_height ~ DataN$death_cause, var.equal = TRUE) #t = 3.3704, df = 26, p-value = 0.002354 
#high nests usually report overall mortrality, not by predation
table(is.na(DataN$Nest_height), DataN$Nest_position) # 30 missing height values are from elevated nests, i.e. only a subset of heights known

# Phylogenetic meta-regression with Egg_number (NOTE: 16 NA) 
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~scale(DataN$Egg_number), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 94.4
plot(DataN$Zr ~ DataN$Egg_number) 
resN <- rbind(resN, data.frame(Model="Egg number (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

# Phylogenetic meta-regression with Study_days (NOTE: 28 NA) 
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~scale(DataN$Study_days), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 98.2
plot(DataN$Zr ~ DataN$Study_days) 
resN <- rbind(resN, data.frame(Model="Study duration (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

# Phylogenetic meta-regression with Median_year as moderator
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~scale(DataN$Median_year), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 90.6
plot(DataN$Zr ~ DataN$Median_year) # all ns
resN <- rbind(resN, data.frame(Model="Median study year (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

hist(DataN$Median_year[DataN$death_cause=="yes"], col="blue")
hist(DataN$Median_year[DataN$death_cause=="no"], col=rgb(0, 0, 0, 0.5), add=TRUE)
table(DataN$death_cause, DataN$Median_year) # recent studies more likely to only measure predation = "yes"
t.test(DataN$Median_year ~ DataN$death_cause, var.equal = FALSE) #t = -2.3125, df = 26.758, p-value = 0.02868

# Phylogenetic meta-regression with Pub_year as moderator
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~scale(DataN$Pub_year), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 90.2
plot(DataN$Zr ~ DataN$Pub_year) 
resN <- rbind(resN, data.frame(Model="Publication year (slope)", M=MR_N_phylo$b[2], CI.lb=MR_N_phylo$ci.lb[2], CI.ub=MR_N_phylo$ci.ub[2], pch=20)) #harvest results (Zr)

# Phylogenetic meta-regression with urbmin as a factor
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(urbmin)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 94.2
plot(DataN$Zr ~ DataN$urbmin) # only one data point at urbmin=4
resN <- rbind(resN, data.frame(Model=c("Min urbanisation score:","   1","   2","   3","   4"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20,20))) #harvest results (Zr)

# Phylogenetic meta-regression with urbmax as a factor
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(urbmax)-1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # all ns
sum(MR_N_phylo$sigma2)/(sum(MR_N_phylo$sigma2)+(sum(1/DataN$VZr)*(MR_N_phylo$k-1)/(sum(1/DataN$VZr)^2-sum((1/DataN$VZr)^2))))*100 # total heterogeneity 90.6
plot(DataN$Zr ~ DataN$urbmax) 
resN <- rbind(resN, data.frame(Model=c("Min urbanisation score:","   3","   4","   5"), M=c(NA,MR_N_phylo$b), CI.lb=c(NA,MR_N_phylo$ci.lb), CI.ub=c(NA,MR_N_phylo$ci.ub), pch=c(20,20,20,20))) #harvest results (Zr)

# meta-regression with urbmin_scaled*urbmax_scaled interaction (as continuous predictors)
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~scale(urbmin) * scale(urbmax), method="REML", data=DataN)
summary(MR_N_phylo) # ns interaction
# meta-regression with Death_cause and Nest_openness, without interaction 
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause) + as.factor(Nest_openness) -1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # (death_cause)yes signif 0.1711  0.0570  0.2851
# meta-regression with Death_cause*Nest_openness interaction
MR_N_phylo <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause) * as.factor(Nest_openness), random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo) # signif interaction
MR_N_phylo1 <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause) * as.factor(Nest_openness) - 1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo1) # signif interaction
resN <- rbind(resN, data.frame(Model=c("Multivariate meta-regression:","   Cup - No","   Cup - Yes *"), M=c(NA,MR_N_phylo1$b[1:2]), CI.lb=c(NA,MR_N_phylo1$ci.lb[1:2]), CI.ub=c(NA,MR_N_phylo1$ci.ub[1:2]), pch=c(20,20,20))) #harvest results (Zr)
MR_N_phylo2 <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause) * as.factor(relevel(Nest_openness, ref="hole")) - 1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo2) # signif interaction
resN <- rbind(resN, data.frame(Model=c("   Hole - No","   Hole - Yes *"), M=c(MR_N_phylo2$b[1:2]), CI.lb=c(MR_N_phylo2$ci.lb[1:2]), CI.ub=c(MR_N_phylo2$ci.ub[1:2]), pch=c(20,20))) #harvest results (Zr)
MR_N_phylo3 <- rma.mv(yi=Zr, V=VZr, mod= ~as.factor(death_cause) * as.factor(relevel(Nest_openness, ref="orb")) - 1, random=list(~1|Paper_ID,~1|Species_latin,~1|ES_ID), R=list(Species_latin=CorMatrix), data=DataN, method="REML") #with phylogeny
summary(MR_N_phylo3) # signif interaction
resN <- rbind(resN, data.frame(Model=c("   Orb - Yes"), M=c(MR_N_phylo3$b[2]), CI.lb=c(MR_N_phylo3$ci.lb[2]), CI.ub=c(MR_N_phylo3$ci.ub[2]), pch=c(20))) #harvest results (Zr)
resN <- rbind(resN, data.frame(Model=c("  Nest openness - Death cause interaction *"), M=c(MR_N_phylo3$b[5]), CI.lb=c(MR_N_phylo3$ci.lb[5]), CI.ub=c(MR_N_phylo3$ci.ub[5]), pch=c(20))) #harvest results (Zr)

table(DataN$death_cause, DataN$Nest_openness) # only 5 data points are hole/yes, 2 are are orb/yes, 0 are orb/no
boxplot(DataN$Zr ~ DataN$death_cause * DataN$Nest_openness, varwidth=TRUE) #most positive values in hole/yes subest - less predation with increasing urbanisation

# tidy up extracted results table
resN$M <- round(transf.ztor(resN$M),3)
resN$CI.lb <- round(transf.ztor(resN$CI.lb),3)
resN$CI.ub <- round(transf.ztor(resN$CI.ub),3)
write.csv(resN,"MA_MR_dataN_resN.csv")
revresN <- resN[rev(rownames(resN)),] #reverse for plotting (from bottom to top)


### PLOT - MA and MR models on all data
opar <- par()      # make a copy of current settings
#par(opar)          # restore original settings

pdf(file="Fig_MA_MR_dataN.pdf",width=4,height=6,pointsize=10)
par(mfrow=c(1,1))
par(mar=c(4,10,2,0))
plot(revresN$M, 1:length(revresN$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.8,0.8), ylim=c(0.25, length(revresN$M)+.5), xlab="effect size [r]", pch=revresN$pch, cex=1.1, cex.axis=.9)
abline(v=0,lty=3)
mtext(revresN$Model, 2, 9, at=1:length(revresN$M), las=2, cex=.8, font=1, adj=0)
segments(revresN$CI.lb, 1:length(revresN$M), revresN$CI.ub, 1:length(revresN$M), lwd=1.25)
dev.off()

## OVERALL for natural nests: tendency for less predation with increasing urbanisation (positive Zr) - small overall effect size, 
# this effect is more pronounced when only nest lost due to predation only included
# and no much effect of urbanisation when other causes of mortality (confounding) are potentially present.
# Tendency for more negative valus in nests higher above ground  - more mortality in urbanised areas, 
# in lower nests more positive values more likely - less predation in urbanised areas, 
# but this result is likely to be related on death_cause variable.



### PLOT - bubble plots for natural nests: r-death_cause, r-nest_height
plot(DataN$Zr ~ DataN$death_cause) # positive values (less predation in urbanised areas) when "yes" = only nests that died from predation are included
plot(DataN$Zr ~ DataN$Nest_height) # K=23

pdf(file="Fig_bubble2_dataN.pdf",width=6,height=4,pointsize=10)
par(mfrow=c(1,2)) 
par(mar=c(4,4,4,2))
#A
symbols(DataN$death_cause, DataN$rho, circles=sqrt(1/DataN$VZr),inches=0.4, xlab="Predation as only source of mortality",ylab="effect size [r]",main="A",xlim=c(0.5,2.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2), labels=c("No","Yes"))
abline(h=0,lty=3)
#B
symbols(DataN$Nest_openness, DataN$rho, circles=sqrt(1/DataN$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="B",xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
dev.off()

pdf(file="Fig_bubble2_dataN_v2.pdf",width=6,height=4,pointsize=10)
par(mfrow=c(1,2)) 
par(mar=c(4,4,4,2))
DataNy <- subset(DataN, DataN$death_cause=="yes", select=c(Nest_openness, rho, VZr))
DataNn <- subset(DataN, DataN$death_cause=="no", select=c(Nest_openness, rho, VZr))
#A
symbols(DataNy$Nest_openness, DataNy$rho, circles=sqrt(1/DataNy$VZr),inches=0.4, xlab="Nest openness",ylab="effect size [r]",main="A. Only mortality from predation",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
#B
symbols(DataNn$Nest_openness, DataNn$rho, circles=sqrt(1/DataNn$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="B. Mortality from all sources",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
dev.off()


pdf(file="Fig_bubble4_dataN.pdf",width=6,height=6,pointsize=10)
par(mfrow=c(2,2)) 
par(mar=c(4,4,4,2))
#A
symbols(DataN$death_cause, DataN$rho, circles=sqrt(1/DataN$VZr),inches=0.4, xlab="Predation as only source of mortality",ylab="effect size [r]",main="A. Mortality sources",xlim=c(0.5,2.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2), labels=c("No","Yes"))
abline(h=0,lty=3)
#B
symbols(DataN$Nest_openness, DataN$rho, circles=sqrt(1/DataN$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="B. Nest openness",xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
#C
DataNy <- subset(DataN, DataN$death_cause=="yes", select=c(Nest_openness, rho, VZr))
symbols(DataNy$Nest_openness, DataNy$rho, circles=sqrt(1/DataNy$VZr),inches=0.4, xlab="Nest openness",ylab="effect size [r]",main="C. Only mortality from predation",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
#D
DataNn <- subset(DataN, DataN$death_cause=="no", select=c(Nest_openness, rho, VZr))
symbols(DataNn$Nest_openness, DataNn$rho, circles=sqrt(1/DataNn$VZr),inches=0.4,xlab="Nest openness",ylab="effect size [r]",main="D. Mortality from all sources",ylim=c(-1,1),xlim=c(0,3.5),xaxp=c(1,2,1), xaxt="n")
axis(1, c(1,2,3), labels=c("Cup","Hole","Orb"))
abline(h=0,lty=3)
dev.off()


### PLOT - Fig_MA_MR_dataA and Fig_MA_MR_dataN subsets in one figure (skip min and max scores)

pdf(file="Fig_MA_MR_dataA_dataN.pdf",width=9,height=5,pointsize=10)
par(mfcol=c(1,2)) 
par(mar=c(4,10,2,0))
revresAA <- revresA[c(15:length(revresA$M)), ]
plot(revresAA$M, 1:length(revresAA$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.6,0.6), ylim=c(0.25, length(revresAA$M)+.5), xlab="effect size [r]", pch=revresAA$pch, cex=1.1, cex.axis=.9, main="A.   Artificial nests")
abline(v=0,lty=3)
mtext(revresAA$Model, 2, 9, at=1:length(revresAA$M), las=2, cex=.8, font=1, adj=0)
segments(revresAA$CI.lb, 1:length(revresAA$M), revresAA$CI.ub, 1:length(revresAA$M), lwd=1.25)
revresNN <- revresN[c(17:length(revresN$M)), ]
plot(revresNN$M, 1:length(revresNN$M), ylab=NA, yaxt="n", bty="n", xlim=c(-0.6,0.6), ylim=c(0.25, length(revresNN$M)+.5), xlab="effect size [r]", pch=revresNN$pch, cex=1.1, cex.axis=.9, main="B.   Natural nests")
abline(v=0,lty=3)
mtext(revresNN$Model, 2, 9, at=1:length(revresNN$M), las=2, cex=.8, font=1, adj=0)
segments(revresNN$CI.lb, 1:length(revresNN$M), revresNN$CI.ub, 1:length(revresNN$M), lwd=1.25)
dev.off()

#mtext("a)",side=2,line=7,at=14,las=2)



###### PUBLICATION BIAS

MR_nest_intc <- rma.mv(yi=Zr, V=VZr, random=list(~1|Paper_ID,~1|ES_ID), mod=~Nest_type-1, method="REML", data=Data)
summary(MR_nest_intc) # in artificial nests signif negative r (more predation with incr urbanization)
# for natural nests ns positive r (less predation with incr urbanization)
Residuals <- residuals(MR_nest_intc)
Precision <- sqrt(1/MR_nest_intc$vi)
plot(Residuals,Precision, xlim=c(-1,1), xlab="Residuals", ylab="Precision [1/SE]")
abline(v=0,lty=3)

model <- rma(yi=Residuals,sei=1/Precision) 
summary(model) # ns est -0.0074 CI -0.0547   0.0400  
funnel(model,yaxi="seinv")
#Trim and fill
TF <- trimfill(model) 
TF #Estimated number of missing studies on the right side: 0 (SE = 6.2815)
funnel(TF)

#Egger's regression test
regtest(model,model="lm") #test for funnel plot asymmetry: t = 0.6458, df = 115, p = 0.5197
ranktest(model) #Kendall's tau = -0.0540, p = 0.3897 (warning about ties, do not use)

### PLOT - funnel plots
pdf(file="Fig_funnels2_alldata.pdf",width=6,height=4,pointsize=10)
par(mfcol=c(1,2)) 
par(mar=c(4,4,2,1))
plot(Data$rho, sqrt(1/Data$variance), xlim=c(-1.2,1.2), xlab = "Effect size [r]", ylab="Precision [1/SE]", main="A")
abline(v=0,lty=3)
plot(Residuals,Precision, xlim=c(-1.2,1.2), xlab = "Residuals", ylab="Precision [1/SE]", main="B")
abline(v=0,lty=3)
dev.off()
