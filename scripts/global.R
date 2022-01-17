############## Global code for the game application
############## Reading in the code for running the stepwise simulation. 
############## 1) I first read in all the code to read data I need to simulate the lice development
############## 2) Initiate the model in the SV object
############## 3) Put the SV object into a mutable R6 object, for easy update of the model
############## 4) The method also includes t
############## ## 5) Method 1 in the R6 object is a stepwise update of R6
############## ## 6) The only parameter in the R6 method is the treatment, for now
############## ## 7) Update t
############## ## 8) Update score

############## Code to simulate salmon lice dynamics in a fish farm with multiple cages
############## based on Aldrin et al. (2017): https://doi.org/10.1016/j.ecolmodel.2017.05.019
############## and Aldrin and Huseby (2019) Norsk Regnesentral report SAMBA 28/19.
############## This code is written by Leif Chr. Stige.


############################################################################

## Clear work space
#rm(list=ls())
#library(R6)
library(ggplot2)
library(tidyverse)
# install.packages("R6")
# library(R6)
## Paths

## Paths
## parameterpath <- "//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/Analyser/Modellfiler fra NR/"
## scriptpath <-"//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/Analyser/Scripts/"
#scriptpath2 <-"//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/AP2/spill/scripts/LusespillTest/"

## Read model data into R
paramSamples <- readRDS("input/paramSamples2020v1.RDS")

## Convert mortality of R and CO to logit scale
logit <- function(p) log(p/(1-p))
paramSamples["lambda0.RCOnat",] <- logit(paramSamples["lambda0.RCOnat",])

## Read environmental data into R
# (these data are real examples, generated in the script FindTempAndExternalLicePressure.R)
load("input/EnvList.Rdata")

## Read typical production data into R
# (these data are averages for three regions based on reported data, 

# calculated in the script FindTempAndExternalLicePressure.R)
load("input/ProdList.Rdata")

# The annual mortality data from
# Norwegian salmon farming from the years 2014 to 2017 ranged between
# 14.2 and 16.2 % (Hjeltnes et al., 2019), corresponding to an average
# monthly mortality rate around 1.3 %.
# Hjeltnes, B., Bang Jensen, B., Bornø, G., Haukaas, A., Walde, C.S., 2019. 
# National Veterinary Institute Fish Health Report 2018. Report 6b/2019. 

## Read model functions
source("ModelFunctions_v2.R")

## Set model settings
Region = c("PO 1-4", "PO 5-7", "PO 8-13")[1]
#Ncages = 4
Ndays = 550 # Maximum 600
POstart = c("Vaar", "Hoest")[1]
ran.ef = c("postsamp","rangen")[1]

## Set function to draw treatment effects
if(ran.ef == "postsamp") trt.sample <- trt.sample.postsamp
if(ran.ef == "rangen") trt.sample <- trt.sample.rangen

## Extract estimated coefficients
Coef.est <- extractCoef(paramSamples, modelVersion="2020v1")

## Set fixed parameters
Coef.fix <- set.Coef.fix()

## Combine estimated and fixed parameters
Coef <- c(Coef.est, Coef.fix)

## Extract posterior samples of random effects
if(ran.ef == "postsamp"){  
  
  # Extract posterior random effects
  RE.all <- extractRanef(paramSamples)
}

set.seed(123456)

## Draw among the posterior samples of the random effects
if(ran.ef == "postsamp"){  
  # Draw farm and cages within farm
  Farm.f <- draw.farm(RE.all, Ndays)
  Cages.f <- draw.cage(RE.all, Farm.f, Ncages)
  
  # Randomly draw one posterior sample ('B')
  B <- draw.b()
  
  # Random effects for posterior sample B of randomly selected cages and farm
  RE <- drawRanef(RE.all, Farm.f, Cages.f, Ndays, B)
}

## Randomly generate random effects
if(ran.ef == "rangen"){ RE <- rangenRanef(Coef, Ncages, Ndays) }

## Find time series for

# Number and weight of salmon
PD <- find.ProdData(ProdList, Region, POstart, Ndays)

# ST, Sea temperature
# N.AF.Ext, Weighted mean adult female lice abundance in nearby cages, and
# A.Ext, Weighted mean adult female lice abundance per salmon in nearby cages
ED <- find.EnvData(EnvList, Region, start.mo = attr(PD,"start.mo"), Ndays)

## Define model variables
MV <- define.model.variables(Ncages, Ndays)

## Combine into one list with state variables
SV <- c(PD, ED, MV)

# Treatment threshold
# thr <- 0.1

## Initialize model
SV2 <- model.initialize(SV, RE, Coef, Ncages, Ndays)
SV <- c(SV, SV2)

# Set days between counting
stepsize <- 7

# Day 1
t <- 1

# Count 20 lice per fish
ncount <- 20
SV$n.SAL[t,] <- ncount

# Update state variables
# lice counts
SV$Y.CH[t,] <- update_Y.CH(t)
SV$Y.AF[t,] <- update_Y.AF(t)
SV$Y.OM[t,] <- update_Y.OM(t)

# Second count of 20 lice per fish
n2count <- 20
SV$n2.SAL[t,] <- n2count

# Update state variables
# lice counts
SV$Y2.CH[t,] <- update_Y2.CH(t)
SV$Y2.AF[t,] <- update_Y2.AF(t)
SV$Y2.OM[t,] <- update_Y2.OM(t)

# true lice abundances
SV$N.R[t, ] <- update_N.R(t)
SV$N.CO[t, ] <- update_N.CO(t)
SV$N.CH[t, , ] <- update_N.CH(t)
SV$N.PA[t, , ] <- update_N.PA(t)
SV$N.AF[t, , ] <- update_N.AF(t)
SV$N.AM[t, , ] <- update_N.AM(t)

# cleaner fish abundances
SV$N.lump[t,] <- update_N.lump(t)
SV$N.wrasse[t,] <- update_N.wrasse(t)   

# cleaner fish-caused mortality of salmon lice
SV$m.clf.PA[t, , ] <- update_m.clf.PA(t)
SV$m.clf.AF[t, , ] <- update_m.clf.AF(t) 

# calculate treatment mortalities in subsequent time-steps 
# if treatments are applied at time t;

# chemical treatments
# (hydrogen peroxide, deltamethrin, azamethiphos, diflubenzuron)
SV[c("ux.HPcht.PA", "ux.HPcht.AF")] <- update_ux.HPcht(t)
SV[c("ux.DMcht.CH", "ux.DMcht.PA", "ux.DMcht.AF")] <- update_ux.DMcht(t)
SV[c("ux.AZcht.PA", "ux.AZcht.AF")] <- update_ux.AZcht(t)
SV[c("ux.EMcht.CH", "ux.EMcht.PA", "ux.EMcht.AF")] <- update_ux.EMcht(t)
SV[c("ux.DBcht.CH", "ux.DBcht.PA")] <- update_ux.DBcht(t)

# physical treatments (thermal, freshwater, mechanical)
SV[c("ux.therm.CH", "ux.therm.PA", "ux.therm.AF")] <- update_ux.therm(t)
SV[c("ux.freshw.CH", "ux.freshw.PA", "ux.freshw.AF")] <- update_ux.freshw(t)
SV[c("ux.mech.CH", "ux.mech.PA", "ux.mech.AF")] <- update_ux.mech(t)

# treatment mortality at time t
SV$m.trt.CH[t, , ] <- update_m.trt.CH(t) 
SV$m.trt.PA[t, , ] <- update_m.trt.PA(t) 
SV$m.trt.AF[t, , ] <- update_m.trt.AF(t) 

## survival rates
SV$s.CH[t, , ] <- update_s.CH(t)
SV$s.PA[t, , ] <- update_s.PA(t)
SV$s.AF[t, , ] <- update_s.AF(t)

## infection rate
SV$d.CO[t, ,] <- update_d.CO(t)

##  reproduction rate for internal recruitment 
SV$r.AF[t, , ] <- update_r.AF(t)

## reproduction rate for external recruitment
SV$r.Ext[t] <- update_r.Ext(t) 

