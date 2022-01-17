############## Code to simulate salmon lice dynamics in a fish farm with multiple cages
############## based on Aldrin et al. (2017): https://doi.org/10.1016/j.ecolmodel.2017.05.019
############## and Aldrin and Huseby (2019) Norsk Regnesentral report SAMBA 28/19.
############## This code is written by Leif Chr. Stige.


############################################################################

## Clear work space
rm(list=ls())

## Paths
parameterpath <- "//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/Analyser/Modellfiler fra NR/"
scriptpath <-"//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/Analyser/Scripts/"

## Read model data into R
paramSamples <- readRDS(paste0(parameterpath,"paramSamples2020v1.RDS"))

## Convert mortality of R and CO to logit scale
logit <- function(p) log(p/(1-p))
paramSamples["lambda0.RCOnat",] <- logit(paramSamples["lambda0.RCOnat",])

## Read environmental data into R
# (these data are real examples, generated in the script FindTempAndExternalLicePressure.R)
load(paste0(scriptpath, "EnvList.Rdata"))

## Read typical production data into R
# (these data are averages for three regions based on reported data, 
# calculated in the script FindTempAndExternalLicePressure.R)
load(paste0(scriptpath, "ProdList.Rdata"))

  # The annual mortality data from
  # Norwegian salmon farming from the years 2014–2017 ranged between
  # 14.2 and 16.2 % (Hjeltnes et al., 2019), corresponding to an average
  # monthly mortality rate around 1.3 %.
  # Hjeltnes, B., Bang Jensen, B., Bornø, G., Haukaas, A., Walde, C.S., 2019. 
  # National Veterinary Institute Fish Health Report 2018. Report 6b/2019. 

## Read model functions
source(paste0(scriptpath,"ModelFunctions_v2.R"))

## Set model settings
Region = c("PO 1-4", "PO 5-7", "PO 8-13")[1]
Ncages = 2
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
thr <- 100#0.1

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



# Set up the plot
ymax.plot <- 20
logoffset <- .01
yat <- c(0,.1,.5,2,5,20,50,100,1000,1e4,1e5,1e6,1e7,1e8,1e9,1e19)
pars <- c('plt','usr')
par(mfrow=c(2, 1), mai=c(.3,.5,.3,.05), omi=c(.5,.2,0,0))
plot(0,0, type="n",xlim=c(1,Ndays),ylim=log10(c(logoffset, ymax.plot)), axes=F, xlab="", ylab="")
box()
axis(1, lab=F)
axis(2, at=log10(yat+logoffset), lab=yat)
mtext(side=2, line=2, "Lus per laks")
mtext(side=3, line=.1, adj=0.03, "Merd 1")
abline(h=log10(logoffset + thr), lty=2, col="red")
legend(x=Ndays*(-.05), y=log10(ymax.plot)*1.05, pch=c(1,1), col=c("blue","red"), legend=c("Andre bevegelige","Voksne hunnlus"), bty="n")
par1 <- c(list(mfg=c(1,1)), par(pars))

plot(0,0, type="n",xlim=c(1,Ndays),ylim=log10(c(logoffset, ymax.plot)), axes=F, xlab="", ylab="")
box()
axis(1)
axis(2, at=log10(yat+logoffset), lab=yat)
mtext(side=1, line=2.5, "Dag")
mtext(side=2, line=2, "Lus per laks")
mtext(side=3, line=.1, adj=0.03, "Merd 2")
abline(h=log10(logoffset + thr), lty=2, col="red")
par2 <- c(list(mfg=c(2,1)), par(pars))

# Plot the counts for the first day
par(par1)
points(t, log10(logoffset + SV$Y.OM[t, 1]/SV$n.SAL[t, 1]), col = "blue")
points(t, log10(logoffset + SV$Y.AF[t, 1]/SV$n.SAL[t, 1]), col = "red")
if(Ncages > 1){
  par(par2)
  points(t, log10(logoffset + SV$Y.OM[t, 2]/SV$n.SAL[t, 2]), col = "blue")
  points(t, log10(logoffset + SV$Y.AF[t, 2]/SV$n.SAL[t, 2]), col = "red")
}

## Iterative parts of the model:
Nsteps <- floor(Ndays/stepsize)
for(stp in 1:Nsteps){ 
  # Run the daily lice dynamics until next counting
  for(i in 1:stepsize){
    t <- t + 1
  #  print(paste(i, "days since last count", t, "days total"))
    # Update state variables
    
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
  }

  # Count lice
  SV$n.SAL[t,] <- ncount

  # lice counts
  SV$Y.CH[t,] <- update_Y.CH(t)
  SV$Y.AF[t,] <- update_Y.AF(t)
  SV$Y.OM[t,] <- update_Y.OM(t)
  
  # A possible second count of lice
  SV$n2.SAL[t,] <- n2count
  
  # lice counts
  SV$Y2.CH[t,] <- update_Y2.CH(t)
  SV$Y2.AF[t,] <- update_Y2.AF(t)
  SV$Y2.OM[t,] <- update_Y2.OM(t)
  
  # Plot
  par(par1)
  points(t, log10(logoffset + SV$Y.OM[t, 1]/SV$n.SAL[t, 1]), col = "blue")
  points(t, log10(logoffset + SV$Y.AF[t, 1]/SV$n.SAL[t, 1]), col = "red")
  if(Ncages > 1){
    par(par2)
    points(t, log10(logoffset + SV$Y.OM[t, 2]/SV$n.SAL[t, 2]), col = "blue")
    points(t, log10(logoffset + SV$Y.AF[t, 2]/SV$n.SAL[t, 2]), col = "red")
  }

  # Decide whether to treat next day
  trt.type <- c("HPcht", "DMcht", "AZcht", "EMcht", "DBcht",
                "therm", "freshw", "mech", 
                "fx")[6]  # Note: "fx" is not implemented in the rest of this script
  strategy <- c("cage-wise","farm-wise")[2]
  natonull <- function(x){x[is.na(x)]<-0; return(x)}
  # Cage-wise treatment:
  if(t<Ndays & strategy == "cage-wise") SV[[paste0("use.",trt.type)]][t+1, ] <- 
    natonull(1*((SV$Y.AF[t, ] + SV$Y.OM[t, ])/SV$n.SAL[t, ] > thr))
  # Farm-wise treatment:
  if(t<Ndays & strategy == "farm-wise") SV[[paste0("use.",trt.type)]][t+1, ] <- 
    natonull(1*(sum(SV$Y.AF[t, ] + SV$Y.OM[t, ])/sum(SV$n.SAL[t, ]) > thr))

  # Plot arrows when treatments are used
  if(SV[[paste0("use.",trt.type)]][t+1, 1]==1){par(par1)
    arrows(x0=t+1,x1=t+1,y0=log10(ymax.plot),y1=log10(ymax.plot*.6), length=.05)}
  if(Ncages > 1) {
    if(SV[[paste0("use.",trt.type)]][t+1, 2]==1){par(par2)
      arrows(x0=t+1,x1=t+1,y0=log10(ymax.plot),y1=log10(ymax.plot*.6), length=.05)}
  }
}

# Plot the "true" lice levels
par(par1)
lines(log10(logoffset + c(apply(SV$N.AM[, ,1], 1, sum) + apply(SV$N.PA[, ,1], 1, sum))/SV$N.SAL[,1]), col = "blue")
lines(log10(logoffset + c(apply(SV$N.AF[, ,1], 1, sum)/SV$N.SAL[,1])), col = "red")
if(Ncages > 1){
  par(par2)
  lines(log10(logoffset + c(apply(SV$N.AM[, ,2], 1, sum) + apply(SV$N.PA[, ,2], 1, sum))/SV$N.SAL[,2]), col = "blue")
  lines(log10(logoffset + c(apply(SV$N.AF[, ,2], 1, sum)/SV$N.SAL[,2])), col = "red")
}

# Treatment effect
# treat1 <- which(use.HPcht[,1]==1); treat1 <- treat1[treat1 < Ndays-7]
# bef1 <- treat1; aft1 <- rep(treat1, 7) + 1:7
# 1 - sum(N.PA[bef1, ,1])/sum(N.PA[aft1, ,1])
# 1 - sum(N.AM[bef1, ,1])/sum(N.AM[aft1, ,1])
# 1 - sum(N.AF[bef1, ,1])/sum(N.AF[aft1, ,1])

# Plot production data
par(mfrow=c(2, 1), mai=c(.3,.5,.3,.05), omi=c(.5,.2,0,0))
plot(apply(SV$N.SAL, 1, sum), type="l", xlab="", ylab = "")
mtext(side=2, line=2, "Antall (mill.)")
mtext(side=3, line=.5, adj=0.5, "Laks", font=2)
plot(apply(SV$W.SAL, 1, sum), type="l", xlab="", ylab = "")
mtext(side=1, line=2.5, "Dag")
mtext(side=2, line=2, "Vekt (kg)")

# Plot environmental data
par(mfrow=c(3, 1), mai=c(.3,.5,.3,.05), omi=c(.5,.2,0,0))
plot(SV$ST, type="l", xlab="", ylab = "")
mtext(side=2, line=2, "Temperatur")
mtext(side=3, line=.5, adj=0.5, "Miljødata", font=2)
plot(SV$SS, type="l", xlab="", ylab = "")
mtext(side=2, line=2, "Salinitet")
plot(SV$N.AF.Ext, type="l", xlab="", ylab = "")
mtext(side=2, line=2, "Eksternt smittepress")
mtext(side=1, line=2.5, "Dag")

# Plot all cages in one panel
par(mfrow=c(1, 1), mai=c(.3,.5,.3,.05), omi=c(.5,.2,0,0))
plot(0,0, type="n",xlim=c(1,Ndays),ylim=log10(c(logoffset, ymax.plot)), axes=F, xlab="", ylab="")
box()
axis(1)
axis(2, at=log10(yat+logoffset), lab=yat)
mtext(side=1, line=2.5, "Dag")
mtext(side=2, line=2, "Lus per laks")
for(i in 1:Ncages) lines(log10(logoffset + c(apply(SV$N.AF[, ,i], 1, sum)/SV$N.SAL[,i])), col = "red")
