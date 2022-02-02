############## Code to simulate salmon lice dynamics in a fish farm with multiple cages
############## based on Aldrin et al. (2017): https://doi.org/10.1016/j.ecolmodel.2017.05.019
############## and Aldrin and Huseby (2019) Norsk Regnesentral report SAMBA 28/19.
############## This code is written by Leif Chr. Stige.


############################################################################

####### 1. Code only run once to save Coef ----

## Clear work space
rm(list=ls())

## Paths
parameterpath <- "//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/Analyser/Modellfiler fra NR/"
scriptpath <-"//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/Analyser/Scripts/"

## Read model data into R

paramSamples <- readRDS(paste0(parameterpath,"paramSamples20211209.RDS"))    # Note jan. 2022: New file

## Convert mortality of R and CO to logit scale
logit <- function(p) log(p/(1-p))
paramSamples["lambda0.RCOnat",] <- logit(paramSamples["lambda0.RCOnat",])

## Read model functions
source(paste0(scriptpath,"ModelFunctions_v3.R"))                             # Note jan. 2022: New file

## Extract estimated coefficients
Coef.est <- extractCoef(paramSamples, modelVersion="2020v1") 
# Note jan. 2022: This is still correct as the parameter names are not changed since that version

## Set fixed parameters
Coef.fix <- set.Coef.fix()

## Combine estimated and fixed parameters
Coef <- c(Coef.est, Coef.fix)
# Note jan. 2022: This file can be saved and later loaded instead of paramSamples

## Write Coef to file
saveRDS(Coef, file=paste0(scriptpath, "Coef", "20211209", ".Rds"))


############################################################################

####### 2. Prepare simulation  ----
####### Assuming Coef is saved and random effects are generated (not drawn from posterior)  

## Clear work space
rm(list=ls())

## Paths
scriptpath <-"//vetinst.no\\dfs-felles/StasjonK/FAG/Akva/FoU-Prosjekter/32080 LuseKontroll/Analyser/Scripts/"

## Read Coef from file.
## Note: Coef should be available in the global environment
Coef <- readRDS(file=paste0(scriptpath, "Coef", "20211209", ".RDS"))

## Read environmental data into R
# (these data are real examples, generated in the script FindTempAndExternalLicePressure.R)
load(paste0(scriptpath, "EnvListShort.Rdata"))                               # Note jan. 2022: New file

## Read model functions
source(paste0(scriptpath,"ModelFunctions_v3b.R"))                             # Note jan. 2022: New file
trt.sample <- trt.sample.rangen   # Set function to generate random treatment effects

## Set default model settings
## All model settings are collected in this list.
## Some settings apply to the whole simulation.
## Others are changed during the simulation, for example to apply a lice treatment.
## Some settings can be potentially changed by the user,
## others should be fixed at default values.
default.model.settings <- list(
  
  # Place and time to simulate:
  Region = c("PO 1-4", "PO 5-7", "PO 8-13")[1], # Region or PO. Note jan. 2022: "PO1", "PO2" etc. also works
  start.mo = 5,       # Start month. Note jan. 2022: Replaces POstart = c("Vaar", "Hoest")[1]
                      # Note: Should be between 4 and 10 to have enough environmental data 
                      # from real production cycles.

  Ncages = 4,  # Number of cages
  Ndays = 550, # Number of days in production cycle. Maximum 600

  # Number and weight of salmon:
  w0 = 0.2, # Initial weight (kg) of salmon,
  nstock = 1, # Total number of salmon (millions) stocked in farm (equally divided between cages),
  dstock = 0, # Cage-to-cage delay (days) in sequential stocking of salmon,
  mnat = 0.005/30, # Baseline daily mortality of salmon,
  tslaught = 600, # Time (days) that slaughter starts,
  dslaught = 0, # Cage-to-cage delay (days) in sequential slaughter of salmon.

  # Weekly counting
  stepsize = 7, # days between lice counts
  ncount = 20, # number of salmon counted per cage in first count

  # Lice skirts
  # Note: Has to be decided at start of simulation (for technical reasons).
  do_applyskirt = 0, # apply lice skirt (0: no, 1: yes)? 
  skirtstartday = 1, # from which day of production are skirts applied?
  skirtduration = 180, # how many days will skirts stay?
  skirteffect = 0.5, # proportion of external lice larvae stopped by skirt
  
  # Second lice count
  do_count2 = 0, # should a second count be performed (0: no, 1: yes)?
  n2count = 20, # number of salmon counted in second count, if performed

  # Lice treatment
  do_treat = 0, # apply treatment (0: no, 1: yes)?
  which_treat = "all", # which cages should be treated? ("all" or a vector with cage numbers)
  trt.type = c("HPcht", "DMcht", "AZcht", "EMcht", "DBcht",
                                "therm", "freshw", "mech", 
                                "fx")[6],
  treat.delay = 4, # days from lice count to treatment, if performed
  M.trt = 0.8, # treatment mortality (if treat.type == "fx")

  # Cleaner fish
  do_addclf = 0, # add cleaner fish (0: no, 1: yes)?
  which_clf = "all", # into which cages should cleaner fish be added?
  clfratio = 0.05 # cleaner fish ratio when added

)


########### 3. Some plotting functions ----
logoffset <- .01
ymax.plot <- 20

plot.t1 <- function(...){
  # Set up the plot
  Ndays <- start.model.settings$Ndays
  Ncages <- start.model.settings$Ncages
  yat <- c(0,.1,.5,2,5,20,50,100,1000,1e4,1e5,1e6,1e7,1e8,1e9,1e19)
  par(mfrow=c(2, 1), mai=c(.3,.5,.3,.05), omi=c(.5,.2,0,0))
  plot(0,0, type="n",xlim=c(1,Ndays),ylim=log10(c(logoffset, ymax.plot)), axes=F, xlab="", ylab="")
  box()
  axis(1, lab=F)
  axis(2, at=log10(yat+logoffset), lab=yat)
  mtext(side=2, line=2, "Lus per laks")
  mtext(side=3, line=.1, adj=0.03, "Merd 1")
  lines(log10(logoffset + SV$Lusegrense[1:Ndays]), lty=2, col="red")           # Note jan. 2022: New code
  legend(x = "topright", inset=c(.01,.01), pch=c(1,1), col=c("blue","red"), 
         legend=c("Andre bevegelige","Voksne hunnlus"),cex=.8)   # Note jan. 2022: changed layout a little
  pars <- c('plt','usr')
  par1 <- c(list(mfg=c(1,1)), par(pars))
  
  plot(0,0, type="n",xlim=c(1,Ndays),ylim=log10(c(logoffset, ymax.plot)), axes=F, xlab="", ylab="")
  box()
  axis(1)
  axis(2, at=log10(yat+logoffset), lab=yat)
  mtext(side=1, line=2.5, "Dag")
  mtext(side=2, line=2, "Lus per laks")
  mtext(side=3, line=.1, adj=0.03, "Merd 2")
  lines(log10(logoffset + SV$Lusegrense[1:Ndays]), lty=2, col="red")            # Note jan. 2022: New code
  par2 <- c(list(mfg=c(2,1)), par(pars))
  return(list(par1, par2))
}

plot.later <- function(...){
  Ndays <- start.model.settings$Ndays
  Ncages <- start.model.settings$Ncages
  treat.delay <- new.model.settings$treat.delay
  trt.type <- new.model.settings$trt.type

  par(par1)
  points(t, log10(logoffset + SV$Y.OM[t, 1]/SV$n.SAL[t, 1]), col = "blue")
  points(t, log10(logoffset + SV$Y.AF[t, 1]/SV$n.SAL[t, 1]), col = "red")
  if(Ncages > 1){
    par(par2)
    points(t, log10(logoffset + SV$Y.OM[t, 2]/SV$n.SAL[t, 2]), col = "blue")
    points(t, log10(logoffset + SV$Y.AF[t, 2]/SV$n.SAL[t, 2]), col = "red")
  }

  # Plot arrows when treatments are used
  if((t+treat.delay)<=Ndays){                                            # Note jan. 2022: new condition
    if(SV[[paste0("use.",trt.type)]][t+treat.delay, 1]==1){par(par1)  # Note jan. 2022: notice "t+treat.delay" not "t+1"
      arrows(x0=t+1,x1=t+1,y0=log10(ymax.plot),y1=log10(ymax.plot*.6), length=.05)}
    if(Ncages > 1) {
      if(SV[[paste0("use.",trt.type)]][t+treat.delay, 2]==1){par(par2) # Note jan. 2022: notice "t+treat.delay" not "t+1"
        arrows(x0=t+1,x1=t+1,y0=log10(ymax.plot),y1=log10(ymax.plot*.6), length=.05)}
    }
  }
}

########### 4. Start simulation ----

## Create SV and RE start conditions
# start with default model settings
start.model.settings <- default.model.settings
SV_RE_start <- create.SV_RE_start(EnvList_local = EnvList, 
                                  model.settings = start.model.settings)

## Code run at time t = 1

# rename model settings 
new.model.settings <- start.model.settings

# extract SV and RE from SV_RE_start
SV = SV_RE_start$SV 
RE = SV_RE_start$RE 

# set time index
t = 1

# Apply lice skirts?
# If so, this is done by temporarily changing model settings
# new.model.settings$do_applyskirt <- 1
SV_updated <- update.skirt(SV_local = SV, t_local = t, 
                           model.settings = new.model.settings)
# new.model.settings$do_applyskirt <- 0
SV <- SV_updated

# Add cleaner fish from start of simulation?
# If so, this is done by temporarily changing model settings
# new.model.settings$do_addclf <- 1
SV_updated <- update.licecontrol(SV_local = SV, t_local = t, 
                                 model.settings = new.model.settings)
# new.model.settings$do_addclf <- 0
SV <- SV_updated


# Set up the plot:
par.list <- plot.t1()
par1 <- par.list[[1]]
par2 <- par.list[[2]]

# Code run for each week:
while(t <= new.model.settings$Ndays){
  # Count lice
  SV_updated <- update.licecount1(SV_local = SV, RE_local = RE, t_local = t, 
                                  model.settings = new.model.settings)
  
  SV <- SV_updated
  
  # A possible second count of lice
  # If so, this is done by temporarily changing model settings
  # new.model.settings$do_count2 <- 1
  SV_updated <- update.licecount2(SV_local = SV, RE_local = RE, t_local = t, 
                                  model.settings = new.model.settings)
  # new.model.settings$do_count2 <- 0
  SV <- SV_updated
  
  # Apply lice control measures
  # Here, this is decided "automatically"; in app, the user will decide
  # If treatment is to be applied, we temporarily change model settings
  thr.t <- SV$Lusegrense[t]
  if((t + new.model.settings$treat.delay) <= new.model.settings$Ndays & 
     sum(SV$Y.AF[t, ])/sum(SV$n.SAL[t, ]) > thr.t){
    new.model.settings$do_treat <- 1
  }
  SV_updated <- update.licecontrol(SV_local = SV, t_local = t, 
                                   model.settings = new.model.settings)
  
  SV <- SV_updated
  new.model.settings$do_treat <- 0 # Afterwards, we reset to do_treat = 0
  
  # Plot the updated state:
  plot.later()

  # Update dynamics until next count:
  SV_t_updated <- update.licedynamics(SV_local = SV, RE_local = RE, t_local = t, 
                                      model.settings = new.model.settings)

  # Extract SV and time from SV_t_updated: 
  SV <- SV_t_updated$SV
  t <- SV_t_updated$t_next
  
}


#### 5. Plot other information ----

# Plot the "true" lice levels
par(par1)
lines(log10(logoffset + c(apply(SV$N.AM[, ,1], 1, sum) + apply(SV$N.PA[, ,1], 1, sum))/SV$N.SAL[,1]), col = "blue")
lines(log10(logoffset + c(apply(SV$N.AF[, ,1], 1, sum)/SV$N.SAL[,1])), col = "red")
par(par2)
lines(log10(logoffset + c(apply(SV$N.AM[, ,2], 1, sum) + apply(SV$N.PA[, ,2], 1, sum))/SV$N.SAL[,2]), col = "blue")
lines(log10(logoffset + c(apply(SV$N.AF[, ,2], 1, sum)/SV$N.SAL[,2])), col = "red")

# Add reported lice counts and treatments for same location (not sure all zeros are true zeros)
points(log10(logoffset + SV$AF.Obs), pch=3, col="green")
arrows(x0=which(SV$Treat.Obs==1), x1=which(SV$Treat.Obs==1),
      y0=log10(ymax.plot),y1=log10(ymax.plot*.6), length=.05, col="green")



# Plot production data
Ncages <- start.model.settings$Ncages
par(mfrow=c(2, 1), mai=c(.3,.5,.3,.05), omi=c(.5,.2,0,0))
plot(apply(SV$N.SAL, 1, sum), type="n", xlab="", ylab = "", ylim=c(0, max(SV$N.SAL, na.rm=T)))
for(i in 1:Ncages) lines(1:(dim(SV$N.SAL)[1]), SV$N.SAL[,i])
mtext(side=2, line=2, "Antall (mill.)")
mtext(side=3, line=.5, adj=0.5, "Laks", font=2)
plot(apply(SV$W.SAL, 1, mean, na.rm=T), type="n", xlab="", ylab = "")
for(i in 1:Ncages) lines(1:(dim(SV$W.SAL)[1]), SV$W.SAL[,i])
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
Ndays <- start.model.settings$Ndays
yat <- c(0,.1,.5,2,5,20,50,100,1000,1e4,1e5,1e6,1e7,1e8,1e9,1e19)
par(mfrow=c(1, 1), mai=c(.3,.5,.3,.05), omi=c(.5,.2,0,0))
plot(0,0, type="n",xlim=c(1,Ndays),ylim=log10(c(logoffset, ymax.plot)), axes=F, xlab="", ylab="")
box()
axis(1)
axis(2, at=log10(yat+logoffset), lab=yat)
mtext(side=1, line=2.5, "Dag")
mtext(side=2, line=2, "Lus per laks")
for(i in 1:Ncages) lines(log10(logoffset + c(apply(SV$N.AF[, ,i], 1, sum)/SV$N.SAL[,i])), col = "red")

