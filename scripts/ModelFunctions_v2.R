############## Functions to simulate salmon lice dynamics in a fish farm with multiple cages
############## based on Aldrin et al. (2017): https://doi.org/10.1016/j.ecolmodel.2017.05.019
############## and Aldrin and Huseby (2019) Norsk Regnesentral report SAMBA 28/19.
############## This code is written by Leif Chr. Stige.


# Inverse logit function
invlogit <- function(x) exp(x)/(1 + exp(x))

# Logit function
logit <- function(p) log(p/(1-p))


# Function to extract coefficients (either mean of posterior samples or for sample no. 'sample')
extractCoef <- function(x,
                        modelVersion="2020v1",
                        ndigits = 10,
                        posterior.sample = "mean"){
  ## Name of selected parameters
  if (modelVersion == "2017v1"){
    outName <-
      c("lambda.RCOnat",
        "lambda.clf",
        "delta.Em10",
        "delta.Nm10",
        "delta.Rs",
        "delta.Rp",
        "delta.CHm10",
        "delta.CHs",
        "delta.CHp",
        "delta.PAm10",
        "delta.PAs",
        "delta.PAp",
        "delta0.CO",
        "delta1.CO",
        "sigma.COdf.2",
        "sigma.COd.2",
        "gamma.r",
        "kappa.clf",
        "rho.CH",
        "rho.OM",
        "rho.AF",
        "beta0.CHcount",
        "sigma.CHCount.2",
        "beta1.CHcount",
        "lambda0.CHnat",
        "phi.CHnat",
        "sigma.CHnat.2",
        "lambda0.PAnat",
        "phi.PAnat",
        "sigma.PAnat.2",
        "lambda0.Anat",
        "phi.Anat",
        "sigma.Anat.2",
        "lambda.DMcht.CH",
        "lambda.DMcht",
        "sigma.DMcht.2",
        "lambda.AZcht",
        "lambda.HPcht",
        "lambda.EMcht",
        "sigma.EMcht.2",
        "lambda.DIcht",
        "mu.Ext",
        "phi.Ext",
        "sigma.Extar.2",
        "sigma.Ext.2")
  } else if (modelVersion == "2020v1"){
    outName <-
      c("lambda0.RCOnat",
        "lambda1.RCOsal",
        "theta.sal",
        "lambda.lump",
        "gamma.lump",
        "lambda.wrasse",
        "gamma.wrasse",
        "delta.Em10",
        "delta.Nm10",
        "delta.Rs",
        "delta.Rp",
        "delta.CHm10",
        "delta.CHs",
        "delta.CHp",
        "delta.PAm10",
        "delta.PAs",
        "delta.PAp",
        "delta0.CO",
        "delta1.CO",
        "delta2.CO",
        "delta3.CO",
        "sigma.COdf.2",
        "sigma.COd.2",
        "gamma.r",
        "kappa.lump",
        "kappa.wrasse",
        "kappa.therm",
        "kappa.mech",
        "kappa.freshw",
        "rho.CH",
        "rho.OM",
        "rho.AF",
        "beta0.CHcount",
        "sigma.CHCount.2",
        "beta1.CHcount",
        "mu0.CHnat",
        "lambda1.CHsal",
        "phi.CHnat",
        "sigma.CHnat.2",
        "mu0.PAnat",
        "lambda1.PAsal",
        "phi.PAnat",
        "sigma.PAnat.2",
        "mu0.Anat",
        "lambda1.Asal",
        "phi.Anat",
        "sigma.Anat.2",
        "lambda.DMcht",
        "sigma.DMcht.2",
        "lambda.AZcht",
        "lambda.HPcht",
        "lambda.EMcht",
        "sigma.EMcht.2",
        "lambda.DBcht",
        "lambda.therm",
        "sigma.therm.2",
        "lambda.mech",
        "lambda.freshw",
        "mu.Ext",
        "phi.Ext",
        "sigma.Extar.2",
        "sigma.Ext.2")
  } else {
    stop(paste("Model", modelVersion, "is not available.\n"))
  }
  rr <- match(outName, row.names(x))
  x <- x[rr,]
  
  ## If sample = "mean", extract mean across posterior samples
  if(posterior.sample == "mean") out <- round(rowMeans(x), ndigits)
  
  ## If not, extract values for posterior sample no. 'sample' 
  if(posterior.sample == "draw") out <- round(x[, B], ndigits)
  
  out
}

# Function to set fixed parameters
set.Coef.fix <- function(...){
  Coef.fix <- c(
    delta.del.HPcht = 0,  # Delay (days) for Hydrogen peroxide effect (Table 1, Aldrin et al. 2017)
    delta.dur.HPcht = 7,  # Duration constant (degree days or days)
    delta.del.DMcht = 2,  # Deltamethrin
    delta.dur.DMcht = 84,
    delta.del.AZcht = 1,  # Azamethiphos
    delta.dur.AZcht = 42, 
    delta.del.EMcht = 5,  # Emamectin
    delta.dur.EMcht = 210,
    delta.del.DBcht = 10, # Diflubenzuron 
    delta.dur.DBcht = 126,
    beta0.R = 172.5,      # Reproduction parameter (number of viable eggs for first extrusion)
    beta1.R = 0.2         # Reproduction parameter (how the number of viable eggs per extrusion increases with stage-age)
  )
  return(Coef.fix)
}

# Function to extract posterior samples of random effects from model output
extractRanef <- function(paramSamples, ...){ 
  ## Random variation between farms in probability of counting chalimus-stage lice
  beta0.CHcount.all <- paramSamples[grep("beta0.CHcount.", row.names(paramSamples)),]
  
  ## Random variation in external infestation pressure
  z.Ext.all <- paramSamples[grep("z.Ext", row.names(paramSamples)),]
  
  ## On arithmetic scale:
  e.Ext.all <- exp(z.Ext.all)
  
  ## Random variation in natural mortality rates of stages CH, PA and A
  z.CHnat.all <- paramSamples[grep("z.CHnat", row.names(paramSamples)),] 
  z.PAnat.all <- paramSamples[grep("z.PAnat", row.names(paramSamples)),] 
  z.Anat.all  <- paramSamples[grep("z.Anat", row.names(paramSamples)),] 
  
  ## Mortality transformed to scale 0,1 
  # (Eq. 15, Aldrin et al. 2017)
  # m.CHnat.all <- invlogit(z.CHnat.all)
  # m.PAnat.all <- invlogit(z.PAnat.all)
  # m.Anat.all <- invlogit(z.Anat.all)
  
  ## Random variation in infection rate
  delta0.CO.all <- paramSamples[grep("delta0.CO.", row.names(paramSamples)),]
  
  ## Random variation in effects of chemical treatments
  uStar.DMcht.all <- paramSamples[grep("uStar.DMcht.", row.names(paramSamples)),]
  uStar.AZcht.all <- paramSamples[grep("uStar.AZcht.", row.names(paramSamples)),]
  uStar.HPcht.all <- paramSamples[grep("uStar.HPcht.", row.names(paramSamples)),]
  uStar.EMcht.all <- paramSamples[grep("uStar.EMcht.", row.names(paramSamples)),]
  uStar.DBcht.all <- paramSamples[grep("uStar.DBcht.", row.names(paramSamples)),] # New

  ## Random variation in effects of physical treatments
  uStar.therm.all <- paramSamples[grep("uStar.therm.", row.names(paramSamples)),] # thermic
  uStar.freshw.all <- paramSamples[grep("uStar.freshw.", row.names(paramSamples)),] # freshwater
  uStar.mech.all <- paramSamples[grep("uStar.mech.", row.names(paramSamples)),]   # mechanical  
  
  # Sampling units:
  Samples.list <- function(uStar.select){
    farm.appl <- unlist(lapply(row.names(uStar.select),function(x){paste(unlist(strsplit(x,"[.]")[[1]][c(6,5)]), collapse=".")}))
    cage <- unlist(lapply(row.names(uStar.select),function(x){strsplit(x,"[.]")[[1]][3]}))
    stage <- unlist(lapply(row.names(uStar.select),function(x){strsplit(x,"[.]")[[1]][4]}))
    out <- list()
    length(out) <- length(unique(farm.appl))
    for(i in 1:length(out)){
      farm.appl.i <- unique(farm.appl)[i]
      out[[i]] <- list()
      length(out[[i]]) <- length(unique(cage[farm.appl==farm.appl.i]))
      for(j in 1:length(out[[i]])){
        cage.j <- unique(cage[farm.appl==farm.appl.i])[j]
        out[[i]][[j]] <- which(farm.appl==farm.appl.i & cage==cage.j)
        names(out[[i]][[j]]) <- stage[farm.appl==farm.appl.i & cage==cage.j]
      }
    }
    out
  }
  uStar.DMcht.samples <- Samples.list(uStar.DMcht.all)
  uStar.AZcht.samples <- Samples.list(uStar.AZcht.all)
  uStar.HPcht.samples <- Samples.list(uStar.HPcht.all)
  uStar.EMcht.samples <- Samples.list(uStar.EMcht.all)
  uStar.DBcht.samples <- Samples.list(uStar.DBcht.all)
  uStar.therm.samples <- Samples.list(uStar.therm.all)
  uStar.freshw.samples <- Samples.list(uStar.freshw.all)
  uStar.mech.samples <- Samples.list(uStar.mech.all)
  
  ## Number of days modelled for each farm
  
  Farmnames <- unlist(lapply(row.names(beta0.CHcount.all), function(x) strsplit(x,"beta0.CHcount.")[[1]][2]))
  FarmMaxTimes <- rep(NA, length(Farmnames))
  names(FarmMaxTimes) <- Farmnames
  for(i in 1:length(Farmnames)){
    FarmMaxTimes[i] <- length(grep(Farmnames[i],row.names(z.Ext.all)))
  }
  
  ## Number of cages per farm
  FarmNcages <- rep(NA, length(Farmnames))
  names(FarmNcages) <- Farmnames
  for(i in 1:length(Farmnames)){
    FarmNcages[i] <- length(grep(Farmnames[i],row.names(delta0.CO.all))) - 1 
  }
  
return(
  list(
    beta0.CHcount.all =  beta0.CHcount.all,
    z.Ext.all = z.Ext.all,
    e.Ext.all = e.Ext.all,
    z.CHnat.all = z.CHnat.all,
    z.PAnat.all = z.PAnat.all,
    z.Anat.all = z.Anat.all,
    delta0.CO.all = delta0.CO.all,
    uStar.DMcht.all = uStar.DMcht.all,
    uStar.AZcht.all = uStar.AZcht.all,
    uStar.HPcht.all = uStar.HPcht.all,
    uStar.EMcht.all = uStar.EMcht.all,
    uStar.DBcht.all = uStar.DBcht.all,
    uStar.therm.all = uStar.therm.all,
    uStar.freshw.all = uStar.freshw.all,
    uStar.mech.all = uStar.mech.all,

    uStar.DMcht.samples = uStar.DMcht.samples,
    uStar.AZcht.samples = uStar.AZcht.samples,
    uStar.HPcht.samples = uStar.HPcht.samples,
    uStar.EMcht.samples = uStar.EMcht.samples,
    uStar.DBcht.samples = uStar.DBcht.samples,
    uStar.therm.samples = uStar.therm.samples,
    uStar.freshw.samples = uStar.freshw.samples,
    uStar.mech.samples = uStar.mech.samples,
    
    Farmnames = Farmnames,
    FarmMaxTimes = FarmMaxTimes,
    FarmNcages = FarmNcages
  )
)  
}

# Function to randomly draw a farm with sufficient number of days
draw.farm <- function(RE.all, Ndays, ....){
  if(Ndays > median(RE.all$FarmMaxTimes)) stop("Too long farm duration")
  Dummy <- 0
  while(Dummy==0){
    Farm.f <- sample(RE.all$Farmnames, 1)
    if(RE.all$FarmMaxTimes[Farm.f]>=Ndays){Dummy <- 1}
  }
  return(Farm.f)
}

# Function to randomly select cages within farms
draw.cage <- function(RE.all, Farm.f, Ncages,...){
  Cages.f <- sample(c(1:RE.all$FarmNcages[Farm.f]), size=Ncages, replace = T)
  return(Cages.f)
}

# Function to randomly select one posterior sample
draw.b <- function(Bmax = 400){sample(1:Bmax, 1)}

# Function to draw random effects for posterior sample B of randomly selected cages and farm
drawRanef <- function(RE.all, Farm.f, Cages.f, Ndays, B,...){
  ## Random variation between farms in probability of counting chalimus-stage lice
  beta0.CHcount.f <- RE.all$beta0.CHcount.all[paste("beta0.CHcount", Farm.f, sep="."), B]
  
  ## Random variation in external recruitment for farm f
  # This is drawn from a single posterior sample to retain the time-series correlation
  e.Ext.f <- RE.all$e.Ext.all[paste("z","Ext",Farm.f,1:Ndays, sep="."), B]
  
  ## Random variation in natural mortality for farm f
  # These are drawn from single posterior samples to retain the time-series correlation
  # m.CHnat.f <- m.CHnat.all[paste("z","CHnat",Farm.f,1:Ndays, sep="."),  B]
  # m.PAnat.f <- m.PAnat.all[paste("z","PAnat",Farm.f,1:Ndays, sep="."), B]
  # m.Anat.f <- m.Anat.all[paste("z","Anat",Farm.f,1:Ndays, sep="."), B]
  z.CHnat.f <- RE.all$z.CHnat.all[paste("z","CHnat",Farm.f,1:Ndays, sep="."),  B]
  z.PAnat.f <- RE.all$z.PAnat.all[paste("z","PAnat",Farm.f,1:Ndays, sep="."), B]
  z.Anat.f <- RE.all$z.Anat.all[paste("z","Anat",Farm.f,1:Ndays, sep="."), B]
  
  ## Random variation in infection rate for farm f
  # This is drawn from the same posterior sample for consistency with the other random effects
  delta0.CO.f <- RE.all$delta0.CO.all[paste("delta0.CO", Cages.f, Farm.f, sep="."), B]
  
  return(
    list(
      beta0.CHcount.f = beta0.CHcount.f,
      e.Ext.f = e.Ext.f,
      z.CHnat.f = z.CHnat.f,
      z.PAnat.f = z.PAnat.f,
      z.Anat.f = z.Anat.f,
      delta0.CO.f = delta0.CO.f
    )
  )
}

# Function to extract random effect of chemical or physical treatment
# for each stage (columns) for Ncages (rows)
trt.sample.postsamp <- function(trt, posterior.sample="mean",...){ # Draw from posterior samples
  
  # The argument 'trt' is the type of chemical or physical treatment

    if(trt=="DMcht") {uStar <- RE.all$uStar.DMcht.all; uStar.samples <- RE.all$uStar.DMcht.samples}
    if(trt=="AZcht") {uStar <- RE.all$uStar.AZcht.all; uStar.samples <- RE.all$uStar.AZcht.samples}
    if(trt=="HPcht") {uStar <- RE.all$uStar.HPcht.all; uStar.samples <- RE.all$uStar.HPcht.samples}
    if(trt=="EMcht") {uStar <- RE.all$uStar.EMcht.all; uStar.samples <- RE.all$uStar.EMcht.samples}
    if(trt=="DBcht") {uStar <- RE.all$uStar.DBcht.all; uStar.samples <- RE.all$uStar.DBcht.samples}
    #  if(trt=="TBcht") {uStar <- RE.all$uStar.TBcht.all; uStar.samples <- RE.all$uStar.TBcht.samples}
    if(trt=="therm") {uStar <- RE.all$uStar.therm.all; uStar.samples <- RE.all$uStar.therm.samples}
    if(trt=="freshw") {uStar <- RE.all$uStar.freshw.all; uStar.samples <- RE.all$uStar.freshw.samples}
    if(trt=="mech") {uStar <- RE.all$uStar.mech.all; uStar.samples <- RE.all$uStar.mech.samples}
    
    # Randomly drawing farm treatment and  cages within treatment
    ran.treat <- sample(sample(uStar.samples,size=1)[[1]], size=Ncages, replace=T)
    
    # Calculating random effects
    # If posterior.sample=="mean", we use posterior mean values for treatment effects
    # If posterior.sample=="draw", we draw sample B from the posterior distribution
    out <- array(dim=c(Ncages,length(ran.treat[[1]])),dimnames=list(1:Ncages, names(ran.treat[[1]])))
    for(i in 1:Ncages){
      if(posterior.sample=="mean") out[i,] <- rowMeans(uStar[ran.treat[[i]], ])
      if(posterior.sample=="draw") out[i,] <- uStar[ran.treat[[i]], B] 
    }

  return(out)
}

# Function to randomly generate random effects
rangenRanef <- function(Coef, Ncages, Ndays, ...){
  
  ## Random variation between farms in probability of counting chalimus-stage lice
  # (Eq. 47, Aldrin et al. 2017)
  beta0.CHcount.f <- rnorm(mean=Coef["beta0.CHcount"], sd=Coef["sigma.CHCount.2"]^.5, n=1)
  
  ## Random variation in external recruitment for farm f
  # (Eqs. 33-35, Aldrin et al. 2017)
  z.Ext.f <- rnorm(mean=Coef["mu.Ext"], sd=Coef["sigma.Ext.2"]^.5, n=1) +     # Farm-specific mean +
    arima.sim(list(ar=Coef["phi.Ext"]), sd=Coef["sigma.Extar.2"]^.5, n=Ndays) # time variation
  
  # On arithmetic scale (Eq. 33, Aldrin et al. 2017:
  e.Ext.f <- exp(z.Ext.f)
  
  ## Random variation in natural mortality for farm f
  
  # (Eqs. 15-17, Aldrin et al. 2017), Eqs. 8-10 Aldrin & Huseby 2019:
  # Note: Salinity is not corrected for here, but in model initialisation
  z.CHnat.f <- log(Coef["mu0.CHnat"]) - Coef["sigma.CHnat.2"]/(2*(1 - Coef["phi.CHnat"]^2)) +
                           arima.sim(list(ar=Coef["phi.CHnat"]), sd=Coef["sigma.CHnat.2"]^.5, n=Ndays)
  z.PAnat.f <- log(Coef["mu0.PAnat"]) - Coef["sigma.PAnat.2"]/(2*(1 - Coef["phi.PAnat"]^2)) +
                           arima.sim(list(ar=Coef["phi.PAnat"]), sd=Coef["sigma.PAnat.2"]^.5, n=Ndays)
  z.Anat.f <- log(Coef["mu0.Anat"]) - Coef["sigma.Anat.2"]/(2*(1 - Coef["phi.Anat"]^2)) +
                          arima.sim(list(ar=Coef["phi.Anat"]), sd=Coef["sigma.Anat.2"]^.5, n=Ndays)

  # Constraints: (0.0006–0.02) for CH, (0.002–0.21) for PA, and(0.0003–0.70) for A.
  constrain <- function(x, min, max){x[x<min] <- min; x[x>max] <- max; x}
  logit <- function(p){log(p/(1-p))}
  z.CHnat.f <- constrain(z.CHnat.f, logit(0.0006), logit(0.02))
  z.PAnat.f <- constrain(z.PAnat.f, logit(0.002), logit(0.21))
  z.Anat.f <- constrain(z.Anat.f, logit(0.0003), logit(0.70))
  
  ## Random variation in infection rate for farm f
  # (Eqs. 29-30, Aldrin et al. 2017)
  delta0.CO.f <- rnorm(n = Ncages, 
                        mean = rnorm(n=1, mean = Coef["delta0.CO"], sd = Coef["sigma.COd.2"]^.5),
                        sd = Coef["sigma.COdf.2"]^.5)
  
  return(
    list(
      beta0.CHcount.f = beta0.CHcount.f,
      e.Ext.f = e.Ext.f,
      z.CHnat.f = z.CHnat.f,
      z.PAnat.f = z.PAnat.f,
      z.Anat.f = z.Anat.f,
      delta0.CO.f = delta0.CO.f
    )
  )
  
}

# Function to randomly generate random effect of chemical or physical treatment
trt.sample.rangen <- function(trt,...){ # Randomly generate cage-specific effects for each stage
    # Note that effect of CM and DM are identical
    # Furthermore, sd of HP, CM, DM an AS are identical, as are sd of EM and DI (Table 5, Aldrin et al. 2017)
    if(trt=="HPcht") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.HPcht"],sd=Coef["sigma.DMcht.2"]^.5),
                                   dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af")))}
    if(trt=="DMcht") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.DMcht"],sd=Coef["sigma.DMcht.2"]^.5),
                                   dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af")))}
    if(trt=="AZcht") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.AZcht"],sd=Coef["sigma.DMcht.2"]^.5),
                                   dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af")))}
    if(trt=="EMcht") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.EMcht"],sd=Coef["sigma.EMcht.2"]^.5),
                                   dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af")))}
    if(trt=="DBcht") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.DBcht"],sd=Coef["sigma.EMcht.2"]^.5),
                                   dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af"))) }
    if(trt=="therm") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.therm"],sd=Coef["sigma.therm.2"]^.5),
                                   dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af"))) }
    if(trt=="freshw") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.freshw"],sd=Coef["sigma.therm.2"]^.5),
                                    dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af"))) }
    if(trt=="mech") {out <- array(rnorm(n=Ncages*3, mean=Coef["lambda.mech"],sd=Coef["sigma.therm.2"]^.5),
                                  dim=c(Ncages,3),dimnames=list(1:Ncages,c("CH","PA","Af"))) }
    
  return(out)
}

# Function to find production data and start month (if not specified)
# Uses Region, POstart, Ncages, Ndays, ProdList; gives N.SAL, W.SAL, start.mo
find.ProdData <- function(ProdList, Region, POstart, Ndays,...){
  
  # Find relevant Region
  if(Region=="PO 1-4" | is.element(Region, paste0("PO", 1:4))) ProdData <- ProdList[[1]]
  if(Region=="PO 5-7" | is.element(Region, paste0("PO", 5:7))) ProdData <- ProdList[[2]]
  if(Region=="PO 8-13"| is.element(Region, paste0("PO", 8:13))) ProdData <- ProdList[[3]]
  
  # Extract monthly n.sal (million individuals)
  if(POstart=="Vaar" | is.element(POstart, 1:6)) n.sal <- ProdData$n.spr / 10^6 
  if(POstart=="Hoest" | is.element(POstart, 7:12)) n.sal <- ProdData$n.aut / 10^6 
  
  # Extract monthly v.sal (kg)
  if(POstart=="Vaar" | is.element(POstart, 1:6)) v.sal <- ProdData$v.spr / 10^3 
  if(POstart=="Hoest" | is.element(POstart, 7:12)) v.sal <- ProdData$v.aut / 10^3 
  
  # Set start month
  if(POstart=="Vaar") start.mo <- ProdData$m0.spr
  if(POstart=="Hoest") start.mo <- ProdData$m0.aut
  if(is.element(POstart, 1:12)) start.mo <- POstart
  
  # Create daily time series
  # Multiplication array for interpolation assuming 30-days months (assume same value for all last month)
  X <- array(0, dim=c(600, 20))
  X[, 1] <- c(c(30:1)/30, rep(0, 30*19))
  for(i in 2:19) X[, i] <- c(rep(0, 30*(i-2)), c(0:30,29:1)/30, rep(0, 30*(20-i)))
  X[, 20] <- c(rep(0,30*18),c(1:30)/30, rep(1, 30))
  
  # The total production per farm is divided on 8, the typical number of cages per farm according to Aldrin et al. (2017)
  n.sal.daily <- c(X %*% n.sal)[1:Ndays] / 8
  v.sal.daily <- c(X %*% v.sal)[1:Ndays]
  
  # Number and weight of salmon
  N.SAL <- array(n.sal.daily, dim = c(Ndays, Ncages)) # Number of Salmon (millions)
  W.SAL <- array(v.sal.daily, dim = c(Ndays, Ncages)) # Mean weight of Salmon (kg)
  
  out <- list(
    N.SAL = N.SAL,
    W.SAL = W.SAL) 
  attr(out, "start.mo") <- start.mo
  
  return(out)
}

# Function to find environmental data (temperature, external lice pressure) 
# Uses Region, EnvList, start.mo, Ndays; gives ST, N.AF.Ext, A.Ext
# Currently sets salinity to 35 psu
find.EnvData <- function(EnvList, Region, start.mo, Ndays, ...){
  # Find environmental data in the right PO(s), month(s) and number of months with data
  POs <- unlist(lapply(EnvList, function(x) x$PO))
  Months <- unlist(lapply(EnvList, function(x) x$Month))
  Nmonths <- unlist(lapply(EnvList, function(x) x$Nmonths))
  if(Region=="PO 1-4") PO.alts <- 1:4
  if(Region=="PO 5-7") PO.alts <- 5:7
  if(Region=="PO 8-13") PO.alts <- 8:13
  if(is.element(Region, paste0("PO", 1:13))) PO.alts <- as.numeric(substr(Region, 3,4))
  env.alts <- which(is.element(POs, PO.alts) & 
                      Months == start.mo &         # This is set in find.ProdData
                      Nmonths >= Ndays/30)
  if(length(env.alts)==0) stop("Not enough environmental data for this area, period and production length")
  
  # Draw a random environmental sample that meets the criteria
  index <- sample(env.alts, size=1)
  EnvData.select <- EnvList[[index]]
  
  ST <- as.numeric(EnvData.select$ST.daily)[1:Ndays]
  N.AF.Ext <- as.numeric(EnvData.select$N.AF.Ext.daily)[1:Ndays]
  A.Ext <- as.numeric(EnvData.select$A.Ext.daily)[1:Ndays]
  
  # Salinity is set at 35 psu
  SS <- ST
  SS[] <- 35
  
  # Stage-average temperatures
  ST.stageav <- EnvData.select$ST.stageav[1:Ndays,]
  
  return(
    list(
      ST = ST,
      N.AF.Ext = N.AF.Ext,
      A.Ext = A.Ext,
      SS = SS,
      ST.stageav = ST.stageav
    )
  )
}


# Function to define model variables
define.model.variables <- function(Ncages, Ndays,...){
  
  ## Set maximum number of days modelled in a stage (before dying)
  a.max.R <- 60
  a.max.CO <- 60
  a.max.CH <- 60
  a.max.PA <- 60
  a.max.AF <- 80
  
  ##################### External influencing variables

  # External recruitment rate
  r.Ext <- array(dim = c(Ndays, 1))

  ##################### State variables
  
  # Number of salmon lice at different stages
  N.R <- array(dim = c(Ndays, a.max.R)) # N recruits (time t_local, stage-age a)
  N.CO <- array(dim = c(Ndays, a.max.CO)) # N copepodids (time t_local, stage-age a)
  N.CH <- array(dim = c(Ndays, a.max.CH, Ncages)) # N chalimus (time t_local, stage-age a, cage c)
  N.PA <- array(dim = c(Ndays, a.max.PA, Ncages)) # N pre-adults (time t_local, stage-age a, cage c)
  N.AF <- array(dim = c(Ndays, a.max.AF, Ncages)) # N adult females (time t_local, stage-age a, cage c)
  N.AM <- N.AF # N adult males (time t_local, stage-age a, cage c)
  
  # Number and weight of salmon, number of cleaner fish and cleaner fish ratio
  # N.SAL <- array(dim = c(Ndays, Ncages)) # Number of Salmon (millions)
  # W.SAL <- array(dim = c(Ndays, Ncages)) # Mean weight of Salmon (kg)
  N.lump <- array(dim = c(Ndays, Ncages)) # Number of lumpsucker cleaner fish
  N.wrasse <- array(dim = c(Ndays, Ncages)) # Number of wrasse cleaner fish
  
  # Number of salmon on which lice were counted in different cages
  n.SAL <- array(0, dim = c(Ndays, Ncages)) # Number 
  
  # Number of salmon on which lice were counted in a potential second count
  n2.SAL <- array(0, dim = c(Ndays, Ncages)) # Number 
  
  # Number of salmon lice counted (total across all fish counted in each cage) 
  Y.CH <- array(0, dim = c(Ndays, Ncages)) # Number of Chalimus
  Y.AF <- array(0, dim = c(Ndays, Ncages)) # Number of adult females
  Y.OM <- array(0, dim = c(Ndays, Ncages)) # Number of other mobiles
  
  # Number of salmon lice counted in a potential second count 
  Y2.CH <- array(0, dim = c(Ndays, Ncages)) # Number of Chalimus
  Y2.AF <- array(0, dim = c(Ndays, Ncages)) # Number of adult females
  Y2.OM <- array(0, dim = c(Ndays, Ncages)) # Number of other mobiles
  
  # Stocking of cleaner fish (lump or wrasse) in each cage
  # This is initially set to zero, but changed in time-steps when supplied
  S.lump <- array(0, dim=c(Ndays, Ncages))
  S.wrasse <- array(0, dim=c(Ndays, Ncages))
  
  ## Initial state at time t_local = 1 is zero lice
  N.R[1,] <- 0
  N.CO[1,] <- 0
  N.CH[1,,] <- 0
  N.PA[1,,] <- 0
  N.AF[1,,] <- 0
  N.AM[1,,] <- 0
  N.lump[1,] <- 0
  N.wrasse[1,] <- 0
  
  ############### Chemical treatment use
  
  # Regression coefficients for chemical treatment effects
  ux.HPcht.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # Hydrogen peroxide effect over time, ages, cages for stage = PA
  ux.HPcht.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # Hydrogen peroxide effect over time, ages, cages for stage = A
  ux.DMcht.CH <- array(0, dim = c(Ndays, a.max.CH, Ncages))  # Deltamethrin
  ux.DMcht.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # Deltamethrin
  ux.DMcht.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # Deltamethrin
  ux.AZcht.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # Azamethiphos
  ux.AZcht.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # Azamethiphos
  ux.EMcht.CH <- array(0, dim = c(Ndays, a.max.CH, Ncages))  # Emamectin
  ux.EMcht.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # Emamectin
  ux.EMcht.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # Emamectin
  ux.DBcht.CH <- array(0, dim = c(Ndays, a.max.CH, Ncages))  # Diflubenzuron 
  ux.DBcht.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # Diflubenzuron 
  
  # Regression coefficients for physical treatment effects
  ux.therm.CH <- array(0, dim = c(Ndays, a.max.CH, Ncages))  # thermal
  ux.therm.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # thermal
  ux.therm.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # thermal
  ux.freshw.CH <- array(0, dim = c(Ndays, a.max.CH, Ncages))  # freshwater
  ux.freshw.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # freshwater
  ux.freshw.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # freshwater
  ux.mech.CH <- array(0, dim = c(Ndays, a.max.CH, Ncages))  # mechanical
  ux.mech.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # mechanical
  ux.mech.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # mechanical

  # Regression coefficients for generic treatment effects with fixed mortality
  ux.fx.CH <- array(0, dim = c(Ndays, a.max.CH, Ncages))  # 
  ux.fx.PA <- array(0, dim = c(Ndays, a.max.PA, Ncages))  # 
  ux.fx.AF <- array(0, dim = c(Ndays, a.max.AF, Ncages))  # 

  # Indicator variables (0, 1) for use of chemical treatments in each cage
  # This is initially set to zero, but set to one in time-steps when applied.
  use.HPcht <- array(0, dim=c(Ndays, Ncages))
  use.DMcht <- array(0, dim=c(Ndays, Ncages))
  use.AZcht <- array(0, dim=c(Ndays, Ncages))
  use.EMcht <- array(0, dim=c(Ndays, Ncages))
  use.DBcht <- array(0, dim=c(Ndays, Ncages))
  
  # Indicator variables (0, 1) for use of physical treatments in each cage
  # This is initially set to zero, but set to one in time-steps when applied.
  use.therm <- array(0, dim=c(Ndays, Ncages))
  use.freshw <- array(0, dim=c(Ndays, Ncages))
  use.mech <- array(0, dim=c(Ndays, Ncages))
  
  # Indicator variables (0, 1) for use of generic treatments in each cage
  # This is initially set to zero, but set to one in time-steps when applied.
  use.fx <- array(0, dim=c(Ndays, Ncages))

  return(
    list(
      r.Ext = r.Ext,
      N.R = N.R,
      N.CO = N.CO,
      N.CH = N.CH,
      N.PA = N.PA,
      N.AF = N.AF,
      N.AM = N.AM,
      N.lump = N.lump,
      N.wrasse = N.wrasse,
      n.SAL = n.SAL,
      n2.SAL = n2.SAL,
      Y.CH = Y.CH,
      Y.AF = Y.AF,
      Y.OM = Y.OM,
      Y2.CH = Y2.CH,
      Y2.AF = Y2.AF,
      Y2.OM = Y2.OM,
      S.lump = S.lump,
      S.wrasse = S.wrasse,
      ux.HPcht.PA = ux.HPcht.PA,
      ux.HPcht.AF = ux.HPcht.AF,
      ux.DMcht.CH = ux.DMcht.CH,
      ux.DMcht.PA = ux.DMcht.PA,
      ux.DMcht.AF = ux.DMcht.AF,
      ux.AZcht.PA = ux.AZcht.PA,
      ux.AZcht.AF = ux.AZcht.AF,
      ux.EMcht.CH = ux.EMcht.CH,
      ux.EMcht.PA = ux.EMcht.PA, 
      ux.EMcht.AF = ux.EMcht.AF,
      ux.DBcht.CH = ux.DBcht.CH, 
      ux.DBcht.PA = ux.DBcht.PA, 
      ux.therm.CH = ux.therm.CH,
      ux.therm.PA = ux.therm.PA,
      ux.therm.AF = ux.therm.AF,
      ux.freshw.CH = ux.freshw.CH,
      ux.freshw.PA = ux.freshw.PA,
      ux.freshw.AF = ux.freshw.AF,
      ux.mech.CH = ux.mech.CH,
      ux.mech.PA = ux.mech.PA,
      ux.mech.AF = ux.mech.AF,
      ux.fx.CH = ux.fx.CH,
      ux.fx.PA = ux.fx.PA,
      ux.fx.AF = ux.fx.AF,
      use.HPcht = use.HPcht,
      use.DMcht = use.DMcht,
      use.AZcht = use.AZcht,
      use.EMcht = use.EMcht,
      use.DBcht = use.DBcht,
      use.therm = use.therm,
      use.freshw = use.freshw,
      use.mech = use.mech,
      use.fx = use.fx
    )
  )
}


# Development rate function (Eq. 24, Aldrin et al. 2017)
dev.rate <- function(d.m, delta.s){ # d.m is array with median development time parameters
  A <- d.m # Array with stage-ages 0, 1, 2, ...
  for(i in 1:ncol(A)) A[,i] <- i-1
  d.raw <- log(2)*d.m^(-delta.s) * delta.s * A^(delta.s - 1) # developmental rate
  d <- pmin(d.raw, 1) # constrain to max. 1
  return(d)
}

# Median development time parameter for different stage-ages (Eq. 25, Aldrin et al. 2017):
med.time <- function(st.av, delta.m10, delta.p){delta.m10*(10/st.av)^delta.p}

# Function to initialize model 
# (define vital rates and fill in values that do not have to be filled in iteratively)
model.initialize <- function(SV, RE, Coef, Ncages, Ndays, ...){
  
  ## Set maximum number of days modelled in a stage (before dying)
  a.max.R <- 60
  a.max.CO <- 60
  a.max.CH <- 60
  a.max.PA <- 60
  a.max.AF <- 80

  #################### Vital rates
  
  ## Transformed salinity (Eq. 7, Aldrin & Huseby 2019)
  x.sal <- pmin(SV$SS, Coef["theta.sal"]) - Coef["theta.sal"]
  
  ## Natural mortality rates of salmon lice
  m.nat.R <- array(dim=c(Ndays, a.max.R))    
  m.nat.CO <- array(dim=c(Ndays, a.max.CO))  
  m.nat.CH <- array(dim=c(Ndays, a.max.CO, Ncages))
  m.nat.PA <- array(dim=c(Ndays, a.max.PA, Ncages))
  m.nat.AF <- array(dim=c(Ndays, a.max.AF, Ncages))
  
  # Fill values into arrays
  # This code requires that time is the first dimension in the arrays

  # Mortality depends here on salinity
  # Eq. 8-10 (Aldrin & Huseby 2019). 
  m.nat.R[,] <- invlogit(Coef["lambda0.RCOnat"] + Coef["lambda1.RCOsal"] * x.sal)  # Eq. 4-6, Aldrin & Huseby 2019
  m.nat.CO[,] <- invlogit(Coef["lambda0.RCOnat"])  # Eq. 4-6, Aldrin & Huseby 2019. NB! No salinity correction here.
  m.nat.CH[,,] <- invlogit(RE$z.CHnat.f + Coef["lambda1.CHsal"] * x.sal)
  m.nat.PA[,,] <- invlogit(RE$z.PAnat.f + Coef["lambda1.PAsal"] * x.sal)
  m.nat.AF[,,] <- invlogit(RE$z.Anat.f + Coef["lambda1.Asal"] * x.sal)
  # Note, in the posterior samples, salinity is already corrected for in the random effects, z.CHnat, z.PAnat, z.Anat
  # The correction here is therefore potentially double
  # With SS set to 35, this means using a random sample of observed salinity effects (but not for R and CO)
  # When randomly generating random effects, this is not an issue
  
  ## Cleaner fish-caused mortality of salmon lice
  m.clf.PA <- array(dim=c(Ndays, a.max.CH, Ncages))
  m.clf.AF <- array(dim=c(Ndays, a.max.AF, Ncages))
  
  ## Chemical or physical treatment-caused mortality of salmon lice
  m.trt.CH <- array(dim=c(Ndays, a.max.CH, Ncages))
  m.trt.PA <- array(dim=c(Ndays, a.max.PA, Ncages))
  m.trt.AF <- array(dim=c(Ndays, a.max.AF, Ncages))
  
  ## Survival rates
  s.R <- (1 - m.nat.R)     # Eq. 12, Aldrin et al. 2017
  s.CO <- (1 - m.nat.CO)   # Eq. 12, Aldrin et al. 2017
  s.CH <- array(dim=c(Ndays, a.max.CO, Ncages))
  s.PA <- array(dim=c(Ndays, a.max.PA, Ncages))
  s.AF <- array(dim=c(Ndays, a.max.AF, Ncages))
  
  ## Developmental rates
  #d.R is calculated below
  d.CO <- array(dim=c(Ndays, a.max.CO, Ncages))
  d.CH <- array(dim=c(Ndays, a.max.CH, Ncages))
  d.PA <- array(dim=c(Ndays, a.max.PA, Ncages))
  
  ## Reproduction rates
  r.AF <- array(dim=c(Ndays, a.max.AF, Ncages))
  
  # Egg and nauplii development time parameter for different time steps and stage-ages (using Eq. 26, Aldrin et al. 2017)
  d.m.R <- med.time(SV$ST.stageav[,1:a.max.R], Coef["delta.Em10"] + Coef["delta.Nm10"], Coef["delta.Rp"])
  
  # Chalimus development time parameter for different time steps and stage-ages
  d.m.CH <- med.time(SV$ST.stageav[,1:a.max.CH], Coef["delta.CHm10"], Coef["delta.CHp"])

  # Pre-adult development time parameter for different time steps and stage-ages
  d.m.PA <- med.time(SV$ST.stageav[,1:a.max.PA], Coef["delta.PAm10"], Coef["delta.PAp"])

  # Arrays with development rates for different stages, time-steps t_local and stage-ages
  # These depend on farm temperatures only, and do not differ between cages
  d.R <- dev.rate(d.m.R, Coef["delta.Rs"])
  d.CH.farm <- dev.rate(d.m.CH, Coef["delta.CHs"])
  d.PA.farm <- dev.rate(d.m.PA, Coef["delta.PAs"])
  for(cage.j in 1:Ncages){ # Fill values into all cages
    d.CH[,,cage.j] <- d.CH.farm
    d.PA[,,cage.j] <- d.PA.farm
  }
  return(
    list(
      m.nat.R = m.nat.R,
      m.nat.CO = m.nat.CO,
      m.nat.CH = m.nat.CH,
      m.nat.PA = m.nat.PA,
      m.nat.AF = m.nat.AF,
      m.clf.PA = m.clf.PA,
      m.clf.AF = m.clf.AF,
      m.trt.CH = m.trt.CH,
      m.trt.PA = m.trt.PA,
      m.trt.AF = m.trt.AF,
      s.R = s.R,
      s.CO = s.CO,
      s.CH = s.CH,
      s.PA = s.PA,
      s.AF = s.AF,
      d.CO = d.CO,
      d.CH = d.CH,
      d.PA = d.PA,
      r.AF = r.AF,
      d.m.R = d.m.R,
      d.m.CH = d.m.CH,
      d.m.PA = d.m.PA,
      d.R = d.R,
      d.CH.farm = d.CH.farm,
      d.PA.farm = d.PA.farm
    )
  )
}

# Functions to count CH lice
update_Y.CH <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV$Y.CH[t_local, , drop = F]
  if(t_local == 1) {
    N.CH.t_local <- SV$N.CH[t_local, , , drop = F]
  } else {
    N.CH.t_local <- SV$N.CH[t_local-1, , , drop = F]
  }
  
  for(cage.j in 1:Ncages){
    if(SV$n.SAL[t_local, cage.j]>0){
      # Probability of sampling chalimus-stage lice (Eqs. 45-47, Aldrin et al. 2017)
      p.CHcount <- invlogit(RE$beta0.CHcount.f + Coef["beta1.CHcount"] * SV$W.SAL[t_local, cage.j])
      
      # Expected number of chalimums-stage lice (Eq. 42, Aldrin et al. 2017)
      E.CH <- SV$n.SAL[t_local, cage.j] * p.CHcount * sum(N.CH.t_local[1, , cage.j])/sum(SV$N.SAL[t_local, cage.j])

      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.CH>0, rnbinom(n=1, mu = E.CH, size = Coef["rho.CH"] * SV$n.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}
update_Y.AF <- function(t_local, ...){
  out <- SV$Y.AF[t_local,  , drop = F]
  if(t_local == 1) {
    N.AF.t_local <- SV$N.AF[t_local, , , drop = F]
  } else {
    N.AF.t_local <- SV$N.AF[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV$n.SAL[t_local, cage.j]>0){
      # Expected number of adult females (Eq. 43, Aldrin et al. 2017)
      E.AF <- SV$n.SAL[t_local, cage.j] * sum(N.AF.t_local[1, , cage.j])/sum(SV$N.SAL[t_local, cage.j])

      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.AF>0, rnbinom(n=1, mu = E.AF, size = Coef["rho.AF"] * SV$n.SAL[t_local, cage.j]), 0) 
    }
  }
  return(out)
}
update_Y.OM <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV$Y.OM[t_local,  , drop = F]
  if(t_local == 1) {
    N.PA.t_local <- SV$N.PA[t_local, , , drop = F]
    N.AM.t_local <- SV$N.AM[t_local, , , drop = F]
  } else {
    N.PA.t_local <- SV$N.PA[t_local-1, , , drop = F]
    N.AM.t_local <- SV$N.AM[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV$n.SAL[t_local, cage.j]>0){
      # Expected number of other mobiles (Eq. 44, Aldrin et al. 2017)
      E.OM <- SV$n.SAL[t_local, cage.j] * (sum(N.PA.t_local[1, , cage.j]) + sum(N.AM.t_local[1, , cage.j]))/sum(SV$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.OM>0, rnbinom(n=1, mu = E.OM, size = Coef["rho.OM"] * SV$n.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}

# # Function to do a second count of lice
update_Y2.CH <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV$Y2.CH[t_local,  , drop = F]
  if(t_local == 1) {
    N.CH.t_local <- SV$N.CH[t_local, , , drop = F]
  } else {
    N.CH.t_local <- SV$N.CH[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV$n2.SAL[t_local, cage.j]>0){
      # Probability of sampling chalimus-stage lice (Eqs. 45-47, Aldrin et al. 2017)
      p.CHcount <- invlogit(RE$beta0.CHcount.f + Coef["beta1.CHcount"] * SV$W.SAL[t_local, cage.j])
      
      # Expected number of chalimums-stage lice (Eq. 42, Aldrin et al. 2017)
      E.CH <- SV$n2.SAL[t_local, cage.j] * p.CHcount * sum(N.CH.t_local[1, , cage.j])/sum(SV$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.CH>0, rnbinom(n=1, mu = E.CH, size = Coef["rho.CH"] * SV$n2.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}
update_Y2.AF <- function(t_local, ...){
  out <- SV$Y2.AF[t_local,  , drop = F]
  if(t_local == 1) {
    N.AF.t_local <- SV$N.AF[t_local, , , drop = F]
  } else {
    N.AF.t_local <- SV$N.AF[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV$n2.SAL[t_local, cage.j]>0){
      # Expected number of adult females (Eq. 43, Aldrin et al. 2017)
      E.AF <- SV$n2.SAL[t_local, cage.j] * sum(N.AF.t_local[1, , cage.j])/sum(SV$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.AF>0, rnbinom(n=1, mu = E.AF, size = Coef["rho.AF"] * SV$n2.SAL[t_local, cage.j]), 0) 
    }
  }
  return(out)
}
update_Y2.OM <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV$Y2.OM[t_local,  , drop = F]
  if(t_local == 1) {
    N.PA.t_local <- SV$N.PA[t_local, , , drop = F]
    N.AM.t_local <- SV$N.AM[t_local, , , drop = F]
  } else {
    N.PA.t_local <- SV$N.PA[t_local-1, , , drop = F]
    N.AM.t_local <- SV$N.AM[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV$n2.SAL[t_local, cage.j]>0){
      # Expected number of other mobiles (Eq. 44, Aldrin et al. 2017)
      E.OM <- SV$n2.SAL[t_local, cage.j] * (sum(N.PA.t_local[1, , cage.j]) + sum(N.AM.t_local[1, , cage.j]))/sum(SV$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.OM>0, rnbinom(n=1, mu = E.OM, size = Coef["rho.OM"] * SV$n2.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}


## Functions to update state variables one time step

## Model for the recruitment stage
update_N.R <- function(t_local, a.max.R = 60,...){
  out <- SV$N.R[t_local, , drop = F] 
    if(t_local > 1){
    # Stage-age 0 (Eq. 1, Aldrin et al. 2017)
    out[1, 1] <- # N recruits at time t_local and stage-age 0
      RE$e.Ext.f[t_local] * # random variation in external recruitment for farm f at time t_local
      SV$N.AF.Ext[t_local-1] * # weighted sum of adult females in surrounding farms
      SV$r.Ext[t_local-1] + # reproductive rate of females
      sum(
        SV$N.AF[t_local-1, , ] * # Number of adult females summed across stage-ages of adult females and cages
          SV$s.AF[t_local-1, , ] * # Survival
          SV$r.AF[t_local-1, , ]   # Development rate
      )
    
    # Stage-age >0 (Eq. 2, Aldrin et al. 2007)
    out[1, 2:a.max.R] <- # N recruits at time t_local and stage-age > 0
      SV$N.R[t_local-1, 1:(a.max.R - 1)] * # Number at same stage at time 1-1 and stage-age a-1
      SV$s.R[t_local-1, 1:(a.max.R - 1)] *                   # Survival
      (1 - SV$d.R[t_local-1, 1:(a.max.R - 1)])               # Development
    
  }
  return(out)
}

## Model for the copepodid stage
update_N.CO <- function(t_local, a.max.CO = 60,...){
  out <- SV$N.CO[t_local, , drop = F]
    if(t_local > 1){
    # Stage-age 0 (Eq. 3, Aldrin et al. 2017)
    out[1, 1] <- sum(
      SV$N.R[t_local-1, ] * # Number of recruits at time t_local-1 at different stage-ages
        SV$s.R[t_local-1, ] * # Survival
        SV$d.R[t_local-1, ]   # Development rate
    )
    
    # Stage-age >0 (Eq. 4, Aldrin et al. 2007)
    if(Ncages == 1){
      out[1, 2:a.max.CO] <- SV$N.CO[t_local-1, 1:(a.max.CO - 1)] * # Number at same stage at time 1-1 and stage-age a-1
        SV$s.CO[t_local-1, 1:(a.max.CO - 1)] *                   # Survival
        (1 - sum(SV$d.CO[t_local-1, 1:(a.max.CO - 1), ]))  # Development
    }
    if(Ncages > 1){
      out[1, 2:a.max.CO] <- SV$N.CO[t_local-1, 1:(a.max.CO - 1)] * # Number at same stage at time 1-1 and stage-age a-1
        SV$s.CO[t_local-1, 1:(a.max.CO - 1)] *                   # Survival
        (1 - rowSums(SV$d.CO[t_local-1, 1:(a.max.CO - 1), ]))  # Development
    }
    
  }
  return(out)
}


## Model for the chalimus stage
update_N.CH <- function(t_local, a.max.CH = 60,...){
  out <- SV$N.CH[t_local, , , drop = F]
  if(t_local > 1){
    # Stage-age 0 (Eq. 5, Aldrin et al. 2017)
    for(cage.j in c(1:Ncages)){
      out[1, 1, cage.j] <- sum( # N chalimus at time t_local and stage-age 0 in cage c
        SV$N.CO[t_local-1, ] * # Number of copepodids at time t_local-1 at different stage-ages
          SV$s.CO[t_local-1, ] * # Survival
          SV$d.CO[t_local-1, , cage.j]   # Development rate
      )
    }
    
    # Stage-age >0 (Eq. 6, Aldrin et al. 2007)
    out[1, 2:a.max.CH, ] <- SV$N.CH[t_local-1, 1:(a.max.CH - 1), ] * # Number at same stage at time 1-1, stage-age a-1, cage c 
      SV$s.CH[t_local-1, 1:(a.max.CH - 1), ] *                   # Survival
      (1 - SV$d.CH[t_local-1, 1:(a.max.CH - 1), ])               # Development
    
  }
  return(out)
}

## Model for the pre-adult stage (assuming no movement of fish between cages)
update_N.PA <- function(t_local, a.max.PA = 60,...){
  out <- SV$N.PA[t_local, , , drop = F]
  if(t_local > 1){
    
    # Stage-age 0 (Eq. 7, Aldrin et al. 2017)
    if(Ncages == 1){
      out[1, 1, ] <- sum( # N pre-adults at time t_local and stage-age 0 in cage c
        SV$N.CH[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV$s.CH[t_local-1, , ] * # Survival
          SV$d.CH[t_local-1, , ]   # Development rate
      )
    }
    if(Ncages > 1){
      out[1, 1, ] <- colSums( # N pre-adults at time t_local and stage-age 0 in cage c
        SV$N.CH[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV$s.CH[t_local-1, , ] * # Survival
          SV$d.CH[t_local-1, , ]   # Development rate
      )
    }
    
    # Stage-age >0 (Eq. 8, Aldrin et al. 2007)
    out[1, 2:a.max.PA, ] <- SV$N.PA[t_local-1, 1:(a.max.PA - 1), ] * # Number at same stage at time 1-1, stage-age a-1, cage c 
      SV$s.PA[t_local-1, 1:(a.max.PA - 1), ] *                   # Survival
      (1 - SV$d.PA[t_local-1, 1:(a.max.PA - 1), ])               # Development
    
  }
  return(out)
}

## Model for the adult stages
update_N.AF <- function(t_local, a.max.AF = 80,...){
  out <- SV$N.AF[t_local, , , drop = F]
  if(t_local > 1){
    
    # Adult females stage-age 0 (Eq. 9, Aldrin et al. 2017)
    if(Ncages == 1){
      out[1, 1, ] <- 0.5 * sum( # N adult females at time t_local and stage-age 0 in cage c
        SV$N.PA[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV$s.PA[t_local-1, , ] * # Survival
          SV$d.PA[t_local-1, , ]   # Development rate
      )
    }
    if(Ncages > 1){
      out[1, 1, ] <- 0.5 * colSums( # N adult females at time t_local and stage-age 0 in cage c
        SV$N.PA[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV$s.PA[t_local-1, , ] * # Survival
          SV$d.PA[t_local-1, , ]   # Development rate
      )
    }
    
    # Adult females stage-age >0 (Eq. 10, Aldrin et al. 2007)
    out[1, 2:a.max.AF, ] <- SV$N.AF[t_local-1, 1:(a.max.AF - 1), ] * # Number at same stage at time 1-1, stage-age a-1, cage c 
      SV$s.AF[t_local-1, 1:(a.max.AF - 1), ]                     # Survival
    
  }
  return(out)
}
    
# Adult males (Eq. 11, Aldrin et al. 2007)
update_N.AM <- function(t_local, ...){
  out <- SV$N.AF[t_local, , , drop=F]
  return(out)
}


## Cleaner fish model (Eq. 40, Aldrin et al. 2017), Eqs. 1-3 Aldrin & Huseby 2019

update_N.lump <- function(t_local,...){
  if(t_local==1) out <- SV$N.lump[t_local,] + SV$S.lump[t_local,]
  if(t_local>1) out <- SV$N.lump[t_local-1,] * (1 - Coef["kappa.lump"]) *
    (1 - Coef["kappa.therm"] * SV$use.therm[t_local-1,]) *
    (1 - Coef["kappa.mech"] * SV$use.mech[t_local-1,]) *
    (1 - Coef["kappa.freshw"] * SV$use.freshw[t_local-1,]) +
    SV$S.lump[t_local,]
  return(out)
}

update_N.wrasse <- function(t_local,...){
  if(t_local==1) out <- SV$N.wrasse[t_local,] + SV$S.wrasse[t_local,]
  if(t_local>1) out <- SV$N.wrasse[t_local-1,] * (1 - Coef["kappa.wrasse"]) *
      (1 - Coef["kappa.therm"] * SV$use.therm[t_local-1,]) *
      (1 - Coef["kappa.mech"] * SV$use.mech[t_local-1,]) *
      (1 - Coef["kappa.freshw"] * SV$use.freshw[t_local-1,]) +
      SV$S.wrasse[t_local,]
  return(out)
}
  
## Cleaner fish-caused mortality of salmon lice (Eq. 18, Aldrin et al. 2017) 
## Eqs. 11-12 Aldrin & Huseby 2019
update_m.clf.PA <- function(t_local,...){
  out <- SV$m.clf.PA[t_local, , , drop = F]
  if(Ncages ==1) N.lice.cages <- sum(SV$N.PA[t_local, , ]) + sum(SV$N.AF[t_local, , ])
  if(Ncages >1) N.lice.cages <- colSums(SV$N.PA[t_local, , ]) + colSums(SV$N.AF[t_local, , ])
  which.lice <- which(N.lice.cages > 0)
  m.lump.cages <- rep(0, Ncages)
  m.lump.cages[which.lice] <- (1 - exp(-Coef["gamma.lump"]) * 
                                 N.lice.cages / SV$N.SAL[t_local, which.lice]) * 
    Coef["lambda.lump"] * SV$N.lump[t_local, which.lice] / N.lice.cages[which.lice]
  m.wrasse.cages <- rep(0, Ncages)
  m.wrasse.cages[which.lice] <- (1 - exp(-Coef["gamma.wrasse"]) * 
                                 N.lice.cages / SV$N.SAL[t_local, which.lice]) * 
    Coef["lambda.wrasse"] * SV$N.wrasse[t_local, which.lice] / N.lice.cages[which.lice]
  m.clf.cages <- pmin(m.lump.cages + m.wrasse.cages, 0.5)
  for(cage.j in 1:Ncages){
    out[1, , cage.j] <- m.clf.cages[cage.j]
  }
  return(out)
}  

update_m.clf.AF <- function(t_local,...){
  out <- SV$m.clf.AF[t_local, , , drop = F]
  if(Ncages ==1) N.lice.cages <- sum(SV$N.PA[t_local, , ]) + sum(SV$N.AF[t_local, , ])
  if(Ncages >1) N.lice.cages <- colSums(SV$N.PA[t_local, , ]) + colSums(SV$N.AF[t_local, , ])
  which.lice <- which(N.lice.cages > 0)
  m.lump.cages <- rep(0, Ncages)
  m.lump.cages[which.lice] <- (1 - exp(-Coef["gamma.lump"]) * 
                                 N.lice.cages / SV$N.SAL[t_local, which.lice]) * 
    Coef["lambda.lump"] * SV$N.lump[t_local, which.lice] / N.lice.cages[which.lice]
  m.wrasse.cages <- rep(0, Ncages)
  m.wrasse.cages[which.lice] <- (1 - exp(-Coef["gamma.wrasse"]) * 
                                   N.lice.cages / SV$N.SAL[t_local, which.lice]) * 
    Coef["lambda.wrasse"] * SV$N.wrasse[t_local, which.lice] / N.lice.cages[which.lice]
  m.clf.cages <- pmin(m.lump.cages + m.wrasse.cages, 0.5)
  for(cage.j in 1:Ncages){
    out[1, , cage.j] <- m.clf.cages[cage.j]
  }
  return(out)
}  

## Calculate mortality in subsequent time-steps 
## caused by chemical treatments applied in time step t_local

# Hydrogen peroxide
update_ux.HPcht <- function(t_local, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV[c("ux.HPcht.PA", "ux.HPcht.AF")]
  if(sum(SV$use.HPcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    HPcht.eff <- trt.sample("HPcht") # Draw a posterior sample with cage-specific effects for each stage
    HPcht.eff <- log(1 + exp(HPcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){             # Loop on cages
      if(SV$use.HPcht[t_local, cage.j]==1) {      # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.HPcht"], Ndays)       # Start time of effect (Eq. 19)     
        t2 <- min(t1 + Coef["delta.dur.HPcht"] - 1, Ndays)  # End time of effect (Eq. 19)
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.PA){
            out$ux.HPcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- HPcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.HPcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- HPcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Deltamethrin
update_ux.DMcht <- function(t_local, a.max.CH = 60, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV[c("ux.DMcht.CH", "ux.DMcht.PA", "ux.DMcht.AF")]
    if(sum(SV$use.DMcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    DMcht.eff <- trt.sample("DMcht") # Draw a posterior sample with cage-specific effects for each stage
    DMcht.eff <- log(1 + exp(DMcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV$use.DMcht[t_local, cage.j]==1) {     # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.DMcht"], Ndays)       # Start time of effect (Eq. 19)    
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.DMcht"] / max(SV$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.CH){
            out$ux.DMcht.CH[t_local.eff, a.min:a.max.CH, cage.j] <- DMcht.eff[cage.j, "CH"] # Effect on CH
          }
          if(a.min <= a.max.PA){
            out$ux.DMcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- DMcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.DMcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- DMcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Azamethiphos
update_ux.AZcht <- function(t_local, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV[c("ux.AZcht.PA", "ux.AZcht.AF")]
    if(sum(SV$use.AZcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    AZcht.eff <- trt.sample("AZcht") # Draw a posterior sample with cage-specific effects for each stage
    AZcht.eff <- log(1 + exp(AZcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      # When azamethiphos is used in combination with deltamethrin or cyper-methrin, 
      # we assume it has the same effect as using deltamethrin or cypermethrin alone.
      use.DMCMcht.t_local <- F # Indicator variable for DM or CM use
      if(SV$use.DMcht[t_local, cage.j]!=0 | SV$use.CMcht[t_local, cage.j]!=0){use.DMCMcht.t_local <- T} # Used the same day
      if(t_local>1){if(SV$use.DMcht[t_local-1, cage.j]!=0 | SV$use.CMcht[t_local-1, cage.j]!=0) use.DMCMcht.t_local <- T} # Used the day before
      if(SV$ux.DMcht.AF[t_local, a.max.AF, cage.j] != 0){use.DMCMcht.t_local <- T}   # Used previously, but still in effect
      if(SV$ux.CMcht.AF[t_local, a.max.AF, cage.j] != 0){use.DMCMcht.t_local <- T}
      if(SV$use.AZcht[t_local, cage.j]==1 & SV$use.DMCMcht.t_local==F) {        # If cage.j was treated and DM or CM not used
        t1 <- min(t_local + Coef["delta.del.AZcht"], Ndays)       # Start time of effect     
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.AZcht"] / max(SV$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.PA){
            out$ux.AZcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- AZcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.AZcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- AZcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Emamectin
update_ux.EMcht <- function(t_local, a.max.CH = 60, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV[c("ux.EMcht.CH", "ux.EMcht.PA", "ux.EMcht.AF")]
    if(sum(SV$use.EMcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    EMcht.eff <- trt.sample("EMcht") # Draw a posterior sample with cage-specific effects for each stage
    EMcht.eff <- log(1 + exp(EMcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV$use.EMcht[t_local, cage.j]==1) {     # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.EMcht"], Ndays)       # Start time of effect     
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.EMcht"] / max(SV$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.CH){
            out$ux.EMcht.CH[t_local.eff, a.min:a.max.CH, cage.j] <- EMcht.eff[cage.j, "CH"] # Effect on CH
          }
          if(a.min <= a.max.PA){
            out$ux.EMcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- EMcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.EMcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- EMcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Diflubenzuron
update_ux.DBcht <- function(t_local, a.max.CH = 60, a.max.PA = 60,...){
  out <- SV[c("ux.DBcht.CH", "ux.DBcht.PA")]
    if(sum(SV$use.DBcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    DBcht.eff <- trt.sample("DBcht") # Draw a posterior sample with cage-specific effects for each stage
    DBcht.eff <- log(1 + exp(DBcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV$use.DBcht[t_local, cage.j]==1) {     # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.DBcht"], Ndays)       # Start time of effect     
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.DBcht"] / max(SV$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.CH){
            out$ux.DBcht.CH[t_local.eff, a.min:a.max.CH, cage.j] <- DBcht.eff[cage.j, "CH"] # Effect on CH
          }
          if(a.min <= a.max.PA){
            out$ux.DBcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- DBcht.eff[cage.j, "PA"] # Effect on PA
          }
        }
      }
    }
  }
  return(out)
}

## Calculate mortality in subsequent time-steps 
## caused by physical treatments applied in time step t_local

# Thermal
update_ux.therm <- function(t_local,...){
  out <- SV[c("ux.therm.CH", "ux.therm.PA", "ux.therm.AF")]
    if(sum(SV$use.therm[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    therm.eff <- trt.sample("therm") # Draw a posterior sample with cage-specific effects for each stage
    therm.eff <- log(1 + exp(therm.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV$use.therm[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.therm.CH[t_local, , cage.j] <- therm.eff[cage.j, "CH"] # Effect on CH
        out$ux.therm.PA[t_local, , cage.j] <- therm.eff[cage.j, "PA"] # Effect on PA
        out$ux.therm.AF[t_local, , cage.j] <- therm.eff[cage.j, "Af"] # Effect on Af
      }
    }
  }
  return(out)
}

# Freshwater
update_ux.freshw <- function(t_local,...){
  out <- SV[c("ux.freshw.CH", "ux.freshw.PA", "ux.freshw.AF")]
    if(sum(SV$use.freshw[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    freshw.eff <- trt.sample("freshw") # Draw a posterior sample with cage-specific effects for each stage
    freshw.eff <- log(1 + exp(freshw.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV$use.freshw[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.freshw.CH[t_local, , cage.j] <- freshw.eff[cage.j, "CH"] # Effect on CH
        out$ux.freshw.PA[t_local, , cage.j] <- freshw.eff[cage.j, "PA"] # Effect on PA
        out$ux.freshw.AF[t_local, , cage.j] <- freshw.eff[cage.j, "Af"] # Effect on Af
      }
    }
  }
  return(out)
}

# Mechanical
update_ux.mech <- function(t_local,...){
  out <- SV[c("ux.mech.CH", "ux.mech.PA", "ux.mech.AF")]
    if(sum(SV$use.mech[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    mech.eff <- trt.sample("mech") # Draw a posterior sample with cage-specific effects for each stage
    mech.eff <- log(1 + exp(mech.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV$use.mech[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.mech.CH[t_local, , cage.j] <- mech.eff[cage.j, "CH"] # Effect on CH
        out$ux.mech.PA[t_local, , cage.j] <- mech.eff[cage.j, "PA"] # Effect on PA
        out$ux.mech.AF[t_local, , cage.j] <- mech.eff[cage.j, "Af"] # Effect on Af
      }
    }
  }
  return(out)
}

# Fixed treatment mortality = M.trt
update_ux.fx <- function(t_local, M.trt = 0.5, ...){
  out <- SV[c("ux.fx.CH", "ux.fx.PA", "ux.fx.AF")]
  if(sum(SV$use.fx[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    # Trt.eff <- log(M.trt / (1 - M.trt))
    # Trt.eff <- log(1 + exp(Trt.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    Trt.eff <- (-log(1 - M.trt)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){             # Loop on cages
      if(SV$use.fx[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.fx.CH[t_local, , cage.j] <- Trt.eff # Effect on CH
        out$ux.fx.PA[t_local, , cage.j] <- Trt.eff # Effect on PA
        out$ux.fx.AF[t_local, , cage.j] <- Trt.eff # Effect on Af
      }
    }
  }
  return(out)
}

# Mortality due to chemical and physical treatments (Eq. 21, Aldrin et al. 2017)
update_m.trt.CH <- function(t_local,...){
  out <- SV$m.trt.CH[t_local, , , drop = F]
  ux.sum <- SV$ux.DMcht.CH[t_local, , ] + 
            SV$ux.EMcht.CH[t_local, , ] +
            SV$ux.DBcht.CH[t_local, , ] +
            SV$ux.therm.CH[t_local, , ] +
            SV$ux.freshw.CH[t_local, , ] +
            SV$ux.mech.CH[t_local, , ] +
            SV$ux.fx.CH[t_local, , ]
  m.sum <- 1 - exp(-ux.sum)
  out[1, , ] <- m.sum  
  return(out)
}

update_m.trt.PA <- function(t_local,...){
  out <- SV$m.trt.PA[t_local, , , drop = F]
  ux.sum <- SV$ux.HPcht.PA[t_local, , ] +
            SV$ux.DMcht.PA[t_local, , ] + 
            SV$ux.AZcht.PA[t_local, , ] +
            SV$ux.EMcht.PA[t_local, , ] +
            SV$ux.DBcht.PA[t_local, , ] +
            SV$ux.therm.PA[t_local, , ] +
            SV$ux.freshw.PA[t_local, , ] +
            SV$ux.mech.PA[t_local, , ] +
            SV$ux.fx.PA[t_local, , ]
  m.sum <- 1 - exp(-ux.sum)
  out[1, , ] <- m.sum  
  return(out)
}

update_m.trt.AF <- function(t_local,...){
  out <- SV$m.trt.AF[t_local, , , drop = F]
  ux.sum <- 
    SV$ux.HPcht.AF[t_local, , ] +
    SV$ux.DMcht.AF[t_local, , ] + 
    SV$ux.AZcht.AF[t_local, , ] +
    SV$ux.EMcht.AF[t_local, , ] +
    SV$ux.therm.AF[t_local, , ] +
    SV$ux.freshw.AF[t_local, , ] +
    SV$ux.mech.AF[t_local, , ] +
    SV$ux.fx.AF[t_local, , ]
  m.sum <- 1 - exp(-ux.sum)
  out[1, , ] <- m.sum  
  return(out)
}

## Survival rates (Eq. 12, Aldrin et al. 2017)
update_s.CH <- function(t_local){
  out <- (1 - SV$m.nat.CH[t_local,,]) * (1 - SV$m.trt.CH[t_local,,])
  return(out)
  }
update_s.PA <- function(t_local){
  out <- (1 - SV$m.nat.PA[t_local,,]) * (1 - SV$m.clf.PA[t_local,,]) * (1 - SV$m.trt.PA[t_local,,])
  return(out)
  }
update_s.AF <- function(t_local){
  out <- (1 - SV$m.nat.AF[t_local,,]) * (1 - SV$m.clf.AF[t_local,,]) * (1 - SV$m.trt.AF[t_local,,])
  return(out)}

## Infection rate
update_d.CO <- function(t_local, a.max.CO = 60,...){
  out <- SV$d.CO[t_local, , , drop = F]
  if(t_local < Ndays){
    # Cage-specific expected rate at logit scale (Eq. 28, Aldrin et al. 2017), Eq. 14 Aldrin & Huseby 2019
    eta.CO.t_local <- RE$delta0.CO.f + log(SV$N.SAL[t_local+1,]) + Coef["delta1.CO"] * (log(SV$W.SAL[t_local+1,]) - 0.55) +
      Coef["delta2.CO"] * (SV$ST[t_local+1] - 9) + Coef["delta3.CO"] * (SV$ST[t_local+1] - 9)^2
    
    # Cage-specific expected rate at proportion scale (Eq. 27, Aldrin et al. 2017), Eq. 13 Aldrin & Huseby 2019
    d.CO.t_local <- exp(eta.CO.t_local) / (1 + sum(exp(eta.CO.t_local)))
    
    # Age-specific values
    out[,,] <- rep(d.CO.t_local, each=a.max.CO)
  } else out[,,] <- 0
  return(out)  
}


## Reproduction

#  Reproduction rate for internal recruitment at time t_local (Eq. 31, Aldrin et al. 2017)
update_r.AF <- function(t_local, a.max.AF = 80, ...){

  # Median hatching time f(Eq. 32, Aldrin et al. 2017)
  d.m.E.t_local <- Coef["delta.Em10"]*(10/SV$ST[t_local])^Coef["delta.Rp"]
  
  out <- SV$r.AF[t_local, , , drop = F]

  r.raw <- Coef["beta0.R"] *          # N eggs first extrusion
    (1:a.max.AF)^Coef["beta1.R"] *    # N eggs different stage-ages of AF
    ((d.m.E.t_local + 1)^-1)                # Extrusion rate

  allee <-                          # Allee effect (differs between cages)
  (1 - exp(-Coef["gamma.r"] * sum(SV$N.AF[t_local, , ]) / SV$N.SAL[t_local, ]))

  for(cage.j in 1:Ncages){
    out[1, , cage.j] <- r.raw * allee[cage.j]
  }
  return(out)
}

## External recruitment (section 2.2.7, Aldrin et al. 2017)
update_r.Ext <- function(t_local, ...){
  
  # Median hatching time f(Eq. 32, Aldrin et al. 2017)
  d.m.E.t <- Coef["delta.Em10"]*(10/SV$ST[t])^Coef["delta.Rp"]
  
  # Allee effect
  allee <-                          
    (1 - exp(-Coef["gamma.r"] * SV$A.Ext[t]))
  
  out <- Coef["beta0.R"] *
    11^Coef["beta1.R"] *
    (1 / (d.m.E.t + 1)) *
    allee
  
  return(out)
}


########################## 


# Functions to count CH lice
update_Y.CH2 <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV_local_local$Y.CH[t_local, , drop = F]
  if(t_local == 1) {
    N.CH.t_local <- SV_local_local$N.CH[t_local, , , drop = F]
  } else {
    N.CH.t_local <- SV_local_local$N.CH[t_local-1, , , drop = F]
  }
  
  for(cage.j in 1:Ncages){
    if(SV_local_local$n.SAL[t_local, cage.j]>0){
      # Probability of sampling chalimus-stage lice (Eqs. 45-47, Aldrin et al. 2017)
      p.CHcount <- invlogit(RE$beta0.CHcount.f + Coef["beta1.CHcount"] * SV_local_local$W.SAL[t_local, cage.j])
      
      # Expected number of chalimums-stage lice (Eq. 42, Aldrin et al. 2017)
      E.CH <- SV_local_local$n.SAL[t_local, cage.j] * p.CHcount * sum(N.CH.t_local[1, , cage.j])/sum(SV_local_local$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.CH>0, rnbinom(n=1, mu = E.CH, size = Coef["rho.CH"] * SV_local_local$n.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}
update_Y.AF2 <- function(t_local, ...){
  out <- SV_local_local$Y.AF[t_local,  , drop = F]
  if(t_local == 1) {
    N.AF.t_local <- SV_local_local$N.AF[t_local, , , drop = F]
  } else {
    N.AF.t_local <- SV_local_local$N.AF[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV_local_local$n.SAL[t_local, cage.j]>0){
      # Expected number of adult females (Eq. 43, Aldrin et al. 2017)
      E.AF <- SV_local$n.SAL[t_local, cage.j] * sum(N.AF.t_local[1, , cage.j])/sum(SV_local$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.AF>0, rnbinom(n=1, mu = E.AF, size = Coef["rho.AF"] * SV_local$n.SAL[t_local, cage.j]), 0) 
    }
  }
  return(out)
}
update_Y.OM2 <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV_local$Y.OM[t_local,  , drop = F]
  if(t_local == 1) {
    N.PA.t_local <- SV_local$N.PA[t_local, , , drop = F]
    N.AM.t_local <- SV_local$N.AM[t_local, , , drop = F]
  } else {
    N.PA.t_local <- SV_local$N.PA[t_local-1, , , drop = F]
    N.AM.t_local <- SV_local$N.AM[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV_local$n.SAL[t_local, cage.j]>0){
      # Expected number of other mobiles (Eq. 44, Aldrin et al. 2017)
      E.OM <- SV_local$n.SAL[t_local, cage.j] * (sum(N.PA.t_local[1, , cage.j]) + sum(N.AM.t_local[1, , cage.j]))/sum(SV_local$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.OM>0, rnbinom(n=1, mu = E.OM, size = Coef["rho.OM"] * SV_local$n.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}

# # Function to do a second count of lice
update_Y2.CH2 <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV_local$Y2.CH[t_local,  , drop = F]
  if(t_local == 1) {
    N.CH.t_local <- SV_local$N.CH[t_local, , , drop = F]
  } else {
    N.CH.t_local <- SV_local$N.CH[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV_local$n2.SAL[t_local, cage.j]>0){
      # Probability of sampling chalimus-stage lice (Eqs. 45-47, Aldrin et al. 2017)
      p.CHcount <- invlogit(RE$beta0.CHcount.f + Coef["beta1.CHcount"] * SV_local$W.SAL[t_local, cage.j])
      
      # Expected number of chalimums-stage lice (Eq. 42, Aldrin et al. 2017)
      E.CH <- SV_local$n2.SAL[t_local, cage.j] * p.CHcount * sum(N.CH.t_local[1, , cage.j])/sum(SV_local$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.CH>0, rnbinom(n=1, mu = E.CH, size = Coef["rho.CH"] * SV_local$n2.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}
update_Y2.AF2 <- function(t_local, ...){
  out <- SV_local$Y2.AF[t_local,  , drop = F]
  if(t_local == 1) {
    N.AF.t_local <- SV_local$N.AF[t_local, , , drop = F]
  } else {
    N.AF.t_local <- SV_local$N.AF[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV_local$n2.SAL[t_local, cage.j]>0){
      # Expected number of adult females (Eq. 43, Aldrin et al. 2017)
      E.AF <- SV_local$n2.SAL[t_local, cage.j] * sum(N.AF.t_local[1, , cage.j])/sum(SV_local$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.AF>0, rnbinom(n=1, mu = E.AF, size = Coef["rho.AF"] * SV_local$n2.SAL[t_local, cage.j]), 0) 
    }
  }
  return(out)
}
update_Y2.OM2 <- function(t_local, ...){
  ## As lice are counted at the start of the day, before mortality and development,
  ## we use the lice numbers from the day before (these are updated later)
  # (section 2.2.1 of Aldrin et al. 2017)
  out <- SV_local$Y2.OM[t_local,  , drop = F]
  if(t_local == 1) {
    N.PA.t_local <- SV_local$N.PA[t_local, , , drop = F]
    N.AM.t_local <- SV_local$N.AM[t_local, , , drop = F]
  } else {
    N.PA.t_local <- SV_local$N.PA[t_local-1, , , drop = F]
    N.AM.t_local <- SV_local$N.AM[t_local-1, , , drop = F]
  }
  for(cage.j in 1:Ncages){
    if(SV_local$n2.SAL[t_local, cage.j]>0){
      # Expected number of other mobiles (Eq. 44, Aldrin et al. 2017)
      E.OM <- SV_local$n2.SAL[t_local, cage.j] * (sum(N.PA.t_local[1, , cage.j]) + sum(N.AM.t_local[1, , cage.j]))/sum(SV_local$N.SAL[t_local, cage.j])
      
      # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
      out[cage.j] <- ifelse(E.OM>0, rnbinom(n=1, mu = E.OM, size = Coef["rho.OM"] * SV_local$n2.SAL[t_local, cage.j]), 0)
    }
  }
  return(out)
}


## Functions to update state variables one time step

## Model for the recruitment stage
update_N.R2 <- function(t_local, a.max.R = 60,...){
  out <- SV_local$N.R[t_local, , drop = F] 
  if(t_local > 1){
    # Stage-age 0 (Eq. 1, Aldrin et al. 2017)
    out[1, 1] <- # N recruits at time t_local and stage-age 0
      RE$e.Ext.f[t_local] * # random variation in external recruitment for farm f at time t_local
      SV_local$N.AF.Ext[t_local-1] * # weighted sum of adult females in surrounding farms
      SV_local$r.Ext[t_local-1] + # reproductive rate of females
      sum(
        SV_local$N.AF[t_local-1, , ] * # Number of adult females summed across stage-ages of adult females and cages
          SV_local$s.AF[t_local-1, , ] * # Survival
          SV_local$r.AF[t_local-1, , ]   # Development rate
      )
    
    # Stage-age >0 (Eq. 2, Aldrin et al. 2007)
    out[1, 2:a.max.R] <- # N recruits at time t_local and stage-age > 0
      SV_local$N.R[t_local-1, 1:(a.max.R - 1)] * # Number at same stage at time 1-1 and stage-age a-1
      SV_local$s.R[t_local-1, 1:(a.max.R - 1)] *                   # Survival
      (1 - SV_local$d.R[t_local-1, 1:(a.max.R - 1)])               # Development
    
  }
  return(out)
}

## Model for the copepodid stage
update_N.CO2 <- function(t_local, a.max.CO = 60,...){
  out <- SV_local$N.CO[t_local, , drop = F]
  if(t_local > 1){
    # Stage-age 0 (Eq. 3, Aldrin et al. 2017)
    out[1, 1] <- sum(
      SV_local$N.R[t_local-1, ] * # Number of recruits at time t_local-1 at different stage-ages
        SV_local$s.R[t_local-1, ] * # Survival
        SV_local$d.R[t_local-1, ]   # Development rate
    )
    
    # Stage-age >0 (Eq. 4, Aldrin et al. 2007)
    if(Ncages == 1){
      out[1, 2:a.max.CO] <- SV_local$N.CO[t_local-1, 1:(a.max.CO - 1)] * # Number at same stage at time 1-1 and stage-age a-1
        SV_local$s.CO[t_local-1, 1:(a.max.CO - 1)] *                   # Survival
        (1 - sum(SV_local$d.CO[t_local-1, 1:(a.max.CO - 1), ]))  # Development
    }
    if(Ncages > 1){
      out[1, 2:a.max.CO] <- SV_local$N.CO[t_local-1, 1:(a.max.CO - 1)] * # Number at same stage at time 1-1 and stage-age a-1
        SV_local$s.CO[t_local-1, 1:(a.max.CO - 1)] *                   # Survival
        (1 - rowSums(SV_local$d.CO[t_local-1, 1:(a.max.CO - 1), ]))  # Development
    }
    
  }
  return(out)
}


## Model for the chalimus stage
update_N.CH2 <- function(t_local, a.max.CH = 60,...){
  out <- SV_local$N.CH[t_local, , , drop = F]
  if(t_local > 1){
    # Stage-age 0 (Eq. 5, Aldrin et al. 2017)
    for(cage.j in c(1:Ncages)){
      out[1, 1, cage.j] <- sum( # N chalimus at time t_local and stage-age 0 in cage c
        SV_local$N.CO[t_local-1, ] * # Number of copepodids at time t_local-1 at different stage-ages
          SV_local$s.CO[t_local-1, ] * # Survival
          SV_local$d.CO[t_local-1, , cage.j]   # Development rate
      )
    }
    
    # Stage-age >0 (Eq. 6, Aldrin et al. 2007)
    out[1, 2:a.max.CH, ] <- SV_local$N.CH[t_local-1, 1:(a.max.CH - 1), ] * # Number at same stage at time 1-1, stage-age a-1, cage c 
      SV_local$s.CH[t_local-1, 1:(a.max.CH - 1), ] *                   # Survival
      (1 - SV_local$d.CH[t_local-1, 1:(a.max.CH - 1), ])               # Development
    
  }
  return(out)
}

## Model for the pre-adult stage (assuming no movement of fish between cages)
update_N.PA2 <- function(t_local, a.max.PA = 60,...){
  out <- SV_local$N.PA[t_local, , , drop = F]
  if(t_local > 1){
    
    # Stage-age 0 (Eq. 7, Aldrin et al. 2017)
    if(Ncages == 1){
      out[1, 1, ] <- sum( # N pre-adults at time t_local and stage-age 0 in cage c
        SV_local$N.CH[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV_local$s.CH[t_local-1, , ] * # Survival
          SV_local$d.CH[t_local-1, , ]   # Development rate
      )
    }
    if(Ncages > 1){
      out[1, 1, ] <- colSums( # N pre-adults at time t_local and stage-age 0 in cage c
        SV_local$N.CH[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV_local$s.CH[t_local-1, , ] * # Survival
          SV_local$d.CH[t_local-1, , ]   # Development rate
      )
    }
    
    # Stage-age >0 (Eq. 8, Aldrin et al. 2007)
    out[1, 2:a.max.PA, ] <- SV_local$N.PA[t_local-1, 1:(a.max.PA - 1), ] * # Number at same stage at time 1-1, stage-age a-1, cage c 
      SV_local$s.PA[t_local-1, 1:(a.max.PA - 1), ] *                   # Survival
      (1 - SV_local$d.PA[t_local-1, 1:(a.max.PA - 1), ])               # Development
    
  }
  return(out)
}

## Model for the adult stages
update_N.AF2 <- function(t_local, a.max.AF = 80,...){
  out <- SV_local$N.AF[t_local, , , drop = F]
  if(t_local > 1){
    
    # Adult females stage-age 0 (Eq. 9, Aldrin et al. 2017)
    if(Ncages == 1){
      out[1, 1, ] <- 0.5 * sum( # N adult females at time t_local and stage-age 0 in cage c
        SV_local$N.PA[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV_local$s.PA[t_local-1, , ] * # Survival
          SV_local$d.PA[t_local-1, , ]   # Development rate
      )
    }
    if(Ncages > 1){
      out[1, 1, ] <- 0.5 * colSums( # N adult females at time t_local and stage-age 0 in cage c
        SV_local$N.PA[t_local-1, , ] * # Number of chalimus at time t_local-1 at different stage-ages
          SV_local$s.PA[t_local-1, , ] * # Survival
          SV_local$d.PA[t_local-1, , ]   # Development rate
      )
    }
    
    # Adult females stage-age >0 (Eq. 10, Aldrin et al. 2007)
    out[1, 2:a.max.AF, ] <- SV_local$N.AF[t_local-1, 1:(a.max.AF - 1), ] * # Number at same stage at time 1-1, stage-age a-1, cage c 
      SV_local$s.AF[t_local-1, 1:(a.max.AF - 1), ]                     # Survival
    
  }
  return(out)
}

# Adult males (Eq. 11, Aldrin et al. 2007)
update_N.AM2 <- function(t_local, ...){
  out <- SV_local$N.AF[t_local, , , drop=F]
  return(out)
}


## Cleaner fish model (Eq. 40, Aldrin et al. 2017), Eqs. 1-3 Aldrin & Huseby 2019

update_N.lump2 <- function(t_local,...){
  if(t_local==1) out <- SV_local$N.lump[t_local,] + SV_local$S.lump[t_local,]
  if(t_local>1) out <- SV_local$N.lump[t_local-1,] * (1 - Coef["kappa.lump"]) *
      (1 - Coef["kappa.therm"] * SV_local$use.therm[t_local-1,]) *
      (1 - Coef["kappa.mech"] * SV_local$use.mech[t_local-1,]) *
      (1 - Coef["kappa.freshw"] * SV_local$use.freshw[t_local-1,]) +
      SV_local$S.lump[t_local,]
  return(out)
}

update_N.wrasse2 <- function(t_local,...){
  if(t_local==1) out <- SV_local$N.wrasse[t_local,] + SV_local$S.wrasse[t_local,]
  if(t_local>1) out <- SV_local$N.wrasse[t_local-1,] * (1 - Coef["kappa.wrasse"]) *
      (1 - Coef["kappa.therm"] * SV_local$use.therm[t_local-1,]) *
      (1 - Coef["kappa.mech"] * SV_local$use.mech[t_local-1,]) *
      (1 - Coef["kappa.freshw"] * SV_local$use.freshw[t_local-1,]) +
      SV_local$S.wrasse[t_local,]
  return(out)
}

## Cleaner fish-caused mortality of salmon lice (Eq. 18, Aldrin et al. 2017) 
## Eqs. 11-12 Aldrin & Huseby 2019
update_m.clf.PA2 <- function(t_local,...){
  out <- SV_local$m.clf.PA[t_local, , , drop = F]
  if(Ncages ==1) N.lice.cages <- sum(SV_local$N.PA[t_local, , ]) + sum(SV_local$N.AF[t_local, , ])
  if(Ncages >1) N.lice.cages <- colSums(SV_local$N.PA[t_local, , ]) + colSums(SV_local$N.AF[t_local, , ])
  which.lice <- which(N.lice.cages > 0)
  m.lump.cages <- rep(0, Ncages)
  m.lump.cages[which.lice] <- (1 - exp(-Coef["gamma.lump"]) * 
                                 N.lice.cages / SV_local$N.SAL[t_local, which.lice]) * 
    Coef["lambda.lump"] * SV_local$N.lump[t_local, which.lice] / N.lice.cages[which.lice]
  m.wrasse.cages <- rep(0, Ncages)
  m.wrasse.cages[which.lice] <- (1 - exp(-Coef["gamma.wrasse"]) * 
                                   N.lice.cages / SV_local$N.SAL[t_local, which.lice]) * 
    Coef["lambda.wrasse"] * SV_local$N.wrasse[t_local, which.lice] / N.lice.cages[which.lice]
  m.clf.cages <- pmin(m.lump.cages + m.wrasse.cages, 0.5)
  for(cage.j in 1:Ncages){
    out[1, , cage.j] <- m.clf.cages[cage.j]
  }
  return(out)
}  

update_m.clf.AF2 <- function(t_local,...){
  out <- SV_local$m.clf.AF[t_local, , , drop = F]
  if(Ncages ==1) N.lice.cages <- sum(SV_local$N.PA[t_local, , ]) + sum(SV_local$N.AF[t_local, , ])
  if(Ncages >1) N.lice.cages <- colSums(SV_local$N.PA[t_local, , ]) + colSums(SV_local$N.AF[t_local, , ])
  which.lice <- which(N.lice.cages > 0)
  m.lump.cages <- rep(0, Ncages)
  m.lump.cages[which.lice] <- (1 - exp(-Coef["gamma.lump"]) * 
                                 N.lice.cages / SV_local$N.SAL[t_local, which.lice]) * 
    Coef["lambda.lump"] * SV_local$N.lump[t_local, which.lice] / N.lice.cages[which.lice]
  m.wrasse.cages <- rep(0, Ncages)
  m.wrasse.cages[which.lice] <- (1 - exp(-Coef["gamma.wrasse"]) * 
                                   N.lice.cages / SV_local$N.SAL[t_local, which.lice]) * 
    Coef["lambda.wrasse"] * SV_local$N.wrasse[t_local, which.lice] / N.lice.cages[which.lice]
  m.clf.cages <- pmin(m.lump.cages + m.wrasse.cages, 0.5)
  for(cage.j in 1:Ncages){
    out[1, , cage.j] <- m.clf.cages[cage.j]
  }
  return(out)
}  

## Calculate mortality in subsequent time-steps 
## caused by chemical treatments applied in time step t_local

# Hydrogen peroxide
update_ux.HPcht2 <- function(t_local, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV_local[c("ux.HPcht.PA", "ux.HPcht.AF")]
  if(sum(SV_local$use.HPcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    HPcht.eff <- trt.sample("HPcht") # Draw a posterior sample with cage-specific effects for each stage
    HPcht.eff <- log(1 + exp(HPcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){             # Loop on cages
      if(SV_local$use.HPcht[t_local, cage.j]==1) {      # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.HPcht"], Ndays)       # Start time of effect (Eq. 19)     
        t2 <- min(t1 + Coef["delta.dur.HPcht"] - 1, Ndays)  # End time of effect (Eq. 19)
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.PA){
            out$ux.HPcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- HPcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.HPcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- HPcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Deltamethrin
update_ux.DMcht2 <- function(t_local, a.max.CH = 60, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV_local[c("ux.DMcht.CH", "ux.DMcht.PA", "ux.DMcht.AF")]
  if(sum(SV_local$use.DMcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    DMcht.eff <- trt.sample("DMcht") # Draw a posterior sample with cage-specific effects for each stage
    DMcht.eff <- log(1 + exp(DMcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV_local$use.DMcht[t_local, cage.j]==1) {     # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.DMcht"], Ndays)       # Start time of effect (Eq. 19)    
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.DMcht"] / max(SV_local$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.CH){
            out$ux.DMcht.CH[t_local.eff, a.min:a.max.CH, cage.j] <- DMcht.eff[cage.j, "CH"] # Effect on CH
          }
          if(a.min <= a.max.PA){
            out$ux.DMcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- DMcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.DMcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- DMcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Azamethiphos
update_ux.AZcht2 <- function(t_local, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV_local[c("ux.AZcht.PA", "ux.AZcht.AF")]
  if(sum(SV_local$use.AZcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    AZcht.eff <- trt.sample("AZcht") # Draw a posterior sample with cage-specific effects for each stage
    AZcht.eff <- log(1 + exp(AZcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      # When azamethiphos is used in combination with deltamethrin or cyper-methrin, 
      # we assume it has the same effect as using deltamethrin or cypermethrin alone.
      use.DMCMcht.t_local <- F # Indicator variable for DM or CM use
      if(SV_local$use.DMcht[t_local, cage.j]!=0 | SV_local$use.CMcht[t_local, cage.j]!=0){use.DMCMcht.t_local <- T} # Used the same day
      if(t_local>1){if(SV_local$use.DMcht[t_local-1, cage.j]!=0 | SV_local$use.CMcht[t_local-1, cage.j]!=0) use.DMCMcht.t_local <- T} # Used the day before
      if(SV_local$ux.DMcht.AF[t_local, a.max.AF, cage.j] != 0){use.DMCMcht.t_local <- T}   # Used previously, but still in effect
      if(SV_local$ux.CMcht.AF[t_local, a.max.AF, cage.j] != 0){use.DMCMcht.t_local <- T}
      if(SV_local$use.AZcht[t_local, cage.j]==1 & SV_local$use.DMCMcht.t_local==F) {        # If cage.j was treated and DM or CM not used
        t1 <- min(t_local + Coef["delta.del.AZcht"], Ndays)       # Start time of effect     
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.AZcht"] / max(SV_local$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.PA){
            out$ux.AZcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- AZcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.AZcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- AZcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Emamectin
update_ux.EMcht2 <- function(t_local, a.max.CH = 60, a.max.PA = 60, a.max.AF = 80,...){
  out <- SV_local[c("ux.EMcht.CH", "ux.EMcht.PA", "ux.EMcht.AF")]
  if(sum(SV_local$use.EMcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    EMcht.eff <- trt.sample("EMcht") # Draw a posterior sample with cage-specific effects for each stage
    EMcht.eff <- log(1 + exp(EMcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV_local$use.EMcht[t_local, cage.j]==1) {     # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.EMcht"], Ndays)       # Start time of effect     
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.EMcht"] / max(SV_local$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.CH){
            out$ux.EMcht.CH[t_local.eff, a.min:a.max.CH, cage.j] <- EMcht.eff[cage.j, "CH"] # Effect on CH
          }
          if(a.min <= a.max.PA){
            out$ux.EMcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- EMcht.eff[cage.j, "PA"] # Effect on PA
          }
          if(a.min <= a.max.AF){
            out$ux.EMcht.AF[t_local.eff, a.min:a.max.AF, cage.j] <- EMcht.eff[cage.j, "Af"]  # Effect on A
          } 
        }
      }
    }
  }
  return(out)
}

# Diflubenzuron
update_ux.DBcht2 <- function(t_local, a.max.CH = 60, a.max.PA = 60,...){
  out <- SV_local[c("ux.DBcht.CH", "ux.DBcht.PA")]
  if(sum(SV_local$use.DBcht[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    DBcht.eff <- trt.sample("DBcht") # Draw a posterior sample with cage-specific effects for each stage
    DBcht.eff <- log(1 + exp(DBcht.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV_local$use.DBcht[t_local, cage.j]==1) {     # If cage.j was treated
        t1 <- min(t_local + Coef["delta.del.DBcht"], Ndays)       # Start time of effect     
        # End time of effect (Eqs. 19, 20, but avoid dividing on ST < 1):
        t2 <- min(t1 + round(Coef["delta.dur.DBcht"] / max(SV_local$ST[t_local], 1), 0) - 1, Ndays)  
        for(t_local.eff in t1:t2){
          a.min <- t_local.eff - t_local + 1        # Minimum age-index (= age + 1) of effect (Eq. 19)
          if(a.min <= a.max.CH){
            out$ux.DBcht.CH[t_local.eff, a.min:a.max.CH, cage.j] <- DBcht.eff[cage.j, "CH"] # Effect on CH
          }
          if(a.min <= a.max.PA){
            out$ux.DBcht.PA[t_local.eff, a.min:a.max.PA, cage.j] <- DBcht.eff[cage.j, "PA"] # Effect on PA
          }
        }
      }
    }
  }
  return(out)
}

## Calculate mortality in subsequent time-steps 
## caused by physical treatments applied in time step t_local

# Thermal
update_ux.therm2 <- function(t_local,...){
  out <- SV_local[c("ux.therm.CH", "ux.therm.PA", "ux.therm.AF")]
  if(sum(SV_local$use.therm[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    therm.eff <- trt.sample("therm") # Draw a posterior sample with cage-specific effects for each stage
    therm.eff <- log(1 + exp(therm.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV_local$use.therm[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.therm.CH[t_local, , cage.j] <- therm.eff[cage.j, "CH"] # Effect on CH
        out$ux.therm.PA[t_local, , cage.j] <- therm.eff[cage.j, "PA"] # Effect on PA
        out$ux.therm.AF[t_local, , cage.j] <- therm.eff[cage.j, "Af"] # Effect on Af
      }
    }
  }
  return(out)
}

# Freshwater
update_ux.freshw2 <- function(t_local,...){
  out <- SV_local[c("ux.freshw.CH", "ux.freshw.PA", "ux.freshw.AF")]
  if(sum(SV_local$use.freshw[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    freshw.eff <- trt.sample("freshw") # Draw a posterior sample with cage-specific effects for each stage
    freshw.eff <- log(1 + exp(freshw.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV_local$use.freshw[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.freshw.CH[t_local, , cage.j] <- freshw.eff[cage.j, "CH"] # Effect on CH
        out$ux.freshw.PA[t_local, , cage.j] <- freshw.eff[cage.j, "PA"] # Effect on PA
        out$ux.freshw.AF[t_local, , cage.j] <- freshw.eff[cage.j, "Af"] # Effect on Af
      }
    }
  }
  return(out)
}

# Mechanical
update_ux.mech2 <- function(t_local,...){
  out <- SV_local[c("ux.mech.CH", "ux.mech.PA", "ux.mech.AF")]
  if(sum(SV_local$use.mech[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    mech.eff <- trt.sample("mech") # Draw a posterior sample with cage-specific effects for each stage
    mech.eff <- log(1 + exp(mech.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){         # Loop on cages
      if(SV_local$use.mech[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.mech.CH[t_local, , cage.j] <- mech.eff[cage.j, "CH"] # Effect on CH
        out$ux.mech.PA[t_local, , cage.j] <- mech.eff[cage.j, "PA"] # Effect on PA
        out$ux.mech.AF[t_local, , cage.j] <- mech.eff[cage.j, "Af"] # Effect on Af
      }
    }
  }
  return(out)
}

# Fixed treatment mortality = M.trt
update_ux.fx2 <- function(t_local, M.trt = 0.5, ...){
  out <- SV_local[c("ux.fx.CH", "ux.fx.PA", "ux.fx.AF")]
  if(sum(SV_local$use.fx[t_local, ]) > 0){ # If treatment was applied in at least one cage at time t_local
    # Trt.eff <- log(M.trt / (1 - M.trt))
    # Trt.eff <- log(1 + exp(Trt.eff)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    Trt.eff <- (-log(1 - M.trt)) # Transformed coefficients (Eq. 22, Aldrin et al. 2017)
    for(cage.j in 1:Ncages){             # Loop on cages
      if(SV_local$use.fx[t_local, cage.j]==1) {     # If cage.j was treated
        out$ux.fx.CH[t_local, , cage.j] <- Trt.eff # Effect on CH
        out$ux.fx.PA[t_local, , cage.j] <- Trt.eff # Effect on PA
        out$ux.fx.AF[t_local, , cage.j] <- Trt.eff # Effect on Af
      }
    }
  }
  return(out)
}

# Mortality due to chemical and physical treatments (Eq. 21, Aldrin et al. 2017)
update_m.trt.CH2 <- function(t_local,...){
  out <- SV_local$m.trt.CH[t_local, , , drop = F]
  ux.sum <- SV_local$ux.DMcht.CH[t_local, , ] + 
    SV_local$ux.EMcht.CH[t_local, , ] +
    SV_local$ux.DBcht.CH[t_local, , ] +
    SV_local$ux.therm.CH[t_local, , ] +
    SV_local$ux.freshw.CH[t_local, , ] +
    SV_local$ux.mech.CH[t_local, , ] +
    SV_local$ux.fx.CH[t_local, , ]
  m.sum <- 1 - exp(-ux.sum)
  out[1, , ] <- m.sum  
  return(out)
}

update_m.trt.PA2 <- function(t_local,...){
  out <- SV_local$m.trt.PA[t_local, , , drop = F]
  ux.sum <- SV_local$ux.HPcht.PA[t_local, , ] +
    SV_local$ux.DMcht.PA[t_local, , ] + 
    SV_local$ux.AZcht.PA[t_local, , ] +
    SV_local$ux.EMcht.PA[t_local, , ] +
    SV_local$ux.DBcht.PA[t_local, , ] +
    SV_local$ux.therm.PA[t_local, , ] +
    SV_local$ux.freshw.PA[t_local, , ] +
    SV_local$ux.mech.PA[t_local, , ] +
    SV_local$ux.fx.PA[t_local, , ]
  m.sum <- 1 - exp(-ux.sum)
  out[1, , ] <- m.sum  
  return(out)
}

update_m.trt.AF2 <- function(t_local,...){
  out <- SV_local$m.trt.AF[t_local, , , drop = F]
  ux.sum <- 
    SV_local$ux.HPcht.AF[t_local, , ] +
    SV_local$ux.DMcht.AF[t_local, , ] + 
    SV_local$ux.AZcht.AF[t_local, , ] +
    SV_local$ux.EMcht.AF[t_local, , ] +
    SV_local$ux.therm.AF[t_local, , ] +
    SV_local$ux.freshw.AF[t_local, , ] +
    SV_local$ux.mech.AF[t_local, , ] +
    SV_local$ux.fx.AF[t_local, , ]
  m.sum <- 1 - exp(-ux.sum)
  out[1, , ] <- m.sum  
  return(out)
}

## Survival rates (Eq. 12, Aldrin et al. 2017)
update_s.CH2 <- function(t_local){
  out <- (1 - SV_local$m.nat.CH[t_local,,]) * (1 - SV_local$m.trt.CH[t_local,,])
  return(out)
}
update_s.PA2 <- function(t_local){
  out <- (1 - SV_local$m.nat.PA[t_local,,]) * (1 - SV_local$m.clf.PA[t_local,,]) * (1 - SV_local$m.trt.PA[t_local,,])
  return(out)
}
update_s.AF2 <- function(t_local){
  out <- (1 - SV_local$m.nat.AF[t_local,,]) * (1 - SV_local$m.clf.AF[t_local,,]) * (1 - SV_local$m.trt.AF[t_local,,])
  return(out)}

## Infection rate
update_d.CO2 <- function(t_local, a.max.CO = 60,...){
  out <- SV_local$d.CO[t_local, , , drop = F]
  if(t_local < Ndays){
    # Cage-specific expected rate at logit scale (Eq. 28, Aldrin et al. 2017), Eq. 14 Aldrin & Huseby 2019
    eta.CO.t_local <- RE$delta0.CO.f + log(SV_local$N.SAL[t_local+1,]) + Coef["delta1.CO"] * (log(SV_local$W.SAL[t_local+1,]) - 0.55) +
      Coef["delta2.CO"] * (SV_local$ST[t_local+1] - 9) + Coef["delta3.CO"] * (SV_local$ST[t_local+1] - 9)^2
    
    # Cage-specific expected rate at proportion scale (Eq. 27, Aldrin et al. 2017), Eq. 13 Aldrin & Huseby 2019
    d.CO.t_local <- exp(eta.CO.t_local) / (1 + sum(exp(eta.CO.t_local)))
    
    # Age-specific values
    out[,,] <- rep(d.CO.t_local, each=a.max.CO)
  } else out[,,] <- 0
  return(out)  
}


## Reproduction

#  Reproduction rate for internal recruitment at time t_local (Eq. 31, Aldrin et al. 2017)
update_r.AF2 <- function(t_local, a.max.AF = 80, ...){
  
  # Median hatching time f(Eq. 32, Aldrin et al. 2017)
  d.m.E.t_local <- Coef["delta.Em10"]*(10/SV_local$ST[t_local])^Coef["delta.Rp"]
  
  out <- SV_local$r.AF[t_local, , , drop = F]
  
  r.raw <- Coef["beta0.R"] *          # N eggs first extrusion
    (1:a.max.AF)^Coef["beta1.R"] *    # N eggs different stage-ages of AF
    ((d.m.E.t_local + 1)^-1)                # Extrusion rate
  
  allee <-                          # Allee effect (differs between cages)
    (1 - exp(-Coef["gamma.r"] * sum(SV_local$N.AF[t_local, , ]) / SV_local$N.SAL[t_local, ]))
  
  for(cage.j in 1:Ncages){
    out[1, , cage.j] <- r.raw * allee[cage.j]
  }
  return(out)
}

## External recruitment (section 2.2.7, Aldrin et al. 2017)
update_r.Ext2 <- function(t_local, ...){
  
  # Median hatching time f(Eq. 32, Aldrin et al. 2017)
  d.m.E.t <- Coef["delta.Em10"]*(10/SV$ST[t])^Coef["delta.Rp"]
  
  # Allee effect
  allee <-                          
    (1 - exp(-Coef["gamma.r"] * SV$A.Ext[t]))
  
  out <- Coef["beta0.R"] *
    11^Coef["beta1.R"] *
    (1 / (d.m.E.t + 1)) *
    allee
  
  return(out)
}
