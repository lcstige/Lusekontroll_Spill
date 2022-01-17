update_SV <- function(SV_local = SV, Nsteps = 1, do_treat, trt.type, Ncages) {
  ## Parameters

  # trt.type <- c("HPcht", "DMcht", "AZcht", "EMcht", "DBcht",
  #               "therm", "freshw", "mech", 
  #               "fx")[6]  # Note: "fx" is not implemented in the rest of this script
  strategy <- c("cage-wise","farm-wise")[2]
  
  # Treatment threshold
  # thr <- 0.1
  Ncages = 4
  
  
  # Update functions --------------------------------------------------------
  
  # Functions to count CH lice
  update_Y.CH2 <- function(t_local, ...){
    ## As lice are counted at the start of the day, before mortality and development,
    ## we use the lice numbers from the day before (these are updated later)
    # (section 2.2.1 of Aldrin et al. 2017)
    out <- SV_local$Y.CH[t_local, , drop = F]
    if(t_local == 1) {
      N.CH.t_local <- SV_local$N.CH[t_local, , , drop = F]
    } else {
      N.CH.t_local <- SV_local$N.CH[t_local-1, , , drop = F]
    }
    
    for(cage.j in 1:Ncages){
      if(SV_local$n.SAL[t_local, cage.j]>0){
        # Probability of sampling chalimus-stage lice (Eqs. 45-47, Aldrin et al. 2017)
        p.CHcount <- invlogit(RE$beta0.CHcount.f + Coef["beta1.CHcount"] * SV_local$W.SAL[t_local, cage.j])
        
        # Expected number of chalimums-stage lice (Eq. 42, Aldrin et al. 2017)
        E.CH <- SV_local$n.SAL[t_local, cage.j] * p.CHcount * sum(N.CH.t_local[1, , cage.j])/sum(SV_local$N.SAL[t_local, cage.j])
        
        # Observed counts (cf. Eq. 41, Aldrin et al. 2017)
        out[cage.j] <- ifelse(E.CH>0, rnbinom(n=1, mu = E.CH, size = Coef["rho.CH"] * SV_local$n.SAL[t_local, cage.j]), 0)
      }
    }
    return(out)
  }
  update_Y.AF2 <- function(t_local, ...){
    out <- SV_local$Y.AF[t_local,  , drop = F]
    if(t_local == 1) {
      N.AF.t_local <- SV_local$N.AF[t_local, , , drop = F]
    } else {
      N.AF.t_local <- SV_local$N.AF[t_local-1, , , drop = F]
    }
    for(cage.j in 1:Ncages){
      if(SV_local$n.SAL[t_local, cage.j]>0){
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
  
  
  
  # Update loop -------------------------------------------------------------
  
  
  #Nsteps <- floor(Ndays/stepsize)
  for(stp in 1:Nsteps){ 
    # Run the daily lice dynamics until next counting
    for(i in 1:stepsize){
      
      t <- t + 1
      #  print(paste(i, "days since last count", t, "days total"))
      # Update state variables
      
      # true lice abundances
      SV_local$N.R[t, ]    <- update_N.R2(t)
      SV_local$N.CO[t, ]   <- update_N.CO2(t)
      SV_local$N.CH[t, , ] <- update_N.CH2(t)
      SV_local$N.PA[t, , ] <- update_N.PA2(t)
      SV_local$N.AF[t, , ] <- update_N.AF2(t)
      SV_local$N.AM[t, , ] <- update_N.AM2(t)
      
      # cleaner fish abundances
      SV_local$N.lump[t,] <- update_N.lump2(t)
      SV_local$N.wrasse[t,] <- update_N.wrasse2(t)   
      
      # cleaner fish-caused mortality of salmon lice
      SV_local$m.clf.PA[t, , ] <- update_m.clf.PA2(t)
      SV_local$m.clf.AF[t, , ] <- update_m.clf.AF2(t) 
      
      # calculate treatment mortalities in subsequent time-steps 
      # if treatments are applied at time t;
      
      # chemical treatments
      # (hydrogen peroxide, deltamethrin, azamethiphos, diflubenzuron)
      SV_local[c("ux.HPcht.PA", "ux.HPcht.AF")] <- update_ux.HPcht2(t)
      SV_local[c("ux.DMcht.CH", "ux.DMcht.PA", "ux.DMcht.AF")] <- update_ux.DMcht2(t)
      SV_local[c("ux.AZcht.PA", "ux.AZcht.AF")] <- update_ux.AZcht2(t)
      SV_local[c("ux.EMcht.CH", "ux.EMcht.PA", "ux.EMcht.AF")] <- update_ux.EMcht2(t)
      SV_local[c("ux.DBcht.CH", "ux.DBcht.PA")] <- update_ux.DBcht2(t)
      
      # physical treatments (thermal, freshwater, mechanical)
      SV_local[c("ux.therm.CH", "ux.therm.PA", "ux.therm.AF")] <- update_ux.therm2(t)
      SV_local[c("ux.freshw.CH", "ux.freshw.PA", "ux.freshw.AF")] <- update_ux.freshw2(t)
      SV_local[c("ux.mech.CH", "ux.mech.PA", "ux.mech.AF")] <- update_ux.mech2(t)
      
      # treatment mortality at time t
      SV_local$m.trt.CH[t, , ] <- update_m.trt.CH2(t) 
      SV_local$m.trt.PA[t, , ] <- update_m.trt.PA2(t) 
      SV_local$m.trt.AF[t, , ] <- update_m.trt.AF2(t) 
      
      ## survival rates
      SV_local$s.CH[t, , ] <- update_s.CH2(t)
      SV_local$s.PA[t, , ] <- update_s.PA2(t)
      SV_local$s.AF[t, , ] <- update_s.AF2(t)
      
      ## infection rate
      SV_local$d.CO[t, ,] <- update_d.CO2(t)
      
      ##  reproduction rate for internal recruitment 
      SV_local$r.AF[t, , ] <- update_r.AF2(t)
      
      ## reproduction rate for external recruitment
      SV_local$r.Ext[t] <- update_r.Ext2(t) 
    }
    
    # Count lice
    SV_local$n.SAL[t,] <- ncount
    
    # lice counts
    SV_local$Y.CH[t,] <- update_Y.CH2(t)
    SV_local$Y.AF[t,] <- update_Y.AF2(t)
    SV_local$Y.OM[t,] <- update_Y.OM2(t)
    sum(SV_local$Y.OM)
    # A possible second count of lice
    SV_local$n2.SAL[t,] <- n2count
    
    # lice counts
    SV_local$Y2.CH[t,] <- update_Y2.CH2(t)
    SV_local$Y2.AF[t,] <- update_Y2.AF2(t)
    SV_local$Y2.OM[t,] <- update_Y2.OM2(t)
    
    # Decide whether to treat next day
    
    if(t<Ndays & strategy == "cage-wise") {
      if(do_treat == 1) SV_local[[paste0("use.",trt.type)]][t+1, ] <- 1
    }
    # Farm-wise treatment:
    if(t<Ndays & strategy == "farm-wise") {
      if(do_treat == 1) SV_local[[paste0("use.",trt.type)]][t+1, ] <- 1
    }
      
    # natonull <- function(x){x[is.na(x)]<-0; return(x)}
    # # Cage-wise treatment:
    # if(t<Ndays & strategy == "cage-wise") SV_local[[paste0("use.",trt.type)]][t+1, ] <- 
    #   natonull(1*((SV_local$Y.AF[t, ] + SV_local$Y.OM[t, ])/SV_local$n.SAL[t, ] > thr))
    # # Farm-wise treatment:
    # if(t<Ndays & strategy == "farm-wise") SV_local[[paste0("use.",trt.type)]][t+1, ] <- 
    #   natonull(1*(sum(SV_local$Y.AF[t, ] + SV_local$Y.OM[t, ])/sum(SV_local$n.SAL[t, ]) > thr))
    ## Return t to global
    
    
  }
  list(t = t, SV = SV_local) %>% 
    return()
}

# SV %>% names()
# SV$Y.AF %>% head


logoffset <- .01
log_transform <- function(lus, cage_no) {
  log10(logoffset + lus[, cage_no]/SV$n.SAL[, cage_no])
}
shift_treatment <- function(treatment, cage_no) { # To fit the treatments to the filtered days
  c(treatment[2:nrow(treatment), cage_no], 0)
}
summarise_data <- function() {
  
  map(1:Ncages, function(x) {
    data.frame(Y.CH  = log_transform(lus = SV$Y.CH, cage_no = x)) %>% 
      mutate(Y.OM  = log_transform(lus = SV$Y.OM, cage_no = x)) %>% 
      mutate(Y.AF  = log_transform(lus = SV$Y.AF, cage_no = x)) %>% 
      mutate(Y2.CH = log_transform(lus = SV$Y2.CH, cage_no = x)) %>% 
      mutate(Y2.OM = log_transform(lus = SV$Y2.OM, cage_no = x)) %>% 
      mutate(Y2.AF = log_transform(lus = SV$Y2.AF, cage_no = x)) %>% 
      # mutate(N.CH  = SV$N.CH[,,x] + logoffset) %>% 
      # mutate(N.OM  = SV$N.OM[,,x] + logoffset) %>% 
      # mutate(N.AF  = SV$N.AF[,,x] + logoffset) #%>% 
      mutate(cage  = as.character(x)) %>% 
      mutate(day   = 1:Ndays) %>% 
      mutate(use.HPcht  = shift_treatment(SV$use.HPcht, x)) %>% 
      mutate(use.DMcht  = shift_treatment(SV$use.DMcht, x)) %>% 
      mutate(use.AZcht  = shift_treatment(SV$use.AZcht, x)) %>% 
      mutate(use.DBcht  = shift_treatment(SV$use.DBcht, x)) %>% 
      mutate(use.therm  = shift_treatment(SV$use.therm, x)) %>% 
      mutate(use.freshw = shift_treatment(SV$use.freshw, x)) %>% 
      mutate(use.mech   = shift_treatment(SV$use.mech, x)) %>% 
      mutate(use.fx     = shift_treatment(SV$use.fx, x)) %>% 
      filter(((day-1)%%7 == 0) | (day == 1))
  }
  ) %>% 
    do.call(rbind, .) %>% 
    mutate(treatment = ((use.HPcht + use.DMcht + use.AZcht + use.DBcht + use.therm + use.freshw + use.mech + use.fx != 0) * day)) %>% 
    dplyr::select(!starts_with("use"))
  
  
}

summarise_data()

