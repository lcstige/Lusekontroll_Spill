library(tidyverse)
library(R6)
source('global.R')


SeaLiceWorld <- R6Class("SeaLiceWorld", list(
  SV = read_SV(),
  do_treat = 1,
  add = function(x = 1) {
    self$SV <- update_SV(do_treat=do_treat, trt.type = "therm")
    invisible(self)
    },
    summarise_data = function(SV = self$SV) {
      logoffset <- .01
      log_transform <- function(lus, cage_no) {
        log10(logoffset + lus[, cage_no]/SV$n.SAL[, cage_no])
      }
      shift_treatment <- function(treatment, cage_no) { # To fit the treatments to the filtered days
        c(treatment[2:nrow(treatment), cage_no], 0)
      }
      
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
  )
)

do_treat = 1

SV <- read_SV()
summarise_data(SV) %>% dim
SV <-update_SV(SV, do_treat = 1, trt.type = 'HPcht')


sim_obj <- SeaLiceWorld$new()

sim_obj$summarise_data() %>% head

atta <- sim_obj$summarise_data()
atta %>% dim
sim_obj$add()


do_treat = TreatmentThreshold = 1

trt.type = "therm"

SeaLiceWorld$SV

dataz <- SeaLiceWorld$new()
dataz$SV %>% dim

dataz$summarised_data %>% dim

dataz$add()
dataz$summarised_data %>% dim

summarised_data = summarise_data()
summarised_data %>% dim


SV_T <- update_SV(do_treat=do_treat, trt.type = "therm") 
SV <- SV_T$SV 

rm(SV)
SV <- read_SV()
summarise_data(SV) %>% head

SV_T <- update_SV(do_treat=do_treat, trt.type = "therm") 
SV <- SV_T$SV 
