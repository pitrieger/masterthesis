library(tidyverse)
library(parallel)
library(lavaan)
library(here)
source(here("Rscripts/simulation", "Simulator_PokropekEtAl.R"))
source(here("Rscripts/simulation", "Detector_Rieger.R"))
source(here("Rscripts/simulation", "Detector_ByrneVandeVijer.R"))
source(here("Rscripts/simulation", "Detector_CheungRensvold.R"))
source(here("Rscripts/simulation", "Detector_MInd.R"))
source(here("Rscripts/simulation", "Detector_Janssens.R"))

# Parallelization setup
RNGkind("L'Ecuyer-CMRG") # multi-core compatible RNG
options(mc.cores = max(detectCores() - 1, 2))
set.seed(9)

# Simulation parameters
sim_param = expand.grid(nsim = 100,
                        n = c(200, 400, 800),
                        p = c(3, 4, 5),
                        g = c(2),
                        h = c(0.5),
                        k = c(1, 2), 
                        interceptbias = c(0, 0.2),
                        loadingbias = c(0, 0.2))
sim_param = sim_param[sim_param$k<sim_param$p,]
sim_param_df = sim_param[sim_param$interceptbias != 0 | sim_param$loadingbias != 0,]
sim_param = split(sim_param_df, 1:nrow(sim_param_df))
sim_out = lapply(sim_param, function(x) set_names(list(x), "sim_param"))

# get entries of confusion matrix
get_confusion = function(pred, true, varnames){
  c("TP" = sum(pred %in% true), # true positive detections
    "TN" = sum(!varnames[!varnames %in% pred] %in% true), # true negative detections = non-detected varnames that are not in noninvariant_true
    "FP" = sum(pred %in% varnames[!varnames %in% true]), # false positive 
    "FN" = sum(varnames[!varnames %in% pred] %in% true)) # false negative
}

# simulation function
run_sim = function(x) {
  pars = x$sim_param
  varnames = paste0("y", 1:pars$p)
  
  D = sim_PMI(n = pars$n,
              g = pars$g,
              p = pars$p,
              h = pars$h,
              k = pars$k)
  # constants
    # varnames of noninvariant items
    noninvariant_true = ifelse(length(D$p_affected)>0, paste0("y", D$p_affected), character())
    # true number of noninvariant items = TP + FN = sensitivity denominator
    
  # run detectors
    # Janssens
    ni_J = tryCatch({
      get_confusion(detect_Janssens(varnames, D$sim_dat)$noninvariant, 
                    noninvariant_true, varnames)
    }, error = function(cond){return(rep(NA, 4))},
    silent = T)
    
    # MInd
    ni_M = tryCatch({
      get_confusion(detect_MInd(varnames, D$sim_dat)$noninvariant, 
                   noninvariant_true, varnames)
    }, error = function(cond){return(rep(NA, 4))},
    silent = T)
    
    # Cheung & Rensveld (1999):
    ni_C = tryCatch({
      get_confusion(detect_CheungRensvold(varnames, D$sim_dat)$noninvariant, 
                   noninvariant_true, varnames)    
      }, error = function(cond){return(rep(NA, 4))},
      silent = T)
    
    # Byrne & van de Vijer (2010):
    ni_B = tryCatch({
      get_confusion(detect_ByrneVandeVijer(varnames, D$sim_dat)$noninvariant, 
                   noninvariant_true, varnames)
    }, error = function(cond){return(rep(NA, 4))},
    silent = T)
    
    # Own Idea v1: simple
    ni_R1 = tryCatch({
      get_confusion(detect_Rieger_v1(varnames, D$sim_dat)$noninvariant, 
                   noninvariant_true, varnames)
    }, error = function(cond){return(rep(NA, 4))},
    silent = T)
    
    # Own Idea v2: stepwise
    ni_R2 = tryCatch({
      get_confusion(detect_Rieger_v2(varnames, D$sim_dat)$noninvariant, 
                   noninvariant_true, varnames)
    }, error = function(cond){return(rep(NA, 4))},
    silent = T)

    rbind(ni_J, ni_M, ni_C, ni_B, ni_R1, ni_R2)
}

# replicate within parallel
system.time(out <- mclapply(sim_out, function(x) replicate(x$sim_param$nsim, run_sim(x))))
save.image(here("data/Prelim_Sim_2021-10-26.RData"))


lapply(out, function(l) apply(l, 1:2, function(x) mean(x, na.rm = T)))


get_sensspec = function(x, alpha = 0.05, type = "sensitivity"){
  stopifnot(type %in% c("sensitivity", "specificity"))
  if(type == "sensitivity"){
    r = x[,"TP"]
    n = r + x[,"FN"]
  } else {
    r = x[,"TN"]
    n = r + x[,"FP"]
  } 
  # proportion
  est = r/n
  # Clopper-Pearson CIs
    # get F-values, but use 0 or 1 if est is 0 or 1 
    lowerF = qf(1-alpha/2, 2*(n-r+1), 2*r)
    lowerF = ifelse(is.nan(lowerF), round(est), lowerF)
    upperF = qf(1-alpha/2, 2*(r+1), 2*(n-r))
    upperF = ifelse(is.nan(upperF), round(est), upperF)
    #CIs
    lowerCI = r / (r + (n-r+1)*lowerF)
    lowerCI = ifelse(is.nan(lowerCI), round(est), lowerCI)
    upperCI = (r + 1) * upperF / ((n-r) + (r+1)*upperF)
  data.frame(method = rownames(x), est, lowerCI, upperCI)
}


ggplot(get_sensspec(x), aes(y = est, ymin = lowerCI, ymax = upperCI, x = method)) + 
  geom_pointrange()
ggplot(get_sensspec(x, type = "specificity"), aes(y = est, ymin = lowerCI, ymax = upperCI, x = method)) + 
  geom_pointrange()


























