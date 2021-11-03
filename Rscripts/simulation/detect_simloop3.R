library(tidyverse)
library(parallel)
library(lavaan)
library(here)
library(RColorBrewer)
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
sim_param_df = expand.grid(nsim = 100,
                        n = c(200, 500, 1000),
                        p = c(3, 4, 5, 6),
                        g = c(2, 4, 8, 16),
                        h = c(0.25, 0.5),
                        k = c(1, 2, 3), 
                        itembias = c(0, 0.2))
#                        interceptbias = c(0, 0.2, 0.4),
#                       loadingbias = c(0, 0.2, 0.4))
sim_param_df = sim_param_df[2*sim_param_df$k<sim_param_df$p & 
                            (sim_param_df$h * sim_param_df$g) %% 1 == 0,]

sim_param = split(sim_param_df, 1:nrow(sim_param_df))
sim_out = lapply(sim_param, function(x) set_names(list(x), "sim_param"))

# get entries of confusion matrix
get_confusion = function(pred, true, varnames){
  c("TP" = sum(pred %in% true), # true positive detections
    "TN" = sum(!varnames[!varnames %in% pred] %in% true), # true negative detections = non-detected varnames that are not in noninvariant_true
    "FP" = sum(pred %in% varnames[!varnames %in% true]), # false positive 
    "FN" = sum(varnames[!varnames %in% pred] %in% true)) # false negative
}

# get sensitivity and specificity with CIs from confusion matrix output
get_sensspec = function(x, alpha = 0.05, type = "sensitivity"){
  stopifnot(type %in% c("sensitivity", "specificity"))
  if(type == "sensitivity"){
    r = as.numeric(as.data.frame(x)[,"TP"])
    n = r + as.numeric(as.data.frame(x)[,"FN"])
  } else {
    r = as.numeric(as.data.frame(x)[,"TN"])
    n = r + as.numeric(as.data.frame(x)[,"FP"])
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

# simulation function
run_sim = function(x) {
  # simulation parameters
  pars = x$sim_param
  pars_list = list()
  for(i in 1:ncol(pars)) pars_list[[i]] = pars[,i]
  names(pars_list) = colnames(pars)
  
  D = do.call(sim_PMI, pars_list)

  # constants
    varnames = paste0("y", 1:pars$p)
    # varnames of noninvariant items
    noninvariant_true = ifelse(length(D$p_affected)>0, paste0("y", D$p_affected), character())
    
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
#save.image(here("data/Prelim_Sim_2021-11-02.1RData"))

# compute sensitivity
out_sensitivity = do.call("rbind", 
                          lapply(lapply(out, function(l) apply(l, 1:2, function(x) sum(x, na.rm = T))), #sum confusion matrices over replications 
                                 function(x) get_sensspec(x))) # compute sensitivity
out_sensitivity = cbind(out_sensitivity, sim_param_df[rep(1:nrow(sim_param_df), each = 6),]) # add simulation parameters

# compute specificity
out_specificity = do.call("rbind", 
                          lapply(lapply(out, function(l) apply(l, 1:2, function(x) sum(x, na.rm = T))), #sum confusion matrices over replications 
                                 function(x) get_sensspec(x, type = "specificity"))) # compute specificity
out_specificity = cbind(out_specificity, sim_param_df[rep(1:nrow(sim_param_df), each = 6),]) # add simulation parameters


## Generate plots
  # relabel method names
  out_sensitivity$method = factor(out_sensitivity$method, 
                                  levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                                  labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  out_specificity$method = factor(out_specificity$method, 
                                  levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                                  labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  # theme
  sensspec_layer = list(geom_pointrange(),
                        geom_line(show.legend = F),
                        facet_grid(rows = vars(k), cols = vars(p),
                                   labeller = function(x) label_both(x, sep = " = ")),
                        scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)),
                        scale_color_brewer(type = "qual", palette = 2),
                        theme_bw(),
                        theme(strip.background = element_rect(color = "black", fill = "white")))

  ## SENSITIVITY
    # Plots for Loadingbias
    df_p1.1 = out_sensitivity %>% filter(loadingbias == 0.2 &
                                         interceptbias == 0)
    ggplot(df_p1.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = n, color = method, shape = method)) + 
      sensspec_layer + 
      labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")
    
    # Plots for Interceptbias
    df_p2.1 = out_sensitivity %>% filter(loadingbias == 0 &
                                       interceptbias == 0.2)
    ggplot(df_p2.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = n, color = method, shape = method)) + 
      sensspec_layer + 
      labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")
    
    # Plots for Intercept & Loadingbias
    df_p3.1 = out_sensitivity %>% filter(loadingbias == 0.2 &
                                         interceptbias == 0.2)
    ggplot(df_p3.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = n, color = method, shape = method)) + 
      sensspec_layer + 
      labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")

  ## SPECIFICITY
    # Plots for Loadingbias
    df_p1.2 = out_specificity %>% filter(loadingbias == 0.2 &
                                         interceptbias == 0)
    ggplot(df_p1.2, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = n, color = method, shape = method)) + 
      sensspec_layer +
      labs(y = "Specificity", x = "n", color = "Method", shape = "Method")
    
    # Plots for Interceptbias
    df_p2.2 = out_specificity %>% filter(loadingbias == 0 &
                                         interceptbias == 0.2)
    ggplot(df_p2.2, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = n, color = method, shape = method)) + 
      sensspec_layer +
      labs(y = "Specificity", x = "n", color = "Method", shape = "Method")
    
    
    # Plots for Intercept & Loadingbias
    df_p3.2 = out_specificity %>% filter(loadingbias == 0.2 &
                                         interceptbias == 0.2)
    ggplot(df_p3.2, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = n, color = method, shape = method)) + 
      sensspec_layer +
      labs(y = "Specificity", x = "n", color = "Method", shape = "Method")
