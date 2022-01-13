run_sim = function(pars) {
  D = sim_PMI(n = pars$n, g = pars$g, p = pars$p, h = 0.5, k = pars$k)
  
  # constants
  # varnames of noninvariant items
  varnames = paste0("y", 1:pars$p)
  noninvariant_true = ifelse(length(D$p_affected)>0 & #only if there are any
                             any(c(pars$loadingbias, pars$interceptbias) > 0),
                             paste0("y", D$p_affected), character())
  
  # run detectors
  # MInd
  ni_M = tryCatch({
    get_confusion(detect_MInd(varnames, D$sim_dat)$noninvariant, 
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
    get_confusion(detect_Rieger(varnames, D$sim_dat)$noninvariant, 
                  noninvariant_true, varnames)
  }, error = function(cond){return(rep(NA, 4))},
  silent = T)
  
  # Own Idea v2: stepwise
  ni_R2 = tryCatch({
    get_confusion(detect_Rieger_step(varnames, D$sim_dat)$noninvariant, 
                  noninvariant_true, varnames)
  }, error = function(cond){return(rep(NA, 4))},
  silent = T)
  
  rbind(ni_M, ni_B, ni_R1, ni_R2)
}

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
  # get F-values, but use 0 or 1 if est is 0 or 1 to fix NaNs
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
