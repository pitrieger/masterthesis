# Detector Byrne & van de Vijer (2010)

detect_ByrneVandeVijer = function(varnames, data, CFI.delta = 0.01, group.constraints = c("intercepts","loadings")){
  # baseline model with all vars
  base_model = paste("eta =~", paste(varnames, collapse = " + "))
  base_fit = cfa(base_model, data = data, group = "grp", group.equal = group.constraints)
  CFI.base = as.numeric(fitmeasures(base_fit, fit.measures = "cfi"))
  
  # var-leave-one-out models
  models = sapply(varnames, function(j) paste("eta =~", paste(varnames[varnames != j], collapse = " + ")))
  CFI = numeric(length(models))
  for(j in 1:length(models)){
    fit = cfa(models[j], data = data, group = "grp", group.equal = group.constraints)
    CFI[j] = as.numeric(fitmeasures(fit, fit.measures = "cfi"))
  }
  
  list(varnames = varnames,
       noninvariant = varnames[which(CFI > CFI.delta + CFI.base)], # interpretation: deletion improves goodness-of-fit
       CFIs = CFI,
       CFI.base = CFI.base)
}