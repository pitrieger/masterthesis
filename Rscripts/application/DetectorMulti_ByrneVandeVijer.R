# Detector Byrne & van de Vijer (2010)
require(tidyverse)

detectMulti_ByrneVandeVijer = function(varnames, base_model, data, CFI.delta = 0.01, group.constraints = c("intercepts","loadings")){
  # baseline model with all vars
  base_fit = cfa(base_model, data = data, group = "grp", group.equal = group.constraints)
  CFI.base = as.numeric(fitmeasures(base_fit, fit.measures = "cfi"))
  
  # var-leave-one-out models
  models = sapply(varnames, function(varremove) gsub(paste0(" ", varremove, " |"), "", base_model) %>%
                    gsub(paste0(" ", varremove, "\n"), "\n", .)) # remove indicator everywhere
  # clean model syntax after removal
  models = models %>% gsub("\\~[[:blank:]]*\\+", "~", .) %>% # replaced var at beginning
    gsub("\\+[[:blank:]]*\\+", "+", .) %>% # replaced var in middle
    gsub("\\+[[:blank:]]*\n|\\+[[:blank:]]*$", "\n", .) # replaced var at end of line or very end
  
  CFI = numeric(length(models))
  for(j in 1:length(models)){
    fit = cfa(models[j], data = data, group = "grp", group.equal = group.constraints)
    CFI.temp = try(as.numeric(fitmeasures(fit, fit.measures = "cfi")))
    if(inherits(CFI.temp, "try-error")){
      CFI[j] = CFI.base # if model couldn't be fitted without item j, identify j as invariant
    } else {
      CFI[j] = CFI.temp
    }
  }
  
  CFI == CFI.base
  
  list(varnames = varnames,
       noninvariant = varnames[which(CFI > CFI.delta + CFI.base)], # interpretation: deletion improves goodness-of-fit
       nofit = varnames[CFI == CFI.base], # couldn't fit leave-one-out model without this variable for whatever reason
       CFIs = CFI,
       CFI.base = CFI.base)
}
