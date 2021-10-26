# Detector MInd
# note, to lift constraints on manifest variable 1, the order of the variables is changed.

detect_MInd = function(varnames, data, alpha = 0.05, group.constraints = c("intercepts","loadings")){
  # fully constrained baseline model
  base_model = paste("eta =~", paste(varnames, collapse = " + "))
  base_fit = cfa(base_model, data = data, group = "grp", 
                 group.equal = group.constraints)
  
  # var-wise lifting of constraints
  rej = numeric(length(varnames))
  for(j in 1:length(varnames)){
    # select params for which constraints will be lifted
    lift.constraints = NULL
    if("intercepts" %in% group.constraints){
      lift.constraints = paste("eta =~", varnames[j])
    } 
    if("loadings" %in% group.constraints){
      lift.constraints = c(lift.constraints, 
                           paste(varnames[j], "~1"))
    }  
    
    # fit model with lifted constraints
    if(j == 1){ # reverse order of model specification s.t. MV 1 isn't reference level with loading = 1
      fit = cfa(paste("eta =~", paste(rev(varnames), collapse = " + ")), data = data, group = "grp", 
                group.equal = group.constraints,
                group.partial = lift.constraints)
    } else {
      fit = cfa(base_model, data = data, group = "grp", 
                group.equal = group.constraints,
                group.partial = lift.constraints)  
    }
    
    # significance test
    cons.test = lavTestLRT(base_fit, fit)
    rej[j] = cons.test$`Pr(>Chisq)`[2]
  }

  list(varnames = varnames,
       noninvariant = varnames[which(rej < alpha)], # interpretation: lifting constraint improves goodness-of-fit
       alpha = alpha)
}
