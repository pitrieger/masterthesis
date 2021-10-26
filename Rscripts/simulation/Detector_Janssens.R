# Janssens 1994
# checking whether parameters are significant in some groups but not others
# this has multiple theoretical issues, see thesis

detect_Janssens = function(varnames, data, alpha = 0.05, group.constraints = c("intercepts","loadings")){
  base_model1 = paste("eta =~", paste(varnames, collapse = " + "))
  base_fit1 = cfa(base_model1, data = data, group = "grp")
  invisible(capture.output({base_fit_sum <- summary(base_fit1)$PE})) # captures list output from summary without print
  
  # check significance of parameters across groups
  term = paste(base_fit_sum$lhs, base_fit_sum$op, base_fit_sum$rhs)
  intercept_ids = which(term %in% paste(varnames, "~1 "))
  loading_ids = which(term %in% paste("eta =~", varnames))
  intercept_sig = matrix(base_fit_sum$pvalue[intercept_ids], 
                         nrow = length(varnames), ncol = length(unique(data$grp)))
  loading_sig = matrix(base_fit_sum$pvalue[loading_ids], 
                         nrow = length(varnames), ncol = length(unique(data$grp)))
  
  # change for last value via reverse s.t. item 1 is not reference item
  base_model2 = paste("eta =~", paste(rev(varnames), collapse = " + "))
  base_fit2 = cfa(base_model2, data = data, group = "grp")
  invisible(capture.output({base_fit_sum <- summary(base_fit2)$PE})) # captures list output from summary without print
  term = paste(base_fit_sum$lhs, base_fit_sum$op, base_fit_sum$rhs)
  intercept_sig[1,] = base_fit_sum$pvalue[which(term %in% paste(varnames[1], "~1 "))]
  loading_sig[1,] = base_fit_sum$pvalue[which(term %in% paste("eta =~", varnames[1]))]

  # difference across groups?
  intercept_dif = apply(intercept_sig, 1, function(x) length(unique(x < alpha))>1)
  loading_dif = apply(loading_sig, 1, function(x) length(unique(x < alpha))>1)
  
  noninvariant_id = rep(F, length(varnames)) 
  if("intercepts" %in% group.constraints){
    noninvariant_id = intercept_dif
  }
  if("loadings" %in% group.constraints){
    noninvariant_id[!noninvariant_id] = loading_dif[!noninvariant_id]
  }
  
  list(varnames = varnames,
       noninvariant = varnames[noninvariant_id],
       alpha = alpha)
}

        