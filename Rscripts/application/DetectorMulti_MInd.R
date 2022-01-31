# Detector MInd
# note, to lift constraints on manifest variable used as marker variable (i.e. if in first position for any latent variables), 
# the order of the variables is changed.

revModLine = function(x){ # reverses the order of variables for a latent variable equation in lavaan syntax
  stopifnot("Equation is not a latent variable equation containing =~" = grepl("\\=\\~", x))
  mainComponents = unlist(strsplit(x, "\\=~"))
  rhs = unlist(strsplit(mainComponents[2], "\\+"))
  rhs = gsub("[[:blank:]]+", "", rhs)
  paste(mainComponents[1], "=~", paste(rev(rhs), collapse = " + "))
}

detectMulti_MInd = function(varnames, base_model, data, alpha = 0.05, group.constraints = c("intercepts","loadings")){
  # fully constrained baseline model
  base_fit = cfa(base_model, data = data, group = "grp", 
                 group.equal = group.constraints)

  # get all parameters relating to latent variables in base_model
  simple_fit = cfa(base_model, data = data, meanstructure = T)
  simple_partab = data.frame(par = paste(simple_fit@ParTable$lhs, simple_fit@ParTable$op, simple_fit@ParTable$rhs),
                             var = simple_fit@ParTable$rhs)
  simple_partab = simple_partab[simple_fit@ParTable$op == "=~",]
  
  #var-wise lifting of constraints
  rej = numeric(length(varnames))
  for(j in 1:length(varnames)){
    # select params for which constraints will be lifted
    lift.constraints = NULL
    if("loadings" %in% group.constraints){
      lift.constraints = simple_partab$par[simple_partab$var == varnames[j]]
    } 
    if("intercepts" %in% group.constraints){
      lift.constraints = c(lift.constraints, 
                           paste(varnames[j], "~1"))
    }  
    
    # fit model with lifted constraints
    if(grepl(paste0("\\=\\~[[:blank:]]*", varnames[j]), base_model)){ # if variable j is a marker variable in any factor equation =~
      # split model into lines
      base_model_split = unlist(strsplit(base_model, "\n"))
      # reverse order for those lines that have variable j in marker variable position
      base_model_split = sapply(base_model_split, function(line) ifelse(grepl("\\=\\~", line) & grepl(varnames[j], line),
                                                                        revModLine(line),
                                                                        line))
      base_model_j = paste0(base_model_split, collapse = "\n")
      fit = cfa(base_model_j, data = data, group = "grp", 
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
       noninvariant_bonferroni = varnames[which(rej < alpha/length(varnames))],
       alpha = alpha)
}
