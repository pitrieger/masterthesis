# Detector Rieger

get_pval = function(X, hat.eta, grp){
  stopifnot(nrow(X) == length(hat.eta) & length(hat.eta) == length(grp))
  g = length(unique(grp))
  p = ncol(X)
  p.vals = numeric(p)
  for(j in 1:p){
    # get residuals
    R = resid(lm(X[,j] ~ hat.eta))
    # group means
    aov.p = summary(aov(R ~ as.factor(grp)))[[1]]$`Pr(>F)`[1]
    # group correlations
    cor.p = numeric(g)
    for(i in 1:g){
      cor.p[i] = cor.test(R[grp == i], hat.eta[grp == i])$p.value # equivalent to coefficient test for lm
    }
    # aggregate p.vals
    p.vals[j] = (g+1)*min(aov.p, cor.p)
  }
  p.vals
}

get_pval_metric = function(X, hat.eta, grp){
  stopifnot(nrow(X) == length(hat.eta) & length(hat.eta) == length(grp))
  g = length(unique(grp))
  p = ncol(X)
  p.vals = numeric(p)
  for(j in 1:p){
    # get residuals
    R = resid(lm(X[,j] ~ hat.eta))
    # equal correlations => nonsignificant interaction of eta & grp
    fit = lm(R ~ hat.eta*as.factor(grp))
    cor.p = summary(fit)$coefficients[grepl("hat\\.eta\\:", names(coef(fit))),4] # p-vals of rows containing interaction

    # aggregate p.vals
    p.vals[j] = (g)*min(cor.p)
  }
  p.vals
}

detect_Rieger = function(varnames, data, alpha = 0.05, detection.type = "both"){
  model = paste("eta =~", paste(varnames, collapse = " + "))
  fit = cfa(model, data = data)
  hat.eta = predict(fit)[,1]
  X = fit@Data@X[[1]]
  grp = data$grp
  if(detection.type == "both"){
    p.vals = get_pval(X, hat.eta, grp)
  } else if(detection.type == "metric"){
    p.vals = get_pval_metric(X, hat.eta, grp)
  }
  list(varnames = varnames,
       noninvariant = varnames[which(p.vals < alpha/(length(p.vals) - order(p.vals) + 1))], # Bonferroni-holm
       p.vals = p.vals)
}

detect_Rieger_step = function(varnames, data, alpha = 0.05, detection.type = "both"){
  p = length(varnames)
  model = paste("eta =~", paste(varnames, collapse = " + "))
  fit = cfa(model, data = data)
  
  detected = detect_Rieger(varnames, data, alpha, detection.type)
  p.val.min = p*min(detected$p.vals) #Bonferroni-Holm == Bonferroni
  
  varnames_it = varnames
  
  while(p.val.min < alpha){
    varnames_it = detected$varnames[-which.min(detected$p.vals)]
    if(length(varnames_it) <= 1){
      break
    }
    detected = detect_Rieger(varnames_it, data, alpha)
    p.val.min = length(varnames_it)*min(detected$p.vals) # Bonferroni-Holm
  }
  list(varnames = varnames,
       noninvariant = varnames[!varnames %in% varnames_it])
}
