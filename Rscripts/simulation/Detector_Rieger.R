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
      cor.p[i] = cor.test(R[grp == i], hat.eta[grp == i])$p.value
    }
    # aggregate p.vals
    p.vals[j] = 2*g*p*min(aov.p, cor.p)
  }
  p.vals
}

detect_Rieger_v1 = function(varnames, data, alpha = 0.05){
  model = paste("eta =~", paste(varnames, collapse = " + "))
  fit = cfa(model, data = data)
  hat.eta = predict(fit)[,1]
  X = fit@Data@X[[1]]
  grp = data$grp
  p.vals = get_pval(X, hat.eta, grp)
  list(varnames = varnames,
       noninvariant = varnames[which(p.vals < alpha)],
       p.vals = p.vals)
}

detect_Rieger_v2 = function(varnames, data, alpha = 0.05){
  model = paste("eta =~", paste(varnames, collapse = " + "))
  fit = cfa(model, data = data)
  
  detected = detect_Rieger_v1(varnames, data, alpha)
  p.val.min = min(detected$p.vals)
  varnames_it = varnames
  
  while(p.val.min < alpha){
    varnames_it = detected$varnames[-which.min(detected$p.vals)]
    detected = detect_Rieger_v1(varnames_it, data, alpha)
    p.val.min = min(detected$p.vals)
  }
  list(varnames = varnames,
       noninvariant = varnames[!varnames %in% varnames_it])
}
