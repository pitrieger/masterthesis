# Detector Rieger
# Changes to V1: enable multiple latent variables & get model as input instead of varnames


# For after vacation:
# problems: order of deletion if multiple latent variables
#           when to stop -> if 1 indicator left => how to interpret


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

update_model = function(model, remove){
  # remove
  model.upd = gsub(remove, "", model)
  # fix
  model.upd = gsub("\\+[[:blank:]]*\n|\\+[[:blank:]]*$", "\n", model.upd) # removed var at end of LV formula
  model.upd = gsub("\\+[[:blank:]]*\\+", "\\+", model.upd) # removed middle
  model.upd = gsub("\\~[[:blank:]]*\\+", "\\~", model.upd) # removed start
  model.upd = gsub("\\~[[:blank:]]*\n|\\~[[:blank:]]*$", "\\~", model.upd) # no vars left
  model.upd
}

detect_Rieger1 = function(model, data, alpha = 0.05, ...){
  # fit user-supplied model
  fit = cfa(model, data = data, ...)
  
  # get lv-predictions
  hat.eta = predict(fit)
  
  # obtain actually used data
  hat.lambda = fit@Model@GLIST$lambda
  hat.lambda.str = (hat.lambda != 0)
  varnames = fit@Model@dimNames[[1]][[1]]
  lvnames = fit@Model@dimNames[[1]][[2]]
  X = fit@Data@X[[1]]
  grp = data$grp
  
  # for each latent variable, test invariance of residuals
  out = list()
  for(j in 1:ncol(hat.eta)){
    # manifest variables used for latent variable j
    select.vars = hat.lambda.str[,j]
    select.varnames = varnames[select.vars]
    
    # test invariance
    p.vals = get_pval(X[, select.vars], hat.eta[,j], grp)
    out[[j]] = list(varnames = select.varnames,
                    noninvariant = select.varnames[which(p.vals < alpha)],
                    p.vals = p.vals)
  } 
  names(out) = lvnames
  out
}

model = 'antiel =~ ow_ae1 + ow_ae2 + ow_ae3 + ow_ae4 + ow_ae5..'
model = 'antiel =~ ow_ae1 + ow_ae2 + ow_ae3 + ow_ae4 + ow_ae5..
 mistexp =~ ow_me1 + ow_me2 + ow_me3 + ow_me4
 nataff =~ ow_na1 + ow_na2 + ow_na3'

detect_Rieger1(model, data)

detect_Rieger2 = function(model, data, alpha = 0.05){
  fit = cfa(model, data = data)
  
  # first round detection
  detected = detect_Rieger1(model, data, alpha)
  p.val.min = unlist(lapply(detected, function(list) min(list$p.vals)))
  varnames_ni = unlist(lapply(detected, function(list) list$varnames[which.min(list$p.vals)]))
  
  # remove noninvariant items from model
  model.upd = update_model(model, paste0(varnames_ni[p.val.min<alpha],collapse = "|"))
  
  varnames = sapply(detected, "[[", 1, simplify = F)
  noninvariant = lapply(detected, function(list) list$varnames[which.min(list$p.vals)])

  while(any(p.val.min < alpha)){
    detected = detect_Rieger1(model.upd, data, alpha)
    p.val.min = unlist(lapply(detected, function(list) min(list$p.vals)))
    varnames_ni = unlist(lapply(detected, function(list) list$varnames[which.min(list$p.vals)]))
  
    for(j in 1:length(varnames)){
      if(p.val.min[j]<alpha){
        noninvariant[[names(p.val.min)[j]]] = c(noninvariant[[names(p.val.min)[j]]], detected[[j]]$varnames[which.min(detected[[j]]$p.vals)]) 
      }
    }
    model.upd = update_model(model.upd, paste0(varnames_ni[p.val.min<alpha],collapse = "|")) 
  }
  out = list()
  for(j in 1:length(varnames)){
    out[[j]] = list(varnames = varnames[[j]],
                    noninvariant = noninvariant[[j]])
  }
  out
}
detect_Rieger2(model, data)
