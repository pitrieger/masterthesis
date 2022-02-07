# Detector Rieger
library(stringr)
library(tidyverse)

get_pvalMulti = function(data, hat.eta, grp, models, varnames){
  g = length(unique(grp))
  
  # identify which latent variables underly each item
  covariates = sapply(varnames, function(vn) which(grepl(vn, models)))

  p.vals = rep(1,length(varnames))
  for(j in 1:length(covariates)){
      # relevant LV for item j
      hat.eta.j = as.matrix(hat.eta[,covariates[[j]]])
      # get residuals
      R = resid(lm(data[,varnames[j]] ~ hat.eta.j))
      # group means
      aov.p = summary(aov(R ~ as.factor(grp)))[[1]]$`Pr(>F)`[1]

    for(i in 1:ncol(hat.eta.j)){
      # group correlations
      cor.p = numeric(g)
      for(l in 1:g){
        cor.p[l] = cor.test(R[grp == l], hat.eta.j[grp == l, i])$p.value # equivalent to coefficient test for lm
      }
      # aggregate p.vals
      p.vals[j] = min(p.vals[j], (ncol(hat.eta.j)*g+1)*min(aov.p, cor.p))
    }  
      
  }
  pmin(1, p.vals)
}


get_pvalMulti_metric = function(data, hat.eta, grp, models, varnames){
  g = length(unique(grp))
  
  # identify which latent variables underly each item
  covariates = sapply(varnames, function(vn) which(grepl(vn, models)))
  
  p.vals = rep(1,length(varnames))
  for(j in 1:length(covariates)){
    # relevant LV for item j
    hat.eta.j = as.matrix(hat.eta[,covariates[[j]]])
    
    # get residuals
    R = resid(lm(data[,varnames[j]] ~ hat.eta.j))
    
    cor.p = 1
    for(i in 1:ncol(hat.eta.j)){
      hat.eta.ji = hat.eta.j[,i]
      fit = lm(R ~ hat.eta.ji*as.factor(grp))
      cor.p = min(cor.p, ncol(hat.eta.j)*anova(fit)["hat.eta.ji:as.factor(grp)", "Pr(>F)"])
    }
    
    # save pval if lower than existing one
    p.vals[j] = min(p.vals[j], length(covariates)*cor.p)
  }
  pmin(1, p.vals)
}

detectMulti_Rieger = function(varnames, base_model, data, alpha = 0.05, detection.type = "both"){
  #model = paste("eta =~", paste(varnames, collapse = " + "))
  fit = cfa(base_model, data = data)
  
  # get model components which contain a latent variable and at least one of the indicators
  base_model_split = unlist(strsplit(base_model, "\n"))
  base_model_split = base_model_split[grepl(paste0(varnames, collapse = "|"), base_model_split) & 
                                        grepl("\\=\\~", base_model_split)]
  
  # predicted latent variables
  hat.eta = predict(fit, newdata = data)
  
  # get p - values
  if(detection.type == "both"){
    p.vals = get_pvalMulti(data = data, hat.eta = hat.eta, grp = data$grp, models = base_model_split, varnames = varnames)
  } else if(detection.type == "metric"){
    p.vals = get_pvalMulti_metric(data = data, hat.eta = hat.eta, grp = data$grp, models = base_model_split, varnames = varnames)
  }
  
  list(varnames = varnames,
       noninvariant = varnames[which(p.vals < alpha/(length(p.vals) - order(order(p.vals)) + 1))], # Bonferroni-holm
       p.vals = p.vals)
}

detectMulti_Rieger_step = function(varnames, base_model, data, alpha = 0.05, detection.type = "both"){
  model_it = base_model
  detected = detectMulti_Rieger(varnames, base_model, data, alpha, detection.type)
  
  p.val.min = min(detected$p.vals)
  varnames_it = varnames
  keep = NULL
  
  while(p.val.min < alpha/length(varnames_it)){ # 
    # identify variable to be removed -> mustn't be in keep
    varname_remove = detected$varnames[which(detected$p.vals == p.val.min)]
    
    # update model & model_split & keep & varnames
    model_it = gsub(varname_remove, "", model_it) %>%
      gsub("\\S*\\*[[:blank:]]*\\+", "", .) %>% # removed var at beginning or middle with coefficient constraints such as lv1 =~ x1 + b1*x2 + b1*x3
      gsub("\\S*\\*[[:blank:]]*\n", "\n", .) %>% # removed var at end with coefficient constraint
      gsub("\\~[[:blank:]]*\\+", "~", .) %>% # removed var at beginning
      gsub("\\+[[:blank:]]*\\+", "+", .) %>% # removed var in middle
      gsub("\\+[[:blank:]]*\n|\\+[[:blank:]]*$", "\n", .) # removed var at end of line or very end
    cat(model_it)
    
    model_split = unlist(strsplit(model_it, "\n"))
    model_split = model_split[grepl(paste0(varnames, collapse = "|"), model_split) & 
                                grepl("\\=\\~", model_split)]
    model_split_rhs = sapply(strsplit(model_split, "\\=\\~"), "[[", 2)
    model_split_vars = sapply(strsplit(model_split_rhs, "\\+"), function(x) gsub("[[:blank:]]*", "", x), simplify = F)
    model_split_nvar = sapply(model_split_vars, length) # number of variables
    if(any(model_split_nvar <= 1)){
      keep = unlist(model_split_vars[which(model_split_nvar <= 1)])
    }
    varnames_it = detected$varnames[-which(detected$p.vals == p.val.min)]
  
    if(all(varnames_it %in% keep)){
      break
    }
    detected = detectMulti_Rieger(varnames_it, model_it, data, alpha, detection.type)
    
    # min p.val among variables not in keep
    p.val.min = min(detected$p.vals[!(detected$varnames %in% keep)]) # Bonferroni-Holm
  }
  list(varnames = varnames,
       noninvariant = varnames[!varnames %in% varnames_it])
}
