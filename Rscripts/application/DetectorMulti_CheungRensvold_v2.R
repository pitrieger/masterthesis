# Detector via CFI
# General idea: Test different pairings of reference variables and 1 other variable which is held constant across 
#               groups while the remaining parameters are estimated freely. Invariant items are those for which 
#               invariance is never rejected in chisqtest of multigroup baseline model and constrained models

# Differs from first version by also implementing constraints for intercepts
require(e1071)
require(stringr)

count_0rows = function(X){ # count number of consecutive rows from head of matrix that contain only zeros
  counter = 0
  if(!any(X[1,] != 0)){
    for(counter in 1:nrow(X)){
      if(any(X[counter,] != 0)){
        counter = counter-1
        break
      }
    }
  }
  counter
}

permuted_noninvariant = function(X, df, alpha){ #obtain noninvariant items
  p = ncol(X)
  n = p*(p-1)/2 # number of triangular entries in X
  crit = qchisq(1-alpha, df) # critical values for rejecting chisq test
  X_bin = matrix(as.numeric(X > crit), ncol = p) # make binary, s.t. significant pairs = 1
  X_bin[upper.tri(X_bin)] = t(X_bin)[upper.tri(X_bin)] # make symmetric

  perms = e1071::permutations(p) # all permutations
  P = sapply(1:nrow(perms), function(i) diag(p)[perms[i,],], simplify = F) # all permutation matrices
  
  perm.0rows = numeric(length(P)) # number of consecutive zero rows from top
  for(j in 1:length(P)){ 
    A = P[[j]] %*% X_bin %*% t(P[[j]]) # permute rows and columns simultaneously
    A[upper.tri(A, diag = T)] = 0 # set upper triangle to 0 
    perm.0rows[j] = count_0rows(A)
  }
  
  perm.0rows.max = which.max(perm.0rows)
  #P_final = P[[perm.0rows.max]]
  #out = P_final %*% X %*% t(P_final)
  #out
  #dimnames(out) = list(rownames(X)[perms[perm.0rows.max,]],
  #                     colnames(X)[perms[perm.0rows.max,]])
  varnames.permuted = rownames(X)[perms[perm.0rows.max,]] # varnames in order of permutation with max number of consecutive 0 rows from head
  if(max(perm.0rows)<p){
    varnames.permuted[(max(perm.0rows)+1):p] # names of noninvariant varnames (for which chisq>crit in at least one pairing, and thus rejected)
  } else {
    NULL
  }

}

detectMulti_CheungRensvold = function(varnames, base_model, data, alpha = 0.05, group.constraints = c("loadings", "intercepts")){
  # baseline model
  base_fit = cfa(base_model, data = data, group = "grp", meanstructure = T)
  
  # split model into components that contain varnames
  model_split = unlist(strsplit(base_model, "\n"))
  model_split_pos = which(grepl(paste0(varnames, collapse = "|"), model_split) &  # contains varnames
                              grepl("\\=\\~", model_split) & # is a latent variable line
                            !grepl("\\*", model_split)) # doesn't have fixed parameter
  model_split_op = str_extract(model_split, "\\=\\~|\\~\\~|\\~1")
  model_split_lhs = sapply(strsplit(model_split, "\\=\\~|\\~\\~|\\~1"), function(x) ifelse(!is.null(x) & length(x)>=1, x[1], NA))
  model_split_rhs = sapply(strsplit(model_split, "\\=\\~|\\~\\~|\\~1"), function(x) ifelse(!is.null(x) & length(x)>=2, x[2], NA))
  model_split_vars = sapply(strsplit(model_split_rhs, "\\+"), function(x) gsub("[[:blank:]]*", "", x), simplify = F)
  #model_split_nvar = sapply(model_split_vars, length)
  
  out = list(varnames = varnames, 
             noninvariant = NULL,
             noninvariant_bonferroni = NULL,
             convergenceissue = NULL,
             alpha = alpha)
  for(k in 1:length(model_split_pos)){
    line = model_split_pos[k]
    varnames_k = unlist(model_split_vars[line])
    p = length(varnames_k)
    
    # pairs of all indicators 
    pairs = matrix(0, ncol = p, nrow = p, 
                   dimnames=list(varnames_k,varnames_k))
    
    # fit constrained models for lower triangle of pairs & test against baseline model
    for(i in 2:p){ # rows = argument with equality constraint across groups
      for(j in 1:p){ # cols = reference item with loading = 1
        if(j<i){
          # generate new model formulation with constraints (implicit by adding coefficients)
          model_it = model_split
          if("loadings" %in% group.constraints){
            # permute variables and add constraint on loading (implicit by giving coefficient)
            model_line = paste(model_split_lhs[line], model_split_op[line], paste(c(varnames_k[j], paste0("COEF0*",varnames_k[i]), varnames_k[-c(i, j)]), collapse = " + "))
            model_it[line] = model_line
          }
          # add constraints on loading
          if("intercepts" %in% group.constraints){
            intercepts = varnames_k[c(i, j)]
            model_it = c(model_it, paste0(intercepts, " ~ COEF_", LETTERS[1:length(intercepts)], "*1"))
          }
          # combine lines to single string
          model_it = paste0(model_it, collapse = "\n")
          
          # fit new model and test against base_fit
          fit = cfa(model_it, data = data, group = "grp")
          cons.test = lavTestLRT(base_fit, fit)
          pairs[i,j] = cons.test$`Chisq diff`[2]
        }
      }
  }
  out$noninvariant = unique(c(out$noninvariant, 
                              permuted_noninvariant(pairs, df = cons.test$`Df diff`[2], alpha = alpha)))
  out$noninvariant_bonferroni = unique(c(out$noninvariant_bonferroni,
                                         permuted_noninvariant(pairs, df = cons.test$`Df diff`[2], alpha = alpha/length(varnames))))
  }
  out
}
