# Detector via CFI
# General idea: Test different pairings of reference variables and 1 other variable which is held constant across 
#               groups while the remaining parameters are estimated freely. Invariant items are those for which 
#               invariance is never rejected in chisqtest of multigroup baseline model and constrained models

## NOTE: 
# this only works for loadings?!
## TODO:  
# Check how to implement this for intercepts as well?! 
# Implement for CFI
# Check if this needs multiple testing correction when computing critical value in permuted_noninvariant

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

detect_CheungRensvold = function(varnames, data, alpha = 0.05, group.constraints = c("loadings")){
  p = length(varnames)
  
  # baseline model
  base_model = paste("eta =~", paste(varnames, collapse = " + "))
  base_fit = cfa(base_model, data = data, group = "grp")
  
  # pairs of all indicators 
  pairs = matrix(0, ncol = p, nrow = p, 
                 dimnames=list(varnames,varnames))
  
  # fit constrained models for lower triangle of pairs & test against baseline model
  for(i in 2:p){ # rows = argument with equality constraint across groups
    for(j in 1:p){ # cols = reference item with loading = 1
      if(j<i){
        model = paste("eta =~", paste(c(varnames[j], varnames[-j]), collapse = " + "))
        
        # NOTE: lavaan doesn't implement selecting individual parameters for constraints, 
        #       only entire families. Thus, group.equal constrains all loadings and group.partial
        #       then lifts the constrains for all but the argument variable
        fit = cfa(model, data = data, group = "grp", 
                  group.equal = group.constraints, # constrain all selected families of parameters across grps
                  group.partial = c(sapply(varnames[-c(i,j)], function(x) paste("eta =~", x)))) # deconstrain all loadings but argument
        cons.test = lavTestLRT(base_fit, fit)
        pairs[i,j] = cons.test$`Chisq diff`[2]
      }
    }
  }
  
  list(varnames = varnames,
       noninvariant = permuted_noninvariant(pairs, df = cons.test$`Df diff`[2], alpha = alpha),
       alpha = alpha)
}
