## Pokropek et al 2019 data generation
library(tidyverse)
## MISC
# discretizer function:
discretize = function(x){
  sign(x) * ifelse(abs(x)>1.3, 2,
                   ifelse(abs(x)>0.47, 1, 0))
}

# simulator
sim_PMI = function(n = 1500, #obs
                   g = 24, #groups
                   p = 5, #items
                   k = 3, # noninvariant items
                   h = 0.5, #affected groups
                   itembias = 0.2, # absolute bias for both intercept and loading in affected groups
                   interceptbias = itembias, 
                   loadingbias = itembias,
                   discrete = T, # return discretized manifest variables
                   ...) { 
  stopifnot(k<=p,
            (h*g)%%1 == 0, 
            h>=0 & h<=1)
  
  # groups
  g_ind = rep(1:g, each = n)
  
  ## Step 1 ====
  # latent variable parameters
  mu_g = c(0, rnorm(g-1, 0, 0.3))
  sd_g = c(1, abs(rnorm(g-1, 1, 0.1)))
  
  # latent variable
  eta = sapply(1:g, function(i) rnorm(n, mu_g[i], sd_g[i]), simplify = F) %>% unlist() 
  
  # criterion variables
  C1 = 0.3*eta + rnorm(n*g)
  C2 = 0.1*eta + rnorm(n*g)
  
  ## Step 2 ====
  # item parameters
  tau = tau_affected = rnorm(p, 0, 0.5)
  lambda = lambda_affected = runif(p, 0.65, 0.85)
  eps = abs(1 - lambda^2) # unique sd
  
  ## Step 3 ====
  # items where MI holds
  Y = sapply(1:p, function(j) rnorm(n*g, tau[j] + lambda[j]*eta, eps[j]))
  
  ## Step 4 & 5 ====
  g_affected = p_affected = NULL
  if(h*p>0){
    ## Step 4 & 5 ====
    # sample affected groups
    g_affected = sample(1:g, size = h*g)
    # sample affected items
    p_affected = sample(1:p, size = k)
    
    # sample bias in each group separately
    for(j in 1:(h*g)){
      # add bias to model parameters
      tau_affected[p_affected] = tau[p_affected] + sample(c(-1, 1), replace = T, k) * interceptbias
      lambda_affected[p_affected] = lambda[p_affected] + sample(c(-1, 1), replace = T, k) * loadingbias
      eps_affected = abs(1 - lambda_affected^2) # unique sd
      
      # generate data where MI doesn't hold ( STEP 5 )
      Y[g_ind == g_affected[j],] = sapply(1:p, function(j) rnorm(n, tau_affected[j] + lambda_affected[j]*eta, eps_affected[j]))
    }
  }
  
  # Discretize 
  if(discrete){
    Y = apply(Y, 2, discretize)
  }
  
  # Output
  dat = as.data.frame(Y)
  colnames(dat) = paste0("y", 1:p)
  dat$grp = g_ind
  dat$C1 = C1
  dat$C2 = C2
  
  list(sim_dat = dat,
       sim_latent = eta,
       g_affected = g_affected,
       p_affected = p_affected,
       sim_params = c("n" = n, 
                  "g" = g,
                  "p" = p,
                  "k" = k,
                  "h" = h, 
                  "interceptbias" = interceptbias,
                  "loadingbias" = loadingbias))
}