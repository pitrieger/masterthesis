### ### ### ### ### ###
## General Problem  ###
### ### ### ### ### ###
library(lavaan)
library(car)
library(tidyverse)

n = 1000 # number of obs
k = 2 # number of latent vars
p = 10 # number of manifest vars
g = 2 # number of groups/populations/environments
grp = rep(1:g, each = floor(n/g)) # group index
gen_error = function(dist = "norm", n, var, ...){
  dist = paste0("r", dist)
  p = length(var)
  sd = sqrt(var)
  sapply(1:p, function(i) do.call(dist, args = list(n = n, sd = sd[i])))
}

eta = matrix(rnorm(n*k), ncol = k) # latent vars
eta = eta + grp - 1
tau = matrix(rep(0, p*g), nrow = p) # identical intercepts across groups
psi = matrix(rep(1, p*g), nrow = p) # identical variances across groups
Lambda = array(c(rep(1, k), rnorm(k*(p-1))), dim = c(k, p, g)) # identical loadings across groups

# violations of invariant loadings
Lambda[1,4,2]= 3 # fix loading of 1st latent variable on 4th indicator for 2nd group
Lambda[1,4,1]

# violations of different intercepts
tau[2, 2] = 2 # change intercept for 2nd indicator for 2nd group

# Generate Data
X = do.call("rbind", sapply(1:g, 
                            function(j) t(tau[,j] + t(Lambda[,,j]) %*% t(eta[grp == j,])) +
                              gen_error(dist = "norm", n = sum(grp == j), var = psi[,j]) , 
                            simplify = F))
X = as.data.frame(X)
colnames(X) = paste0("x", 1:10)
X$grp = grp

# fit CFA
mod = "eta1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10"
fit = cfa(mod, data = X)
summary(fit)

# performance
eta_pred = predict(fit)
plot(eta[,1], eta_pred)

# item with DIF
plot(eta_pred, X[,4], col = grp)
plot(eta_pred, X[,2], col = grp)

# test residuals
X_resid = matrix(nrow = n, ncol = p)
p_val = numeric(p)
for(i in 1:p){
  X_resid[,i] = resid(lm(X[,i] ~ eta_pred))
  
  # test equal mean, equal variance and zero correlation with latent variable
  p.mean = summary(aov(X_resid[,i] ~ as.factor(grp)))[[1]][1,5]
  p.var = bartlett.test(X_resid[,i], grp)$p.value
  #p.cor = numeric(g)
  #for(j in 1:g){
  #  p.cor[j] = cor.test(X_resid[grp == j,i], eta_pred[grp == j,])$p.value
  #}
  #p_val[i] = min(3*c(p.mean, p.var, g*min(p.cor)), 1)
  p_val[i] = min(2*c(p.mean, p.var), 1)
}
plot(p_val)

plot(eta_pred, X_resid[,2], col = grp)


## Testing via lavTestLRT
fit.grp <- cfa(mod, data = X, group = "grp")

# weak invariance
fit.grp_load <- cfa(mod, data = X, group = "grp",
                    group.equal = "loadings")

# strong invariance
fit.grp_load.inter <- cfa(mod, data = X, group = "grp",
                          group.equal = c("intercepts", "loadings"))

# model comparison tests
lavTestLRT(fit.grp, fit.grp_load, fit.grp_load.inter)

X = X %>% group_by(grp) %>%
  mutate(across(1:p, function(x) scale(x, scale = F)))
mean(X$x2[X$grp == 2])

#
fit.grp <- cfa(mod, data = X, group = "grp")

# weak invariance
fit.grp_load <- cfa(mod, data = X, group = "grp",
                    group.equal = "loadings")

# strong invariance
fit.grp_load.inter <- cfa(mod, data = X, group = "grp",
                          group.equal = c("intercepts", "loadings"))
lavTestLRT(fit.grp, fit.grp_load, fit.grp_load.inter)


fit = cfa(mod, data = X)
eta_pred = predict(fit)
plot(eta[,1], eta_pred, col = grp)
mean(eta[grp == 1,1])
mean(eta[grp == 2,1])
mean(eta_pred[grp == 1,1])
mean(eta_pred[grp == 2,1])

# item with DIF
plot(eta_pred, X$x4, col = grp)
plot(eta_pred, X$x2, col = grp)


## Solutions: 
# Step 1: Demean => I think this is the wrong intuition. There may be valid differences that would 
#        be eradicated due to this...
# Get rid of full violations
# Step 2: 



# Stuff to do tomorrow: 
# - write down the implication of demeaning (see above). If this really is just a measurement issue 
#   then this should not eradicate valid differences between groups. Does this assume there are no 
#   differences on LV? 
# - implement simulation for when there actually is a difference between groups on LV 
# - 


library(sirt)
mod = "eta1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
       eta1 ~~ 1*eta1"
fit.grp <- cfa(mod, data = X, group = "grp")
summary(fit.grp)
fit.grp@ParTable
fit.grp@ParTable$est[fit.grp@ParTable$group == 1 & fit.grp@ParTable$op == "=~"]
lel = summary(fit.grp)[[1]]
lel[1:10,]
lambda = rbind(fit.grp@ParTable$est[fit.grp@ParTable$group == 1 & fit.grp@ParTable$op == "=~"], 
               fit.grp@ParTable$est[fit.grp@ParTable$group == 2 & fit.grp@ParTable$op == "=~"])
nu = rbind(fit.grp@ParTable$est[fit.grp@ParTable$group == 1 & fit.grp@ParTable$op == "~1" & fit.grp@ParTable$lhs != "eta1"],
           fit.grp@ParTable$est[fit.grp@ParTable$group == 2 & fit.grp@ParTable$op == "~1" & fit.grp@ParTable$lhs != "eta1"])
n = 1000
wgt <- matrix(sqrt(500), 2, 10)
mod1 <- sirt::invariance.alignment(lambda, nu, wgt)
par = summary(mod1)
class(mod1)
lambda[,4]
mod1$lambda.aligned[,4]

nu
mod1$nu.aligned
