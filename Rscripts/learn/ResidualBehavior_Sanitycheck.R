

# METRIC VIOLATION
n = 1000
nsim = 1000
alpha_bias = alpha_est = numeric(nsim)
beta_bias  = beta_est = numeric(nsim)
phi_est = numeric(nsim)
mu_est = numeric(nsim)
p = 0.8
tau = c(3, 5)
lam = c(0.7, 0.7)

# eta parameters
mu = c(0.5, 2.5)
phi = c(1, 9)

for(j in 1:nsim){
  G = rbinom(n, 1, p)
  G[G == 0] = 2 # for more efficient coding
  eta = sapply(G, function(l) rnorm(1, mu[l], sqrt(phi[l])))
  Y = tau[G] + lam[G]*eta + rnorm(n, 0, 0.25)
  #plot(Y~eta, col = G)
  
  fit = lm(Y ~ eta)
  alpha_est[j] = coef(fit)[1]
  beta_est[j] = coef(fit)[2]
  phi_est[j] = var(eta)
  mu_est[j] = mean(eta)
}

# pooled variance
mean(phi_est)
(phi_true = p * (phi[1] + (1-p) * mu[1]^2) + 
          (1-p) * (phi[2] + p * mu[2]^2) -
          2 * p * (1-p) * mu[1] * mu[2])
hist(phi_est)
abline(v = phi_true, lwd = 2)

# pooled latent mean
mean(mu_est)
(mu_true = p*mu[1] + (1-p) * mu[2])
hist(mu_est)
abline(v = mu_true, lwd = 2)

# pooled beta
mean(beta_est)
(beta_true = p * (tau[1] * (mu[1] * (1-p) + mu[2] * (p-1)) +
                 tau[2] * (mu[1] * (p-1) + mu[2] * (1-p)))/phi_true + lam[1])
hist(beta_est)
abline(v = beta_true, lwd = 2)

# pooled alpha
mean(alpha_est)
(alpha_true = (p * tau[1] + (1-p) * tau[2]) + lam[1]*mu_true - beta_true * mu_true)
hist(alpha_est)
abline(v = alpha_true, lwd = 2)





# check if alpha is correctly derived
nsim = 1000
alpha_bias = numeric(nsim)
beta_bias = numeric(nsim)
n_k = c(50, 100, 20)
g = c( rep(1, each = n_k[1]), rep(2, each = n_k[2]), rep(3, each = n_k[3]))
alpha_k = c(0, 2, 1)
beta_k = c(2, 2, 2)
mu_k = rnorm(3, 0, 3)

beta = numeric(nsim)
for(j in 1:nsim){
  #mu_k = c(-1, 0, 2)
  X = sapply(g, function(l) rnorm(1, mu_k[l]))
  #X = rnorm(sum(n_k))
  #X = scale(X, scale = T)
  Y = alpha_k[g] + beta_k[g]*X + rnorm(sum(n_k), 0, 0.1)
  #plot(Y ~ X, col = g)
  
  fit = lm(Y ~ X)
  alpha_est = coef(fit)[1]
  beta[j] = coef(fit)[2]
  (alpha_true = sum(n_k/sum(n_k) * (alpha_k)))
  #(alpha_true = sum(n_k/sum(n_k) * (alpha_k + beta_k*mu_k)))
  alpha_bias[j] = alpha_est - alpha_true
}

hist(beta)
hist(alpha_bias)

G = rbinom(1000, 1, 0.5)
X = sapply(G, function(x) rnorm(1, 2*x))
cov(G*X, (1-G)*X)
plot(G*X, (1-G)*X)
plot(X, G)
G*X
summary(alpha_bias)
plot(Y ~ X)


  (beta_est = coef(fit)[2])
  
  N = sum(n_k)
  a = sum(n_k/(N-1) * (mu_k * alpha_k + beta_k + beta_k * mu_k^2))
  b.1 = sum(n_k/(N-1) * mu_k * sum(n_k/N * alpha_k))
  b.2 = sum(n_k/(N-1) * x)
  
  b.2 = numeric(1)
  for(k in 1:3){
    inner = numeric(3)
    for(l in 1:3){
      inner[l] = n_k[k]*n_k[l]/N * beta_k[l]*ifelse(l == k, 1 + mu_k[k]^2, mu_k[k]*mu_k[l])
    }
    b.2 = b.2 + 1/(N-1) * sum(inner[l]) 
  }
  b.2
  a - (b.1 + sum(b.2))
  
  
  (beta_true = sum(n_k/sum(n_k) * (mu_k * alpha_k + beta_k - n_k/sum(n_k) * (mu_k * alpha_k + beta_k))))
  beta_bias[j] = beta_est - beta_true
}
hist(alpha_bias)
summary(alpha_bias)

hist(beta_bias)
summary(beta_bias)

n = 1000
Z = rnorm(n)
sum((Z - mean(Z))^2)/(n-1)
var(Z)




n = 1000
G = rbinom(n, 1, 0.3)
X = rnorm(1000)
cov(G*X, X)
sum(G)/n*cov(X, X)

nsim = 10000
alpha_bias = beta_bias = numeric(nsim)

for(j in 1:nsim){
  # pop params
  alpha1 = rnorm(1, 1)
  alpha2 = rnorm(1)
  beta1 = rnorm(1, 2)
  beta2 = rnorm(1, 0, 5)
  p = rbeta(1, 2, 2)
  
  # temp
  alpha1 = 0
  alpha2 = 0
  beta1 = 1
  beta2 = -2
  
  # RVs
  G = rbinom(n, 1, p)
  X = rnorm(n, G*2)
  Y = alpha1 + G*alpha2 + X*beta1 + G*X*beta2 + rnorm(n)
  #plot(Y ~ X, col = G+1)
  
  # regress
  fit = lm(Y ~ X)

  # estimated in reg
  alpha_est = coef(fit)[1]
  beta_est = coef(fit)[2]
  
  # by hand from pop values
  beta_true = beta1 + p*beta2
  alpha_true = alpha1 + p*alpha2 - 
  
  # bias
  alpha_bias[j] = alpha_est - alpha_true
  beta_bias[j] = beta_est - beta_true
}
hist(alpha_bias)
hist(beta_bias)
summary(alpha_bias)
summary(beta_bias)

mean(X[G == 0])
mean(X[G == 1])



g = 2
(tau_k = 5*runif(g))
(lambda_k = rnorm(g, 2))
mu_k = rnorm(2)
n = 1000
grp = rep(1:g, each = n)
alpha_bias = numeric(1000)
for(j in 1:1000){
  X = sapply(1:g, function(l) rnorm(n, mu_k[l]), simplify = F) %>% unlist()
  X = scale(X)
  Y = tau_k[grp] + lambda_k[grp]*X + rnorm(g*n, 0, 0.5)
  #plot(Y ~ X, col = grp)
  
  (fit = lm(Y ~ X))
  alpha_bias[j] = coef(fit)[1] - mean(tau_k + lambda_k * mu_k)
}
hist(alpha_bias)

X = sapply(1:g, function(l) rnorm(n, mu_k[l]), simplify = F) %>% unlist()
X = scale(X)
Y = tau_k[grp] + lambda_k[grp]*X + rnorm(g*n, 0, 0.5)
plot(Y ~ X, col = grp)
abline(lm(Y ~ X))

cov(X, Y)

(fit = lm(Y ~ X))

mean(tau_k) + mean(lambda_k * mu_k)

Y_tilde = Y - mean(tau_k + lambda_k * mu_k)
plot(Y_tilde ~ X, col = grp)
mean(Y_tilde)
lambda_k 





# 8 cases
library(tidyverse)
par = expand.grid(tau_eq = c(T,F),
                  lam_eq = c(T,F),
                  eta_eq = c(T,F))
g = 2
n = 200
plots = list()
for(j in 1:nrow(par)){
  if(par$eta_eq[j]){
    mu_k = c(0, 0)
  } else {
    mu_k = c(0, 2)
  }
  if(par$tau_eq[j]){
    tau_k = c(2, 2)
  } else {
    tau_k = c(2, 5)
  } 
  if(par$lam_eq[j]){
    lam_k = c(1, 1)
  } else {
    lam_k = c(1, -1)
  }
  G = rep(c(1, 2), each = n)
  X = c(rnorm(n, mu_k[1]), rnorm(n, mu_k[2]))
  Y = rep(tau_k, each = n) + rep(lam_k, each = n)*X + rnorm(2*n)
  
  df = data.frame(Y, X, G = as.factor(G))
  plots[[j]] = ggplot(df, aes(y = Y, x = X, color = G)) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    geom_smooth(inherit.aes = F, aes(y = Y, x = X), method = "lm") +
    ggtitle(paste("Tau", ifelse(par$tau_eq[j], "equal", "diff"),
                  "Lambda", ifelse(par$lam_eq[j], "equal", "diff"),
                  "Eta", ifelse(par$eta_eq[j], "equal", "diff")
    ))
}
library(egg)
ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2)
ggarrange(plots[[5]], plots[[6]], plots[[7]], plots[[8]], nrow = 2)


plots



# SCALAR VIOLATION
n = 1000
nsim = 1000
alpha_bias = alpha_est = numeric(nsim)
beta_bias  = beta_est = numeric(nsim)
phi_est = numeric(nsim)
mu_est = numeric(nsim)
p = 0.6
tau = c(4, 4)
lam = c(0, 10)

# eta parameters
mu = c(0, 1)
phi = c(1, 1)

tau[1] - tau[2]
(-lam[1]) * (mu[1] - mu[2])
r1_est = r2_est = numeric(nsim)
for(j in 1:nsim){
  G = rbinom(n, 1, p)
  G[G == 0] = 2 # for more efficient coding
  eta = sapply(G, function(l) rnorm(1, mu[l], sqrt(phi[l])))
  Y = tau[G] + lam[G]*eta + rnorm(n, 0, 0.25)
  #plot(Y~eta, col = G)
  
  fit = lm(Y ~ eta)
  r1_est[j] = mean(Y[G == 1] - predict(fit)[G==1])
  r2_est[j] = mean(Y[G == 2] - predict(fit)[G==2])
}
hist(r1_est, freq = F, xlim = range(c(r1_est, r2_est)))
hist(r2_est, freq = F, add = T)

tau = c(4, 5)
lam = c(1, 1)

# eta parameters
mu = c(1.5, 0.5)
phi = c(1, 2.5)

G = rbinom(n, 1, p)
G[G == 0] = 2 # for more efficient coding
eta = sapply(G, function(l) rnorm(1, mu[l], sqrt(phi[l])))
Y = tau[G] + lam[G]*eta + rnorm(n, 0, 0.25)
plot(Y~eta, col = G)





x = seq(-10, 10, l = 500)
plot(x, 0.5 - 0.5^3 * x^2 / (1+ (x/2)^2), type = "l")

x = 1000
max(0.5^3 * x^2 / (1+ (x/2)^2))







install.packages("rootSolve")
f = function(tau_1, tau_2, mu_1, mu_2, lam_1, lam_2){
  (tau_1 - tau_2)/2 + 0.5*lam_1*mu_1 - 0.5*lam_2*mu_2 - (0.5^2 * (tau_1*(mu_1 - mu_2) - tau_2*(mu_1-mu_2)) +
                                                           lam_1 * 0.5*(1 + 0.5 * mu_1^2 - 0.5*mu_1*mu_2) +
                                                           lam_2 * 0.5*(1 + 0.5 * mu_2^2 - 0.5*mu_1*mu_2))/(1 + ((mu_1 - mu_2)/2)^2) * 
    (mu_1 - mu_2)/2}
library(rootSolve)
multiroot(f, c(tau_1 = 1, tau_2 = 1, 
               mu_1 = 1, mus_2 = 1, 
               lam_1 = 1, lam_2 = 1))

f(1, 1, 1, 1, 1, 1)

G = rmultinom(1000, 1, c(1, 1, 1))
G[,1]

eta = apply(G,2, function(x) rnorm(1, which(x == 1)))
dim(G)











