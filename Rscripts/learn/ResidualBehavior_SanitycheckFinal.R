# SCALAR VIOLATION
n = 1000
nsim = 1000
alpha_bias = alpha_est = numeric(nsim)
beta_bias  = beta_est = numeric(nsim)
phi_est = numeric(nsim)
mu_est = numeric(nsim)
p = 0.6
tau = c(-4, 5)
lam = c(0.7, 0.7)

# eta parameters
mu = c(0.5, 2.5)
phi = c(1, 2.5)

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
(beta_true = p * (1-p) * (tau[1] - tau[2]) * (mu[1] - mu[2])/phi_true + lam[1])
hist(beta_est)
abline(v = beta_true, lwd = 2)

# pooled alpha
mean(alpha_est)
(alpha_true = (p * tau[1] + (1-p) * tau[2]) + lam[1]*mu_true - beta_true * mu_true)
hist(alpha_est)
abline(v = alpha_true, lwd = 2)


# METRIC VIOLATION
n = 1000
nsim = 1000
alpha_bias = alpha_est = numeric(nsim)
beta_bias  = beta_est = numeric(nsim)
phi_est = numeric(nsim)
mu_est = numeric(nsim)
p = 0.6
tau = c(2, 2)
lam = c(1.5, 0.6)

# eta parameters
mu = c(0.5, 2.5)
phi = c(1, 2.5)

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
(beta_true = (lam[1] * p * (phi[1] + (1 - p) * mu[1]^2 + (p-1)*mu[1]*mu[2]) +
             lam[2] * ((1-p) * (phi[2] + p * mu[2]^2) + (p-1)*p*mu[1]*mu[2]))/phi_true)
hist(beta_est)
abline(v = beta_true, lwd = 2)

# pooled alpha
mean(alpha_est)
(alpha_true = tau[1] + p * lam[1] * mu[1] + (1-p) * lam[2]*mu[2] - beta_true * mu_true)
hist(alpha_est)
abline(v = alpha_true, lwd = 2)





# SIMULTANEOUS VIOLATION
n = 1000
nsim = 1000
alpha_bias = alpha_est = numeric(nsim)
beta_bias  = beta_est = numeric(nsim)
phi_est = numeric(nsim)
mu_est = numeric(nsim)
p = 0.6
tau = c(1, 2.45)
lam = c(1.5, 0.6)

# eta parameters
mu = c(0.5, 2.5)
phi = c(1, 2.5)

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
(beta_true = (p*(1-p)*(tau[1] * (mu[1] - mu[2]) + 
          tau[2] * (mu[2] - mu[1])) +
 lam[1] * p * (phi[1] + (1-p) * mu[1]^2 + mu[1]*mu[2] * (p-1)) +
 lam[2] * ((1-p) * (phi[2] + p * mu[2]^2) + mu[1]*mu[2]*p*(p-1)))/phi_true)
hist(beta_est)
abline(v = beta_true, lwd = 2)

# pooled alpha
mean(alpha_est)
(alpha_true = p*tau[1] + (1-p)*tau[2] + p*lam[1]*mu[1] + (1-p)*lam[2]*mu[2] - beta_true * mu_true)
hist(alpha_est)
abline(v = alpha_true, lwd = 2)

