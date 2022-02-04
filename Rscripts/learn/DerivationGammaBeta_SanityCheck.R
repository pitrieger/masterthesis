nsim = 100
n = 1000
alpha_bias = alpha_est = numeric(nsim)
beta_bias  = beta_est = numeric(nsim)
phi_est = numeric(nsim)
mu_est = numeric(nsim)
g = 4
p = c(0.2, 0.4, 0.3, 0.1)
tau = runif(g, -2, 2)
lam = rnorm(g, 0, 2)

# eta parameters
mu = rnorm(g, 0, 2)
phi = rgamma(g, 2, 1)

for(j in 1:nsim){
  G = rmultinom(n, 1, p)
  G = apply(G, 2, function(l) which(l == 1))
  eta = sapply(G, function(l) rnorm(1, mu[l], sqrt(phi[l])))
  Y = tau[G] + lam[G]*eta + rnorm(n, 0, 0.25)
  #plot(Y~eta, col = G)
  
  fit = lm(Y ~ eta)
  alpha_est[j] = coef(fit)[1]
  beta_est[j] = coef(fit)[2]
  phi_est[j] = var(eta)
  mu_est[j] = mean(eta)
}

# pooled latent mean
mean(mu_est)
(mu_true = sum(p*mu))
#hist(mu_est)
#abline(v = mu_true, lwd = 2)

# pooled variance
print("phi")
mean(phi_est)
ndterm <- 0
for(l in 1:g){
  for(k in 1:g){
    if(l != k){
      ndterm = ndterm + (p[l]*p[k]*mu[l]*mu[k])
    }
  }
}
(phi_true = sum(p * (phi + (1 - p) * mu^2)) - ndterm)
(phi_true = sum(p * (phi + mu * (mu - mu_true))))

# pooled beta
mean(beta_est)
nd_termA = nd_termB = 0
for(l in 1:g){
  for(k in 1:g){
    if(l != k){
      nd_termA = nd_termA + tau[l] * p[l] * p[k] * mu[k]
      nd_termB = nd_termB + lam[l] * p[l] * p[k] * mu[l] * mu[k]
    }
  }
}
(beta_true = (sum(tau * p * (mu-mu_true)) +
                sum(lam * p * (phi + (1 - p) * mu^2)) - nd_termB)/phi_true)
(beta_true = (sum(tau * p * (mu-mu_true)) +
                sum(lam * p * (phi + mu * (mu - mu_true))))/phi_true)
#hist(beta_est)
#abline(v = beta_true, lwd = 2)

#pooled alpha
mean(alpha_est)
(alpha_true = sum(p * (tau + lam*mu)) - beta_true * mu_true)
(alpha_true = sum(p * (tau + lam*mu - mu_true/phi_true * (tau * (mu - mu_true) + lam * (phi + mu * (mu - mu_true))))))

sum(p * (tau * (1 - ((mu_true * (mu - mu_true)) / phi_true)) + lam * (mu - (mu_true * (phi + mu * (mu - mu_true)))/phi_true)))

hist(alpha_est)
abline(v = alpha_true, lwd = 2)