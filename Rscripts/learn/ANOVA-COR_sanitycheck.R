# ANOVA & Correlation component sanity check
library(tidyverse)
set.seed(51)

# constants
g = 4
n_l = rep(c(150, 200), g/2)
grp = sapply(1:g, function(x) rep(x, n_l[x])) %>% unlist()
N = sum(n_l)
eta_mean_l = rnorm(g)

# sim data - perfect MI
eta = rnorm(N) + eta_mean_l[grp]
boxplot(eta ~ grp)
Y = 0 + 2.2 * eta + rnorm(N, 0.5)
eps = residuals(lm(Y ~ eta))

# visualize
plot(eps ~ eta, col = grp)
boxplot(eps ~ grp)
par(mfrow = c(g/2, 2))
for(i in 1:g){
  plot(eps[grp == i] ~ eta[grp == i], col = i)
  abline(h = 0, lty = 2, col = 2)
  abline(lm(eps[grp == i] ~ eta[grp == i]), lwd = 2)
}
par(mfrow = c(1, 1))

# ANOVA -> subscipt i dropped b/c only one example item
summary(aov(eps ~ as.factor(grp)))
hat_mu = mean(eps)
hat_mu_l = sapply(1:g, function(l) mean(eps[grp == l]))
(MST = 1/(g - 1) * sum(n_l * (hat_mu_l - hat_mu)^2))
(MSE = 1/(N - g) * sum(sapply(1:g, function(l) sum((eps[grp == l] - hat_mu_l[l])^2))))

Fstat = MST/MSE
pf(Fstat, g-1, N-g, lower.tail = F)

# Correlation
regs = list()
for(l in 1:g){
  regs[[l]] = lm(eps[grp == l] ~ eta[grp == l])
}

hat_alpha_l = lapply(regs, function(fit) fit$coefficient[1]) %>% unlist
hat_beta_l = lapply(regs, function(fit) fit$coefficient[2]) %>% unlist
hat_hat_eps = lapply(regs, fitted) %>% unlist
hat_hat_eps_check = sapply(1:g, function(l) hat_alpha_l[l] + hat_beta_l[l]*eta[grp == l]) %>% unlist()
plot(hat_hat_eps, hat_hat_eps_check)

hat_se_hat_beta_l = sqrt((sapply(1:g, function(l) sum((eps[grp == l] - hat_hat_eps[grp == l])^2)) %>% unlist())/
  ((n_l - 2)*(sapply(1:g, function(l) sum((eta[grp == l] - eta_mean_l[l])^2)) %>% unlist())))

Tstat = hat_beta_l / hat_se_hat_beta_l

2 * (1 - pt(abs(Tstat), n_l-2, lower.tail = T))
lapply(regs, function(fit) summary(fit)$coefficients[2,4]) %>% unlist
