library(ggplot2)
library(ggpubr)
library(lavaan)
library(here)
set.seed(123)
n = 400
envir = rep(1:2, each = n/2)
eta = c(rnorm(n/2), rnorm(n/2, 0))
mod = "eta =~ Y1 + Y2 + Y3 + Y4"

# full measurement invariance
df_MI = data.frame(Y1 = eta + rnorm(n, 0, 0.2),
                   Y2 = 0.6 * eta + rnorm(n, 0, 0.2),
                   Y3 = 0.9 * eta + rnorm(n, 0, 0.2),
                   Y4 = 2 + 0.4 * eta + rnorm(n, 0, 0.2))
fit_MI = cfa(mod, df_MI)
eta_MI = predict(fit_MI)[,1]

df1 = data.frame(x = predict(lm(df_MI$Y4 ~ eta_MI)),
                 y = resid(lm(df_MI$Y4 ~ eta_MI)),
                 envir = as.factor(envir))
p1 = ggplot(df1, aes(x = x, y = y, color = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  labs(x = expression(hat(eta)),
       y = expression(hat(epsilon)[i])) + 
  theme_minimal()
p1

# unequal taus
df_tau = data.frame(Y1 = eta + rnorm(n, 0, 0.2),
                   Y2 = 0.6 * eta + rnorm(n, 0, 0.2),
                   Y3 = 0.9 * eta + rnorm(n,  0, 0.2),
                   Y4 = (envir == 2)*0.9 + 0.4 * eta + rnorm(n, 0, 0.2))
fit_tau = cfa(mod, df_tau)
eta_tau = predict(fit_tau)[,1]
summary(lm(df_tau$Y4 ~ eta_MI))


df2 = data.frame(x = predict(lm(df_tau$Y4 ~ eta_tau)),
                 y = resid(lm(df_tau$Y4 ~ eta_tau)),
                 envir = as.factor(envir))
p2 = ggplot(df2, aes(x = x, y = y, color = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  labs(x = expression(hat(eta)),
       y = expression(hat(epsilon)[i])) + 
  theme_minimal()
p2

# unequal lambdas
df_lam = data.frame(Y1 = eta + rnorm(n, 0, 0.2),
                    Y2 = 0.6 * eta + rnorm(n, 0, 0.2),
                    Y3 = 0.9 * eta + rnorm(n, 0, 0.2),
                    Y4 = (0.4+(envir == 2)*0.5) * eta + rnorm(n, 0, 0.2),
                    grp = envir)
fit_lam = cfa(mod, df_lam)
eta_lam = predict(fit_lam)[,1]

df3 = data.frame(x = predict(lm(df_lam$Y4 ~ eta_lam)),
                 y = resid(lm(df_lam$Y4 ~ eta_lam)),
                 envir = as.factor(envir))
p3 = ggplot(df3, aes(x = x, y = y, color = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  labs(x = expression(hat(eta)),
       y = expression(hat(epsilon)[i])) + 
  theme_minimal()
p3
detect_Rieger(c("Y1", "Y2", "Y3", "Y4"), df_lam)

# unequal lambdas & taus
df_taulam = data.frame(Y1 = eta + rnorm(n, 0, 0.2),
                    Y2 = 0.6 * eta + rnorm(n, 0, 0.2),
                    Y3 = 0.9 * eta + rnorm(n, 0, 0.2),
                    Y4 = (envir == 2)*0.9 + (0.4+(envir == 2)*0.5) * eta + rnorm(n, 0, 0.2),
                    grp = envir)
fit_taulam = cfa(mod, df_taulam)
eta_taulam = predict(fit_taulam)[,1]

df4 = data.frame(x = predict(lm(df_taulam$Y4 ~ eta_taulam)),
                 y = resid(lm(df$taulam$Y4 ~ eta_taulam)),
                 envir = as.factor(envir))
p4 = ggplot(df4, aes(x = x, y = y, color = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  labs(x = expression(hat(eta)),
       y = expression(hat(epsilon)[i])) + 
  theme_minimal()
p4

#
ggarrange(p1, p2, p3, p4,
          labels = c("A)", "B)", "C)", "D)"),
          ncol = 2, nrow = 2)

#ggsave(here("figures", "Proposal_Motivating2.pdf"), width = 7, height = 5)

# actual effect on false latent mean diff:
summary(lm(eta_MI ~ as.factor(envir)))
summary(lm(eta_lam ~ as.factor(envir)))
summary(lm(eta_tau ~ as.factor(envir)))
summary(lm(eta_taulam ~ as.factor(envir)))



# check lavaan resid function
library(ggplot2)
library(ggpubr)
library(lavaan)
library(here)
set.seed(123)
n = 400
envir = rep(1:2, each = n/2)
eta = c(rnorm(n/2), rnorm(n/2, 0))
mod = "eta =~ Y1 + Y2 + Y3 + Y4"

# full measurement invariance with correlated errors
errorlel = rnorm(n, 0, 1)
df_MI = data.frame(Y1 = eta + rnorm(n, 0, 0.2),
                   Y2 = 0.6 * eta + rnorm(n, 0, 0.2),
                   Y3 = 0.9 * eta + errorlel + rnorm(n, 0, 0.2),
                   Y4 = (envir == 2)*0.9 + (0.4+(envir == 2)*0.5)*eta + errorlel + rnorm(n, 0, 0.2))
fit_MI = cfa(mod, df_MI)
residuals(fit_MI)
eta_hat = predict(fit_MI)[,1]
res1 = matrix(c(residuals(lm(df_MI$Y1 ~ eta_hat)),
                residuals(lm(df_MI$Y2 ~ eta_hat)),
                residuals(lm(df_MI$Y3 ~ eta_hat)),
                residuals(lm(df_MI$Y4 ~ eta_hat))), 
              nrow = n)
cor(res1)
plot(residuals(fit_MI, "obs")[,1], res1[,1], col = envir)
plot(residuals(fit_MI, "obs")[,2], res1[,2], col = envir)
plot(residuals(fit_MI, "obs")[,3], res1[,3], col = envir)
plot(residuals(fit_MI, "obs")[,4], res1[,4], col = envir)
plot(eta_hat, residuals(fit_MI, "obs")[,3], col = envir)
plot(eta_hat, res1[,3], col = envir)
plot(eta_hat, residuals(fit_MI, "obs")[,4], col = envir)
plot(eta_hat, res1[,4], col = envir)

fit_MI = cfa("eta =~ Y1 + Y2 + Y3 + Y4
             Y3 ~~ Y4", df_MI)
residuals(fit_MI)
eta_hat = predict(fit_MI)[,1]
res1 = matrix(c(residuals(lm(df_MI$Y1 ~ eta_hat)),
                residuals(lm(df_MI$Y2 ~ eta_hat)),
                residuals(lm(df_MI$Y3 ~ eta_hat)),
                residuals(lm(df_MI$Y4 ~ eta_hat))), 
              nrow = n)
cor(res1)
plot(residuals(fit_MI, "obs")[,1], res1[,1], col = envir)
plot(residuals(fit_MI, "obs")[,2], res1[,2], col = envir)
plot(residuals(fit_MI, "obs")[,3], res1[,3], col = envir)
plot(residuals(fit_MI, "obs")[,4], res1[,4], col = envir)

  
  
eta_MI = predict(fit_MI)[,1]

residuals(fit_MI, type = "raw")
lavResiduals(fit_MI, "obs")
lel = resid(lm(df_MI$Y3 ~ eta_MI))
R = resid(fit_MI, "obs")
plot(R[,3], lel)



df1 = data.frame(x = predict(lm(df_MI$Y4 ~ eta_MI)),
                 y = resid(lm(df_MI$Y4 ~ eta_MI)),
                 envir = as.factor(envir))
p1 = ggplot(df1, aes(x = x, y = y, color = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  labs(x = expression(hat(eta)),
       y = expression(hat(epsilon)[i])) + 
  theme_minimal()
p1





