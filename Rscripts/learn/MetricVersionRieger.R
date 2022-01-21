


ibrary(ggplot2)
library(ggpubr)
library(here)
set.seed(123)
n = 600
g = 3
envir = rep(1:3, each = n/3)
(eta_mu = rnorm(g-1, 0, 4))
eta = c(rnorm(n/g), rnorm(n/g, eta_mu[1]), rnorm(n/g, eta_mu[2]))

# unequal lambdas
Y3 = (0.3+(envir == 2)*0.5) * eta + rnorm(n, 0.5)
#Y3 = (envir == 2)*2 + (0.3) * eta + rnorm(n, 0.5)
#plot(eta, Y3, col = envir, xlab = expression(eta), ylab = expression(Y[3]))
df3 = data.frame(x = predict(lm(Y3 ~ eta)),
                 y = resid(lm(Y3 ~ eta)),
                 envir = as.factor(envir))
ggplot(df3, aes(x = x, y = y, color = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3", "deepskyblue3")) +
  labs(x = expression(hat(eta)),
       y = expression(hat(epsilon))) + 
  theme_minimal()
p3

summary(lm(Y3 ~ eta*as.factor(envir)))
