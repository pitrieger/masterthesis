library(ggplot2)
library(ggpubr)
library(here)
set.seed(123)
n = 400
envir = rep(1:2, each = n/2)
eta = c(rnorm(n/2), rnorm(n/2, 1))

# full measurement invariance
Y1 = 0.6 * eta + rnorm(n, 0.5)
#plot(eta, Y1, col = envir, xlab = expression(eta), ylab = expression(Y[1]))
#plot(lm(Y1 ~ eta), which = 1, col = envir)
df1 = data.frame(x = predict(lm(Y1 ~ eta)),
                 y = resid(lm(Y1 ~ eta)),
                 envir = as.factor(envir))
p1 = ggplot(df1, aes(x = x, y = y, color = envir, fill = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.4, show.legend = F) +
  geom_hline(yintercept = sapply(unique(df1$envir), function(i) mean(df1$y[df1$envir == i])),
             color = c("black", "orangered3")) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  scale_fill_manual(values = c("black", "orangered3")) +
  scale_x_continuous(breaks = seq(-1, 4, 1)) + 
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1.5)) + 
  labs(x = expression(eta),
       y = expression(r)) + 
  theme_minimal()

# unequal taus
Y2 = (envir == 2)*0.9 + 0.6 * eta + rnorm(n, 0.5)
#plot(eta, Y2, col = envir, xlab = expression(eta), ylab = expression(Y[1]))
df2 = data.frame(x = predict(lm(Y2 ~ eta)),
                 y = resid(lm(Y2 ~ eta)),
                 envir = as.factor(envir))
p2 = ggplot(df2, aes(x = x, y = y, color = envir, fill = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.4, show.legend = F) +
  geom_hline(yintercept = sapply(unique(df2$envir), function(i) mean(df2$y[df2$envir == i])),
             color = c("black", "orangered3")) +  
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  scale_fill_manual(values = c("black", "orangered3")) +
  scale_x_continuous(breaks = seq(-1, 4, 1)) + 
  scale_y_continuous(breaks = seq(-3, 3, 1.5)) + 
  labs(x = expression(eta),
       y = expression(r)) + 
  theme_minimal()
p2

# unequal lambdas
Y3 = (0.3+(envir == 2)*0.5) * eta + rnorm(n, 0.5)
#plot(eta, Y3, col = envir, xlab = expression(eta), ylab = expression(Y[3]))
df3 = data.frame(x = predict(lm(Y3 ~ eta)),
                 y = resid(lm(Y3 ~ eta)),
                 envir = as.factor(envir))
p3 = ggplot(df3, aes(x = x, y = y, color = envir, fill = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.4, show.legend = F) +
  geom_hline(yintercept = sapply(unique(df3$envir), function(i) mean(df3$y[df3$envir == i])),
             color = c("black", "orangered3")) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  scale_fill_manual(values = c("black", "orangered3")) +
  scale_x_continuous(breaks = seq(-1, 4, 1)) + 
  scale_y_continuous(breaks = seq(-3, 3, 1.5)) + 
  labs(x = expression(eta),
       y = expression(r)) + 
  theme_minimal()
p3

# unequal lambdas & taus
Y4 = (envir == 2)*0.9 + (0.3+(envir == 2)*0.5) * eta + rnorm(n, 0.5)
#plot(eta, Y4, col = envir, xlab = expression(eta), ylab = expression(Y[4]))
df4 = data.frame(x = predict(lm(Y4 ~ eta)),
                 y = resid(lm(Y4 ~ eta)),
                 envir = as.factor(envir))
p4 = ggplot(df4, aes(x = x, y = y, color = envir, fill = envir)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.4, show.legend = F) +
  geom_hline(yintercept = sapply(unique(df1$envir), function(i) mean(df1$y[df1$envir == i])),
             color = c("black", "orangered3")) +
  geom_smooth(method = "lm", show.legend = F) + 
  scale_color_manual(values = c("black", "orangered3")) +
  scale_fill_manual(values = c("black", "orangered3")) +
  scale_x_continuous(breaks = seq(-1, 4, 1)) + 
  scale_y_continuous(breaks = seq(-3, 3, 1.5)) + 
  labs(x = expression(eta),
       y = expression(r)) + 
  theme_minimal()

# save
ggarrange(p1, p2, p3, p4,
          labels = c("A)", "B)", "C)", "D)"),
          ncol = 2, nrow = 2, align = "hv")

ggsave(here("figures", "Proposal_Motivating2.pdf"), width = 7, height = 5)
ggsave(paste0("~/Dropbox/Apps/Overleaf/MA_MeasurementEquivalence/figures/Proposal_Motivating2.pdf"), width = 7, height = 5)
