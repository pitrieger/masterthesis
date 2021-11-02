# Preliminatory analysis
library(tidyverse)
library(here)
library(stargazer)
load(here("data/Prelim_Sim_2021-10-28.1RData"))

get_sensspec = function(x, alpha = 0.05, type = "sensitivity"){
  stopifnot(type %in% c("sensitivity", "specificity"))
  if(type == "sensitivity"){
    r = as.numeric(as.data.frame(x)[,"TP"])
    n = r + as.numeric(as.data.frame(x)[,"FN"])
  } else {
    r = as.numeric(as.data.frame(x)[,"TN"])
    n = r + as.numeric(as.data.frame(x)[,"FP"])
  } 
  # proportion
  est = r/n
  # Clopper-Pearson CIs
  # get F-values, but use 0 or 1 if est is 0 or 1 
  lowerF = qf(1-alpha/2, 2*(n-r+1), 2*r)
  lowerF = ifelse(is.nan(lowerF), round(est), lowerF)
  upperF = qf(1-alpha/2, 2*(r+1), 2*(n-r))
  upperF = ifelse(is.nan(upperF), round(est), upperF)
  #CIs
  lowerCI = r / (r + (n-r+1)*lowerF)
  lowerCI = ifelse(is.nan(lowerCI), round(est), lowerCI)
  upperCI = (r + 1) * upperF / ((n-r) + (r+1)*upperF)
  data.frame(method = rownames(x), est, lowerCI, upperCI)
}

sensspec_layer = list(geom_pointrange(alpha = 0.7),
                      geom_line(show.legend = F, alpha = 0.7),
                      scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)),
                      scale_color_brewer(type = "qual", palette = 2),
                      scale_shape_manual(values = c(0, 1, 2, 15, 16, 17)),
                      theme_bw(),
                      theme(strip.background = element_rect(color = "black", fill = "white")))

# bind all confusion matrix info together in dataframe
out_df = lapply(out, function(l) apply(l, 1:2, function(x) sum(x, na.rm = T)))
out_df = do.call("rbind", out_df) %>% 
  as.data.frame
out_df$method = rep(rownames(out$`1`), length(out))
out_df$id = rep(1:length(out), each = nrow(out$`1`))
out_df = left_join(out_df, sim_param_df %>% mutate(id = 1:nrow(sim_param_df)))




#########################
## Function of n, p, k ##
#########################
out_df_npk = out_df %>% filter(itembias > 0) %>%
  group_by(method, n, p, k) %>%
  summarize(TP = sum(TP, na.rm = T),
            TN = sum(TN, na.rm = T),
            FP = sum(FP, na.rm = T),
            FN = sum(FN, na.rm = T))

out_df_npk = cbind(out_df_npk, get_sensspec(out_df_npk)[,2:4])
out_df_npk$method = factor(out_df_npk$method, 
                                levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                                labels = c("J", "MInd", "CR", "BV", "R1", "R2"))

ggplot(out_df_npk, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer + 
  facet_grid(rows = vars(k), cols = vars(p),
             labeller = function(x) label_both(x, sep = " = ")) +
  labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")


#########################
## Function of g, p, k ##
#########################
out_df_gpk = out_df %>% filter(itembias > 0) %>%
  group_by(method, g, p, k) %>%
  summarize(TP = sum(TP, na.rm = T),
            TN = sum(TN, na.rm = T),
            FP = sum(FP, na.rm = T),
            FN = sum(FN, na.rm = T))

out_df_gpk = cbind(out_df_gpk, get_sensspec(out_df_gpk)[,2:4])
out_df_gpk$method = factor(out_df_gpk$method, 
                           levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                           labels = c("J", "MInd", "CR", "BV", "R1", "R2"))

ggplot(out_df_gpk, aes(y = est, ymin = lowerCI, ymax = upperCI,
                       x = g, color = method, shape = method)) + 
  sensspec_layer + 
  facet_grid(rows = vars(k), cols = vars(p),
             labeller = function(x) label_both(x, sep = " = ")) +
  labs(y = "Sensitivity", x = "g", color = "Method", shape = "Method")












# compute sensitivity
out_sensitivity = do.call("rbind", 
                          lapply(lapply(out, function(l) apply(l, 1:2, function(x) sum(x, na.rm = T))), #sum confusion matrices over replications 
                                 function(x) get_sensspec(x))) # compute sensitivity
out_sensitivity = cbind(out_sensitivity, sim_param_df[rep(1:nrow(sim_param_df), each = 6),]) # add simulation parameters
# relabel method names
out_sensitivity$method = factor(out_sensitivity$method, 
                                levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                                labels = c("J", "MInd", "CR", "BV", "R1", "R2"))

## SENSITIVITY
# theme
sensspec_layer = list(geom_pointrange(alpha = 0.7),
                      geom_line(show.legend = F, alpha = 0.7),
                      facet_grid(rows = vars(k), cols = vars(p),
                                 labeller = function(x) label_both(x, sep = " = ")),
                      scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)),
                      scale_color_brewer(type = "qual", palette = 2),
                      scale_shape_manual(values = c(0, 1, 2, 15, 16, 17)),
                      theme_bw(),
                      theme(strip.background = element_rect(color = "black", fill = "white")))

## SENSITIVITY
# Plots for itembias
# as function of n
df_p1.1 = out_sensitivity %>% filter(itembias == 0.5 & g == 16 & h == 0.5)
ggplot(df_p1.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer + 
  labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")

# as function of g
df_p1.1 = out_sensitivity %>% filter(itembias == 0.25 & n == 500 & h == 0.5)
ggplot(df_p1.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = g, color = method, shape = method)) + 
  sensspec_layer + 
  labs(y = "Sensitivity", x = "g", color = "Method", shape = "Method")






fitJ = lm(est ~ as.factor(n) + as.factor(p) + as.factor(k) + as.factor(h) + as.factor(g), data = out_sensitivity, subset = (out_sensitivity$method == "J"))
fitMInd = lm(est ~ as.factor(n) + as.factor(p) + as.factor(k) + as.factor(h) + as.factor(g), data = out_sensitivity, subset = (out_sensitivity$method == "MInd"))
fitCR = lm(est ~ as.factor(n) + as.factor(p) + as.factor(k) + as.factor(h) + as.factor(g), data = out_sensitivity, subset = (out_sensitivity$method == "CR"))
fitBV = lm(est ~ as.factor(n) + as.factor(p) + as.factor(k) + as.factor(h) + as.factor(g), data = out_sensitivity, subset = (out_sensitivity$method == "BV"))
fitR1 = lm(est ~ as.factor(n) + as.factor(p) + as.factor(k) + as.factor(h) + as.factor(g), data = out_sensitivity, subset = (out_sensitivity$method == "R1"))
fitR2 = lm(est ~ as.factor(n) + as.factor(p) + as.factor(k) + as.factor(h) + as.factor(g), data = out_sensitivity, subset = (out_sensitivity$method == "R2"))
summary(fitJ)
stargazer(fitJ, fitMInd, fitCR, fitBV, fitR1, fitR2, column.labels = c("J", "MInd", "CR", "BV", "R1", "R2"))

