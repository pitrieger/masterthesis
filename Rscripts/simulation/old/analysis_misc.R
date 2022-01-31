#ANALYSIS MISC
library(tidyverse)
library(parallel)
library(lavaan)
library(here)
library(RColorBrewer)

load(here("data/Prelim_Sim_2021-11-03.RData"))
lel
# compute sensitivity
out_sensitivity = do.call("rbind", 
                          lapply(lapply(out, function(l) apply(l, 1:2, function(x) sum(x, na.rm = T))), #sum confusion matrices over replications 
                                 function(x) get_sensspec(x))) # compute sensitivity
out_sensitivity = cbind(out_sensitivity, sim_param_df[rep(1:nrow(sim_param_df), each = 6),]) # add simulation parameters

# compute specificity
out_specificity = do.call("rbind", 
                          lapply(lapply(out, function(l) apply(l, 1:2, function(x) sum(x, na.rm = T))), #sum confusion matrices over replications 
                                 function(x) get_sensspec(x, type = "specificity"))) # compute specificity
out_specificity = cbind(out_specificity, sim_param_df[rep(1:nrow(sim_param_df), each = 6),]) # add simulation parameters


## Generate plots
# relabel method names
out_sensitivity$method = factor(out_sensitivity$method, 
                                levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                                labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
out_specificity$method = factor(out_specificity$method, 
                                levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                                labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
# theme
sensspec_layer = list(geom_pointrange(),
                      geom_line(show.legend = F),
                      facet_grid(rows = vars(k), cols = vars(p),
                                 labeller = function(x) label_both(x, sep = " = ")),
                      scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)),
                      scale_color_brewer(type = "qual", palette = 2),
                      theme_bw(),
                      theme(strip.background = element_rect(color = "black", fill = "white")))

## SENSITIVITY
# Plots for Loadingbias
df_p1.1 = out_sensitivity %>% filter(loadingbias == 0.2 &
                                       interceptbias == 0)
ggplot(df_p1.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer + 
  labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")

# Plots for Interceptbias
df_p2.1 = out_sensitivity %>% filter(loadingbias == 0 &
                                       interceptbias == 0.2)
ggplot(df_p2.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer + 
  labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")

# Plots for Intercept & Loadingbias
df_p3.1 = out_sensitivity %>% filter(loadingbias == 0.2 &
                                       interceptbias == 0.2)
ggplot(df_p3.1, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer + 
  labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")

## SPECIFICITY
# Plots for Loadingbias
df_p1.2 = out_specificity %>% filter(loadingbias == 0.2 &
                                       interceptbias == 0)
ggplot(df_p1.2, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer +
  labs(y = "Specificity", x = "n", color = "Method", shape = "Method")

# Plots for Interceptbias
df_p2.2 = out_specificity %>% filter(loadingbias == 0 &
                                       interceptbias == 0.2)
ggplot(df_p2.2, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer +
  labs(y = "Specificity", x = "n", color = "Method", shape = "Method")


# Plots for Intercept & Loadingbias
df_p3.2 = out_specificity %>% filter(loadingbias == 0.2 &
                                       interceptbias == 0.2)
ggplot(df_p3.2, aes(y = est, ymin = lowerCI, ymax = upperCI,
                    x = n, color = method, shape = method)) + 
  sensspec_layer +
  labs(y = "Specificity", x = "n", color = "Method", shape = "Method")
