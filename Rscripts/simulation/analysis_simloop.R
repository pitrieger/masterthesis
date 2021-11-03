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

#####################
#### SENSITIVITY ####
#####################
  #########################
  ## Across all settings ##
  #########################
  out_df_sum = out_df %>% filter(itembias > 0) %>% 
    group_by(method) %>% 
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  out_df_sum = data.frame(method = out_df_sum$method,
                          n = out_df_sum$TP + out_df_sum$TN + out_df_sum$FP + out_df_sum$FN,
                           sensitivity = get_sensspec(out_df_sum)$est,
                           specificity = get_sensspec(out_df_sum, type = "specificity")$est)
  stargazer(out_df_sum, summary = F)
  
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
  
  
  #########################
  ## Function of h, p, k ##
  #########################
  out_df_hpk = out_df %>% filter(itembias > 0) %>%
    group_by(method, h, p, k) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  
  out_df_hpk = cbind(out_df_hpk, get_sensspec(out_df_hpk)[,2:4])
  out_df_hpk$method = factor(out_df_hpk$method, 
                             levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  
  ggplot(out_df_hpk, aes(y = est, ymin = lowerCI, ymax = upperCI,
                         x = h, color = method, shape = method)) + 
    sensspec_layer + 
    facet_grid(rows = vars(k), cols = vars(p),
               labeller = function(x) label_both(x, sep = " = ")) +
    labs(y = "Sensitivity", x = "h", color = "Method", shape = "Method")

#####################
#### SPECIFICITY ####
#####################
  #########################
  ## Function of n, p, k ##
  #########################
  out_df_npk = out_df %>% filter(itembias > 0) %>%
    group_by(method, n, p, k) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  
  out_df_npk = cbind(out_df_npk, get_sensspec(out_df_npk, type = "specificity")[,2:4])
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
  
  out_df_gpk = cbind(out_df_gpk, get_sensspec(out_df_gpk, type = "specificity")[,2:4])
  out_df_gpk$method = factor(out_df_gpk$method, 
                             levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  
  ggplot(out_df_gpk, aes(y = est, ymin = lowerCI, ymax = upperCI,
                         x = g, color = method, shape = method)) + 
    sensspec_layer + 
    facet_grid(rows = vars(k), cols = vars(p),
               labeller = function(x) label_both(x, sep = " = ")) +
    labs(y = "Sensitivity", x = "g", color = "Method", shape = "Method")
  
  
  #########################
  ## Function of h, p, k ##
  #########################
  out_df_hpk = out_df %>% filter(itembias > 0) %>%
    group_by(method, h, p, k) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  
  out_df_hpk = cbind(out_df_hpk, get_sensspec(out_df_hpk, type = "specificity")[,2:4])
  out_df_hpk$method = factor(out_df_hpk$method, 
                             levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  
  ggplot(out_df_hpk, aes(y = est, ymin = lowerCI, ymax = upperCI,
                         x = h, color = method, shape = method)) + 
    sensspec_layer + 
    facet_grid(rows = vars(k), cols = vars(p),
               labeller = function(x) label_both(x, sep = " = ")) +
    labs(y = "Sensitivity", x = "h", color = "Method", shape = "Method")

  
#######################
#### Zero itembias ####
#######################
# Note that k and h are irrelevant in this case
# Also, any positive classification is wrong when itembias == 0
# but the way it is coded, there still are TP & FP classifications, thus 
# compute accuracy as (TN + FN) / (TN + FN + TP + FP), i.e. how many of all
#  were classified as not non-invariant out of all, because all are invariant

################################
#### As function of n, p, g ####
################################
out_df_0bias = out_df %>% filter(itembias == 0) %>%
    group_by(method, n, p, g) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  r = (out_df_0bias$TN + out_df_0bias$FN)
  n = (out_df_0bias$TN + out_df_0bias$FN + out_df_0bias$TP + out_df_0bias$FP)
  out_df_0bias$est = r/n
  lowerF = qf(1-0.05/2, 2*(n-r+1), 2*r)
  lowerF = ifelse(is.nan(lowerF), round(est), lowerF)
  upperF = qf(1-0.05/2, 2*(r+1), 2*(n-r))
  upperF = ifelse(is.nan(upperF), round(est), upperF)
  #CIs
  lowerCI = r / (r + (n-r+1)*lowerF)
  out_df_0bias$lowerCI = ifelse(is.nan(lowerCI), round(est), lowerCI)
  out_df_0bias$upperCI = (r + 1) * upperF / ((n-r) + (r+1)*upperF)
  
      
  out_df_0bias$method = factor(out_df_0bias$method, 
                             levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  
  ggplot(out_df_0bias, aes(y = est, ymin = lowerCI, ymax = upperCI,
                         x = n, color = method, shape = method)) + 
    sensspec_layer + 
    facet_grid(rows = vars(g), cols = vars(p),
               labeller = function(x) label_both(x, sep = " = ")) +
    labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")



  
  
  
  
  
  
  
  
  
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

