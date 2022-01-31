# Preliminatory analysis
library(tidyverse)
library(here)
library(stargazer)

mostrecent = list.files(here("data"))[list.files(here("data")) %>% 
  str_extract(., "\\d{4}-\\d{2}-\\d{2}_\\d{4}") %>%
  as.Date(., format = "%Y-%m-%d_%H%M") %>% which.max]
load(here("data", mostrecent))

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

ggsavewrap = function(x, ...){
  ggsave(here("figures", x), ...)
  ggsave(paste0("~/Dropbox/Apps/Overleaf/MA_MeasurementEquivalence/figures/", x), ...)
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
  # with itembias
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
  out_df_sum$method = factor(out_df_sum$method, 
                             levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  out_df_sum = out_df_sum %>% arrange(method)
  stargazer(out_df_sum, summary = F)
  
  # without itembias
  out_df_sum0 = out_df %>% filter(itembias == 0) %>% 
    group_by(method) %>% 
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  out_df_sum0 = data.frame(method = out_df_sum0$method,
                          n = out_df_sum0$TP + out_df_sum0$TN + out_df_sum0$FP + out_df_sum0$FN,
                          specificity = get_sensspec(out_df_sum0, type = "specificity")$est)
  out_df_sum0$method = factor(out_df_sum0$method, 
                             levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  out_df_sum0 = out_df_sum0 %>% arrange(method)
  stargazer(out_df_sum0, summary = F)
  
  # identify which parameter settings lead BV to yield NAs
  out_df[which(rowSums(out_df[,1:4]) != out_df$nsim*out_df$p),]
  
  
  #cbind and print as tex file
#  out_df_sumtex = left_join(out_df_sum, out_df_sum0, by = "method") 
#  out_df_sumtex$method = factor(out_df_sumtex$method, 
#                                            levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
#                                            labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
#  body = cat(paste(apply(out_df_sumtex, 1, function(x) paste(x, collapse = " & ")), collapse = "\\\\ \n"))
#  paste(c("\\\begin{table}[!ht]\n
#    \\\centering\n
#    \\\begin{tabular}{@{\\\extracolsep{5pt}}llllll}\n
#\\\hline\n
# & \\\multicolumn{3}{c}{$\\\delta > 0$} & \\\multicolumn{2}{c}{$\\\delta = 0$}\\\\ \n
#\\\cline{2-4} \\\cline{5-6} \\\\[-1.8ex] method & n & Sensitivity & Specificity & n & Specificity \\\\ \\\hline \\\hline", body, "x"))
#  colnames(out_df_sumtex) = c("Method", "n", "Sensitivity", "Specificity", "n", "Specificity")

  
  #########################
  ## Function of n, p, k ##
  #########################
  out_df_npk = out_df %>% filter(itembias > 0) %>%
    group_by(method, k) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  
  out_df_npk = cbind(out_df_npk, get_sensspec(out_df_npk)[,2:4])
  out_df_npk$method = factor(out_df_npk$method, 
                                  levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                                  labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  
  {ggplot(out_df_npk, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = k, color = method, shape = method)) + 
    sensspec_layer + 
    labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")} #%>%
    #ggplot2::ggplot_build() %>% ggplot2::ggplot_gtable() %>%
    #gtable::gtable_filter(pattern = "panel-2-1$", invert = TRUE) %>%
    #gtable::gtable_filter(pattern = "panel-2-2$", invert = TRUE) %>%
    #{grid::grid.newpage(); grid::grid.draw(.)}
  ggsavewrap("Sensitivity_linegrid_npk.pdf", width = 8, height = 5)
  
  
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
  ggsavewrap("Sensitivity_linegrid_gpk.pdf", width = 8, height = 5)
  
  #########################
  ## Function of h, p, k ##
  #########################
  #ut_df_hpk = out_df %>% filter(itembias > 0) %>%
  # group_by(method, h, p, k) %>%
  # summarize(TP = sum(TP, na.rm = T),
  #           TN = sum(TN, na.rm = T),
  #           FP = sum(FP, na.rm = T),
  #           FN = sum(FN, na.rm = T))
  
  #ut_df_hpk = cbind(out_df_hpk, get_sensspec(out_df_hpk)[,2:4])
  #ut_df_hpk$method = factor(out_df_hpk$method, 
  #                          levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
  #                          labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  
  #gplot(out_df_hpk, aes(y = est, ymin = lowerCI, ymax = upperCI,
  #                      x = h, color = method, shape = method)) + 
  # sensspec_layer + 
  # facet_grid(rows = vars(k), cols = vars(p),
  #            labeller = function(x) label_both(x, sep = " = ")) +
  # labs(y = "Sensitivity", x = "h", color = "Method", shape = "Method")

  
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
    labs(y = "Specificity", x = "n", color = "Method", shape = "Method")
  ggsavewrap("Specificity_linegrid_npk.pdf", width = 8, height = 5)
  
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
    labs(y = "Specificity", x = "g", color = "Method", shape = "Method")
  ggsavewrap("Specificity_linegrid_gpk.pdf", width = 8, height = 5)
  
  #########################
  ## Function of h, p, k ##
  #########################
  #out_df_hpk = out_df %>% filter(itembias > 0) %>%
  #  group_by(method, h, p, k) %>%
  #  summarize(TP = sum(TP, na.rm = T),
  #            TN = sum(TN, na.rm = T),
  #            FP = sum(FP, na.rm = T),
  #            FN = sum(FN, na.rm = T))
  #
  #out_df_hpk = cbind(out_df_hpk, get_sensspec(out_df_hpk, type = "specificity")[,2:4])
  #out_df_hpk$method = factor(out_df_hpk$method, 
  #                           levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
  #                           labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
  #
  #ggplot(out_df_hpk, aes(y = est, ymin = lowerCI, ymax = upperCI,
  #                       x = h, color = method, shape = method)) + 
  #  sensspec_layer + 
  #  facet_grid(rows = vars(k), cols = vars(p),
  #             labeller = function(x) label_both(x, sep = " = ")) +
  #  labs(y = "Specificity", x = "h", color = "Method", shape = "Method")

  
#######################
#### Zero itembias ####
#######################
# Note that k and h are irrelevant in this case
# Note that sensitivity is 
  
## SPECIFICITY
################################
#### As function of n, p, g ####
################################
  out_df_0bias = out_df %>% filter(itembias == 0) %>%
    group_by(method, n, p, g) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  out_df_0bias = cbind(out_df_0bias, get_sensspec(out_df_0bias, type = "specificity")[,2:4])  
  out_df_0bias$method = factor(out_df_0bias$method, 
                               levels = c("ni_J", "ni_M", "ni_C", "ni_B", "ni_R1", "ni_R2"),
                               labels = c("J", "MInd", "CR", "BV", "R1", "R2"))
 
  ggplot(out_df_0bias, aes(y = est, ymin = lowerCI, ymax = upperCI,
                         x = n, color = method, shape = method)) + 
    sensspec_layer + 
    facet_grid(cols = vars(p), rows = vars(g),
               labeller = function(x) label_both(x, sep = " = ")) +
    labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method")
  ggsavewrap("Specificity_linegrid_0bias.pdf", width = 8, height = 10)
  
  
  
  
  
  
  
  
  
  
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

