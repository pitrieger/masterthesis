library(tidyverse)
library(here)
library(stargazer)
library(grid)
library(ggpubr)

mostrecent = list.files(here("publicdata"))[list.files(here("publicdata")) %>% 
  str_extract(., "\\d{4}-\\d{2}-\\d{2}_\\d{4}") %>%
  as.Date(., format = "%Y-%m-%d_%H%M") %>% which.max]

load(here("publicdata", mostrecent))

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
  # get F-values, but use 0 or 1 if est is 0 or 1 to fix NaNs
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
                      #scale_color_discrete(na.translate = F),
                      scale_color_brewer(na.translate = F, type = "qual", palette = 2),
                      #scale_color_manual(na.translate = F, values = c("#1B9E77","#7570B3","#66A61E","#E6AB02","#A6761D","#666666")),
                      scale_shape_manual(na.translate = F, values = c(15, 17, 19, 5, 4, 8)),
                      theme_bw(),
                      theme(strip.background = element_rect(color = "black", fill = "white")))

# bind all confusion matrix info together in dataframe
out_df = lapply(out, function(l) apply(l, 1:2, function(x) sum(x, na.rm = T)))
out_df = do.call("rbind", out_df) %>% 
  as.data.frame
out_df$method = rep(rownames(out$`1`), length(out))
out_df$id = rep(1:length(out), each = nrow(out$`1`))
out_df = left_join(out_df, sim_param_df %>% mutate(id = 1:nrow(sim_param_df)))

# check number where trycatch was actually used
out_df %>% 
  mutate(numba = TP + FP + TN + FN,
         diff = numba - p*nsim) %>% 
  group_by(method) %>% 
  summarize(t = min(diff)) # 0

#####################
#### SENSITIVITY ####
#####################
  #########################
  ## Across all settings ##
  #########################
  # with itembias
  out_df_sum = out_df %>% filter(interceptbias > 0 & loadingbias > 0) %>% 
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
                             levels = c("ni_J", "ni_M", "ni_M_bonf", "ni_C", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "MInd (Bonf.)", "CR", "CR (Bonf.)", "BV", "R1", "R2"))
  out_df_sum_both = out_df_sum %>% arrange(method)
  stargazer(out_df_sum_both, summary = F)
  
  # with just interceptbias
  out_df_sum = out_df %>% filter(interceptbias > 0 & loadingbias == 0) %>% 
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
                             levels = c("ni_J", "ni_M", "ni_M_bonf", "ni_C", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "MInd (Bonf.)", "CR", "CR (Bonf.)", "BV", "R1", "R2"))
  out_df_sum_intercept = out_df_sum %>% arrange(method)
  stargazer(out_df_sum_intercept, summary = F)
  
  # with just loadingbias
  out_df_sum = out_df %>% filter(interceptbias == 0 & loadingbias > 0) %>% 
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
                             levels = c("ni_J", "ni_M", "ni_M_bonf", "ni_C", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "MInd (Bonf.)", "CR", "CR (Bonf.)", "BV", "R1", "R2"))
  out_df_sum_loading = out_df_sum %>% arrange(method)
  stargazer(out_df_sum_loading, summary = F)
  
  # without itembias
  out_df_sum = out_df %>% filter(interceptbias == 0 & loadingbias == 0) %>% 
    group_by(method) %>% 
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  out_df_sum = data.frame(method = out_df_sum$method,
                          n = out_df_sum$TP + out_df_sum$TN + out_df_sum$FP + out_df_sum$FN,
                          specificity = get_sensspec(out_df_sum, type = "specificity")$est)
  out_df_sum$method = factor(out_df_sum$method, 
                             levels = c("ni_J", "ni_M", "ni_M_bonf", "ni_C", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                             labels = c("J", "MInd", "MInd (Bonf.)", "CR", "CR (Bonf.)", "BV", "R1", "R2"))
  out_df_sum0 = out_df_sum %>% arrange(method)
  stargazer(out_df_sum0, summary = F)
  
  # combined 
    # first 2 small checks
    sapply(list(out_df_sum_both$n, out_df_sum_loading$n, out_df_sum_intercept$n, out_df_sum0$n), 
           function(x) identical(x, out_df_sum_both$n)) 
    sapply(list(out_df_sum_both$method, out_df_sum_loading$method, out_df_sum_intercept$method, out_df_sum0$method), 
           function(x) identical(x, out_df_sum_both$method))
    #bind
    out_df_sum_all = data.frame(method = out_df_sum_both$method,
                                n = out_df_sum_intercept$n,
                                Sens_both = out_df_sum_both$sensitivity,
                                Spec_both = out_df_sum_both$specificity,
                                Sens_intercept = out_df_sum_intercept$sensitivity,
                                Spec_intercept = out_df_sum_intercept$specificity,
                                Sens_loading = out_df_sum_loading$sensitivity,
                                Spec_loading = out_df_sum_loading$specificity,
                                Spec_0 = out_df_sum0$specificity)
    stargazer(out_df_sum_all, summary = F)

    #missing classifications due to errors:
    min(c(out_df_sum_both$n, out_df_sum_intercept$n, out_df_sum_loading$n, out_df_sum0$n))
    c(out_df_sum_both$n, out_df_sum_intercept$n, out_df_sum_loading$n, out_df_sum0$n) - 109200
    8*sum(sim_param_df$p * sim_param_df$nsim) - sum(c(out_df_sum_both$n, out_df_sum_intercept$n, out_df_sum_loading$n, out_df_sum0$n))
  
    
    # ROC plot version of table
    p1 = ggplot(out_df_sum_both, aes(y = sensitivity, x = 1-specificity, color = method, shape = method)) + 
      geom_point(size = 2.5) +
      geom_abline(intercept = 0, slope = 1, color = "grey25") + 
      coord_fixed() +
      scale_shape_manual(values = c(15, 2, 17, 1, 19, 5, 4, 8)) + 
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) + 
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) + 
      scale_color_manual(values = c(c("#1B9E77", "#D95F02", "#D95F02", "#7570B3", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"))) +
      theme_bw() +
      labs(x = "1 - Specificity", y = "Sensitivity",
           title = "Simultaneous Violation") +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    p2 = ggplot(out_df_sum_intercept, aes(y = sensitivity, x = 1-specificity, color = method, shape = method)) + 
      geom_point(size = 2.5) +
      geom_abline(intercept = 0, slope = 1, color = "grey25") + 
      coord_fixed() +
      scale_shape_manual(values = c(15, 2, 17, 1, 19, 5, 4, 8)) + 
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) + 
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) + 
      scale_color_manual(values = c(c("#1B9E77", "#D95F02", "#D95F02", "#7570B3", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"))) +
      theme_bw() +
      labs(x = "1 - Specificity", y = "Sensitivity",
           title = "Violation of Scalar MI") +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    p3 = ggplot(out_df_sum_loading, aes(y = sensitivity, x = 1-specificity, color = method, shape = method)) + 
      geom_point(size = 2.5) +
      geom_abline(intercept = 0, slope = 1, color = "grey25") + 
      coord_fixed() +
      scale_shape_manual(values = c(15, 2, 17, 1, 19, 5, 4, 8)) + 
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) + 
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) + 
      scale_color_manual(values = c(c("#1B9E77", "#D95F02", "#D95F02", "#7570B3", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"))) +
      theme_bw() +
      labs(x = "1 - Specificity", y = "Sensitivity",
           title = "Violation of Metric MI") +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    out_df_sum0$method = factor(out_df_sum0$method, levels = rev(levels(out_df_sum0$method)))
    p4 = ggplot(out_df_sum0, aes(y = method, x = 1-specificity, fill = method)) + 
      geom_bar(stat = "identity", show.legend = F) + 
      coord_fixed(ratio = 1/8) + 
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) + 
      scale_fill_manual(values = rev(c("#1B9E77", "#D95F02", "#D95F02", "#7570B3", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"))) +
      theme_bw() +
      labs(y = "1 - Specificity",
           title = "Perfect MI") +
      theme(legend.title = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))
      
    p1234 = ggarrange(p1 + theme(axis.title.y = element_text(vjust = -13)), 
                      p2 + theme(axis.title.y = element_text(vjust = -13)), 
                      p3 + theme(axis.title.y = element_text(vjust = -13)), 
                      p4,
                      labels = c("A)", "B)", "C)", "D)"),
              ncol = 2, nrow = 2, 
              common.legend = T, legend = "right",
              align = "hv")
    p1234
    ggsavewrap("ROC_simultaneousMetricScalar.pdf", width = 8, height = 5)
    
      
  #########################
  ## Function of n, p, k ##
  #########################
  out_df_npk = out_df %>% filter(interceptbias > 0 & loadingbias > 0) %>%
    group_by(method, n, p, k) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  
  out_df_npk_sens = cbind(out_df_npk, get_sensspec(out_df_npk)[,2:4])
  out_df_npk_sens$method = factor(out_df_npk_sens$method, 
                                  levels = c("ni_J", "ni_M_bonf", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                                  labels = c("J", "MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"))
  out_df_npk_sens$m = out_df_npk_sens$k
  
  ggplot(out_df_npk_sens, aes(y = est, ymin = lowerCI, ymax = upperCI,
                      x = n, color = method, shape = method)) + 
    sensspec_layer + 
    facet_grid(rows = vars(m), cols = vars(p),
               labeller = function(x) label_both(x, sep = " = ")) +
    scale_x_continuous(breaks = c(100, 200, 500, 1000)) +
    labs(y = "Sensitivity", x = "n", color = "Method", shape = "Method") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1.1, vjust = 1.2))
  
  ggsavewrap("Sensitivity_linegrid_npk.pdf", width = 8, height = 5)
  
  #######################
  ## Function of h & g ##
  #######################
  out_df_gh = out_df %>% 
    filter(interceptbias > 0 & loadingbias > 0) %>% 
    group_by(method, h, g) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  
  out_df_gh_sens = cbind(out_df_gh, get_sensspec(out_df_gh)[,2:4])
  out_df_gh_sens$method = factor(out_df_gh_sens$method, 
                                 levels = c("ni_J", "ni_M_bonf", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                                 labels = c("J", "MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"))
  
  ggplot(out_df_gh_sens, aes(y = est, ymin = lowerCI, ymax = upperCI,
                       x = g, color = method, shape = method)) + 
    facet_wrap(~h, 
               labeller = function(x) label_both(x, sep = " = ")) + 
    sensspec_layer + 
    scale_x_continuous(breaks = 2^(1:4)) +
    labs(y = "Sensitivity", x = "g", color = "Method", shape = "Method")
  ggsavewrap("Sensitivity_linegrid_gh.pdf", width = 8, height = 2)
  
#####################
#### SPECIFICITY ####
#####################
  #########################
  ## Function of n, p, k ##
  #########################
  out_df_npk_spec = cbind(out_df_npk, get_sensspec(out_df_npk, type = "specificity")[,2:4])
  out_df_npk_spec$method = factor(out_df_npk_spec$method, 
                                  levels = c("ni_J", "ni_M_bonf", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                                  labels = c("J", "MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"))
  out_df_npk_spec$m = out_df_npk_spec$k
  
  ggplot(out_df_npk_spec, aes(y = est, ymin = lowerCI, ymax = upperCI,
                         x = n, color = method, shape = method)) + 
    sensspec_layer + 
    facet_grid(rows = vars(m), cols = vars(p),
               labeller = function(x) label_both(x, sep = " = ")) +
    scale_x_continuous(breaks = c(100, 200, 500, 1000)) +
    labs(y = "Specificity", x = "n", color = "Method", shape = "Method") + 
    theme(axis.text.x = element_text(angle = 40, hjust = 1.1, vjust = 1.2))
  
  ggsavewrap("Specificity_linegrid_npk.pdf", width = 8, height = 5)
  
  #######################
  ## Function of h & g ##
  #######################
  out_df_gh_spec = cbind(out_df_gh, get_sensspec(out_df_gh, type = "specificity")[,2:4])
  out_df_gh_spec$method = factor(out_df_gh_spec$method, 
                                 levels = c("ni_J", "ni_M_bonf", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                                 labels = c("J", "MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"))
  
  ggplot(out_df_gh_spec, aes(y = est, ymin = lowerCI, ymax = upperCI,
                        x = g, color = method, shape = method)) + 
    facet_wrap(~h, 
               labeller = function(x) label_both(x, sep = " = ")) + 
    sensspec_layer + 
    scale_x_continuous(breaks = 2^(1:4)) +
    labs(y = "Specificity", x = "g", color = "Method", shape = "Method")
  ggsavewrap("Specificity_linegrid_gh.pdf", width = 8, height = 2)
  

#######################
#### Zero itembias ####
#######################
# Note that k and h are irrelevant in this case
# Note that sensitivity is 
  
## SPECIFICITY
################################
#### As function of n, p, g ####
################################
  out_df_0bias = out_df %>% filter(interceptbias == 0 & loadingbias == 0) %>%
    group_by(method, n) %>%
    summarize(TP = sum(TP, na.rm = T),
              TN = sum(TN, na.rm = T),
              FP = sum(FP, na.rm = T),
              FN = sum(FN, na.rm = T))
  out_df_0bias = cbind(out_df_0bias, get_sensspec(out_df_0bias, type = "specificity")[,2:4])  
  out_df_0bias$method = factor(out_df_0bias$method, 
                               levels = c("ni_J", "ni_M_bonf", "ni_C_bonf", "ni_B", "ni_R1", "ni_R2"),
                               labels = c("J", "MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"))
  
  ggplot(out_df_0bias, aes(y = est, ymin = lowerCI, ymax = upperCI,
                         x = n, color = method, shape = method)) + 
    sensspec_layer + 
  #  facet_grid(cols = vars(p)) + 
    scale_x_continuous(breaks = c(100, 200, 500, 1000)) +
    labs(y = "Specificity", x = "n", color = "Method", shape = "Method")
  ggsavewrap("Specificity_lineplot_0bias.pdf", width = 8, height = 3)
