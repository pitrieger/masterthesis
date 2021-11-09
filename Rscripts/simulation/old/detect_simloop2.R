library(tidyverse)
library(parallel)
library(lavaan)
library(here)
source(here("Rscripts/simulation", "Simulator_PokropekEtAl.R"))
source(here("Rscripts/simulation", "Detector_Rieger.R"))
source(here("Rscripts/simulation", "Detector_ByrneVandeVijer.R"))
source(here("Rscripts/simulation", "Detector_CheungRensvold.R"))
#MMGFA_files = list.files(here("Rscripts", "KimDeRoover_MixtureMG_FA"))
#MMGFA_files = MMGFA_files[grepl("R$", MMGFA_files)]
#for(j in 1:length(MMGFA_files)){
#  source(here("Rscripts", "KimDeRoover_MixtureMG_FA", MMGFA_files[j]))
#}

confusion_metric = function()

# parallel RNG seed
RNGkind("L'Ecuyer-CMRG")
set.seed(9)

options(mc.cores = max(detectCores() - 1, 2))

# Grid
sim_param = expand.grid(nsim = 100,
                        n = c(200, 400, 850),
                        p = c(3, 4, 5),
                        g = c(2),
                        h = c(0, 0.5),
                        k = c(1, 2), 
                        interceptbias = c(0, 0.2),
                        loadingbias = c(0, 0.2))
sim_param = sim_param[sim_param$k<sim_param$p,]
sim_param_df = sim_param[sim_param$interceptbias != 0 | sim_param$loadingbias != 0,]
sim_param = split(sim_param_df, 1:nrow(sim_param_df))
sim_out = lapply(sim_param, function(x) set_names(list(x), "sim_param"))

run_sim = function(x) {
  pars = x$sim_param
  varnames = paste0("y", 1:pars$p)
  contain_correct = count_identified = matrix(nrow = pars$nsim, ncol = 4)
  colnames(contain_correct) = colnames(count_identified) = c("Rieger1", "Rieger2", "Byrne", "Cheung")
  
  for(s in 1:pars$nsim){
    # Simulate data
    D = sim_PMI(n = pars$n,
                g = pars$g,
                p = pars$p,
                h = pars$h,
                k = pars$k)
    noninvariant_true = ifelse(length(D$p_affected)>0, paste0("y", D$p_affected), character())
    
    # Own Idea v1: simple
    noninvariant_R1 = try(detect_Rieger_v1(varnames = varnames,
                                           data = D$sim_dat)$noninvariant)
    if(class(noninvariant_R1)!="try-error"){
      count_identified[s,1] = length(noninvariant_R1)
      contain_correct[s,1] = sum(noninvariant_true %in% noninvariant_R1)
    }
    
    # Own Idea v2: stepwise
    noninvariant_R2 = try(detect_Rieger_v2(varnames = varnames,
                                           D$sim_dat)$noninvariant)
    if(class(noninvariant_R2)!="try-error"){
      count_identified[s,2] = length(noninvariant_R2)
      contain_correct[s,2] = sum(noninvariant_true %in% noninvariant_R2)
    }
    
    ## Via CFI
    # Byrne & van de Vijer (2010):
    noninvariant_B = try(detect_ByrneVandeVijer(varnames = varnames,
                                                D$sim_dat)$noninvariant)
    if(class(noninvariant_B)!="try-error"){
      count_identified[s,3] = length(noninvariant_B)
      contain_correct[s,3] = sum(noninvariant_true %in% noninvariant_B)
    }
    
    # Cheung & Rensveld (1999):
    noninvariant_C = try(detect_CheungRensvold(varnames = varnames,
                                               D$sim_dat)$noninvariant)
    if(class(noninvariant_C)!="try-error"){
      count_identified[s,4] = length(noninvariant_C)
      contain_correct[s,4] = sum(noninvariant_true %in% noninvariant_C)
    }
  }
  #c(x,
  #  list(count_correct = count_correct,
  #       contain_correct = contain_correct))
  list(sim_pars = pars, 
       count_identified = count_identified,
       contain_correct = contain_correct)
}

system.time(out <- mclapply(sim_out, function(x) run_sim(x)))
#save.image("temp.RData")
load("temp.RData")


sim_param_df = lapply(out, function(x) apply(x$contain_correct, 2, function(y) prop.test(sum(y, na.rm = T), x$sim_pars$k*(x$sim_pars$n - sum(is.na(y))), correct = F)$estimate)) %>%
  unlist() %>%
  matrix(., nrow = 108, ncol = 4, byrow = T) %>% 
  as.data.frame() %>% 
  rename(., R1_est = V1, R2_est = V2, BvdV_est = V3, CR_est = V4) %>%
  cbind(sim_param_df, .)

sim_param_df = lapply(out, function(x) apply(x$contain_correct, 2, function(y) prop.test(sum(y, na.rm = T), x$sim_pars$k*(x$sim_pars$n - sum(is.na(y))), correct = F)$conf.int[1])) %>%
  unlist() %>%
  matrix(., nrow = 108, ncol = 4, byrow = T) %>% 
  as.data.frame() %>% 
  rename(., R1_lower = V1, R2_lower = V2, BvdV_lower = V3, CR_lower = V4) %>%
  cbind(sim_param_df, .)

sim_param_df = lapply(out, function(x) apply(x$contain_correct, 2, function(y) prop.test(sum(y, na.rm = T), x$sim_pars$k*(x$sim_pars$n - sum(is.na(y))), correct = F)$conf.int[2])) %>%
  unlist() %>%
  matrix(., nrow = 108, ncol = 4, byrow = T) %>% 
  as.data.frame() %>% 
  rename(., R1_upper = V1, R2_upper = V2, BvdV_upper = V3, CR_upper = V4) %>%
  cbind(sim_param_df, .)

lel = sim_param_df %>% 
  filter(h == 0.5 &  k == 1 & interceptbias == 0.2 & loadingbias == 0.2) %>%
  pivot_longer(contains("_est"), names_to = "detector", values_to = "est")
lel$lower = lel$upper =  NA
lel$lower[lel$detector == "R1_est"] = lel$R1_lower[lel$detector == "R1_est"]
lel$lower[lel$detector == "R2_est"] = lel$R2_lower[lel$detector == "R2_est"]
lel$lower[lel$detector == "BvdV_est"] = lel$BvdV_lower[lel$detector == "BvdV_est"]
lel$lower[lel$detector == "CR_est"] = lel$CR_lower[lel$detector == "CR_est"]

lel$upper[lel$detector == "R1_est"] = lel$R1_upper[lel$detector == "R1_est"]
lel$upper[lel$detector == "R2_est"] = lel$R2_upper[lel$detector == "R2_est"]
lel$upper[lel$detector == "BvdV_est"] = lel$BvdV_upper[lel$detector == "BvdV_est"]
lel$upper[lel$detector == "CR_est"] = lel$CR_upper[lel$detector == "CR_est"]

ggplot(lel,
       aes(x = n, y = est, color = detector, shape = detector, ymin = lower, ymax = upper)) +
  geom_pointrange() + 
  geom_line() +
  facet_wrap(~p, labeller = function(x) label_both(x, sep = " = ")) + 
  scale_x_continuous(breaks = c(200, 400, 800))
  

colnames(sim_param_df)
#lapply(out, function(x) apply(x$count_identified, 2, function(y) mean(y == x$sim_pars$k, na.rm = T)))
out$`10`$sim_pars$k
out$`10`$contain_correct

contain_correct
colMeans(contain_correct, na.rm = T)
colMeans(count_correct, na.rm = T)

n_denom = apply(count_correct, 2, function(x) sum(!is.na(x)))
n_num = apply(count_correct, 2, function(x) sum(x[!is.na(x)] == 1))
confint(prop.test(n_num, n_denom))
prop.test(n_num[2], n_denom[2])

prop.test(apply(count_correct, 2, function(x) sum(x[] == 1)))
