library(tidyverse)
library(parallel)
library(lavaan)
library(here)
source(here("Rscripts/simulation", "Simulator_PokropekEtAl.R"))
source(here("Rscripts/simulation", "Detector_Rieger.R"))
source(here("Rscripts/simulation", "Detector_ByrneVandeVijer.R"))
source(here("Rscripts/simulation", "Detector_CheungRensvold.R"))
MMGFA_files = list.files(here("Rscripts", "KimDeRoover_MixtureMG_FA"))
MMGFA_files = MMGFA_files[grepl("R$", MMGFA_files)]
for(j in 1:length(MMGFA_files)){
  source(here("Rscripts", "KimDeRoover_MixtureMG_FA", MMGFA_files[j]))
}
RNGkind("L'Ecuyer-CMRG")
set.seed(9)

detectCores()

# Grid
sim_param = expand.grid(nsim = 100,
                        n = 1000,
                        p = c(3, 4, 5, 10, 20),
                        g = c(4, 8, 16, 32),
                        h = c(0.25, 0.5),
                        k = c(1, 2, 3, 4, 5, 10))
sim_param = sim_param[sim_param$k<sim_param$p,]
sim_param = split(sim_param, 1:nrow(sim_param))
sim_out = lapply(sim_param, function(x) set_names(list(x), "sim_param"))

mclapply(sim_out, )

getOption("mc.cores", 2L)

for(i in 1:nrow(sim_param)){
  cat(i, "\n")
  varnames = paste0("y", 1:sim_param$p[i])
  contain_correct = count_correct = matrix(nrow = sim_param$nsim[i], ncol = 4)
  colnames(contain_correct) = colnames(count_correct) = c("Rieger1", "Rieger2", "Byrne", "Cheung")
  
  for(s in 1:sim_param$nsim[i]){
    # Simulate data
    D = sim_PMI(n = sim_param$n[i], 
                g = sim_param$g[i], 
                p = sim_param$p[i], 
                h = sim_param$h[i], 
                k = sim_param$k[i])
    noninvariant_true = ifelse(length(D$p_affected)>0, paste0("y", D$p_affected), character())
    
    # Own Idea v1: simple
    noninvariant_R1 = try(detect_Rieger_v1(varnames = varnames,
                                           data = D$sim_dat)$noninvariant)
    if(class(noninvariant_R1)!="try-error"){
      count_correct[s,1] = as.numeric(length(noninvariant_R1) == sim_param$k[i])
      contain_correct[s,1] = sum(noninvariant_true %in% noninvariant_R1)
    }
    
    # Own Idea v2: stepwise
    noninvariant_R2 = try(detect_Rieger_v2(varnames = varnames,
                                           D$sim_dat)$noninvariant)
    if(class(noninvariant_R2)!="try-error"){
      count_correct[s,2] = as.numeric(length(noninvariant_R2) == sim_param$k[i])
      contain_correct[s,2] = sum(noninvariant_true %in% noninvariant_R2)
    }
    
    ## Via CFI
    # Byrne & van de Vijer (2010):
    noninvariant_B = try(detect_ByrneVandeVijer(varnames = varnames,
                                                D$sim_dat)$noninvariant)
    if(class(noninvariant_B)!="try-error"){
      count_correct[s,3] = as.numeric(length(noninvariant_B) == sim_param$k[i])
      contain_correct[s,3] = sum(noninvariant_true %in% noninvariant_B)
    }
    
    # Cheung & Rensveld (1999):
    noninvariant_C = try(detect_CheungRensvold(varnames = varnames,
                                               D$sim_dat)$noninvariant)
    if(class(noninvariant_C)!="try-error"){
      count_correct[s,4] = as.numeric(length(noninvariant_C) == sim_param$k[i])
      contain_correct[s,4] = sum(noninvariant_true %in% noninvariant_C)
    }
    
    res = list()
    res[1] = list(sim_param = as.numeric(sim_param[i,]),
                  count_correct,
                  contain_correct)
    cat("*")
  }
}


contain_correct
colMeans(contain_correct, na.rm = T)
colMeans(count_correct, na.rm = T)

n_denom = apply(count_correct, 2, function(x) sum(!is.na(x)))
n_num = apply(count_correct, 2, function(x) sum(x[!is.na(x)] == 1))
confint(prop.test(n_num, n_denom))
prop.test(n_num[2], n_denom[2])

prop.test(apply(count_correct, 2, function(x) sum(x[] == 1)))
