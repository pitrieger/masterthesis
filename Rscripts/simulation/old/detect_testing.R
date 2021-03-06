library(tidyverse)
library(lavaan)
library(here)
source(here("Rscripts/simulation", "Simulator_PokropekEtAl.R"))
source(here("Rscripts/simulation", "Detector_Rieger_v3.R"))
source(here("Rscripts/simulation", "Detector_ByrneVandeVijer.R"))
#source(here("Rscripts/simulation", "Detector_CheungRensvold.R"))
source(here("Rscripts/simulation", "Detector_CheungRensvold_v2.R"))
source(here("Rscripts/simulation", "Detector_MInd.R"))
source(here("Rscripts/simulation", "Detector_Janssens.R"))
MMGFA_files = list.files(here("Rscripts", "KimDeRoover_MixtureMG_FA"))
MMGFA_files = MMGFA_files[grepl("R$", MMGFA_files)]
for(j in 1:length(MMGFA_files)){
  source(here("Rscripts", "KimDeRoover_MixtureMG_FA", MMGFA_files[j]))
}

# Simulate data
n = 800
p = 3
g = 4
h = 0.5
k = 1
D = sim_PMI(n = n, g = g, p = p, h = h, k = k, interceptbias = 0.5, loadingbias = 0)
varnames = paste0("y", 1:p)

detect_Rieger(varnames, D$sim_dat)
detect_Rieger(varnames, D$sim_dat, detection.type = "metric")
detect_Rieger_step(varnames, D$sim_dat)
detect_Rieger_step(varnames, D$sim_dat, detection.type = "metric")

data = D$sim_dat

D$p_affected


fit = cfa("eta =~ y1 + y2 + y3", D$sim_dat)

summary(fit)
fit = cfa("eta =~ y1 + y2 + y3", D$sim_dat, group = "grp")
summary(fit)
fit = cfa("eta =~ y1 + y2 + y3", D$sim_dat, group = "grp", group.equal = "residuals")
summary(fit)

detect_CheungRensvold(varnames, D$sim_dat)

detect_ByrneVandeVijer(varnames = paste0("y", 1:p),
                       D$sim_dat)
D$p_affected


# Fit full model
  #mod = paste("eta =~", paste(paste0("y", 1:p), collapse = " + "))
  #fit = cfa(mod, data = D$sim_dat)
  #summary(fit)
  #fit@ParTable
  #fit = x

# Own Idea v1: simple
detect_Rieger_v1(varnames = varnames,
                 data = D$sim_dat)
D$p_affected

# Own Idea v2: stepwise
detect_Rieger_v2(varnames = paste0("y", 1:p),
                 D$sim_dat)
D$p_affected

## Via CFI
# Byrne & van de Vijer (2010):
detect_ByrneVandeVijer(varnames = paste0("y", 1:p),
                       D$sim_dat)
D$p_affected

# Cheung & Rensveld (1999):
detect_CheungRensvold(varnames = paste0("y", 1:p),
                      D$sim_dat)

# MInd 
detect_MInd(varnames = varnames,
            D$sim_dat)
D$p_affected

# Janssens
detect_Janssens(varnames, D$sim_dat)
D$p_affected

# Via Roover SCA-P
Y = D$sim_dat[,1:(p+1)]

# groupwise center & aggregate scale
Y = Y %>% group_by(grp) %>%
  mutate(across(starts_with("y"), function(x) scale(x, center = T, scale = F)))
Y = as.matrix(Y[,1:p])
Y = Y %*% solve(cov(Y))


# Simulate data
D = sim_PMI(n = out$`3`$sim_pars$n, 
            g = out$`3`$sim_pars$g,
            p = out$`3`$sim_pars$p, 
            h = out$`3`$sim_pars$h, 
            k = out$`3`$sim_pars$k)
varnames = paste0("y", 1:out$`3`$sim_pars$p)
detect_ByrneVandeVijer(varnames, D$sim_dat)
D$p_affected
D$p_affected

detect_MInd(varnames, D$sim_dat)

detect_Janssens(varnames, D$sim_dat)




n = 500
p = 5
g = 8
h = 0.5
k = 2
D = sim_PMI(n = n, g = g, p = p, h = h, k = k, interceptbias = 0.25, loadingbias = 0)
varnames = paste0("y", 1:p)
D$p_affected

detect_CheungRensvold(varnames, D$sim_dat)
detect_CheungRensvold2(varnames, D$sim_dat, group.constraints = c("loadings"))
