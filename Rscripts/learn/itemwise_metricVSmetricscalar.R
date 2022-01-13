# Checking whether metric implementation makes sense

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

# Simulate data
n = 800
p = 5
g = 8
h = 0.5
k = 2
D = sim_PMI(n = n, g = g, p = p, h = h, k = k,
            itembias = 0, loadingbias = 0.25)
varnames = paste0("y", 1:p)

D$p_affected
detect_ByrneVandeVijer(varnames, D$sim_dat, group.constraints = "loadings")
detect_ByrneVandeVijer(varnames, D$sim_dat)


get_pval2 = function(X, hat.eta, grp){
  stopifnot(nrow(X) == length(hat.eta) & length(hat.eta) == length(grp))
  g = length(unique(grp))
  p = ncol(X)
  p.vals = numeric(p)
  for(j in 1:p){
    # get residuals
    R = resid(lm(X[,j] ~ hat.eta))
    # group correlations
    cor.p = numeric(g)
    for(i in 1:g){
      cor.p[i] = cor.test(R[grp == i], hat.eta[grp == i])$p.value # equivalent to coefficient test for lm
    }
    # aggregate p.vals
    p.vals[j] = (g)*min(cor.p)
  }
  p.vals
}

detect_Rieger = function(varnames, data, alpha = 0.05, type = "both"){
  model = paste("eta =~", paste(varnames, collapse = " + "))
  fit = cfa(model, data = data)
  hat.eta = predict(fit)[,1]
  X = fit@Data@X[[1]]
  grp = data$grp
  if(type == "both"){
    p.vals = get_pval(X, hat.eta, grp)
  } else if(type == "loadings"){
    p.vals = get_pval2(X, hat.eta, grp)
  }
  list(varnames = varnames,
       noninvariant = varnames[which(p.vals < alpha/(length(p.vals) - order(p.vals) + 1))], # Bonferroni-holm
       p.vals = p.vals)
}

D = sim_PMI(n = 800, g = 4, p = 6, h = 0.5, k = 2,
            itembias = 0, loadingbias = 0.25)
varnames = paste0("y", 1:p)
detect_Rieger(varnames, D$sim_dat, type = "both")
detect_Rieger(varnames, D$sim_dat, type = "loadings")

fit = cfa(paste("eta =~", paste(varnames, collapse = " + ")), D$sim_dat)
par(mfrow = c(2, 3))
eta = predict(fit)[,1]
for(j in 1:D$sim_params["p"]){
  reg = lm(D$sim_dat[,varnames[j]] ~ eta)
  plot(eta, resid(reg), col = rgb(0, 0, 0, 0.2), main = varnames[j])
  for(l in 1:D$sim_params["g"]){
    abline(lm(resid(reg) ~ eta, subset = D$sim_dat$grp == l))
  }

}
D$p_affected

plot(predict(fit)[,1], resid(lm(D$sim_dat$y1 ~ predict(fit)[,1])), col = D$sim_dat$grp)

par(mfrow = c(1,1))


