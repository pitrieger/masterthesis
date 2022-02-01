## Application - Populism (Castanho et al. 2020)
library(lavaan)
library(here)
library(tidyverse)
library(stringr)
library(e1071)
library(stargazer)
set.seed(629)

# load detectors for single-factor models
detectors = list.files(here("Rscripts/simulation")) 
detectors = detectors[startsWith(detectors, "Detector")]
sapply(detectors, function(i) source(here("Rscripts/simulation", i)))
# load detectors for multi factor models
detectorsMulti = list.files(here("Rscripts/application")) 
detectorsMulti = detectorsMulti[startsWith(detectorsMulti, "DetectorMulti")]
sapply(detectorsMulti, function(i) source(here("Rscripts/application", i)))

# read data
data = read.csv(here('data/CastanhoEtAl2019/', 'data_replication_prq1.csv'),header=F, sep=' ')
varnames = unlist(strsplit(readLines(here('data/CastanhoEtAl2019/','varnames_mar17.txt')), '\t'))
names(data) = varnames
data[data==-999] = NA

# Recode trust variables to make it more intuitive:
data$t_party.r = 6-data$t_party
data$t_parl.r = 6-data$t_parl
data$t_gov.r = 6-data$t_gov
data$cses2.r = 6-data$cses2

# for Oliver & Rahn: recode manich1 to ow_me4 from 1-7 scale to 1-5
range1 <- function(x){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (5 - (1)) + 1 }
data$ow_me4 <- range1(data$manich1)

# for Castanho: 
data$simple8.r<-8-data$simple8
data$rwpop8.r<-8-data$rwpop8
data$manich13.r<-8-data$manich13

##########################
## Single-factor models ## 
##########################
# Akkerman model ----
Akker.varnames = paste0("akker", 1:6)
Akker.data = data[,Akker.varnames]
Akker.data$grp = data$country
Akker.data = na.omit(Akker.data)
Akker.model = '
akk =~ akker1 + akker2 + akker3 + akker4 + akker5 + akker6
'

# Akkerman detection
(Akker.R1 = detect_Rieger(Akker.varnames, Akker.data))
(Akker.R2 = detect_Rieger_step(Akker.varnames, Akker.data))
(Akker.BV = detect_ByrneVandeVijer(Akker.varnames, Akker.data))
(Akker.CR_bonf = detect_CheungRensvold(Akker.varnames, Akker.data))
(Akker.MInd_bonf = detect_MInd(Akker.varnames, Akker.data))
(Akker.R1_metric = detect_Rieger(Akker.varnames, Akker.data, detection.type = "metric"))
(Akker.R2_metric = detect_Rieger_step(Akker.varnames, Akker.data, detection.type = "metric"))
(Akker.BV_metric = detect_ByrneVandeVijer(Akker.varnames, Akker.data, group.constraints = "loadings"))
(Akker.CR_bonf_metric = detect_CheungRensvold(Akker.varnames, Akker.data, group.constraints = "loadings"))
(Akker.MInd_bonf_metric = detect_MInd(Akker.varnames, Akker.data, group.constraints = "loadings"))

# save results
Akker.out = matrix("", ncol = 10, nrow = length(Akker.varnames), 
                   dimnames = list(Akker.varnames, rep(c("MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"), 2)))
Akker.out[rownames(Akker.out) %in% Akker.MInd_bonf$noninvariant_bonferroni,1] = "X"
Akker.out[rownames(Akker.out) %in% Akker.CR_bonf$noninvariant_bonferroni,  2] = "X"
Akker.out[rownames(Akker.out) %in% Akker.BV$noninvariant,  3] = "X"
Akker.out[rownames(Akker.out) %in% Akker.R1$noninvariant,  4] = "X"
Akker.out[rownames(Akker.out) %in% Akker.R2$noninvariant,  5] = "X"
Akker.out[rownames(Akker.out) %in% Akker.MInd_bonf_metric$noninvariant_bonferroni,6] = "X"
Akker.out[rownames(Akker.out) %in% Akker.CR_bonf_metric$noninvariant_bonferroni,  7] = "X"
Akker.out[rownames(Akker.out) %in% Akker.BV_metric$noninvariant,  8] = "X"
Akker.out[rownames(Akker.out) %in% Akker.R1_metric$noninvariant,  9] = "X"
Akker.out[rownames(Akker.out) %in% Akker.R2_metric$noninvariant,  10] = "X"
Akker.out
stargazer(Akker.out, summary = F)


# CSES ----
CSES.varnames = c('akker6', 'cses1', 'cses2.r', 'cses3', 'cses4', 'cses5', 'akker2')
CSES.data = data[, CSES.varnames]
CSES.data$grp = data$country
CSES.data = na.omit(CSES.data)
CSES.model = '
cses =~ akker6 + cses1 + cses2.r + cses3 + cses4 + cses5 + akker2
'

# CSES detection
(CSES.R1 = detect_Rieger(CSES.varnames, CSES.data))
(CSES.R2 = detect_Rieger_step(CSES.varnames, CSES.data))
(CSES.BV = detect_ByrneVandeVijer(CSES.varnames, CSES.data))
(CSES.CR_bonf = detect_CheungRensvold(CSES.varnames, CSES.data))
(CSES.MInd_bonf = detect_MInd(CSES.varnames, CSES.data))
(CSES.R1_metric = detect_Rieger(CSES.varnames, CSES.data, detection.type = "metric"))
(CSES.R2_metric = detect_Rieger_step(CSES.varnames, CSES.data, detection.type = "metric"))
(CSES.BV_metric = detect_ByrneVandeVijer(CSES.varnames, CSES.data, group.constraints = "loadings"))
(CSES.CR_bonf_metric = detect_CheungRensvold(CSES.varnames, CSES.data, group.constraints = "loadings"))
(CSES.MInd_bonf_metric = detect_MInd(CSES.varnames, CSES.data, group.constraints = "loadings"))

# save results
CSES.out = matrix("", ncol = 10, nrow = length(CSES.varnames), 
                   dimnames = list(CSES.varnames, rep(c("MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"), 2)))
CSES.out[rownames(CSES.out) %in% CSES.MInd_bonf$noninvariant_bonferroni,1] = "X"
CSES.out[rownames(CSES.out) %in% CSES.CR_bonf$noninvariant_bonferroni,  2] = "X"
CSES.out[rownames(CSES.out) %in% CSES.BV$noninvariant,  3] = "X"
CSES.out[rownames(CSES.out) %in% CSES.R1$noninvariant,  4] = "X"
CSES.out[rownames(CSES.out) %in% CSES.R2$noninvariant,  5] = "X"
CSES.out[rownames(CSES.out) %in% CSES.MInd_bonf_metric$noninvariant_bonferroni,6] = "X"
CSES.out[rownames(CSES.out) %in% CSES.CR_bonf_metric$noninvariant_bonferroni,  7] = "X"
CSES.out[rownames(CSES.out) %in% CSES.BV_metric$noninvariant,  8] = "X"
CSES.out[rownames(CSES.out) %in% CSES.R1_metric$noninvariant,  9] = "X"
CSES.out[rownames(CSES.out) %in% CSES.R2_metric$noninvariant,  10] = "X"
CSES.out
stargazer(CSES.out, summary = F)

# Elchardus and Spruyt ----
ES.varnames  = c('es1', 'es2', 'es3', 'es4')
ES.data = data[, ES.varnames]
ES.data$grp = data$country
ES.data = na.omit(ES.data)
ES.model = '
es =~ es1 + es2 + es3 + es4
'

# ES detection
(ES.R1 = detect_Rieger(ES.varnames, ES.data))
(ES.R2 = detect_Rieger_step(ES.varnames, ES.data))
(ES.BV = detect_ByrneVandeVijer(ES.varnames, ES.data))
(ES.CR_bonf = detect_CheungRensvold(ES.varnames, ES.data))
(ES.MInd_bonf = detect_MInd(ES.varnames, ES.data))
(ES.R1_metric = detect_Rieger(ES.varnames, ES.data, detection.type = "metric"))
(ES.R2_metric = detect_Rieger_step(ES.varnames, ES.data, detection.type = "metric"))
(ES.BV_metric = detect_ByrneVandeVijer(ES.varnames, ES.data, group.constraints = "loadings"))
(ES.CR_bonf_metric = detect_CheungRensvold(ES.varnames, ES.data, group.constraints = "loadings"))
(ES.MInd_bonf_metric = detect_MInd(ES.varnames, ES.data, group.constraints = "loadings"))

# save results
ES.out = matrix("", ncol = 10, nrow = length(ES.varnames), 
                  dimnames = list(ES.varnames, rep(c("MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"), 2)))
ES.out[rownames(ES.out) %in% ES.MInd_bonf$noninvariant_bonferroni,1] = "X"
ES.out[rownames(ES.out) %in% ES.CR_bonf$noninvariant_bonferroni,  2] = "X"
ES.out[rownames(ES.out) %in% ES.BV$noninvariant,  3] = "X"
ES.out[rownames(ES.out) %in% ES.R1$noninvariant,  4] = "X"
ES.out[rownames(ES.out) %in% ES.R2$noninvariant,  5] = "X"
ES.out[rownames(ES.out) %in% ES.MInd_bonf_metric$noninvariant_bonferroni,6] = "X"
ES.out[rownames(ES.out) %in% ES.CR_bonf_metric$noninvariant_bonferroni,  7] = "X"
ES.out[rownames(ES.out) %in% ES.BV_metric$noninvariant,  8] = "X"
ES.out[rownames(ES.out) %in% ES.R1_metric$noninvariant,  9] = "X"
ES.out[rownames(ES.out) %in% ES.R2_metric$noninvariant,  10] = "X"
ES.out
stargazer(ES.out, summary = F)

#########################
## Multi-factor models ##
#########################
# Oliver & Rahn ----
colnames(data)[colnames(data) == "ow_ae5.."] = "ow_ae5" 
OR.varnames = c("ow_ae1", "ow_ae2", "ow_ae3", "ow_ae4", "ow_ae5",
                "ow_me1", "ow_me2", "ow_me3", "ow_me4",
                "ow_na1", "ow_na2", "ow_na3")
OR.data = data[, OR.varnames]
OR.data$grp = data$country
OR.data = na.omit(OR.data)
OR.model = '
antiel =~ ow_ae1 + ow_ae2 + ow_ae3 + ow_ae4 + ow_ae5
mistexp =~ ow_me1 + ow_me2 + ow_me3 + ow_me4
nataff =~ ow_na1 + ow_na2 + ow_na3
'

# OR detection
(OR.R1 = detectMulti_Rieger(OR.varnames, OR.model, OR.data))
(OR.R2 = detectMulti_Rieger_step(OR.varnames, OR.model, OR.data))
(OR.BV = detectMulti_ByrneVandeVijer(OR.varnames, OR.model, OR.data))
#(OR.CR_bonf = detectMulti_CheungRensvold(OR.varnames, OR.model, OR.data)) # CONVERGENCE ISSUES
(OR.MInd_bonf = detectMulti_MInd(OR.varnames, OR.model, OR.data))
(OR.R1_metric = detectMulti_Rieger(OR.varnames, OR.model, OR.data, detection.type = "metric"))
(OR.R2_metric = detectMulti_Rieger_step(OR.varnames, OR.model, OR.data, detection.type = "metric"))
(OR.BV_metric = detectMulti_ByrneVandeVijer(OR.varnames, OR.model, OR.data, group.constraints = "loadings"))
#(OR.CR_bonf_metric = detectMulti_CheungRensvold(OR.varnames, OR.model, OR.data, group.constraints = "loadings"))
(OR.MInd_bonf_metric = detectMulti_MInd(OR.varnames, OR.model, OR.data, group.constraints = "loadings"))

# save results 
OR.out = matrix("", ncol = 10, nrow = length(OR.varnames), 
                dimnames = list(OR.varnames, rep(c("MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"), 2)))
OR.out[rownames(OR.out) %in% OR.MInd_bonf$noninvariant_bonferroni,1] = "X"
OR.out[, 2] = "-"
OR.out[rownames(OR.out) %in% OR.BV$noninvariant,  3] = "X"
OR.out[rownames(OR.out) %in% OR.R1$noninvariant,  4] = "X"
OR.out[rownames(OR.out) %in% OR.R2$noninvariant,  5] = "X"
OR.out[rownames(OR.out) %in% OR.MInd_bonf_metric$noninvariant_bonferroni,6] = "X"
OR.out[, 7] = "-"
OR.out[rownames(OR.out) %in% OR.BV_metric$noninvariant,  8] = "X"
OR.out[rownames(OR.out) %in% OR.R1_metric$noninvariant,  9] = "X"
OR.out[rownames(OR.out) %in% OR.R2_metric$noninvariant,  10] = "X"
OR.out
stargazer(OR.out, summary = F)


# Stanley ----
Stan.varnames = c("stanley1", "stanley2", "stanley3", "stanley4", "stanley5", "stanley6", "stanley7", "stanley8")
Stan.data = data[,Stan.varnames]
Stan.data$grp = data$country
Stan.data = na.omit(Stan.data)
Stan.model<-'
sta =~ stanley1 + stanley2 + stanley3 + stanley4 + stanley5 + stanley6 + stanley7 + stanley8
method =~ stanley3 + stanley8
method ~~ 0*sta
'

# OR detection
(Stan.R1 = detectMulti_Rieger(Stan.varnames, Stan.model, Stan.data))
(Stan.R2 = detectMulti_Rieger_step(Stan.varnames, Stan.model, Stan.data))
(Stan.BV = detectMulti_ByrneVandeVijer(Stan.varnames, Stan.model, Stan.data))
#(Stan.CR_bonf = detectMulti_CheungRensvold(Stan.varnames, Stan.model, Stan.data)) # Convergence issues 
(Stan.MInd_bonf = detectMulti_MInd(Stan.varnames, Stan.model, Stan.data))
(Stan.R1_metric = detectMulti_Rieger(Stan.varnames, Stan.model, Stan.data, detection.type = "metric"))
(Stan.R2_metric = detectMulti_Rieger_step(Stan.varnames, Stan.model, Stan.data, detection.type = "metric"))
(Stan.BV_metric = detectMulti_ByrneVandeVijer(Stan.varnames, Stan.model, Stan.data, group.constraints = "loadings"))
(Stan.CR_bonf_metric = detectMulti_CheungRensvold(Stan.varnames, Stan.model, Stan.data, group.constraints = "loadings"))
(Stan.MInd_bonf_metric = detectMulti_MInd(Stan.varnames, Stan.model, Stan.data, group.constraints = "loadings"))

# save results 
Stan.out = matrix("", ncol = 10, nrow = length(Stan.varnames), 
                dimnames = list(Stan.varnames, rep(c("MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"), 2)))
Stan.out[rownames(Stan.out) %in% Stan.MInd_bonf$noninvariant_bonferroni,1] = "X"
Stan.out[, 2] = "-"
Stan.out[rownames(Stan.out) %in% Stan.BV$noninvariant,  3] = "X"
Stan.out[rownames(Stan.out) %in% Stan.R1$noninvariant,  4] = "X"
Stan.out[rownames(Stan.out) %in% Stan.R2$noninvariant,  5] = "X"
Stan.out[rownames(Stan.out) %in% Stan.MInd_bonf_metric$noninvariant_bonferroni,6] = "X"
Stan.out[, 7] = "-"
Stan.out[rownames(Stan.out) %in% Stan.BV_metric$noninvariant,  8] = "X"
Stan.out[rownames(Stan.out) %in% Stan.R1_metric$noninvariant,  9] = "X"
Stan.out[rownames(Stan.out) %in% Stan.R2_metric$noninvariant,  10] = "X"
Stan.out
stargazer(Stan.out, summary = F)

# Schulz et al ----
Schulz.varnames = c("nccr_ant1", "nccr_ant2", "nccr_ant3", 
                    "nccr_sov1", "nccr_sov2", "akker2", 
                    "nccr_hom1", "nccr_hom2", "nccr_hom3")
Schulz.data = data[,Schulz.varnames]
Schulz.data$grp = data$country
Schulz.data = na.omit(Schulz.data)
Schulz.model<-'
antiel =~ nccr_ant1 + nccr_ant2 + nccr_ant3
sov =~ nccr_sov1 + nccr_sov2 + akker2
hom =~ nccr_hom1 + nccr_hom2 + nccr_hom3
pop =~ antiel + sov + hom
'

# OR detection
(Schulz.R1 = detectMulti_Rieger(Schulz.varnames, Schulz.model, Schulz.data))
(Schulz.R2 = detectMulti_Rieger_step(Schulz.varnames, Schulz.model, Schulz.data))
(Schulz.BV = detectMulti_ByrneVandeVijer(Schulz.varnames, Schulz.model, Schulz.data))
(Schulz.CR_bonf = detectMulti_CheungRensvold(Schulz.varnames, Schulz.model, Schulz.data))
(Schulz.MInd_bonf = detectMulti_MInd(Schulz.varnames, Schulz.model, Schulz.data))
(Schulz.R1_metric = detectMulti_Rieger(Schulz.varnames, Schulz.model, Schulz.data, detection.type = "metric"))
(Schulz.R2_metric = detectMulti_Rieger_step(Schulz.varnames, Schulz.model, Schulz.data, detection.type = "metric"))
(Schulz.BV_metric = detectMulti_ByrneVandeVijer(Schulz.varnames, Schulz.model, Schulz.data, group.constraints = "loadings"))
(Schulz.CR_bonf_metric = detectMulti_CheungRensvold(Schulz.varnames, Schulz.model, Schulz.data, group.constraints = "loadings"))
(Schulz.MInd_bonf_metric = detectMulti_MInd(Schulz.varnames, Schulz.model, Schulz.data, group.constraints = "loadings"))

# save results 
Schulz.out = matrix("", ncol = 10, nrow = length(Schulz.varnames), 
                dimnames = list(Schulz.varnames, rep(c("MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"), 2)))
Schulz.out[rownames(Schulz.out) %in% Schulz.MInd_bonf$noninvariant_bonferroni,1] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.CR_bonf$noninvariant_bonferroni, 2] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.BV$noninvariant,  3] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.R1$noninvariant,  4] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.R2$noninvariant,  5] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.MInd_bonf_metric$noninvariant_bonferroni,6] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.CR_bonf_metric$noninvariant_bonferroni, 7] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.BV_metric$noninvariant,  8] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.R1_metric$noninvariant,  9] = "X"
Schulz.out[rownames(Schulz.out) %in% Schulz.R2_metric$noninvariant,  10] = "X"
Schulz.out
stargazer(Schulz.out, summary = F)


# Castanho et al ----
Castanho.varnames = c("antiel23", "rwpop8.r", "antiel21",
                      "gewill17", "simple8.r", "gewill3",
                      "manich15", "manich13.r", "manich14")
Castanho.data = data[, Castanho.varnames]
Castanho.data$grp = data$country
Castanho.data = na.omit(Castanho.data)

Castanho.model = '
antiel =~ antiel23 + rwpop8.r + antiel21
people =~ gewill17 + simple8.r + gewill3
manich =~ manich15 + manich13.r + manich14
method =~ 1*gewill3 + b1*antiel23 +  b1*antiel21 + b1*gewill17 + b1*manich15 + b1*manich14
method ~~ 0*antiel + 0*people + 0*manich
'

# OR detection
(Castanho.R1 = detectMulti_Rieger(Castanho.varnames, Castanho.model, Castanho.data))
(Castanho.R2 = detectMulti_Rieger_step(Castanho.varnames, Castanho.model, Castanho.data))
(Castanho.BV = detectMulti_ByrneVandeVijer(Castanho.varnames, Castanho.model, Castanho.data))
(Castanho.CR_bonf = detectMulti_CheungRensvold(Castanho.varnames, Castanho.model, Castanho.data))
(Castanho.MInd_bonf = detectMulti_MInd(Castanho.varnames, Castanho.model, Castanho.data))
(Castanho.R1_metric = detectMulti_Rieger(Castanho.varnames, Castanho.model, Castanho.data, detection.type = "metric"))
(Castanho.R2_metric = detectMulti_Rieger_step(Castanho.varnames, Castanho.model, Castanho.data, detection.type = "metric"))
(Castanho.BV_metric = detectMulti_ByrneVandeVijer(Castanho.varnames, Castanho.model, Castanho.data, group.constraints = "loadings"))
(Castanho.CR_bonf_metric = detectMulti_CheungRensvold(Castanho.varnames, Castanho.model, Castanho.data, group.constraints = "loadings"))
(Castanho.MInd_bonf_metric = detectMulti_MInd(Castanho.varnames, Castanho.model, Castanho.data, group.constraints = "loadings"))

# save results 
Castanho.out = matrix("", ncol = 10, nrow = length(Castanho.varnames), 
                dimnames = list(Castanho.varnames, rep(c("MInd (Bonf.)", "CR (Bonf.)", "BV", "R1", "R2"), 2)))
Castanho.out[rownames(Castanho.out) %in% Castanho.MInd_bonf$noninvariant_bonferroni,1] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.CR_bonf$noninvariant_bonferroni, 2] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.BV$noninvariant,  3] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.R1$noninvariant,  4] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.R2$noninvariant,  5] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.MInd_bonf_metric$noninvariant_bonferroni,6] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.CR_bonf_metric$noninvariant_bonferroni, 7] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.BV_metric$noninvariant,  8] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.R1_metric$noninvariant,  9] = "X"
Castanho.out[rownames(Castanho.out) %in% Castanho.R2_metric$noninvariant,  10] = "X"
Castanho.out
stargazer(Castanho.out, summary = F)


# ALL RESULTS ----
Akker.out
CSES.out
ES.out
OR.out
Stan.out
Schulz.out
Castanho.out

stargazer(Akker.out, summary = F)
stargazer(CSES.out, summary = F)
stargazer(ES.out, summary = F)
stargazer(OR.out, summary = F)
stargazer(Stan.out, summary = F)
stargazer(Schulz.out, summary = F)
stargazer(Castanho.out, summary = F)
