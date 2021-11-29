## Application - Populism (Castanho et al. 2020)
library(lavaan)
library(here)
library(tidyverse)
library(stringr)

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
# Akkerman model 
Akker.varnames = paste0("akker", 1:6)
Akker.data = data[,Akker.varnames]
Akker.data$grp = data$country
Akker.data = na.omit(Akker.data)
Akker.model = '
akk =~ akker1 + akker2 + akker3 + akker4 + akker5 + akker6
'
detect_Rieger(Akker.varnames, Akker.data)
detect_Rieger_step(Akker.varnames, Akker.data)
detect_ByrneVandeVijer(Akker.varnames, Akker.data)
detectMulti_ByrneVandeVijer(Akker.varnames, Akker.model, Akker.data) # should be same as single factor BV function
detect_CheungRensvold(Akker.varnames, Akker.data)
detect_MInd(Akker.varnames, Akker.data)
detectMulti_MInd(Akker.varnames, Akker.model, Akker.data) # should be same as single factor MInd function

# CSES
CSES.varnames = c('akker6', 'cses1', 'cses2.r', 'cses3', 'cses4', 'cses5', 'akker2')
CSES.data = data[, CSES.varnames]
CSES.data$grp = data$country
CSES.data = na.omit(CSES.data)
CSES.model = '
cses =~ akker6 + cses1 + cses2.r + cses3 + cses4 + cses5 + akker2
'
detect_Rieger(CSES.varnames, CSES.data)
detect_Rieger_step(CSES.varnames, CSES.data)
detect_ByrneVandeVijer(CSES.varnames, CSES.data)
detectMulti_ByrneVandeVijer(CSES.varnames, CSES.model, CSES.data) # should be same as single factor BV function
detect_CheungRensvold(CSES.varnames, CSES.data)
detect_MInd(CSES.varnames, CSES.data)
detectMulti_MInd(CSES.varnames, CSES.model, CSES.data) # should be same as single factor MInd function

# Elchardus and Spruyt
ES.varnames  = c('es1', 'es2', 'es3', 'es4')
ES.data = data[, ES.varnames]
ES.data$grp = data$country
ES.data = na.omit(ES.data)
ES.model = '
es =~ es1 + es2 + es3 + es4
'
detect_Rieger(ES.varnames, ES.data)
detect_Rieger_step(ES.varnames, ES.data)
detect_ByrneVandeVijer(ES.varnames, ES.data)
detectMulti_ByrneVandeVijer(ES.varnames, ES.model, ES.data) # should be same as single factor BV function
detect_CheungRensvold(ES.varnames, ES.data)
detect_MInd(ES.varnames, ES.data)
detectMulti_MInd(ES.varnames, ES.model, ES.data) # should be same as single factor MInd function

#########################
## Multi-factor models ##
#########################
# Oliver & Rahn
OR.varnames = c("ow_ae1", "ow_ae2", "ow_ae3", "ow_ae4", "ow_ae5..",
             "ow_me1", "ow_me2", "ow_me3", "ow_me4",
             "ow_na1", "ow_na2", "ow_na3")
OR.data = data[, OR.varnames]
OR.data$grp = data$country
OR.data = na.omit(OR.data)

OR.model = '
antiel =~ ow_ae1 + ow_ae2 + ow_ae3 + ow_ae4 + ow_ae5..
mistexp =~ ow_me1 + ow_me2 + ow_me3 + ow_me4
nataff =~ ow_na1 + ow_na2 + ow_na3
'

detectMulti_ByrneVandeVijer(OR.varnames, OR.model, OR.data)
detectMulti_MInd(OR.varnames, OR.model, OR.data)
#base_model = OR.model
#data = OR.data
#varnames = OR.varnames

# Stanley
Stan.varnames = c("stanley1", "stanley2", "stanley3", "stanley4", "stanley5", "stanley6", "stanley7", "stanley8")
Stan.data = data[,Stan.varnames]
Stan.data$grp = data$country
Stan.data = na.omit(Stan.data)

Stan.model<-'
sta =~ stanley1 + stanley2 + stanley3 + stanley4 + stanley5 + stanley6 + stanley7 + stanley8
method =~ stanley3 + stanley8
method ~~ 0*sta
'

detectMulti_ByrneVandeVijer(Stan.varnames, Stan.model, Stan.data)
detectMulti_MInd(Stan.varnames, Stan.model, Stan.data)

# Schulz et al
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

detectMulti_ByrneVandeVijer(Schulz.varnames, Schulz.model, Schulz.data)
detectMulti_MInd(Schulz.varnames, Schulz.model, Schulz.data)

# Castanho et al
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
detectMulti_ByrneVandeVijer(Castanho.varnames, Castanho.model, Castanho.data)
detectMulti_MInd(Castanho.varnames, Castanho.model, Castanho.data)


