## Application - Populism (Castanho et al. 2019)
library(lavaan)
library(here)
source(here("Rscripts", "Simulator_PokropekEtAl.R"))
source(here("Rscripts", "Detector_Rieger.R"))
source(here("Rscripts", "Detector_ByrneVandeVijer.R"))
source(here("Rscripts", "Detector_CheungRensvold.R"))

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

# Akkerman model
varnames = paste0("akker", 1:6)
akker.data = data[, varnames]
akker.data$grp = data$country
akker.data = na.omit(akker.data)
detect_Rieger_v1(varnames, akker.data)
detect_Rieger_v2(varnames, akker.data)
detect_ByrneVandeVijer(varnames, akker.data)
detect_CheungRensvold(varnames, akker.data)

# CSES
varnames = c('akker6', 'cses1', 'cses2.r', 'cses3', 'cses4', 'cses5', 'akker2')
CSES.data = data[, varnames]
CSES.data$grp = data$country
CSES.data = na.omit(CSES.data)
detect_Rieger_v1(varnames, CSES.data)
detect_Rieger_v2(varnames, CSES.data)
detect_ByrneVandeVijer(varnames, CSES.data)
detect_CheungRensvold(varnames, CSES.data)

# Oliver & Rahn
range1 <- function(x){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (5 - (1)) + 1 }
data$ow_me4<-range1(data$manich1)
model = 'antiel =~ ow_ae1 + ow_ae2 + ow_ae3 + ow_ae4 + ow_ae5..
 mistexp =~ ow_me1 + ow_me2 + ow_me3 + ow_me4
 nataff =~ ow_na1 + ow_na2 + ow_na3'
OR.data = data[,c('ow_ae1', 'ow_ae2', 'ow_ae3', 'ow_ae4', 'ow_ae5..', 'ow_me1', 'ow_me2', 'ow_me4',
               'ow_me3', 'ow_na1', 'ow_na2', 'ow_na3')]
OR.data$grp = data$country
data = na.omit(OR.data)

## Elchardus and Spruyt
varnames  = c('es1', 'es2', 'es3', 'es4')
ES.data = data[, varnames]
ES.data$grp = data$country
ES.data = na.omit(ES.data)
detect_Rieger_v1(varnames, ES.data)
detect_Rieger_v2(varnames, ES.data)
detect_ByrneVandeVijer(varnames, ES.data)
detect_CheungRensvold(varnames, ES.data) # TODO: depends on order of varnames => coding error? 




