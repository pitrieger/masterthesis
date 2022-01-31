## Application - Populism (Castanho et al. 2020)
library(lavaan)
library(here)
# load detectors for single-factor models
detectors = list.files(here("Rscripts/simulation")) 
detectors = detectors[startsWith(detectors, "Detector")]
sapply(detectors, function(i) source(here("Rscripts/simulation", i)))

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

##########################
## Single-factor models ## 
##########################

# Akkerman model 
varnames = paste0("akker", 1:6)
akker.data = data[, varnames]
akker.data$grp = data$country
akker.data = na.omit(akker.data)
detect_Rieger(varnames, akker.data)
detect_Rieger_step(varnames, akker.data)
detect_ByrneVandeVijer(varnames, akker.data)
detect_CheungRensvold(varnames, akker.data)
detect_MInd(varnames, akker.data)

# CSES
varnames = c('akker6', 'cses1', 'cses2.r', 'cses3', 'cses4', 'cses5', 'akker2')
CSES.data = data[, varnames]
CSES.data$grp = data$country
CSES.data = na.omit(CSES.data)
detect_Rieger(varnames, CSES.data)
detect_Rieger_step(varnames, CSES.data)
detect_ByrneVandeVijer(varnames, CSES.data)
detect_CheungRensvold(varnames, CSES.data)
detect_MInd(varnames, CSES.data)

# Elchardus and Spruyt
varnames  = c('es1', 'es2', 'es3', 'es4')
ES.data = data[, varnames]
ES.data$grp = data$country
ES.data = na.omit(ES.data)
detect_Rieger(varnames, ES.data)
detect_Rieger_step(varnames, ES.data)
detect_ByrneVandeVijer(varnames, ES.data)
detect_CheungRensvold(varnames, ES.data)
detect_MInd(varnames, ES.data)


