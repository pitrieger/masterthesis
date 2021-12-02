# load the lavaan package (only needed once per session)
library(lavaan)
# specify the model
HS.model <- ' visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9 
x1 ~ 1
x2 ~ 1
x3 ~ 1
x4 ~ 1
x5 ~ 1
x6 ~ 1
x7 ~ 1
x8 ~ 1
x9 ~ 1
'
# fit the model
fit <- cfa(HS.model, data = HolzingerSwineford1939)
# display summary output
summary(fit, fit.measures = TRUE)
fitmeasures(fit, "cfi")
1 - (fitmeasures(fit, "chisq")- 24)/(fitmeasures(fit, "baseline.chisq") - 36)

fitbase = cfa(mod_baseline, data = HolzingerSwineford1939, meanstructure = T)
summary(fitbase)

apply(HolzingerSwineford1939[, paste0("x", 1:9)], 2, mean)
fit@ParTable$est[10:19]
(ui = 9*10/2 + 9) # = 54 unique infos
(pars = 9 + # intercepts
3*(3-1) + # factor loadings with first one per latent variable fixed to 1
3*4/2 + # factor covariances
9) # unique variances
(dfs = ui - pars)
  
sigma.M = fit@implied$cov[[1]]
sigma.Mbase = fitbase@implied$cov[[1]]
S = fit@SampleStats@cov[[1]]
n = nrow(HolzingerSwineford1939)
tr = function(x) sum(diag(x))
F.M = log(det(sigma.M)) - log(det(S)) + tr(S%*%solve(sigma.M)) - 9
n*F.M
F.Mbase = log(det(sigma.Mbase)) - log(det(S)) + tr(S%*%solve(sigma.Mbase)) - 9
n*F.Mbase


summary(fit)

diag(9)
sum(1:9)
lel = 6
sum(lower.tri(matrix(1, ncol = lel, nrow = lel), diag = T))
lel*(lel + 3) / 2

?lower.tri
p = predict(fit)
p


data(HolzingerSwineford1939)

par(mfrow = c(3, 3))
plot(HolzingerSwineford1939$x1, p[,1])
abline(lm(p[,1] ~ HolzingerSwineford1939$x1))

plot(HolzingerSwineford1939$x2, p[,1])
abline(lm(p[,1] ~ HolzingerSwineford1939$x2))

plot(HolzingerSwineford1939$x3, p[,1])
abline(lm(p[,1] ~ HolzingerSwineford1939$x3))


plot(HolzingerSwineford1939$x4, p[,2])
abline(lm(p[,2] ~ HolzingerSwineford1939$x4))

plot(HolzingerSwineford1939$x5, p[,2])
abline(lm(p[,2] ~ HolzingerSwineford1939$x5))

plot(HolzingerSwineford1939$x6, p[,2])
abline(lm(p[,2] ~ HolzingerSwineford1939$x6))


plot(HolzingerSwineford1939$x7, p[,3])
abline(lm(p[,3] ~ HolzingerSwineford1939$x7))

plot(HolzingerSwineford1939$x8, p[,3])
abline(lm(p[,3] ~ HolzingerSwineford1939$x8))

plot(HolzingerSwineford1939$x9, p[,3])
abline(lm(p[,3] ~ HolzingerSwineford1939$x9))
par(mfrow = c(1, 1))




n = 500
lat1 = rnorm(n)
lat2 = 0.5*rnorm(n)
envir = rep(c(0, 1), each = n/2)

ind1 = 1*lat1 + rnorm(n)
ind2 = 4*lat1 + rnorm(n)
ind3 = 0.5*lat1 + rnorm(n)
ind4 = (envir==0)*2*lat1 + (envir==1)*lat2 + rnorm(n)
ind4 = (envir==0)*2*lat1 + (envir==1)*rnorm(n, 0) + rnorm(n)
ind5 = 4*lat2 + rnorm(n)
ind6 = 0.5*lat2 + rnorm(n)
ind7 = 2*lat2 + rnorm(n)


df = data.frame(lat1, lat2,
                ind1, ind2, ind3, ind4, ind5, ind6, ind7)

mod <- ' 
l1 =~ ind1 + ind2 + ind3 + ind4
l2 =~ ind5 + ind6 + ind7
'
# fit the model
fit <- cfa(mod, data = df)
summary(fit)

preds = predict(fit)
lat1_pred = preds[,1]
lat2_pred = preds[,2]
plot(lat1_pred, ind1, col = envir + 1)
plot(lat1_pred, ind2, col = envir + 1)
plot(lat1_pred, ind3, col = envir + 1)
plot(lat1_pred, ind4, col = envir + 1)

fit.r4 = lm(ind4 ~ lat1_pred)
plot(fit.r4, which = 1, col = envir + 1)
R4 = resid(fit.r4)
mean(R4[envir == 0])
mean(R4[envir == 1])
var(R4[envir == 0])
var(R4[envir == 1])
bartlett.test(R4, envir)
t.test(R4[envir == 0], R4[envir == 1])
cor.test(lat1_pred, R4)
cor.test(lat1_pred[envir == 0], R4[envir == 0])
cor.test(lat1_pred[envir == 1], R4[envir == 1])

plot(lat1_pred, R4)

summary(lm(ind4 ~ preds*envir))
R4 = resid(lm(ind4 ~ preds))
plot(R4, col = envir + 1)

?dhsic.test
plot(lm(ind4 ~ preds), col = envir + 1, which = 1)
lm(R4 ~ envir)
mean.p = oneway.test(R4 ~ envir)$p.val
var.p = bartlett.test(R4 ~ envir)$p.val








library(here) # v0.1
library(tidyverse) # v1.3.0
library(MASS) # v7.3-51.6
library(emmeans) # v1.5.2-1
library(ggplot2) # v3.3.2
library(openxlsx) # v4.2.2
library(survey) # v4.0
library(ggpubr) # v0.4.0
library(kableExtra) # v1.3.1
library(stargazer) # v5.2.2

## Multiple testing correction ====
adjust = 'holm' # argument for emmeans::contrast

## read data ====
load("~/Dropbox/Far from the eye, far from the heart/KSFR21_Replication/data/KSFR21_data_clean.Rdata")


df = dtaDE
dtaDE$su

df = data.frame(ind1 = dtaDE$SuppCollab, 
                ind2 = dtaDE$PartyAttach, 
                ind3 = dtaDE$turnout, 
                ind4 = dtaDE$sticky,
                ind5 = dtaDE$EUAttOut,
                treat = dtaDE$treat, 
                level = dtaDE$level,
                TreatIND = dtaDE$TreatIND)
df$Ti = as.numeric(df$TreatIND)
df = na.exclude(df)
summary(df)
library(lavaan)

cor(df[,1:5])
mod = "f1 =~ ind1 + ind2 + ind3 + ind4 + ind5"

fit = cfa(model = mod, data = df)
summary(fit)

df$preds = predict(fit)

library(ggplot2)

ggplot(df, aes(x = preds)) +
  geom_histogram() +
  facet_wrap(~TreatIND)

plot(df$ind1, df$preds)
abline(lm(preds ~ ind1, data = df))



plot(df)
summary(df)

?na.omit()





8.192       6     0.2244 
pchisq(40.059, 6, lower.tail = F)






HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9'

# configural invariance
fit1 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school")

# weak invariance
fit2 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
            group.equal = "loadings")

#fit2 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
#            group.equal = "intercepts")

# strong invariance
fit3 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
            group.equal = c("intercepts", "loadings"))

# model comparison tests
lavTestLRT(fit1, fit2, fit3)
2 * (logLik(fit1) - logLik(fit2))
2 * (logLik(fit2) - logLik(fit3))

fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
lavTestLRT(fit1)
summary(fit1, fit.measures = T)
fitmeasures(fit1)
2 * (fit1@h1$loglik - fit1@loglik$loglik)
fit1@baseline
qchisq(2 * (fit1@h1$loglik - fit1@loglik$loglik), )
fit1@loglik$npar
fit1@h1$implied$cov

modsat = '
  x1 ~~ x1
  x2 ~~ x1 + x2
  x3 ~~ x1 + x2 + x3
  x4 ~~ x1 + x2 + x3 + x4
  x5 ~~ x1 + x2 + x3 + x4 + x5 
  x6 ~~ x1 + x2 + x3 + x4 + x5 + x6 
  x7 ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7 
  x8 ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 
  x9 ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 
'
fit.sat = cfa(modsat, data = HolzingerSwineford1939)
fit.sat
fit.sat@implied$cov
fit.sat@implied$cov

X = as.matrix(HolzingerSwineford1939[, 7:15])
X = scale(X, scale = F)
t(X) %*% X /301

cov(HolzingerSwineford1939[, 7:15])

HolzingerSwineford1939

mod0 <-'
  x1 ~~ x1
  x2 ~~ x2
  x3 ~~ x3
  x4 ~~ x4
  x5 ~~ x5
  x6 ~~ x6
  x7 ~~ x7
  x8 ~~ x8
  x9 ~~ x9
  
'
fit0 = cfa(mod0, data = HolzingerSwineford1939)
fit0@ParTable
fit0@loglik
fit0@h1$implied

?cfa
logLik(fit1)
getAnywhere(lavTestLRT())
sqrt(0.142/21)




## Countering DIF with multiple items
# note that % of items with DIF cannot be to large
# initially tried with 2/4 DIF on LV1, but this 
# messed with the overall results too much.

n = 500
lat1 = rnorm(n)
lat2 = 0.5*rnorm(n)
envir = rep(c(0, 1), each = n/2)

N = rnorm(n, 0)
ind1 = 1*lat1 + rnorm(n)
ind2 = 4*lat1 + rnorm(n)
#ind3 = (envir==0)*2*lat1 + (envir==1)*N + rnorm(n)
ind3 = -(envir==0)*2*lat1 - (envir==1)*N +  rnorm(n)
ind4 = 2*lat1 + rnorm(n)
ind4 = (envir==0)*2*lat1 + 2*(envir==1)*N + rnorm(n)
ind5 = 4*lat1 + rnorm(n)
ind6 = 0.5*lat1 + rnorm(n)
ind7 = 2*lat1 + rnorm(n)


df = data.frame(lat1, lat2,
                ind1, ind2, ind3, ind4, ind5, ind6, ind7)

mod <- ' 
l1 =~ ind1 + ind2 + ind3 + ind4 + ind5 + ind6 + ind7
'
# fit the model
fit <- cfa(mod, data = df)
summary(fit)

preds = predict(fit)
lat1_pred = preds[,1]

plot(lat1_pred ~ lat1, col = envir + 1)
abline(lm(lat1_pred ~ lat1))
abline(lm(lat1_pred[envir==0] ~ lat1[envir == 0]), col = 1)
abline(lm(lat1_pred[envir==1] ~ lat1[envir == 1]), col = 2)

plot(lat1_pred, ind1, col = envir + 1)
plot(lat1_pred, ind2, col = envir + 1)
plot(lat1_pred, ind3, col = envir + 1)
plot(lat1_pred, ind4, col = envir + 1)

fit.r4 = lm(ind4 ~ lat1_pred)
plot(fit.r4, which = 1, col = envir + 1)
R4 = resid(fit.r4)
mean(R4[envir == 0])
mean(R4[envir == 1])
var(R4[envir == 0])
var(R4[envir == 1])
bartlett.test(R4, envir)
t.test(R4[envir == 0], R4[envir == 1])
cor.test(lat1_pred, R4)
cor.test(lat1_pred[envir == 0], R4[envir == 0])
cor.test(lat1_pred[envir == 1], R4[envir == 1])

df$envir = envir
fit.grouped = cfa(mod, data = df, group = "envir")
fit.grouped_restrictedloadings = cfa(mod, data = df, group = "envir", group.equal = "loadings")
fit.grouped_restricted = cfa(mod, data = df, group = "envir", group.equal = c("intercepts","loadings"))

lavTestLRT(fit.grouped, fit.grouped_restrictedloadings, fit.grouped_restricted)

library(semPlot)
semPaths(fit, "par", weighted = FALSE, nCharNodes = 7, shapeMan = "rectangle",
         sizeMan = 8, sizeMan2 = 5)









library(lavaan)
library(car)
library(tidyverse)

n = 1000 # number of obs
k = 2 # number of latent vars
p = 10 # number of manifest vars
g = 2 # number of groups/populations/environments
grp = rep(1:g, each = floor(n/g)) # group index
gen_error = function(dist = "norm", n, var, ...){
  dist = paste0("r", dist)
  p = length(var)
  sd = sqrt(var)
  sapply(1:p, function(i) do.call(dist, args = list(n = n, sd = sd[i])))
}

eta = matrix(rnorm(n*k), ncol = k) # latent vars
eta = eta + grp - 1
tau = matrix(rep(0, p*g), nrow = p) # identical intercepts across groups
psi = matrix(rep(1, p*g), nrow = p) # identical variances across groups
Lambda = array(c(rep(1, k), rnorm(k*(p-1))), dim = c(k, p, g)) # identical loadings across groups

# violations of invariant loadings
Lambda[1,4,2]= 3 # fix loading of 1st latent variable on 4th indicator for 2nd group
Lambda[1,4,1]

# violations of different intercepts
tau[2, 2] = 2 # change intercept for 2nd indicator for 2nd group

# Generate Data
X = do.call("rbind", sapply(1:g, 
                            function(j) t(tau[,j] + t(Lambda[,,j]) %*% t(eta[grp == j,])) +
                              gen_error(dist = "norm", n = sum(grp == j), var = psi[,j]) , 
                            simplify = F))
X = as.data.frame(X)
colnames(X) = paste0("x", 1:10)
X$grp = grp

# fit CFA
mod = "eta1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10"
fit = cfa(mod, data = X)
summary(fit)

library(sirt)
mod = "eta1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
       eta1 ~~ 1*eta1"
fit.grp <- cfa(mod, data = X, group = "grp")
summary(fit.grp)
fit.grp@ParTable
fit.grp@ParTable$est[fit.grp@ParTable$group == 1 & fit.grp@ParTable$op == "=~"]
lel = summary(fit.grp)[[1]]
lel[1:10,]
lambda = rbind(fit.grp@ParTable$est[fit.grp@ParTable$group == 1 & fit.grp@ParTable$op == "=~"], 
               fit.grp@ParTable$est[fit.grp@ParTable$group == 2 & fit.grp@ParTable$op == "=~"])
nu = rbind(fit.grp@ParTable$est[fit.grp@ParTable$group == 1 & fit.grp@ParTable$op == "~1" & fit.grp@ParTable$lhs != "eta1"],
           fit.grp@ParTable$est[fit.grp@ParTable$group == 2 & fit.grp@ParTable$op == "~1" & fit.grp@ParTable$lhs != "eta1"])
n = 1000
wgt <- matrix(sqrt(500), 2, 10)
mod1 <- sirt::invariance.alignment(lambda, nu, wgt)
par = summary(mod1)
class(mod1)
lambda[,4]
mod1$lambda.aligned[,4]

nu
mod1$nu.aligned

mean(eta[grp == 1,1])


library(rrcov3way)

congruence(c(1, 2, 1), c(1, 1, 1))
X = matrix(rnorm(9), 3, 3)
Y =  matrix(rnorm(9), 3, 3)

sum(X * Y) / sqrt(sum(X^2) * sum(Y^2))
congruence(X)


library(tidyverse)
library(lavaan)
library(here)
source(here("Rscripts", "Simulator_PokropekEtAl.R"))
MMGFA_files = list.files(here("Rscripts", "KimDeRoover_MixtureMG_FA"))
MMGFA_files = MMGFA_files[grepl("R$", MMGFA_files)]
for(j in 1:length(MMGFA_files)){
  source(here("Rscripts", "KimDeRoover_MixtureMG_FA", MMGFA_files[j]))
}

# Simulate data
n = 50
p = 5
g = 24
h = 0.5
k = 2
D = sim_PMI(n = n, g = g, p = p, h = h, k = k, discrete = F)
Lambda_structure = matrix(rep(1, p), ncol = 1)
X = D$sim_dat[, grepl("^y", colnames(D$sim_dat))]
X = cbind(D$sim_dat$grp, X)

rm(fit)
fit = MixtureMG_FA(X, "loadings", nsclust = c(1, 3), nfactors = 1, Maxiter = 10000, nruns = 5, 
                   design = Lambda_structure, preselect = 50)

Cload = fit$MMGFA_solutions$`2.clusters`$clusterspecific.loadings[1,]
fit$MMGFA_solutions$`2.clusters`$clustermemberships
cbind(Cload[[1]], Cload[[2]])
Cload[[1]] - Cload[[2]]
D$p_affected
fit$MMGFA_solutions$`2.clusters`$clustermemberships %>% round
D$g_affected

Cload = fit$MMGFA_solutions$`3.clusters`$clusterspecific.loadings[1,]
fit$MMGFA_solutions$`3.clusters`$clustermemberships %>% round
D$g_affected
cbind(Cload[[1]], Cload[[2]], Cload[[3]])
Cload[[1]] - Cload[[2]] - Cload[[3]]
D$p_affected
fit$MMGFA_solutions$`2.clusters`$clustermemberships %>% round
D$g_affected


library(rrcov3way)
congruence(cbind(Cload[[1]], Cload[[2]], Cload[[3]]))
?congruence
