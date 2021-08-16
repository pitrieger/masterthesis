# load the lavaan package (only needed once per session)
library(lavaan)
# specify the model
HS.model <- ' visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9 '
# fit the model
fit <- cfa(HS.model, data = HolzingerSwineford1939)
# display summary output
summary(fit, fit.measures = TRUE)

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
