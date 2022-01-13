# FACTOR EXTRACTION MANUALLY AND VIA LAVAAN PREDICT

# generate data
library(lavaan)
n = 1000
eta1 = rnorm(n)
eta2 = rnorm(n)

Y1 = eta1 + rnorm(n, 0.2)
Y2 = 2 + 0.5*eta1 + rnorm(n, 0.2)
Y3 = 1.5*eta1 + rnorm(n, 0.2)
Y4 = 0.5*eta1 + 0.5*eta2 + rnorm(n, 0.2)
Y5 = eta2 + rnorm(n, 0.2)
Y6 = 2*eta2 + rnorm(n, 0.2)
df = data.frame(Y1, Y2, Y3, Y4, Y5, Y6)

# model
m = "
eta1 =~ Y1 + Y2 + Y3 + Y4
eta2 =~ Y5 + Y6 + Y4
"
summary(fit <- cfa(m, df))

# Regression method by hand: (unstandardized outputs a bit difficult to extract, so code is pretty ugly)
(identifier = paste0(fit@ParTable$lhs, fit@ParTable$op, fit@ParTable$rhs))
phi = matrix(fit@ParTable$est[c(which(identifier == "eta1~~eta1"),
                                which(identifier == "eta1~~eta2"),
                                which(identifier == "eta1~~eta2"),
                                which(identifier == "eta2~~eta2"))], ncol = 2)
phi
lambda = matrix(c(fit@ParTable$est[c(which(identifier == "eta1=~Y1"),
                                     which(identifier == "eta1=~Y2"),
                                     which(identifier == "eta1=~Y3"),
                                     which(identifier == "eta1=~Y4"))],
                  0, 0, 0, 0, 0,
                  fit@ParTable$est[c(which(identifier == "eta2=~Y4"),
                                     which(identifier == "eta2=~Y5"),
                                     which(identifier == "eta2=~Y6"))]), ncol = 2)
lambda
sigma = fit@implied$cov[[1]]
#lambda %*% phi %*% t(lambda) + psi
B = solve(sigma) %*% lambda %*% phi
Y = as.matrix(df)
Ytilde = scale(Y, scale = F)
eta = Ytilde %*% B

# compare with regression method via lavaan
all.equal(eta[,1], predict(fit)[,1])
plot(eta[,1], predict(fit)[,1])
plot(eta[,2], predict(fit)[,2])


#By hand Bartlett method
psi = diag(fit@ParTable$est[c(which(identifier == "Y1~~Y1"),
                              which(identifier == "Y2~~Y2"),
                              which(identifier == "Y3~~Y3"),
                              which(identifier == "Y4~~Y4"),
                              which(identifier == "Y5~~Y5"),
                              which(identifier == "Y6~~Y6"))])
Bt = solve(expm::sqrtm(t(lambda) %*% solve(psi) %*% sigma %*% solve(psi) %*% lambda)) %*% t(lambda) %*% solve(psi)
eta_bartlett = as.matrix(df) %*% t(Bt)

  # some difference remains, perhaps because of use of sqrtm
plot(eta_bartlett[,1], lavPredict(fit, method = "bartlett")[,1])

















