# Factor rotations
rot = function(angle){
  matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol = 2)
}
n = 1000
eta1 = rnorm(n)
eta2 = rnorm(n)
Y1 = 1*eta1 + 1*eta2 + rnorm(n)
Y2 = -1*eta1 - 0.4 * eta2 +  rnorm(n)
Y3 = 0.9*eta1 + 0.6*eta2 + rnorm(n)
Y4 = -0.3*eta1 + 0.3*eta2 + rnorm(n)
Y5 = 0.5*eta1 -0.2*eta2 + rnorm(n)
Y = cbind(Y1, Y2, Y3, Y4, Y5)
eta = cbind(eta1, eta2)

# Rotation
fit.fa = factanal(Y, 2)
plot(fit.fa$loadings %*% rot(0), asp = 1, xlim = c(-1, 1), ylim = c(-1, 1))
fit.fa = factanal(Y, 2, rotation = "none")
angle = 12.2
angle = 2*pi*angle / 360 
segments(x0 = rep(0, 5), y0 = rep(0, 5), 
         x1 = (fit.fa$loadings%*%rot(angle))[,1], y1 = (fit.fa$loadings%*%rot(angle))[,2])
points((fit.fa$loadings%*%rot(angle)), col = 2)

# scale indeterminacy
M = cov(matrix(rnorm(1000*5), ncol = 5))
M = diag(c(1, 2))
M %*% t(M)
solve(M)%*% M
M %*% solve(M)

fit.fa = factanal(Y, 2)
plot(fit.fa$loadings %*% rot(0), asp = 1, xlim = c(-1, 1), ylim = c(-1, 1))
fit.fa = factanal(Y, 2, rotation = "none")
angle = 12.2
angle = 2*pi*angle / 360 
segments(x0 = rep(0, 5), y0 = rep(0, 5), 
         x1 = (fit.fa$loadings%*%rot(angle))[,1], y1 = (fit.fa$loadings%*%rot(angle))[,2])
points((fit.fa$loadings%*%rot(angle)), col = 2)

# factor scores
fit.fa = factanal(Y, 2, rotation = "none", scores = "regression")
pred = Y %*% solve(fit.fa$correlation) %*% fit.fa$loadings
par(mfrow = c(1, 2))
plot(pred[,1], fit.fa$scores[,1])
plot(pred[,2], fit.fa$scores[,2])
par(mfrow = c(1, 1))

