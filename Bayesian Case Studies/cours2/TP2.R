rm(list=ls())
d = read.csv("deathrate2.csv")
y = d[ , 2]
X = as.matrix(d[ , -c(1, 2)])

n = length(y)
p = ncol(X)

#Q1
summary(lm(y~X))
betahat = (lm(y~X))$coefficients
residuals = lm(y~X)$residuals
s2 = summary(lm(y~X))$sigma^2

X = cbind(1, X) # add a column of 1s for the intercept

#Q2a
g = c(.1, 1, 10, 100, 1000)
#mean beta0:
betahat[1] * g/(g+1)
#mean for sigma^2
a = n/2
b = s2/2 + 1/(2*g+2)*t(betahat)%*%t(X)%*%X%*%betahat
b/(a-1)

#Q2b
g = n
q = 2
X0 = X[, -(8:9)]
BF = (g+1)^(-q/2) * 
  ((t(y)%*%y-g/(g+1)*t(y)%*%X0%*%solve(t(X0)%*%X0)%*%t(X0)%*%y)/
     (t(y)%*%y-g/(g+1)*t(y)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%y))^(n/2)
log10(BF) # Strong evidence in favour of H1

# Q3
marglkd = function(gamma, mat=X, g=nrow(mat)){
  q = sum(gamma)
  X1 = mat[, c(T, gamma)]
  if(q==0){return(-(q+1)/2*log(g+1) -n/2*log(t(y)%*%y))}
  m = -q/2*log(g+1) -
    n/2*log(t(y)%*%y - g/(g+1)* t(y)%*% X1 %*%solve(t(X1)%*%X1) %*%
              t(X1)%*%y)
  return(m)
}


marglkd(c(F, F, F), X0)
marglkd(c(F, F, T), X0)
marglkd(c(F, T, F), X0)
marglkd(c(F, T, T), X0)
marglkd(c(T, F, F), X0)
marglkd(c(T, F, T), X0)
marglkd(c(T, T, F), X0)
marglkd(c(T, T, T), X0)


#Q6

niter = 1e5 # number of iterations
gamma = matrix(F, nrow=niter, ncol=p)
gamma0 = sample(c(T, F), size=p, replace=TRUE) #random initial value of the Gibbs sampler
lkd = rep(0, niter)
modelnumber = rep(0, niter)

oldgamma = gamma0
for(i in 1:niter){
  newgamma = oldgamma
  for(j in 1:p){
    g1 = newgamma; g1[j] = TRUE
    g2 = newgamma; g2[j] = FALSE
    ml1 = marglkd(g1)
    ml2 = marglkd(g2)
    alpha = c(ml1, ml2) - min(ml1, ml2)
    # We want to draw from a Bernoulli, with probability of drawing TRUE equal to exp(p[1])/(exp(p[1])+exp(p[2])).
    # The next line does this. Note that the function sample takes care of the normalization.
    newgamma[j] = sample(c(T, F), size=1, prob=exp(alpha)) 
  }
  gamma[i, ] = newgamma
  lkd[i] = marglkd(newgamma)
  modelnumber[i] = sum(newgamma*2^(0:(p-1)))
  oldgamma = newgamma
}

apply(gamma, 2, "mean")

#Check mixing with autocorrelations
par(mfrow=c(4, 2))
for(i in 2:9) acf(as.numeric(gamma[, i]))
# Autocorrelation decreases quickly. No need to subsample.

#Check convergence+mixing with the trace (sliding mean since the variables are binary)
require(zoo)
for(i in 2:9) plot(rollapply(gamma[, i], width=100, FUN=mean), type="l")

burnin = 500 #500 burn-in steps
gammab = modelnumber[burnin+1:niter] 
res = as.data.frame(table(gammab))
odo = order(res$Freq, decreasing=T)[1:10]
modcho = res$gammab[odo]
probtop10 = res$Freq[odo]/(niter-burnin)

indices = match(modcho, modelnumber)
cbind(probtop10, gamma[indices, ])