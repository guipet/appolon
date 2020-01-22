rm(list=ls())
d = read.csv("titanic.csv")

# Q1
summary(d)

y = d$Survived

X = data.frame(intercept=1, class1=(d$Pclass==1), class2=(d$Pclass==2), 
               sex=(d$Sex=="male"), age=d$Age, sibsp=d$SibSp, 
               parch=d$Parch, fare=d$Fare, EmbarkedS=(d$Embarked=="S"), 
               EmbarkedC=(d$Embarked=="C"))

summary(X)
# For now, remove the passengers whose age is unknown
X = X[!is.na(d$Age), ]
y = y[!is.na(d$Age)]

X = as.matrix(X)
q = dim(X)[2] # in this code, q corresponds to the value q+1 in the question sheet
n = dim(X)[1]

# install.packages("mvtnorm") # For the multivariate normal
require(mvtnorm)

# Q2
post_flat = function(beta, y=y, X=X){
  #For a flat prior, the posterior is proportional to the likelihood
  prob = pnorm(X%*%matrix(beta, ncol=1))
  return(prod(prob^y * (1-prob)^(1-y)))
}

# Q3
mle = summary(glm(y~X-1, family=binomial(link="probit")))
betahat = mle$coefficients[ , 1]
sigma.asymp = mle$cov.unscaled

# Q4a: Metropolis-Hastings
MH = function(niter, tau, y, X, mlestart=T){
  q = dim(X)[2]
  beta = matrix(0, nrow=niter, ncol=q)
  mle = summary(glm(y~X-1, family=binomial(link="probit")))
  if(mlestart){
    betahat = mle$coefficients[, 1]
    sigma.asymp = mle$cov.unscaled
    beta[1, ] = betahat
  }
  else{
    beta[1, ] = 0
  }
  
  for(k in 2:niter){
    proposal = rmvnorm(1, beta[k-1, ], tau^2*sigma.asymp)
    alpha = min(1, post_flat(proposal, y, X) / 
                  post_flat(beta[k-1, ], y, X))
    if(runif(1)<alpha){
      beta[k, ] = proposal
    } 
    else{
      beta[k, ] = beta[k-1, ]
    }
  }
  return(beta)
}

b1 = MH(1e4, 0.1, y, X)
b2 = MH(1e4, 0.7, y, X)
b3 = MH(1e4, 2, y, X)

# Q4b: Plot the trajectory of one parameter, the autocorrelation, and a histogram
par(mfcol=c(3, 3))
i=1 # Change this to look at another parameter
plot(b1[, i], type="l")
plot(b2[, i], type="l")
plot(b3[, i], type="l")

# need to play with the lag to see the auto-correlation go to 0
acf(b1[, i], lag.max=500)
acf(b2[, i], lag.max=50)
acf(b3[, i], lag.max=1000)

hist(b1[1e3:1e4, i], breaks=50)
hist(b2[1e3:1e4, i], breaks=50)
hist(b3[1e3:1e4, i], breaks=50)

# Remove 1000 iterations for burn-in, and thin to every 10th iteration
post_sample = b2[seq(1e3, 1e4, by=10), ]
par(mfcol=c(3, 1))
plot(post_sample[ , i], type="l")
acf(post_sample[ , i])
hist(post_sample[ , i], breaks=50)

# Effective sample size
1e4/(2*sum(acf(b1[, i], lag.max=400, plot=F)$acf)-1)
1e4/(2*sum(acf(b2[, i], lag.max=50, plot=F)$acf)-1)
1e4/(2*sum(acf(b3[, i], lag.max=1000, plot=F)$acf)-1)
