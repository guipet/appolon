rm(list=ls())
deputes=read.csv2('deputes2019.csv', fileEncoding = "UTF-8")
View(deputes)
attach(deputes)    #permet d'éviter de faire liste$colonne

#Q1
#Data exploration
summary(questions_orales)
head(questions_orales, 20)

hist(questions_orales, breaks=(0:16-.5))
barplot(table(questions_orales))
(n = length(questions_orales))
(nh = sum(sexe=="H"))
(nf = sum(sexe=="F"))
(qtot = sum(questions_orales))           #la somme des observations Yi
(qh = sum(questions_orales[sexe=="H"]))  #la somme des observations sachant que c'est des hommes
(qf = sum(questions_orales[sexe=="F"]))  #la somme des observations sachant que c'est des femmes


#Is the Poisson model suitable?
lambdahat = qtot/n # MLE
plot(dnorm(0:15, lambdahat))
table(questions_orales) # observed counts
round(n*dpois(0:15, lambdahat)) #expected counts
chisq.test(table(factor(questions_orales, 
                        levels=0:5)), 
           round(n*dpois(0:5, lambdahat))) # large p-value: accept H0 that the data come from a Poisson

hist(questions_orales, freq=FALSE, breaks=0:16-.5)
points(0:15, dpois(0:15, lambdahat), col=2)

#Q5
#For a prior Gamma(2, 2)
a = 0.5
b = 0
par(mfrow=c(3, 1))
curve(dgamma(x, a, b), xlim=c(0, 5), main="Prior", ylab="density")
curve(dgamma(x, a+qtot, b+n), xlim=c(0, 5), main="Posterior model 1", ylab="density")
curve(dgamma(x, a+qh, b+nh), xlim=c(0, 5), main= "Posterior model 2", ylab="density")
curve(dgamma(x, a+qf, b+nf), col=2, add=T)
legend("topright", c("H", "F"), col=1:2, lty=1)

# Q6: credibility intervals
qgamma(c(.025, .975), a, b) # prior
qgamma(c(.025, .975), a+qtot, b+n) # model 1
qgamma(c(.025, .975), a+qh, b+nh) # model 2: lambda1
qgamma(c(.025, .975), a+qf, b+nf) # model 2: lambda2 

# Q7: Monte Carlo estimation of r_lambda
niter = 1e5

# prior
a = 2
b = 2
rlambda = rgamma(niter, a, b) / rgamma(niter, a, b)
mean(rlambda)
var(rlambda)
quantile(rlambda, c(.025, .975))
hist(rlambda, breaks=50)
plot(density(rlambda), xlim=c(0, 10))

# posterior
rlambdab = rgamma(niter, a+qh, b+nh) / 
  rgamma(niter, a+qf, b+nf)
mean(rlambdab)
var(rlambdab)
quantile(rlambdab, c(.025, .975))

# HPD intervals
library(HDInterval)
hdi(rlambda, .95)
hdi(rlambdab, .95)

#Q8: Vanilla Monte Carlo
# Create functions for the likelihood
lkd.model1 = function(y, n, lambda){
  return(exp(-n*lambda + y*log(lambda)))      #tout mettre sur l'exponentielle pour avoir des résultats
}

lkd.model2 = function(y1, n1, y2, n2, lambda1, lambda2){
  return(exp(-n1*lambda1 + y1*log(lambda1) -
               n2*lambda2 + y2*log(lambda2)))
  #return(lambda1^y1*exp(-n1*lambda1)*lambda2^y2*exp(-n2*lambda2))
}

BF_MC = function(a, b, y1, n1, y2, n2, M){
  lambda = rgamma(M, a, b)
  m1 = cumsum(lkd.model1(y1+y2, n1+n2, lambda))/(1:M)
  lambda1 = rgamma(M, a, b)
  lambda2 = rgamma(M, a, b)
  m2 = cumsum(lkd.model2(y1, n1, y2, n2, lambda1, lambda2))/(1:M)
  return(m2/m1)
}

M = 1e6
#a=2; b=2; y1=qh; n1=nh; y2=qf; n2=nf
MCest = BF_MC(2, 2, qh, nh, qf, nf, M)
MCest[M]
par(mfrow=c(1, 1))
plot(200000:M, MCest[200000:M], type="l")
#abline(h=trueBF, col=2)


#Q9: Importance sampling
BF_IS = function(a, b, y1, n1, y2, n2, M){
  mean.mod1 = (a+y1+y2)/(b+n1+n2)
  mean.mod2.lambda1 = (a+y1)/(b+n1)
  mean.mod2.lambda2 = (a+y2)/(b+n2)
  
  sigma.mod1 = sqrt((a+y1+y2)/(b+n1+n2)^2)
  sigma.mod2.lambda1 = sqrt((a+y1)/(b+n1)^2)
  sigma.mod2.lambda2 = sqrt((a+y2)/(b+n2)^2)
  lambda = rnorm(M, mean.mod1, sigma.mod1)          #sample de notre nouvelle loi
  m1 = cumsum( lkd.model1(y1+y2, n1+n2, lambda) * 
                 dgamma(lambda, a, b) / 
                 dnorm(lambda, mean.mod1, sigma.mod1))/(1:M)
  
  lambda1 = rnorm(M, mean.mod2.lambda1, sigma.mod2.lambda1)
  lambda2 = rnorm(M, mean.mod2.lambda2, sigma.mod2.lambda2)
  m2 = cumsum( lkd.model2(y1, n1, y2, n2, lambda1, lambda2) * 
                 dgamma(lambda1, a, b) * 
                 dgamma(lambda2, a, b) / 
                 (dnorm(lambda1, mean.mod2.lambda1, sigma.mod2.lambda1)*
                    dnorm(lambda2, mean.mod2.lambda2, sigma.mod2.lambda2)))/(1:M)
  
  return(m2/m1)
}

ISest = BF_IS(2, 2, qh, nh, qf, nf, M)
ISest[M]
plot(100:1000, ISest[100:1000], type="l")
#abline(h=trueBF, col=2)

# Note that the Importance Sampling estimator converges much faster than the vanilla Monte Carlo. 


#Q10: Harmonic mean estimator
BF_HM = function(a, b, y1, n1, y2, n2, M){
  lambda1 = rgamma(M, a+y1+y2, b+n1+n2)
  m1 = (1:M)/(cumsum(1/lkd.model1(y1+y2, n1+n2, lambda1)))
  lambda2.1 = rgamma(M, a+y1, b+n1)
  lambda2.2 = rgamma(M, a+y2, b+n2)
  m2 = (1:M)/cumsum(1/lkd.model2(y1, n1, y2, n2, lambda2.1, lambda2.2))
  return(m2/m1)  
}

HMest = BF_HM(2, 2, qh, nh, qf, nf, M)
HMest[M]
plot(100:M, HMest[100:M], type="l")
#abline(h=trueBF, col=2)
# This estimator does not seem to converge. Indeed, there is not guarantee that it has finite variance.


#Q11: exact computation
BF_analytical = function(a, b, y1, n1, y2, n2){
  m1 = b^a/gamma(a)*gamma(a+y1+y2)/(b+n1+n2)^(a+y1+y2)
  m2 = b^(2*a)/gamma(a)^2*gamma(a+y1)/(b+n1)^(a+y1)*
    gamma(a+y2)/(b+n2)^(a+y2)
  return(m2/m1)
}

# This does not work because n is too large
BF_analytical(2, 2, qh, nh, qf, nf)

# It works on the log scale
BF_analytical2 = function(a, b, y1, n1, y2, n2){
  m1 = a*log(b) - lgamma(a) + lgamma(a+y1+y2) - (a+y1+y2)*log(b+n1+n2)
  m2 = 2*a*log(b) - 2*lgamma(a) + lgamma(a+y1) - (a+y1)*log(b+n1) + 
    lgamma(a+y2) - (a+y2)*log(b+n2)
  return(exp(m2-m1))
}
(trueBF = BF_analytical2(2, 2, qh, nh, qf, nf))

# Q12
log10(trueBF) #"strong" evidence in favor of model 1

# Q13 posterior probability of model 1 and model 2
p1 = 1/(trueBF+1)
p2 = trueBF/(trueBF+1)

#Q12 Draw of the posterior sample for lambda
nsamp = 1e5
lambdasampa = rgamma(nsamp, 2+qtot, 2+n)
lambdasampb = rgamma(nsamp, 2+qf, 2+nf)

lambdasampc = rep(NA, nsamp)
for(i in 1:nsamp){
  if(runif(1)<p2){
    lambdasampc[i] = rgamma(1, 2+qf, 2+nf)
  }
  else{
    lambdasampc[i] = rgamma(1, 2+qtot, 2+n)
  }
}

par(mfrow=c(3, 1))
plot(density(lambdasampa), xlim=c(.5, 4), main="Model 1")
plot(density(lambdasampb), xlim=c(.5, 4), main="Model 2")
plot(density(lambdasampc), xlim=c(.5, 4), main="Weighted models")

par(mfrow=c(1, 1))
plot(density(lambdasampa), xlim=c(2, 4))
points(density(lambdasampb), type="l", col=3)
points(density(lambdasampc), type="l", col=2)

curve(dgamma(x, 2+qtot, 2+n), xlim=c(2, 4))
curve(dgamma(x, 2+qf, 2+nf), add=T, col=3)
curve(p2*dgamma(x, 2+qf, 2+nf) + (1-p2)*dgamma(x, 2+qtot, 2+n), col=2, add=T)


