rm(list = ls())

# Q3: Basic Capture-recapture
n1 = 22
n2 = 60
m2 = 11
D = c(n1, n2, m2)

S = 5000

nc = n1 + n2
ni = n1 + n2 - m2

#log probability of each value
valN = (ni:S)
lprob = lfactorial(valN) - lfactorial(valN - ni) +
  lfactorial(2 * valN - nc) - lfactorial(2 * valN + 1)

prob = c(rep(0, ni - 1), exp(lprob - max(lprob)))
prob = prob / sum(prob)

#Posterior mean
(mu = sum((1:S) * prob))

#Posterior standard deviation
sqrt(sum((1:S)^2 * prob) - mu ^ 2)

#Credible interval
cp = cumsum(prob)
ci1 = match(T, cp >= 0.025) - 1
ci2 = match(T, cp >= 0.975)
print(paste0("95% credible interval for N: [", ci1, " , ", ci2, "]"))

credint = function(D, S){
  n1 = D[1]
  n2 = D[2]
  m2 = D[3]
  
  nc = n1 + n2
  ni = n1 + n2 - m2
  
  #log probability of each value
  valN = (ni:S)
  lprob = lfactorial(valN) - lfactorial(valN - ni) +
    lfactorial(2 * valN - nc) - lfactorial(2 * valN + 1)
  
  prob = c(rep(0, ni - 1), exp(lprob - max(lprob)))
  prob = prob / sum(prob)
  
  #Posterior mean
  mu = sum((1:S) * prob)
  
  #Posterior variance
  sum((1:S) ^ 2 * prob) - mu ^ 2
  
  #Confidence interval
  cp = cumsum(prob)
  ci1 = match(T, cp >= 0.025) - 1
  ci2 = match(T, cp >= 0.975)
  print(paste0("95% credible interval for N: [", ci1, " , ", ci2, "]"))
}

credint(D, S)

# impact of S
credint(D, 200)
credint(D, 500)
credint(D, 5000)


#In-model validation
# Synthesize data and check that we recover the true parameter

SynthData = function(N, p = 0.5) {
  n1 = rbinom(1, N, p)
  m2 = rbinom(1, n1, p)
  n2 = m2 + rbinom(1, N - n1, p)
  #   Alternate version which samples at the individual level
  #   caught1 = rbinom(N, 1, p)
  #   caught2 = rbinom(N, 1, p)
  #   n1 = sum(caught1)
  #   n2 = sum(caught2)
  #   m2 = sum(caught1 * caught2)
  return(c(n1, n2, m2))
}

D = SynthData(200, 0.2)

credint(D, 500)

# Out-of-model validation for part 1: 
# what if individuals are not identically distributed?
# Simulate synthetic data with different probabilities per individual
SynthData2 = function(N, alpha = 1, beta = 1) {
  p = rbeta(N, alpha, beta)
  caught1 = rbinom(N, 1, p)
  caught2 = rbinom(N, 1, p)
  n1 = sum(caught1)
  n2 = sum(caught2)
  m2 = sum(caught1 * caught2)
  return(c(n1, n2, m2))
}

a = 20
b = 30
D = SynthData2(200, a, b)
curve(dbeta(x, a, b))

credint(D, 1000)

# Simulate synthetic data with different probabilities
# for different days
SynthData3 = function(N, p1, p2) {
  caught1 = rbinom(N, 1, p1)
  caught2 = rbinom(N, 1, p2)
  n1 = sum(caught1)
  n2 = sum(caught2)
  m2 = sum(caught1 * caught2)
  return(c(n1, n2, m2))
}

D = SynthData3(200, .5, .3)
credint(D, 1000)


# Q9: function to sample r1
accrej = function(n1, r3, m2, m3, p, q) {
  s1 = min(n1 - m2, n1 - r3 - m3)
  s2 = max(n1 - m2, n1 - r3 - m3)
  M = exp(lfactorial(n1) - lfactorial(s2))
  q2 = q / (q + (1 - p) ^ 2 * (1 - q) ^ 2)
  
  found = F
  while (!found) {
    x = rbinom(1, s1, q2)
    u = runif(1)
    if (u * M < exp(lfactorial(n1 - x) - lfactorial(s2 - x))) {
      found = T
    }
  }
  return(x)
}

# Gibbs' sampler using accept-reject

gibbsAR = function(niter, n1, m2, m3) {
  p = rep(NA, niter)
  q = rep(NA, niter)
  N = rep(NA, niter)
  r2 = rep(NA, niter)
  r3 = rep(NA, niter)
  
  #Initial values
  p[1] = 0.5
  q[1] = 0.5
  N[1] = 2 * n1
  r2[1] = 0
  r3[1] = 0
  pb = txtProgressBar(min = 0, max = niter, style = 3)
  
  for (i in 2:niter) {
    p[i] = rbeta(1, n1 + m2 + m3 + 1, N[i - 1] + n1 - 2 * r2[i - 1] - r3[i - 1] - m2 - m3 + 1)
    q[i] = rbeta(1, r2[i - 1] + r3[i - 1] + 1, 2 * n1 - 2 * r2[i - 1] -
                   r3[i - 1] + 1)
    N[i] = rnbinom(1, n1, p[i]) + n1
    r3[i] = rbinom(1, n1 - r2[i - 1] - m3, q[i] / (q[i] + (1 - q[i]) * (1 - p[i])))
    r2[i] = accrej(n1, r3[i], m2, m3, p[i], q[i])
    setTxtProgressBar(pb, i)
  }
  return(list(p, q, N, r2, r3))
}

run1 = gibbsAR(1e3, 22, 11, 6)


#Q10: Metropolis within Gibbs

#log density for r2
lpr1 = function(r2, n1, r3, m2, m3, p, q) {
  return(
    -lfactorial(r2) + lfactorial(n1 - r2) - lfactorial(n1 - r2 - m2)
    - lfactorial(n1 - r2 - r3 - m3) - 2 * r2 * log(1 - p) - 2 * r2 *
      log(1 - q) + r2 * log(q)
  )
}

#Symmetric proposal for r1: current value, plus or minus 1
gibbsMH = function(niter, n1, m2, m3) {
  p = rep(NA, niter)
  q = rep(NA, niter)
  N = rep(NA, niter)
  r2 = rep(NA, niter)
  r3 = rep(NA, niter)
  
  pb = txtProgressBar(min = 0, max = niter, style = 3)
  
  #Initial values
  p[1] = 0.5
  q[1] = 0.1
  N[1] = 2 * n1
  r2[1] = 1
  r3[1] = 1
  
  
  
  for (i in 2:niter) {
    p[i] = rbeta(1, n1 + m2 + m3 + 1, N[i - 1] + n1 - 2 * r2[i - 1] - r3[i -
                                                                           1] - m2 - m3 + 1)
    q[i] = rbeta(1, r2[i - 1] + r3[i - 1] + 1, 2 * n1 - 2 * r2[i - 1] -
                   r3[i - 1] + 1)
    N[i] = rnbinom(1, n1, p[i]) + n1
    r3[i] = rbinom(1, n1 - r2[i - 1] - m3, q[i] / (q[i] + (1 - q[i]) * (1 -
                                                                          p[i])))
    
    r2prop = r2[i - 1] + sample(c(-1, 1), 1)
    if (r2prop < 0 | n1 - r2prop - r3[i] - m3 < 0 | n1 - r2prop - m2 < 0) {
      r2[i] = r2[i - 1]
    }
    else{
      alpha = lpr1(r2prop, n1, r3[i], m2, m3, p[i], q[i]) - lpr1(r2[i - 1], n1, r3[i], m2, m3, p[i], q[i])
      if (runif(1) < exp(alpha)) {
        r2[i] = r2prop
      } else {
        r2[i] = r2[i - 1]
      }
    }
    
    setTxtProgressBar(pb, i)
  }
  return(list(p, q, N, r2, r3))
}

run2 = gibbsMH(1e5, 22, 11, 6)

par(mfrow = c(3, 1))
plot(run2[[3]], type = "l")
acf(run2[[3]])
hist(run2[[3]])

#ESS
1e5 / (2 * sum(acf(run2[[3]], plot = F)$acf))
mean(run2[[3]])
sd(run2[[3]])
library(HDInterval)
hdi(run2[[3]])

# Simulate synthetic data
SynthDataOpen = function(N, p, q){
  n1 = rbinom(1, N, p)
  r2 = rbinom(1, n1, q)
  m2 = rbinom(1, n1 - r2, p)
  r3 = rbinom(1, n1 - r2, q)
  m3 = rbinom(1, n1 - r2 - r3, p)
  return(c(n1, m2, m3))
}

D = SynthDataOpen(200, .3, .6)
run3 = gibbsMH(1e5, D[1], D[2], D[3])
hdi(run3[[3]], .95) # 95% credible interval

# Q14
# install.packages("gtools")
library(gtools)

eurodip = read.table("eurodip.txt")

gibbsAS = function(niter, z) {
  #entries in z are coded as 0 when not observed
  x = (z == 0)
  zc = z # "completed" version of z. Coded as R+1 if dead.
  n = dim(z)[1] #number of individuals
  K = dim(z)[2] #number of samplings
  R = max(max(z)) #number of locations
  
  pb = txtProgressBar(min = 0, max = niter, style = 3)
  
  p = array(NA, c(R + 1, niter))
  q = array(NA, c(R + 1, R + 1, niter))
  zc[z == 0] = sample.int(R)
  
  #Initial values
  p[, 1] = 0.1
  q[-(R + 1), , 1] = rdirichlet(R, rep(1, R + 1))
  
  #Once dead, you stay dead
  q[R + 1, R + 1, ] = 1
  q[R + 1, -(R + 1), ] = 0
  #Never capture dead individuals
  p[R + 1, ] = 0
  
  for (i in 2:niter) {
    setTxtProgressBar(pb, i)
    for (j in 1:R) {
      a1 = sum(sum(zc == j))
      a2 = sum(sum(zc == j & x == 1))
      p[j, i] = rbeta(1, a1 + 1, a2 + 1)
      
      b = rep(NA, R + 1)
      for (r in 1:(R + 1)) {
        b = sum(sum(zc[, -K] == j & zc[, -1] == r))
      }
      q[j, , i] = rdirichlet(1, b + 1)
    }
    
    for (l in 1:n) {
      if (x[l, 1] == 0) {
        zc[l, 1] = sample(R + 1, 1, prob = q[, zc[l, 2], i] * (1 - p[, i]))
      }
      for (k in 2:(K - 1)) {
        if (x[l, k] == 0) {
          zc[l, k] = sample(R + 1, 1, prob = q[zc[l, k - 1], , i] * (1 - p[, i]) *
                              q[, zc[l, k + 1], i])
        }
      }
      if (x[l, K] == 0) {
        zc[l, K] = sample(R + 1, 1, prob = q[zc[l, K - 1], , i] * (1 - p[, i]))
      }
    }
  }
  return(list(p, q))
}


run3 = gibbsAS(1e3, eurodip)
