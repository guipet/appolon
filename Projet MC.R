install.packages("plot3D")
library(plot3D)

#exercice 1
n=10000
par(mfrow=c(2,2))
ftilde=function(x){
  ifelse(x[1]<=4 & x[1]>=0 & x[2]>=0 & x[2]<=2, (cos(x[1])**2+2*(sin(x[2])**2)*cos(x[1])**4)/(1+4*(x[2]-1)**2)*exp(-0.5*(x[1]-2)**2),0)
}


#1er rejet: on maximise par rapport à normale 2,1 en x et par une uniforme 0,2 en y
#constante de normalisation
M=6*sqrt(2*pi) 

#densité instrumentale
g=function(x){
  return(dnorm(x[1],2,1)*dunif(x[2],0,2))
}

#algo du rejet
rejet1=function(n){
  solution=matrix(numeric(2*n), ncol=2) #création d'une matrice à n lignes et 2 colonnes que de 0 qui va à chaque ligne être composé de x,y à la fin du rejet
  for (i in 1:n){
    u=runif(1,0,1)
    m=c(rnorm(1,2,1),runif(1,0,2))
    beta=ftilde(m)/(M*g(m))
    while (u>beta){
      u=runif(1,0,1)
      m=c(rnorm(1,2,1),runif(1,0,2))
      beta=ftilde(m)/(M*g(m))
    }
    solution[i,]=m
  }
  return(solution)
}

#on plot les vectuer sen 2D car impossible d'installer package 3D
testrejet1=rejet1(n)
plot(testrejet1, col="red")

#proportion de x entre 0 et 4 (car on a vu que les 2 sont plutot bien répartis) hist ou barplot??
hist(testrejet1[,1], col='red', main='1er rejet', xlab="réalisation des X", ylab='')




#2ème rejet: normale 2,1 en x et par une cauchy 1,0.5 en y
#constante de normalisation
M2=1.5*pi*sqrt(2*pi)

#densité instrumentale
gbis=function(x){
  return(dnorm(x[1],2,1)*dcauchy(x[2],1,0.5))
}

#algo du rejet
rejet_bis=function(n){
  solution_bis=matrix(numeric(2*n), ncol=2)
  for (i in 1:n){
    u=runif(1,0,1)
    m=c(rnorm(1,2,1),rcauchy(1,1,0.5))
    beta=ftilde(m)/(M2*gbis(m))
    while (u>beta){
      u=runif(1,0,1)
      m=c(rnorm(1,2,1),rcauchy(1,1,0.5))
      beta=ftilde(m)/(M2*gbis(m))
    }
    solution_bis[i,]=m
  }
  return(solution_bis)
}

#on plot les vectuer sen 2D car impossible d'installer package 3D
testrejet_bis=rejet_bis(n)


#test rejet 3: uniforme 0,4 pour x et uniforme 0,2 pour y
#constante de normalisation
M3=24

#densite instrumentale
gterce=function(x){
  return(dunif(x[1],0,4)*dunif(x[2],0,2))
}

#fonction du rejet
rejet_terce=function(n){
  solution_terce=matrix(numeric(2*n), ncol=2)
  for (i in 1:n){
    u=runif(1,0,1)
    m=c(runif(1,0,4),runif(1,0,2))
    beta=ftilde(m)/(M2*gterce(m))
    while (u>beta){
      u=runif(1,0,1)
      m=c(runif(1,0,4),runif(1,0,2))
      beta=ftilde(m)/(M2*gterce(m))
    }
    solution_terce[i,]=m
  }
  return(solution_terce)
}

#on plot les vectuer sen 2D car impossible d'installer package 3D
testrejet_terce=rejet_terce(n)
points(testrejet_terce, col="blue")

##proportion de x entre 0 et 4 (car on a vu que les 2 sont plutot bien répartis)
hist(testrejet_terce[,1], col='blue', main='rejet 3', xlab="réalisation des X", ylab='')



#rejet numéro 4 

#constante de normalisation
M4=6*pi

#densité instrumentale
gquatre=function(x){
  return(dunif(x[1],0,4)*dcauchy(x[2],1,0.5))
}

#algorithme du rejet
rejetquatre=function(n){
  solution_quatre=matrix(numeric(2*n), ncol=2)
  for (i in 1:n){
    u=runif(1,0,1)
    m=c(runif(1,0,4),rcauchy(1,1,0.5))
    beta=ftilde(m)/(M4*gquatre(m))
    while (u>beta){
      u=runif(1,0,1)
      m=c(runif(1,0,4),rcauchy(1,1,0.5))
      beta=ftilde(m)/(M4*gquatre(m))
    }
    solution_quatre[i,]=m
  }
  return(solution_quatre)
}

#on plot les vectuer sen 2D car impossible d'installer package 3D
testrejet_quatre=rejetquatre(n)
points(testrejet_quatre, col="grey")

##proportion de x entre 0 et 4 (car on a vu que les 2 sont plutot bien répartis)
hist(testrejet_quatre[,1], col='grey', main='rejet 4', xlab="réalisation des X", ylab='')




#Algorithme de Metropolis-Hastings question 4
#expliquer que le denominateur s'annule jamais, que le support de f est inclus dans g.


alpha=function(xt,eps){
  bob=ftilde(eps)*g(xt)/(ftilde(xt)*g(eps))
  return(min(bob,1))
}

#quelles valeurs pour x0 et justifier
#espilon varie-t-il en fonction du temps?

MH1=function(n,x0=c(0,0)){
  x0=as.vector(x0)
  X=matrix(numeric(2*n), ncol=2)
  for (i in 1:n){
    eps=c(rnorm(1,2,1),runif(1,0,2))
    p=alpha(x0,eps)
    u=rbinom(1,1,p)
    X[i,]=eps*u+(1-u)*x0
    x0=X[i,]
  }
  return(X)
}


##############################

alphabis=function(xt,eps){
  bob=ftilde(eps)*gbis(xt)/(ftilde(xt)*gbis(eps))  #par le calcul, simuler f ou ftilde revient au même car f(xt/f(eps)=ftilde(xt)/ftilde(eps)
  return(min(bob,1))
}
#quelles valeurs pour x0 et justifier
#espilon varie-t-il en fonction du temps?
MHbis=function(n,x0=c(0,0)){
  x0=as.vector(x0)
  X=matrix(numeric(2*n), ncol=2)
  for (i in 1:n){
    eps=c(rnorm(1,2,1),rcauchy(1,1,0.5))
    p=alphabis(x0,eps)
    u=rbinom(1,1,p)
    X[i,]=eps*u+(1-u)*x0
    x0=X[i,]
  }
  return(X)
}


##############################


alphaterce=function(xt,eps){
  bob=ftilde(eps)*gterce(xt)/(ftilde(xt)*gterce(eps))  #par le calcul, simuler f ou ftilde revient au même car f(xt/f(eps)=ftilde(xt)/ftilde(eps)
  return(min(bob,1))
}
#quelles valeurs pour x0 et justifier
#espilon varie-t-il en fonction du temps?
MHterce=function(n,x0=c(0,0)){
  x0=as.vector(x0)
  X=matrix(numeric(2*n), ncol=2)
  for (i in 1:n){
    eps=c(runif(1,0,4),runif(1,0,2))
    p=alphaterce(x0,eps)
    u=rbinom(1,1,p)
    X[i,]=eps*u+(1-u)*x0
    x0=X[i,]
  }
  return(X)
}


###########################


alphaquatre=function(xt,eps){
  bob=ftilde(eps)*gquatre(xt)/(ftilde(xt)*gquatre(eps))  #par le calcul, simuler f ou ftilde revient au même car f(xt/f(eps)=ftilde(xt)/ftilde(eps)
  return(min(bob,1))
}
#quelles valeurs pour x0 et justifier
#espilon varie-t-il en fonction du temps?
MHquatre=function(n,x0=c(0,0)){
  x0=as.vector(x0)
  X=matrix(numeric(2*n), ncol=2)
  for (i in 1:n){
    eps=c(runif(1,0,4),rcauchy(1,1,0.5))
    p=alphaquatre(x0,eps)
    u=rbinom(1,1,p)
    X[i,]=eps*u+(1-u)*x0
    x0=X[i,]
  }
  return(X)
}








##################################################################################
#réinitialiser la page plot.
#Exercice 2, partie 1, question 1,a)
#estimation de delta par MC classique
n=10000
#fonction pour l'estimateur de MC
indicatrice=function(x,t=2){
  ifelse(x>=t,1,0)
}

#estimateur
W=indicatrice(rweibull(n,2))
deltaMC=mean(W)

#graphe de la convergence de l'estimateur
plot(cumsum(W)/(1:length(W)),type='l', col='blue', xlab="n", ylab="y", ylim=c(0,0.022), main="convergence de l'estimateur")
abline(h=1-pweibull(2,2), col='red')
legend("bottomright", legend=c("estimateur", "valeur exacte"), col=c('blue','red'), lty=c(1), lwd=c(3), title="courbes", bg='aliceblue')


#intervalle de confiance
ICdeltaMC=c(deltaMC-qnorm(0.975)*sqrt(var(W)/n),deltaMC+qnorm(0.975)*sqrt(var(W)/n))


#stratification
#on crée une fonction qui renvoie une matrice. 1ère ligne, les sommes de h(xl) et 2ème ligne la variance de chaque échantillon
Strat=function(n,L){
  if (n%%L!=0){
    return("Impossible, L doit diviser n")
  }
  a=seq(0,1,length.out=L)
  b=qweibull(a,2)
  premieresomme=c()
  variance=c()
  mat=matrix(numeric(2*L), nrow=2)
  u=runif(n/L,0,1)
  for (l in 1:(L-1)){
    xl=qweibull((l-1)/L+u/L,2)
    hxl=indicatrice(xl)
    mat[,l]=c(sum(hxl), var(hxl))
  }
  return(mat)
}

#application avec n=10000 et L=1000
echantillonstrat=Strat(10000,1000)

#on prend les sommes et on crée des échantillons
ff=echantillonstrat[1,]
estistrat=(1/10000)*sum(ff)

#on prend la variance et on crée l'intervalle de confiance
variance=mean(echantillonstrat[2,])
ICstrat=c(estistrat-qnorm(0.975)*sqrt(variance/10000),estistrat+qnorm(0.975)*sqrt(variance/10000))

#plot la convergence

#si j'essaie de garder le même intervalle que dans le td, je me retrouve avec des puissances trop petites encore.



#Partie 1 question 2,a
#on simule x3 par la méthode de l'inverse
F.inv=function(x){
  ifelse(0<x & x<0.25, 4*x, ifelse(x<=1 & x>=0.75, 4*x-2, 1))
}


#fonction pour estimer l'intégrale (la renommer du même nom que dans le rapport)
indicatrice1=function(x,t=1){
  ifelse(x>=t, 1, 0)
}

#X1+X2 suit une loi gamma(2,1). On fera appel directement à rgamma
#on crée une fonction qui revoie l'échantillon
MC.esti=function(n){
  u=runif(n)
  x3=F.inv(u)
  t=rgamma(n,2,1)
  m=colSums(matrix(c(t,x3), nrow=2, byrow=TRUE)) #matrice à 2 lignes. 1ere ligne, les n gamma et 2ème ligne les n x3. ON aditionne les colonnes, pour obtenir 100 réalisations de X1+X2+X3
  return(indicatrice1(m))
}

#estimateur et convergence
h=MC.esti(n)
f=mean(h)
plot(cumsum(h)/(1:length(h)),type='l', col='chartreuse', xlab="n", ylab="y", main="convergences des estimateurs", ylim=c(0.90,1))
legend("bottomright", legend=c("estimateur monte carlo"), col=c('chartreuse'), lty=c(1), lwd=c(3), bg='aliceblue')

#intervalle de confiance
IC1=c(f-qnorm(0.975)*sqrt(var(h)/n),f+qnorm(0.975)*sqrt(var(h)/n))






#exo2, partie1, question c)
#prend en entrée des vecteurs et renvoie le vecteur transformé
F.=function(x){
  ifelse( x>=0 & x<1, x/4, ifelse( x<2 & x>=1, x/4+0.5, ifelse(x>=2, 1, 0)))
}

#fonction qui crée h1 pour l'estimateur de MC
estiF=function(n,t=1){
  m=rgamma(n,2,1)
  a=F.(t-m)            # on crée un vecteur de F(t-X1-X2)
  return(1-a)
}

#estimateur
h1=estiF(n)
f1=mean(h1)

#intervalle de confiance
IC1=c(f1-qnorm(0.975)*sqrt(var(h1)/n), f1+qnorm(0.975)*sqrt(var(h1)/n))
points(cumsum(h1)/(1:length(h1)),type='l', col='blue')
legend("bottomright", legend=c("estimateur MC classique", "estimateur MC avec h1"), col=c('chartreuse', 'blue'), lty=c(1), lwd=c(3), bg='aliceblue')


#convergence de l'estimateur 

#on fait avec G
#fonction qui crée la fonction h2
estiG=function(n, t=1){
  x=runif(n)
  x3=F.inv(x)     #on crée un vecteur de simulation de X3
  echantillon=1-pgamma(t-x3,2,1)
  return(echantillon)
}

#estimateur
h2=estiG(n)
f2=mean(h2)
points(cumsum(h2)/(1:length(h2)),type='l' , col="red")
legend("bottomright", legend=c("estimateur MC classique", "estimateur MC avec h1", "estimateur MC avec h2"), col=c('chartreuse', 'blue', 'red'), lty=c(1), lwd=c(3), bg='aliceblue')

#intervalle de confiance
IC2=c(f2-qnorm(0.975)*sqrt(var(h2)/n), f2+qnorm(0.975)*sqrt(var(h2)/n))

#performance: calcul? Ou juste interprétation des résultats obtenus? peut-on prendre à partir de n=4000?
plot(cumsum(h)/(1:length(h)),type='l', col='chartreuse', xlab="n", ylab="y", main="zoom convergence des estimateurs", ylim=c(0.970,.975))
points(cumsum(h1)/(1:length(h1)),type='l', col='blue')
points(cumsum(h2)/(1:length(h2)),type='l' , col="red")

################################################################################
#x = rnorm(1e4)
#plot(cumsum(x)/(1:1e4),type='l')
#plot(Vectorize(MC.esti)(1:1000)[1,],type='l')
##################################################################################



#PARTIE 2, question 1
n=10000

#pour obtenir ce qu'on veut, on prend t=1.5 sachant que 5 est trop grand
t=1.5


#on définit la fonction longueur
d<-function(x){
  return(min(x[1]+x[4], x[1]+x[3]+x[5], x[2]+x[5], x[2]+x[3]+x[4]))
}

#on définit la fonction pour l'estimateur de montecarlo
indicatrice_2=function(x,t=1.5){
  ifelse(d(x)>=t,1,0)
}

#fonction qui crée un 
MC.sample=function(n,t=1.5){
  x1=rexp(n,6)  #verifier qu'il faut bien donne 1/lambda en arg sur R
  x2=rexp(n,7) 
  x3=rexp(n,3) 
  x4=rexp(n,2) 
  x5=rexp(n,1) 
  P=matrix(c(x1,x2,x3,x4,x5), ncol=5, byrow=FALSE)
  sample=c()
  for (i in 1:n){
    sample=c(sample,indicatrice_2(P[i,],t)) #reflechir avec colSums (sous entend redéfinir d(x) alors)
  }
  return(sample)
}

#estimateur et convergence
b=MC.sample(n)
MC.esti2=mean(b)
plot(cumsum(b)/(1:length(b)), type='l', xlab='n', ylab='estimateur', main='convergence estimateur de monte carlo' )







#Echantillonage préférentiel, question 2b)
#on simule des lois suivant la densité g par rejet
#on crée la densité f(x,lambda)
fexp<-function(x){
  return(dexp(x[1],6)*dexp(x[2],7)*dexp(x[3],3)*dexp(x[4],2)*dexp(x[5],1))
}

g2<-function(x,t){
  return(indicatrice_2(x,t)*fexp(x))
}


###########    manque d'optimalité totale , meilleure beta? faire méthode du prof en 5dimensions? #######################
rejetg=function(n,t=1.5){
  solution=matrix(numeric(5*n), ncol=5)
  for (i in 1:n){
    u=runif(1)
    exp=c(rexp(1,6), rexp(1,7), rexp(1,3), rexp(1,2), rexp(1,1))
    alpha=indicatrice_2(exp,t)  #expliquer pourquoi celui la
    while (u>alpha){
      u=runif(1)
      exp=c(rexp(1,6), rexp(1,7), rexp(1,3), rexp(1,2), rexp(1,1))
      alpha=indicatrice_2(exp,t)
    }
    solution[i,]=exp
  }
  return(solution)
}

#comencons par calculer les ai

a0=c(1,2,3,4,5)

simuler_h=function(n){
  solution=matrix(numeric(5*n), ncol=5)
  for (i in 1:n){
    exp=c(rexp(1,a[1]), rexp(1,a0[2]), rexp(1,a0[3]), rexp(1,a0[4]), rexp(1,a0[5]))
    solution[i,]=exp
  }
  return(solution)
}
solution=simuler_h(n)

a_opti=c()

for (j in 1:5){
  sum1=0
  sum2=0
  for (i in 1:n){
    sum1=sum1+(indicatrice_2(solution[i,],t)*(6/a0[1])*exp(-(6-a0[1])*solution[i,1])*(7/a0[2])*exp(-(7-a0[2])*solution[i,2])*(3/a0[3])*exp(-(3-a0[3])*solution[i,3])*(2/a0[4])*exp(-(2-a0[4])*solution[i,4])*(1/a0[5])*exp(-(1-a0[5])*solution[i,5]))
    sum2=sum2+(indicatrice_2(solution[i,],t)*solution[i,j]*(6/a0[1])*exp(-(6-a0[1])*solution[i,1])*(7/a0[2])*exp(-(7-a0[2])*solution[i,2])*(3/a0[3])*exp(-(3-a0[3])*solution[i,3])*(2/a0[4])*exp(-(2-a0[4])*solution[i,4])*(1/a0[5])*exp(-(1-a0[5])*solution[i,5]))
  }
  a_opti=c(a_opti,sum1/sum2)
}


h_opti<-function(x){
  return(dexp(x[1],a_opti[1])*dexp(x[2],a_opti[2])*dexp(x[3],a_opti[3])*dexp(x[4],a_opti[4])*dexp(x[5],a_opti[5]))
}

#echantillonage preferentiel opti

#rejet pour simuler selon h_opti


simuler_h_opti=function(n){
  solution=matrix(numeric(5*n), ncol=5)
  for (i in 1:n){
    exp=c(rexp(1,a_opti[1]), rexp(1,a_opti[2]), rexp(1,a_opti[3]), rexp(1,a_opti[4]), rexp(1,a_opti[5]))
    solution[i,]=exp
  }
  return(solution)
}

Y=simuler_h_opti(n)

le_truc=c()

for (i in 1:n){
  le_truc=c(le_truc,(fexp(Y[i,])/h_opti(Y[i,])*indicatrice_2(Y[i,],t)))
}

plot(cumsum(le_truc)/(1:n), type='l')

pn_opti=mean(le_truc)

IC_borne_sup=pn_opti+qnorm(0.975)*var(h_opti)/n
IC_borne_inf=pn_opti-qnorm(0.975)*var(h_opti)/n
IC=c(IC_borne_inf,IC_borne_sup)



###############ex8
###FAUX####
rho<-function(x){
  return(min((-1/a_opti[1])*log(x[1])-(1/a_opti[4])*log(x[4]), (-1/a_opti[1])*log(x[1])-(1/a_opti[3])*log(x[3])-(1/a_opti[5])*log(x[5]), (-1/a_opti[2])*log(x[2])-(1/a_opti[5])*log(x[5]), (-1/a_opti[2])*log(x[2])-(1/a_opti[3])*log(x[3])-(1/a_opti[4])*log(x[4]))) 
}


indicatrice_3<-function(x,t=1.5){
  ifelse(rho(x)>=t, 1, 0)
}

fct_8<-function(u){
  return((6*7*3*2*1)/(a_opti[1]*a_opti[2]*a_opti[3]*a_opti[4]*a_opti[5])*u[1]^((6/a_opti[1])-1)*u[2]^((7/a_opti[2])-1)*u[3]^((3/a_opti[3])-1)*u[4]^((2/a_opti[4])-1)*u[5]^((1/a_opti[5])-1))
}

simuler_uniforme=function(n){
  solution=matrix(numeric(5*n), ncol=5)
  for (i in 1:n){
    unif=c(runif(5,0,1))
    solution[i,]=unif
  }
  return(solution)
}


question8=simuler_uniforme(n)
question8_1=1-question8

v8=c()
v8_1=c()

for (i in 1:n){
  v8=c(v8, fct_8(question8[i,])*indicatrice_3(question8[i,]))
  v8_1=c(v8_1, fct_8(question8_1[i,])*indicatrice_3(question8_1[i,]))
}

pn_unif=mean(v8+v8_1)/2














