
### Ejercicio 2 ###
### SOLUCIÓN CON CÓPULA DE CLAYTON ###

C = function(u,v, theta = 1) (u^(-theta) + v^(-theta) - 1)^(-1/theta)

c_ = function(u,v,theta=1) (theta + 1)*(u*v)^(-(theta + 1))*(u^(-theta) + v^(-theta) - 1)^(-(2*theta + 1)/theta)

log_c_ = function(u,v,theta=1) log(theta + 1) -(theta + 1)*(log(u) + log(v)) -(2*theta + 1)/theta*log(u^(-theta) + v^(-theta)-1)

f = function(x1,x2,theta=1) c_(pbeta(x1,1/2,1/2),pgamma(x2,20,4),theta)*dbeta(x1,1/2,1/2)*dgamma(x2,20,4)
log_f = function(x1,x2,theta=1) log_c_(pbeta(x1,1/2,1/2),pgamma(x2,20,4),theta) + dbeta(x1,1/2,1/2,log = T) + dgamma(x2,20,4,log = T)

#log_ratio = function(X,Y,alpha=6,beta=4,a=20,b=4,theta=1) -(theta + 1)*(log(pbeta(Y[1],alpha,beta)/pbeta(X[1],alpha,beta)) + pgamma(Y[2],a,b,log=T) - pgamma(X[2],a,b,log=T)) -(2*theta + 1)/theta*(log(pbeta(Y[1],alpha,beta)^(-theta) + pgamma(Y[2],a,b)^(-theta) - 1 ) - log(pbeta(X[1],alpha,beta) + pgamma(X[2],a,b) - 1 )) + alpha*(log(Y[1]) - log(X[1])) + beta*(log(1-Y[1]) - log(1-X[1])) + a*(log(Y[2]) - log(X[2])) -beta*(Y[2]-X[2])

library(mvtnorm)

M_H = function(Nsim=1000,theta0=1,seed=1234){
  
  covariance_mat = matrix(c(1,0,0,1),2,2)
  set.seed(seed)
  g1 = function(x1) exp(x1)/(exp(x1)+1)
  g1_inv = function(x1) log(x1/(1-x1))
  
  g2 = function(x2) exp(x2)
  g2_inv = function(x2) log(x2)
  
  X = matrix(NA, Nsim+1,2)
  X[1,] = c(rbeta(1,1/2,1/2),rgamma(1,20,4))
  counter = c()
  
  #waste = 0
  for(i in 1:Nsim){
    U <- runif(1)
    Z_n <- rmvnorm(1,mean = c(g1_inv(X[i,1]),g2_inv(X[i,2])),sigma = covariance_mat)
    Y_n <- c(g1(Z_n[1]),g2(Z_n[2]))
    
    ratio = min(0,
                (log_f(Y_n[1],Y_n[2],theta0) + log(Y_n[2]*Y_n[1]*(1-Y_n[1]))) - (log_f(X[i,1],X[i,2],theta0) + log(X[i,2]*X[i,1]*(1-X[i,1]))))
    #ratio = min(0,log_ratio(X=X[i,],Y=Y_n,theta=theta0))
    #print(log_ratio(X=X[i,],Y=Y_n,theta=theta0))
    if(is.nan(ratio)){
      print("Partes de log_ratio")
      print("Parte 1")
      print(-(theta0 + 1)*(log(pbeta(Y_n[1],6,4)*pgamma(Y_n[2],20,4))))
      print("Parte 2")
      print((theta0 + 1)*(log(pbeta(X[i,1],6,4)) + log(pgamma(X[i,2],20,4))))
      print("Parte 3")
      print(6*(log(Y_n[1]) - log(X[i,1])) + 4*(log(1-Y_n[1]) - log(1-X[i,1])) + 20*(log(Y_n[2]) - log(X[i,2])) -4*(Y_n[2]-X[i,2]))
      print("numerador log_f")
      print(log_f(Y_n[1],Y_n[2],theta0))
      print("numerador log derivadas inversas")
      print(log(Y_n[2]*Y_n[1]*(1-Y_n[1])))
      print("denominador")
      print(log_f(X[i,1],X[i,2],theta0))
      print("denominador log derivadas inversas")
      print(log(X[i,2]*X[i,1]*(1-X[i,1])))
      print("Y_n")
      print(Y_n)
      print(paste("X_n al tiempo i=",i))
      print(X[i,])
    }
    #print(ratio)
    if(log(U) <= ratio){
      X[i+1,] = Y_n
      counter <- c(counter,1)
    } else {
      X[i+1,] = X[i,]
      counter <- c(counter,0)
    }
  }
  X = X[-1,]
  #print(waste)
  return(list(X=X,counter=counter))
  
}

test1 = M_H(Nsim=2000,theta0=,seed=423)

Nsim=10000
X = M_H(Nsim,theta0=50)$X

final_beta = curve(dbeta(x,1/2,1/2),0,1)
plot(X[,1],type="l")
hist(X[,1],freq = F)
#hist(X[,1],breaks = 80,freq = F)
lines(final_beta)

final_gamma = curve(dgamma(x,20,4),0,10)
plot(X[,2],type="l")
hist(X[,2],freq = F)
#hist(X[,2],breaks = 60,freq = F)
lines(final_gamma)
plot(X)




Copula = matrix(NA,Nsim,2)
for(i in 1:Nsim){
  Copula[i,] =  c(pbeta(X[i,1],.5,.5),pgamma(X[i,2],20,4))
}

plot(Copula)

### Ejercicio 3 ###

### Generación de datos X para el muestreador de Gibbs e implementación de dicho muestreador ###

library(gtools)

set.seed(45)
pi_ = rdirichlet(1,c(2,2)) #(0.719,0.281)
set.seed(1)
#lambda_ = rgamma(2,4.5,1.5)
lambda_ = rgamma(2,4.5,1)
#curve(dgamma(x,4.5,1),0,10)
X = c()
set.seed(1234)
Nsam = 200
for(i in 1:Nsam){
  U <- runif(1)
  if(U <= pi_[1]){
    X = c(X,rpois(1,lambda_[1]))
    
  } else{
    X = c(X,rpois(1,lambda_[2]))
  }
}


cond_gamma_j = function(alpha=1,beta=1,Y) rgamma(1,alpha + sum(Y), beta + length(Y))

cond_dirichlet = function(vec_alpha=c(1,1),C) rdirichlet(1,vec_alpha + C )

cond_Z_j = function(x,Pi,lambda)  c(1:d)[rmultinom(1,1,Pi*lambda^x*exp(-lambda)/sum(Pi*lambda^x*exp(-lambda))) == 1]

Nsim=2000
d=2

lambda = matrix(NA,Nsim+1,d)
Z = matrix(NA,Nsim+1,Nsam)
Pi = matrix(NA,Nsim+1,d)

lambda[1,] = rgamma(2,1,1)
Z[1,] = rbinom(Nsam,1,rdirichlet(1,c(1,1)))+1
Pi[1,] = rdirichlet(1,c(1,1))

for(i in 1:Nsim){
  for(j in 1:d){
    Y = X[Z[i,] == j]
    lambda[i+1,j] = cond_gamma_j(Y=Y)
  }
  
  Pi[i+1,] = cond_dirichlet(C=as.vector(table(Z[i,])))
  
  for(j in 1:Nsam){
    Z[i+1,j] = cond_Z_j(x=X[j],Pi = Pi[i+1,],lambda = lambda[i+1,])
  }
}

lambda = lambda[-1,]
Pi = Pi[-1,]

lambda[Nsim,]

apply(lambda,2,mean)
apply(Pi,2,mean)

lambda_
pi_

plot(lambda[,1],type = "l")
hist(lambda[,1])

plot(cumsum(lambda[,1])/1:Nsim,type = "l",col="red")
abline(h=lambda_[1],lty="dashed")
plot(cumsum(lambda[,2])/1:Nsim,type = "l",col="blue")
abline(h=lambda_[2],lty="dashed")


#library(tidyverse)
#write_csv(as.data.frame(X),"datos_Gibbs.txt",col_names = F)
#X_test = read_csv("datos_Gibbs.txt",col_names = F)
