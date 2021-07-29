#Cópula Gumbel con parámetro theta > 1

h = function(u,v,theta = 1) ((-log(u))^theta + (-log(v))^theta )^(1/theta) 
C = function(u,v, theta = 1) exp(-h(u,v,theta))

c_ = function(u,v,theta=1) ((log(u)*log(v))^(theta-1))/(u*v)*C(u,v,theta)*(h(u,v,theta)^(1-2*theta))*(h(u,v,theta) + (theta-1))

#Haremos M-H para las distribuciones marginales y LogNormal(0,.25) y exponencial truncada por la derecha Exp(1,a), a = 10
a = 10
log_c_ = function(u,v,theta=1) (theta-1)*log(log(u)*log(v)) - log(u*v) + log(C(u,v,theta)) + (1-(2*theta))*log(h(u,v,theta)) + log(h(u,v,theta) + (theta-1))

f = function(x1,x2,theta=1) c_(plnorm(x1,0,.25),pexp(x2-10,1),theta)*dlnorm(x1,0,.25)*dexp(x2-10,1)
log_f = function(x1,x2,theta=1) log_c_(plnorm(x1,0,.25),pexp(x2-10,1),theta) + dlnorm(x1,0,.25,log = T) + dexp(x2-10,1,log = T)

#Necesitamos esta biblioteca para extraer funciones de la Normal Multivariada
library(mvtnorm)

M_H = function(Nsim=10000,theta0=1,seed=1234){
  
  set.seed(seed)
  
  g1 = function(x1) exp(x1)
  g1_inv = function(x1) log(x1)
  
  g2 = function(x2) exp(x2)+a
  g2_inv = function(x2) log(x2-a) 
  
  X = matrix(NA, Nsim+1,2)
  X[1,] = c(rlnorm(1,0,.25),rexp(1)+10)
  counter = c()
  
  for(i in 1:Nsim){
    U <- runif(1)
    Z_n <- rmvnorm(1,c(g1_inv(X[i,1]),g2_inv(X[i,2])))
    Y_n <- c(g1(Z_n[1]),g2(Z_n[2]))
    #print(Z_n)
    print("Y_n")
    print(Y_n)
    #print("f(Y_n)")
    #print(f(Y_n[1],Y_n[2],theta = theta0))
    #print("f(X_n)")
    #print(f(X[i,1],X[i,2]))
    print("numerador")
    print(log_f(Y_n[1],Y_n[2],theta =theta0))
    #print(f(Y_n[1],Y_n[2],theta =theta0)*Y_n[1]*(Y_n[2]-a))
    #print("denominador")
    #print(f(X[i,1],X[i,2],theta0)*X[i,1]*(X[i,2]-a))
    #ratio = min(0,
    #            log(f(Y_n[1],Y_n[2],theta =theta0)*Y_n[1]*(Y_n[2]-a))-log(f(X[i,1],X[i,2],theta0)*X[i,1]*(X[i,2]-a)))
    ratio = min(0,
                log_f(Y_n[1],Y_n[2],theta =theta0)+log(Y_n[1]*(Y_n[2]-a))-log_f(X[i,1],X[i,2],theta0)-log(X[i,1]*(X[i,2]-a)))
    if(log(U) <= ratio){
      X[i+1,] = Y_n
      counter <- c(counter,1)
    } else {
      X[i+1,] = X[i,]
      counter <- c(counter,0)
    }
  }
  
  X = X[-1,]
  
  return(list(X=X,counter=counter))
}

samples = M_H()

X = samples$X

plot(X[,1],type="l")
final_logNormal = curve(dlnorm(x,0,.25),0,3)
hist(X[,1],freq = F)
lines(final_logNormal)


plot(X[,2],type="l")
final_truncExp = curve(dexp(x-10,1),9,15)
hist(X[,2],freq=F,breaks=40)
lines(final_truncExp)

samples = M_H(theta0 = 10)
X = samples$X
