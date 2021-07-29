# Target distribution : Gamma(4.85,1)

#Usando el método Aceptación-Rechazo
n=5000
set.seed(10)
X_AR= c()
counter_AR = c()
alpha = 4.85
beta = 1
a = floor(alpha)
b= floor(alpha)/alpha
max = (alpha - a)/(1-b)
M=b^(-a)*max^(alpha-a)*exp(-(1-b)*max)

Y = c()

for(i in 1:n){
  U <- runif(1)
  Y_t <- rgamma(1,a,b)
  Y = c(Y,Y_t)
  ratio = dgamma(Y_t,alpha,beta)/(M*dgamma(Y_t,a,b))
  if(U <= ratio){
    X_AR = c(X_AR,Y_t)
    counter_AR <- c(counter_AR,1)
  } else {
    #X[i+1] = X[i]
    counter_AR <- c(counter_AR,0)
  }
}

mean(counter_AR)
mean(X_AR)
x = seq(0,15,length.out = 1000)
y = dgamma(x,4.85,1)
hist(X_AR,freq=F,breaks=40)
lines(x,y)


set.seed(10)
X_MH = c(rgamma(1,a,b))
counter_MH = c()
for(i in 1:n){
  U <- runif(1)
  Y_n <- Y[i]
  ratio = min(1,
              dgamma(Y_n,alpha,beta)*dgamma(X_MH[i],a,b)/
                (dgamma(X_MH[i],alpha,beta)*dgamma(Y_n,a,b)))
  if(U <= ratio){
    X_MH[i+1] = Y_n
    counter_MH <- c(counter_MH,1)
  } else {
    X_MH[i+1] = X_MH[i]
    counter_MH <- c(counter_MH,0)
  }
}

#Quitamos el primer término para tener exactamente n muestras
X_MH=X_MH[-1]

mean(counter_MH)
mean(X_MH)
x = seq(0,15,length.out = 1000)
y = dgamma(x,4.85,1)
hist(X_MH,freq=F,breaks=40)
lines(x,y)

plot(cumsum(X_AR)/1:sum(counter_AR),type = "l")
lines(cumsum(X_MH)/1:n,type = "l",col="blue")
abline(h=4.85,lty="dashed",col="red")

#Con esto podemos observar que se obtienen muestras igual de buenas
#que el Aceptación-Rechazo y además es mucho más eficiente en la
#generación de muestras
acf(X_AR,ylim=c(0,1),xlim=c(0,24))
acf(X_MH,ylim=c(0,1),xlim=c(0,24))



#############################################################
### Comparación de candidatos dentro de M-H independiente ###
#############################################################

Cauchy_samples_MH <- function(candidate = "normal",Nsim = 10^5){
  
  set.seed(15)
  X = c(rt(1,1))
  
  if(candidate == "normal"){
    
    for (t in 2:Nsim){
      Y=rnorm(1) # candidate normal
      rho=dt(Y,1)*dnorm(X[t-1])/(dt(X[t-1],1)*dnorm(Y))
      X[t]=X[t-1] + (Y-X[t-1])*(runif(1)<rho)
    }
    
    return(X)
    
  } else if( candidate == "student"){
    
    for (t in 2:Nsim){
      Y=rt(1,df= .5) # candidate student
      rho=dt(Y,df=1)*dt(X[t-1],df=.5)/(dt(X[t-1],df=1)*dt(Y,df=.5))
      X[t]=X[t-1] + (Y-X[t-1])*(runif(1)<rho)
    }
    
    return(X)
    
  } else {
    warning("El parámetro candidate sólo admite los valores \"normal\" o \"student\". 
      Elija una de estos por favor.")
  }
}

Nsim=10^5
X_normal  = Cauchy_samples_MH("normal")
length(unique(X_normal))/Nsim
plot(X_normal,type="l")

X_student = Cauchy_samples_MH("student")
length(unique(X_student))/Nsim
plot(X_student,type="l")

#Histograma de las muestras obtenidas para cada distribución candidata
Cauchy_target = curve(dt(x,1),-10,10)
hist(X_normal,breaks = 40,freq = F,xlim = c(-10,10),ylim=c(0,.4))
hist(X_student,breaks = 30^4,freq = F,xlim = c(-10,10),ylim=c(0,.35))
lines(Cauchy_target,col="blue")

#lag de las muestras correspondientes a cada distribución candidata
acf(X_normal)
acf(X_student)

#Convergencia ergódica de la estimación P(X < 3), X ~ Cauchy(0,1)
plot(cumsum(X_student<3)/1:Nsim,type="l",ylim=c(.85,1),col="red")
lines(cumsum(X_normal<3)/1:Nsim,type="l")
abline(h=pt(3,1),lty="dashed", col="blue")



########################
### Acceptance rates ###
########################

# Target distribution: N(0,1)

#Esta librería nos permitirá obtener muestras y la distribución
# de la distribución doble exponencial (distribución de Laplace)
library(nimble)
curve(ddexp(x,),-10,10)
MH_target_Normal = function(alpha,n=1000,seed=34){
  
  #alpha=1
  set.seed(seed)
  X = c(rdexp(1,rate=alpha))
  counter = c()
  
  for(i in 1:n){
    U <- runif(1)
    Y_n <- rdexp(1,rate=alpha)
    ratio = min(1,
                dnorm(Y_n)*ddexp(X[i],rate=alpha)/
                  (dnorm(X[i])*ddexp(Y_n,rate=alpha)))
    if(U <= ratio){
      X[i+1] = Y_n
      counter <- c(counter,1)
    } else {
      X[i+1] = X[i]
      counter <- c(counter,0)
    }
  }
  X = X[-1]
  return(list(X=X,counter=counter))
}

L1 = MH_target_Normal(1,seed=150)
mean(L1$counter)
hist(L1$X,freq=F,breaks=40,xlim=c(-4,4))
mean(L1$X)
length(L1$X)

L3 = MH_target_Normal(3,seed=150)
mean(L3$counter)
hist(L3$X,freq=F,breaks=40,xlim=c(-4,4))
mean(L3$X)

x = seq(-5,5,length.out = 1000)
y = dnorm(x)
lines(x,y)

#Convergencia ergódica de las muestras
plot(cumsum(L1$X)/1:1000,type="l",col="green")
lines(cumsum(L3$X)/1:1000,type="l")
abline(h=0,col="red",lty="dashed")

#Función de autocorrelación (autocovarianza)
acf(L1$X)
acf(L3$X)

plot(L1$X,type="l")
plot(L3$X,type="l")
