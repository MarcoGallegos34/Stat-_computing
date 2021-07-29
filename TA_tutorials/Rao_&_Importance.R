###############################################
#                                             #
# Continuación : Estimador Rao-blackwellizado #
#                                             #
###############################################

# Parámetros y valores iniciales

n <- 10000
set.seed(14)
#Para alpha = beta = 1

alpha <- 1
beta <- alpha
v <- 2*alpha


#Generemos las muestras para cada estimador de E[X^2] = Infinito

#X ~ t(v), v= 2
student <- rt(n,v)

# 1/tau ~ InvGamma(1,1)
invGa <- invgamma::rinvgamma(n,alpha,beta)

# tau ~ Gamma(1,1) (forma alterna de generarlo)
tau <- rgamma(n,alpha,beta)

# Estimadores para E[X^2] = Infinito

plot(cumsum(student^2)/1:n,type="l")
lines(cumsum(invGa)/1:n,type="l",col="blue")
lines(cumsum(1/tau)/1:n,type="l",col="red")


# Para alpha = beta = 3
alpha <- 3
beta <- alpha
v <- 2*alpha

#X ~ t(v), v= 6
student <- rt(n,v)

# 1/tau ~ InvGamma(3,3)
invGa <- invgamma::rinvgamma(n,alpha,beta)

# tau ~ Gamma(3,3) (forma alterna de generarlo)
tau <- rgamma(n,alpha,beta)

#Calculemos las varianzas de cada estimador (esto aplica cuando alpha = beta > 2)
k <- 4

variance_X2 <- (gamma((k+1)/2)*gamma((v-k)/2)*v^(k/2))/(sqrt(pi)*gamma(v/2)) - (v/(v-2))^2
variance_Rao <- beta^2/((alpha-1)^2*(alpha-2))

#Proporción entre varianzas
variance_X2/variance_Rao


plot(cumsum(student^2)/1:n,type="l")
lines(cumsum(invGa)/1:n,type="l",col="blue")
lines(cumsum(1/tau)/1:n,type="l",col="red")
abline(h=3/2,lty="dashed",col="brown")

mean(student^2)
mean(invGa)
mean(1/tau)


##################################
#                                #
#    Muestreo por importancia    #
#                                #
##################################

pnorm(-4.5)
#pnorm(-4.5,log=T)


vec <- c()
for(i in 1:200){
  vec <- c(vec,sum(rnorm(1000000) < -4.5))
}

mean(vec)

Nsim=10^4
y=rexp(Nsim)+4.5
weit=dnorm(y)/dexp(y-4.5)
plot(cumsum(weit)/1:Nsim,type="l")
abline(a=pnorm(-4.5),b=0,col="red")
mean(weit)

sum(vec)/(1000000*200)
