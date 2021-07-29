#Función a integrar
h = function(x) return((exp(-(x-3)^2/2) + exp(-(x-6)^2/2)))

#Aquí observamos la función que queremos integrar
curve(h(x)*dnorm(x),-8,5)


n=10000

#El valor de E[h(X)], X ~ N(0,1)
Expected=(exp(-9/4)+exp(-9))/sqrt(2)

set.seed(1)
#Generamos n muestras de una normal estándar
normales = rnorm(n)

#Generamos n muestras para nuestro estimador antitético
antic <- rnorm(n/2)
antic_samples <- as.vector(rbind(antic,-antic))

#Correlación estimada de (h(X),h(-X))
cor(h(antic),h(-antic))

#Aquí podemos ver la convergencia ergódica
plot(cumsum(h(normales))/1:n,type="l")
lines(cumsum(h(antic_samples))/1:n,col="blue")
abline(h=Expected,lty="dashed",col="red")

#Error absoluto de las estimaciones
abs(Expected - mean(h(normales)))
abs(Expected - mean(h(antic_samples)))


########################################################
### Calculo de la integral de 0 a 1 de sqrt(-log(x)) ###
########################################################

raiz_log = function(x) return(sqrt(-log(x)))
n=5000

#Aquí tenemos a la función a integrar
curve(raiz_log(x),from=0,to=1)

set.seed(2)
#CMC
U <- runif(n)

#Varianza aproximada de este estimador
var(raiz_log(U))

#Antithec samples
antic <- runif(n/2)
antic_samples <- as.vector(rbind(antic,1- antic))


#Importance sampling

#Aquí simulamos Z ~ Beta(1/2,1/2)
alpha = 1/2
beta =1/2
set.seed(3)
importBeta <- rbeta(n,alpha,beta)

#Varianza estimada de este estimador
var(raiz_log(importBeta)/dbeta(importBeta,alpha,beta))


#Valor estimado de E[sqrt(-log(U))]
Expected=sqrt(pi)/2

#Convergencia ergódica de cada uno de los estimadores
plot(cumsum(raiz_log(U))/1:n,type="l")
lines(cumsum(raiz_log(antic_samples))/1:n,col="blue")
lines(cumsum(raiz_log(importBeta)/dbeta(importBeta,alpha,beta))/1:n,col="green")
abline(h=Expected,col="red",lty="dashed")


#Errores absolutos de nuestras estimaciones después de n muestras
abs(Expected - mean(raiz_log(U)))
abs(Expected - mean(raiz_log(antic_samples)))
abs(Expected - mean(raiz_log(importBeta)/dbeta(importBeta,alpha,beta)))


###Variables de control###

#Para calcular E[exp(U)], g(x)=exp(x)
m <- 1000
a <- - 12 + 6 * (exp(1) - 1)
U <- runif(m)
T1 <- exp(U) #simple  MC
T2 <- exp(U) + a * (U - 1/2) #controlled

plot(cumsum(T1)/1:m,type="l")
lines(cumsum(T2)/1:m,type="l",col="blue")
abline(h=exp(1)-1,lty="dashed",col="red")

#Aquí podemos observar el porcentaje de reducción de varianza con
#respecto a E[g(X)]
(var(T1) - var(T2)) / var(T1)

#Error absoluto de los estimadores
abs(exp(1)-1 - mean(T1))
abs(exp(1)-1 - mean(T2))

#Para calcular la integral de 0 a 1 de g(x)=exp(-x)/(1 + x^2)

f <- function(u) exp(-.5)/(1+u^2)
g <- function(u) exp(-u)/(1+u^2)
set.seed (510)
m <- 1000
u <- runif(m)
B <- f(u)
A <- g(u)
cor(A,B)
a <- -cov(A,B)/var(B)

T1 <- g(u)
T2 <- T1 + a*(f(u) - exp(-.5)*pi/4)

#Valor esperado de E[g(X)]
Expected = integrate(g,lower=0,upper=1)$value

plot(cumsum(T1)/1:m,type="l")
lines(cumsum(T2)/1:m,type="l",col="blue")
abline(h=Expected,col="red",lty="dashed")

#Error absoluto de los estimadores
abs(Expected - mean(T1))
abs(Expected - mean(T2))
