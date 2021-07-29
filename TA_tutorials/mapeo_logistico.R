#Creación mapeo logístico
mapeo_log <- function(x){
  return(4*x*(1-x))
}

#Aplicación del mapeo logístico de forma iterativa
iter_map_log <- function(x,k){
  iteraciones <- c(x)
  for(i in 1:k){
    iteraciones <- c(iteraciones,mapeo_log(iteraciones[i]))
  }
  return(iteraciones)
}

set.seed(1234)
x_0 <- runif(1)

muestras <- iter_map_log(x_0,1000)
hist(muestras,freq = F)
dist_beta <- curve(expr = dbeta(x,.5,.5) ,from=0,to=1)
lines(dist_beta)

paraUnif <- function(x){
  return((2/pi)*asin(sqrt(x)))
}

hist(paraUnif(muestras))

plot(muestras)

plot(paraUnif(muestras))
