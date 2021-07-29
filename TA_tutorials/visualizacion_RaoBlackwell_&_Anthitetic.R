#Ejemplo estimador antitético

exponent <- 2
no_samples <- 1000
cmc_samples <- runif(no_samples)
antic_samples <- runif(no_samples/2)
antic_samples <- as.vector(rbind(antic_samples,1-antic_samples))

plot(cumsum(cmc_samples^exponent)/1:no_samples,type="l")
abline(h=1/(exponent+1),col="red",lty="dashed")

lines(cumsum(antic_samples^exponent)/1:no_samples,type="l",col="blue")

#Muestras para estimar E[X], X ~ t-student con 2 grados de libertad
muestras_t <- rt(no_samples,2)

#Convergencia ergódica de Monte Carlo Crudo
plot(cumsum(muestras_t)/1:no_samples,type="l",col="blue")

#Convergencia ergódica del estimador Rao-blackwellizado
abline(h=0,col="red")
