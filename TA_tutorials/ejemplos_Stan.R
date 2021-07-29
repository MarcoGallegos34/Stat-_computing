library(rstan)
setwd("C:/Users/marco/Documents")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### bernoulli.stan ###

set.seed(1234)
y <- c()
n <- 100
for(i in 1:n){
  candidate <- runif(1)
  if(candidate < .6){
    y <- c(y,1)
  } else{
    y <- c(y,0)
  }
}



bernoulli_dat <- list(y=y,N=-10)

fit1 <- stan(file = 'bernoulli.stan', data = bernoulli_dat)

print(fit1)


la <- extract(fit1,permuted = TRUE) # return a list of arrays 

alpha = 1
beta = 1

alpha_post = alpha + sum(y)
beta_post = beta + n - sum(y)

final_post = curve(dbeta(x,alpha_post,beta_post),.5,1)
hist(la$theta,freq = F,breaks=40)
lines(final_post)

### Schools.stan ###

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

separate_estimates <- function(){
  
  int_A = 28 + 2*c(-15,15)
  int_B = 8 + 2*c(-10,10)
  int_C = -3 + 2*c(-16,16)
  int_D = 7 + 2*c(-11,11)
  int_E = -1 + 2*c(-9,9)
  int_F = 1 + 2*c(-11,11)
  int_G = 18 + 2*c(-10,10)
  int_H = 12 + 2*c(-18,18)
  plot(int_A,c(10,10),type="l",xlim=c(-37,60),ylim=c(2,10),ylab="Escuelas",xlab="95% CI")
  lines(28,10,type="o",col="red",pch=16)
  lines(int_B,c(9,9))
  lines(8,9,type="o",col="red",pch=16)
  lines(int_C,c(8,8))
  lines(-3,8,type="o",col="red",pch=16)
  lines(int_D,c(7,7))
  lines(7,7,type="o",col="red",pch=16)
  lines(int_E,c(6,6))
  lines(-1,6,type="o",col="red",pch=16)
  lines(int_F,c(5,5))
  lines(1,5,type="o",col="red",pch=16)
  lines(int_G,c(4,4))
  lines(18,4,type="o",col="red",pch=16)
  lines(int_H,c(3,3))
  lines(12,3,type="o",col="red",pch=16)
}

separate_estimates()

fit2 <- stan(file = 'schools.stan', data = schools_dat)


print(fit2)
plot(fit2)
pairs(fit2, pars = c("mu", "tau", "lp__"))

la <- extract(fit2) # return a list of arrays 
mu <- la$mu

hist(mu)
### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit2, permuted = FALSE) 


### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)

dim(a2)
dim(m)


#Para utilizar más modelos
model <- stan_demo()
