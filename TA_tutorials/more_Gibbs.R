x=c(91,504,557,609,693,727,764,803,857,929,970,1043,1089,1195,1384,1713)

set.seed(4321)
xbar=mean(x)
n = length(x)
Nsim = 5000
a=b=3
sh1=(n/2)+a
tau2 =5
theta0 = 10

sigma2=theta=rep(0,Nsim) #init arrays

sigma2[1]=1/rgamma(1,shape=a,rate=b) #init chains

B=sigma2[1]/(sigma2[1] + n*tau2)

theta[1]=rnorm(1,B*theta0+(1-B)*xbar, sd=sqrt(tau2*B))

for (i in 2:Nsim){
  
  B=sigma2[i-1]/(sigma2[i-1]+n*tau2) # 
  
  theta[i]=rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B)) # Aquí tenemos X_{t+1} | Y_{t}
  ra1=(1/2)*(sum((x-theta[i])^2))+b                       # Aquí tenemos Y_{t+1} | X_{t+1}
  sigma2[i]=1/rgamma(1,shape=sh1,rate=ra1)
}

#sigma2 = sigma2[-1]

sigma = sqrt(sigma2)

hist(sqrt(sigma2),breaks = 40)
hist(log(sigma),breaks = 80)
hist(log(sigma2)/2,breaks = 80)
quantile(log(sigma2)/2,probs=c(.05,.95))
mean(log(sigma))

#theta = theta[-1]

quantile(log(theta),probs=c(.05,.95))
hist(theta, breaks =1000, xlim=c(0,20))
hist(log(theta) , breaks = 80)


#####################################################
### Ahora, hagamos muestreo e inferencia con Stan ###
#####################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

data <- list(N=n, x=x, a=a, b=b, tau2=tau2, theta0 = theta0)

fit1 <- stan(file ="simpleGibbs.stan",data = data)

print(fit1)
plot(fit1)
pairs(fit1, pars = c("theta", "sigma2", "lp__"))

la <- rstan::extract(fit1, permuted = TRUE)

hist(la$theta)
hist(log(la$theta),breaks=40)

mean(log(la$sigma2)/2)
hist(log(la$sigma2)/2,breaks=40)

quantile(log(la$theta)/2,probs=c(.05,.95))
quantile(log(la$sigma2)/2,probs=c(.05,.95))



library(mcsm)
data("Energy")
Gibss_hierarchical = function(Nsim = 5000,
                              data = Energy,
                              seed = 1234,
                              mu_0 = 1000,
                              a_1 = 2,
                              a_2 = 2,
                              a_3 = 2,
                              b_1 = 5,
                              b_2 = 5,
                              b_3 = 5,
                              k = 2){
  
  set.seed(seed)
  
  n = dim(data)[1]
  
  sigma2_mu = c(1/rgamma(1,a_3,b_3))
  mu = c(rnorm(1,mu_0,sqrt(sigma2_mu[1])))
  sigma2 = c(1/rgamma(1,a_1,b_1))
  theta = matrix(NA,Nsim+1,k)
  tau2 = c(1/rgamma(1,a_2,b_2))
  theta[1,] = rnorm(2,mu[1],sqrt(tau2[1]))
  
  theta_i = function(mu,B_i,tau2) rnorm(1,B_i*mu + (1-B_i)*mean(cluster_j), sqrt(tau2*B_i))
  
  mu_ = function(theta,C_i,sigma2_mu) rnorm(1,C_i*mu_0 + (1-C_i)*mean(theta), sqrt(sigma2_mu*C_i))
  
  sigma2_ = function(theta) 1/rgamma(1,n*k/2 + a_1, (1/2)*sum((data[,1] - theta[1])^2 + (data[,2] - theta[2])^2) + b_1)
  
  tau2_ = function(theta,mu) 1/rgamma(1,k/2 + a_2, (1/2)*sum((theta - mu)^2) + b_2)
  
  sigma2_mu_ = function(mu) 1/rgamma(1, 1/2 + a_3, (1/2)*(mu-mu_0)^2 + b_3)
  
  
  for(i in 2:(Nsim+1)){
    
    #Se actualiza theta
    for(j in 1:k){
      cluster_j = Energy[,j]
      B_i = sigma2[i-1]/(sigma2[i-1] + length(cluster_j)*tau2[i-1])
      theta[i,j] = theta_i(mu[i-1],B_i,tau2[i-1])
    }
    
    #Se actualiza mu
    C_i = tau2[i-1]/(tau2[i-1] + k*sigma2_mu[i-1])
    mu[i] = mu_(theta[i,],C_i,sigma2_mu[i-1])
    
    #Se actualiza sigma2
    sigma2[i] = sigma2_(theta[i,])
    
    #Se actualiza tau2
    tau2[i] = tau2_(theta[i,],mu[i])
    
    #Se actualiza sigma2_mu
    sigma2_mu[i] = sigma2_mu_(mu[i])
  }
  
  sigma2_mu = sigma2_mu[-1]
  mu = mu[-1]
  sigma2 = sigma2[-1]
  theta = theta[-1,]
  tau2 = tau2[-1]
  
  return(list(sigma2_mu = sigma2_mu, mu = mu, sigma2 = sigma2, theta = theta, tau2 = tau2))
  
}

results = Gibss_hierarchical()

#Extraemos resultados generales

mu = results$mu
theta =results$theta

#Obtenemos desviaciones estándar

sigma_mu = sqrt(results$sigma2_mu)
sigma = sqrt(results$sigma2)
tau = sqrt(results$tau2)

#Histograma de resultados (aplicando log)

hist(log(sigma_mu),breaks=80)
hist(log(mu),breaks=100,xlim=c(6.9,6.920))
hist(log(sigma),breaks=40)
hist(log(theta[,1]),breaks=40)
hist(log(theta[,2]),breaks=40)
hist(log(tau),breaks=40)
