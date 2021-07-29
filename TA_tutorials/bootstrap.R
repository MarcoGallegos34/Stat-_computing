set.seed(12345)
n = 1000
X <- rgamma(n,4,1)
samp_curve = curve(dgamma(x,4,1),0,12)
hist(X,breaks=40,freq = F,border="grey"#,density = 2
     )
lines(samp_curve)


visualization_bootstrap = function(n=10000,B=n,seed=12345){
  
  set.seed(seed)
  X <- rgamma(n,4,1)
  samp_boot = c()
  
  for(i in 1:B){
    samp_boot = c(samp_boot,mean(sample(X,n,replace = T)))
  }
  
  samp_boot_curve = curve(dgamma(x,4*n,n),3.3,4.7)

  hist(samp_boot,breaks = 40,
       freq = F,
       border = "grey",
       main= expression(paste("Distribución de ",bar(X)," generada por medio de Bootstrap")))
  lines(samp_boot_curve)
  abline(v=mean(X),col="red",lty="dashed")
  abline(v=4,col="orange")

}

visualization_bootstrap(n=600)


#######################################################################
### Algoritmo Bootstrap sin definir función visualization_bootstrap ###
#######################################################################

B = n
samp_boot = c()

for(i in 1:B){
  samp_boot = c(samp_boot,mean(sample(X,n,replace = T)))
}

if(n==1000){
  
  samp_boot_curve = curve(dgamma(x,4*n,n),3.3,4.7)
  
} else if(n==10000){
  
  samp_boot_curve = curve(dgamma(x,4*n,n),3.8,4.2)
}

hist(samp_boot,breaks = 40,
     freq = F,
     border = "grey",
     main= expression(paste("Distribución de ",bar(X)," generada por medio de Bootstrap")),
     xlab= expression(bar(X)))
lines(samp_boot_curve)
abline(v=mean(X),col="red",lty="dashed")
abline(v=4,col="orange")

if(n==1000){
  
  legend(4.15,5.5,legend = c(expression(mu),expression(bar(x))),col=c("orange","red"),lty=1:2)
  
} else if(n==10000){
  
  legend(4.06,15,legend = c(expression(mu),expression(bar(x))),col=c("orange","red"),lty=1:2)
}