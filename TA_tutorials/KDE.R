set.seed(1234)
n = 400
X = c()
curve(.7*dnorm(x,2) + .1*dnorm(x,7) + .2*dnorm(x,-2),-5,10)

for(i in 1:n){
  U = runif(1)
  if(U <= .1){
    X = c(X, rnorm(1,7))
  } else if (U <= .2){
    X = c(X, rnorm(1,-2))
  } else{
    X = c(X, rnorm(1,2))
  }
}

original_density = curve(.7*dnorm(x,2) + .1*dnorm(x,7) + .2*dnorm(x,-2),-5,10)
hist(X,freq = F)
sample_ker = density(X)
lines(sample_ker)
lines(original_density,col="blue")

KDE = function(x,h) mean(dnorm(x-X,0,h))
h_opt = 1.06*sd(X)*length(X)^(-1/5)
KDE.fit <- Vectorize( function(x) KDE(x,h_opt) )

plot(KDE.fit,-5,10)
lines(original_density,col="blue")
