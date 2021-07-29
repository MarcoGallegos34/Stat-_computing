set.seed(1234)
Nsim=5000 #initial values
n=15
a=3
b=7
X=t=array(0,dim=c(Nsim,1)) #init arrays

t[1]=rbeta(1,a,b) #init chains

X[1]=rbinom(1,n,t[1])

for (i in 2:Nsim){ #sampling loop
  X[i]=rbinom(1,n,t[i-1])
  t[i]=rbeta(1,a+X[i],n-X[i]+b)
}


betabi = function(x,a,b,n) gamma(n+1)/(gamma(x+1)*gamma(n-x+1))*gamma(x+a)*gamma(n-x+b)/gamma(n+a+b)*gamma(a+b)/(gamma(a)*gamma(b))


target_X = curve(betabi(x,a,b,n),-1,20)
hist(X,freq=F)
lines(target_X)

target_theta = curve(dbeta(x,a,b),0,1)
hist(t,freq=F)
lines(target_theta)
