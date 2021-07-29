library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(dplyr)

#Aquí introduzcan en donde guardaron el archivo LungCapData.txt
df <- read.csv("C:/Users/marco/Desktop/Simulación Estocástica/Ayudantias/LungCapData.txt",sep="")



x_mat1 <- df %>% select(Age,Height)

info_list1 <- list(N=nrow(df),K=2,x=x_mat1,y=df$LungCap)

model1 <- stan("linear_regression.stan",data=info_list1)

x_mat2 <- df[,-c(1)]

info_list2 <- list(N=nrow(df),K=5,x=x_mat2,y=df$LungCap)

model2 <- stan("linear_regression.stan",data=info_list2)


fit1 <- rstan::extract(model1)

plot(density(fit1$alpha,type="l"))


print(model1)


fit2 <- rstan::extract(model2)

print(model2)

hist(fit2$alpha)
hist(fit2$beta[,1])
hist(fit2$beta[,2])
hist(fit2$beta[,3])
hist(fit2$beta[,4])
hist(fit2$sigma)
