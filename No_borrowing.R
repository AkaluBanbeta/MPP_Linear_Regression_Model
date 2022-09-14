#Bayesian analysis based on Current data (ignoring historical data)

# Preparing data for analysis
source("../Data_preparation",local=TRUE)

Lambda_Current = t(X)%*%X + Lambda0
mu_Current=  solve(Lambda_Current)%*%(t(X)%*%Y + Lambda0 %*%mu0)
a_Current=a0+(n)/2
b_Current = b0 + 1/2*( t(Y)%*%Y + t(mu0)%*% Lambda0%*%mu0 - 
                         t(mu_Current)%*% Lambda_Current%*%mu_Current)
sigma2_Current<-rinvgamma(n=niter,shape=a_Current,b_Current)
beta_Current<-matrix(0,nrow=niter,ncol=length(mu_Current))
for(i in 1:niter){
  beta_Current[i,]<-rmnorm(n=1,mu_Current,sigma2_Current[i]*solve(Lambda_Current))
}