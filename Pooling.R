#Bayesian Pooled data analysis (including all data)


#Preparing data for analysis
source("../Data_preparation",local=TRUE)


#analysis
N=sum(n0)+n
Y_Pooled=response
Lambda_Pooled = t(X_Pooled)%*%X_Pooled + Lambda0
mu_Pooled=  solve(Lambda_Pooled)%*%(t(X_Pooled)%*%Y_Pooled + Lambda0 %*%mu0)
a_Pooled=a0+(N)/2
b_Pooled = b0 + 1/2*( t(Y_Pooled)%*%Y_Pooled + t(mu0)%*% Lambda0%*%mu0 - 
                        t(mu_Pooled)%*% Lambda_Pooled%*%mu_Pooled)
sigma2_Pooled<-rinvgamma(n=niter,shape=a_Pooled,b_Pooled)
beta_Pooled<-matrix(0,nrow=niter,ncol=length(mu_Pooled))
for(i in 1:niter){
  beta_Pooled[i,]<-rmnorm(n=1,mu_Pooled,sigma2_Pooled[i]*solve(Lambda_Pooled))
}