
#Metropolis-Hastings (MH) algorithm for the DMPP
  MH_DMPP = function(n00,weight0,X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n,post_mean,post_sd,n_hist_trials){  
  s_weight=matrix(0,nrow=(n00+1),ncol=n_hist_trials)  ### the logit(weights) to be stored 
  weights=matrix(0,nrow=(n00+1),ncol=n_hist_trials)
  yesno=matrix(0,nrow=(n00+1),ncol=1)                 ### for acceptance rate  # weight[n00+1] to be removed later
  return_f03=list(n00)
  ML_SC = rep(NA,n00)
  posterior = rep(NA,n00)
  
  mu_weight= rep(NA,n00)
  sigmasq_weight= rep(NA,n00)
  alpha_DMPP= rep(NA,n00)
  beta_DMPP= rep(NA,n00)
  
  weights[1,]=weight0
  s_weight[1,] = log(weight0/(1-weight0))
  
  return_f03[[1]]<- f0(weight=expit(s_weight[1,]),X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
  ML_SC[1]<- return_f03[[1]]$ML0_SC0
  mu_weight[1]<-rtnorm(1,0.3,0.3,0,1)
  sigmasq_weight[1]<-rtgamma(n=1,shape=1,scale=5,truncation=mu_weight[1]*(1-mu_weight[1]))
  a<-estBetaParams(mu_weight[1],sigmasq_weight[1])
  alpha_DMPP[1]<-a$alpha_weight
  beta_DMPP[1]<-a$beta_weight
  posterior[1] = ML_SC[1] + 
    sum(dlogitbeta(s_weight[1,],alpha_DMPP[1],beta_DMPP[1],log=TRUE)) + 
    dmvnorm(log(c(alpha_DMPP[1],beta_DMPP[1])),mean=c(1,1),sigma=matrix(c(1,0.5,0.5,1),nrow=2),log=TRUE)
  
  for (i in 1:(n00)){  
    
    mu_weight[i+1]<-rtnorm(1,0.3,0.3,0,1)
    sigmasq_weight[i+1]<-rtgamma(n=1,shape=1,scale=5,truncation=mu_weight[i+1]*(1-mu_weight[i+1]))
    
    a<-estBetaParams(mu_weight[i+1],sigmasq_weight[i+1])
    alpha_DMPP[i+1]<-a$alpha_weight
    beta_DMPP[i+1]<-a$beta_weight  
    
    s_weight[(i+1),]<-rnorm(n=n_hist_trials,mean=digamma(alpha_DMPP[i+1])-digamma(beta_DMPP[i+1]),
                            sd=sqrt(trigamma(alpha_DMPP[i+1])+trigamma(beta_DMPP[i+1])))
    
    return_f03[[i+1]] <- f0(weight=expit(s_weight[(i+1),]),X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
    ML_SC[i+1] <- return_f03[[i+1]]$ML0_SC0
    posterior[(i+1)] = ML_SC[(i+1)] + 
      sum(dlogitbeta(s_weight[(i+1),],alpha_DMPP[i+1],beta_DMPP[i+1],log=TRUE)) + 
      dmvnorm(log(c(alpha_DMPP[i+1],beta_DMPP[i+1])),mean=c(1,1),sigma=matrix(c(1,0.5,0.5,1),nrow=2),log=TRUE)
    
    lr = posterior[i+1]-posterior[i]-
      sum(dnorm(s_weight[i+1,],mean=digamma(alpha_DMPP[i+1])-digamma(beta_DMPP[i+1]),
                sd=sqrt(trigamma(alpha_DMPP[i+1])+trigamma(beta_DMPP[i+1])),log=TRUE))+
      sum(dnorm(s_weight[i,],mean=digamma(alpha_DMPP[i])-digamma(beta_DMPP[i]),
                sd=sqrt(trigamma(alpha_DMPP[i])+trigamma(beta_DMPP[i])),log=TRUE)) -
      dtgamma(sigmasq_weight[i+1],shape=1,scale=5,truncation=mu_weight[i+1]*(1-mu_weight[i+1]),log=TRUE)-
      msm::dtnorm(mu_weight[i+1],0.3,0.3,0,1,log=TRUE) +
      dtgamma(sigmasq_weight[i],shape=1,scale=5,truncation=mu_weight[i]*(1-mu_weight[i]),log=TRUE)+
      msm::dtnorm(mu_weight[i],0.3,0.3,0,1,log=TRUE) 
    
    if(lr>700){
      lr<-700
    }
    if(lr< -700){
      lr<- -700
    }
    A = min(exp(lr),1)
    
    if(runif(1)<A){
      yesno[i]<-1 
    }
    else{    
      s_weight[(i+1),]<-s_weight[i,]   
      mu_weight[i+1]<-mu_weight[i]
      sigmasq_weight[i+1]<-sigmasq_weight[i]
      alpha_DMPP[i+1]<-alpha_DMPP[i]
      beta_DMPP[i+1]<-beta_DMPP[i]
      posterior[i+1]<-posterior[i]
      yesno[i]<-0         
    }   
  }
  return(list(s_weight=expit(s_weight[-(n00+1),]),yesno=yesno,posterior=posterior,
              mu_weight=mu_weight,sigmasq_weight=sigmasq_weight))
  
}

