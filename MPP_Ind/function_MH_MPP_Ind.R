#Metropolis-Hastings (MH) algorithm for the MPP Ind
  MH_MPP_Ind = function(n00,weight0,X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n,post_mean,post_sd,n_hist_trials){  
  s_weight=matrix(0,nrow=(n00+1),ncol=n_hist_trials)  ### the logit(weights) to be stored 
  weights=matrix(0,nrow=(n00+1),ncol=n_hist_trials)
  yesno=matrix(0,nrow=(n00+1),ncol=1)                ### for acceptance rate  # weight[n00+1] to be removed later
  return_f03=list(n00)
  ML_SC = rep(NA,n00)
  posterior = rep(NA,n00)
  
  weights[1,]=weight0
  s_weight[1,] = log(weight0/(1-weight0))
  
  return_f03[[1]]<- f0(weight=expit(s_weight[1,]),X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
  ML_SC[1]<- return_f03[[1]]$ML0_SC0
  posterior[1] = ML_SC[1] + 
    sum(dlogitbeta(s_weight[1,],1,1,log=TRUE))
  
  for (i in 1:(n00)){  
    
    s_weight[(i+1),]<-rnorm(n=n_hist_trials,mean=digamma(1)-digamma(1),
                            sd=sqrt(trigamma(1)+trigamma(1)))
    
    return_f03[[i+1]] <- f0(weight=expit(s_weight[(i+1),]),X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
    ML_SC[i+1] <- return_f03[[i+1]]$ML0_SC0
    posterior[(i+1)] = ML_SC[(i+1)] + 
      sum(dlogitbeta(s_weight[(i+1),],1,1,log=TRUE)) 
    
    lr = posterior[i+1]-posterior[i]-
      sum(dnorm(s_weight[i+1,],mean=digamma(1)-digamma(1),
                sd=sqrt(trigamma(1)+trigamma(1)),log=TRUE))+
      sum(dnorm(s_weight[i,],mean=digamma(1)-digamma(1),
                sd=sqrt(trigamma(1)+trigamma(1)),log=TRUE)) 
    
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
      posterior[i+1]<-posterior[i]
      yesno[i]<-0         
    }   
  }
  return(list(s_weight=expit(s_weight[-(n00+1),]),yesno=yesno,posterior=posterior))
}

