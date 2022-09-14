    #Making analysis based on eight methods; four versions of the modified power prior (MPP),
    #two versions of Meta-analytic predictive (MAP) prior, pooling all data and ignoring the historical data
    

    #Data preparation
    source("../Data_preparation",local=TRUE)
    
  
    #Function which returns ML, SC,lambda,mu,a,b with: given data and weight values
    source("../The MPP/function0_ML.R",local=TRUE)
 
    
    #Estimating the weight parameters without prior
    source("../The MPP/weightwithoutprior.R",local=TRUE)

    #Function that gives alpha and beta given mu and sigmasq of beta distribution
    source("../The MPP/function_estBetaParams.R",local=TRUE)
    
    
    #Function used to sample the weights in logit scale 
    source("../The MPP/function_logitbeta.R",local=TRUE)
    
    
    #Metropolist-Hastings (MH) algorithm for four versions of MPP 
    source("../The MPP/function_MH_MPPs.R",local=TRUE)
    
    
    #MH with independent power paarameter (MPP_Ind) 
    return_MH_MPP_Ind=MH_MPP_Ind(niter,post_mean,X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n,post_mean,post_sd,n_hist_trials)
    postweight_MPP_Ind=return_MH_MPP_Ind$s_weight
    n_MPP_Ind_post <- length(postweight_MPP_Ind[,1]) 
    ##  Computation of other model parameters MPP_Ind
    Lambda_MPP_Ind = list(n_MPP_Ind_post)
    mu_MPP_Ind = list(n_MPP_Ind_post)
    a_MPP_Ind = list(n_MPP_Ind_post)
    b_MPP_Ind = list(n_MPP_Ind_post)
    return_f04 = list(n_MPP_Ind_post)
    for (i in 1: n_MPP_Ind_post){
      return_f04[[i]] <- f0(weight=postweight_MPP_Ind[i,],X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
      Lambda_MPP_Ind[[i]] =  return_f04[[i]]$Lambda_star
      mu_MPP_Ind[[i]]=  return_f04[[i]]$mu_star
      a_MPP_Ind[[i]] = return_f04[[i]]$a_star
      b_MPP_Ind[[i]] = return_f04[[i]]$b_star
    }
    sigma2_MPP_Ind <- rinvgamma(n=length(unlist(a_MPP_Ind)),shape=unlist(a_MPP_Ind),unlist(b_MPP_Ind))
    beta_MPP_Ind <- matrix(0,nrow=length(unlist(a_MPP_Ind)),ncol=length(mu_MPP_Ind[[1]]))
    for(i in 1:length(unlist(a_MPP_Ind))){
      beta_MPP_Ind[i,] <- rmnorm(n=1,mu_MPP_Ind[[i]],sigma2_MPP_Ind[i]*solve(Lambda_MPP_Ind[[i]]))
    }
    result_MPP_Ind<-as.mcmc(as.data.frame((beta_MPP_Ind))) 
    summary(result_MPP_Ind)
    
    
    #MH with dependent weights (DMPP)
    return_MH_DMPP=MH_DMPP(niter,post_mean,X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n,post_mean,post_sd,n_hist_trials)
    #posterior of power parameters Robust_DMPP_1
    postweight_DMPP=return_MH_DMPP$s_weight
    n_DMPP_post <- length(postweight_DMPP[,1])
    Lambda_DMPP = list(n_DMPP_post)
    mu_DMPP = list(n_DMPP_post)
    a_DMPP = list(n_DMPP_post)
    b_DMPP = list(n_DMPP_post)
    return_f04 = list(n_DMPP_post)
    for (i in 1: n_DMPP_post){
      return_f04[[i]] <- f0(weight=postweight_DMPP[i,],X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
      Lambda_DMPP[[i]] =  return_f04[[i]]$Lambda_star
      mu_DMPP[[i]]=  return_f04[[i]]$mu_star
      a_DMPP[[i]] = return_f04[[i]]$a_star
      b_DMPP[[i]] = return_f04[[i]]$b_star
    }
    sigma2_DMPP <- rinvgamma(n=length(unlist(a_DMPP)),shape=unlist(a_DMPP),unlist(b_DMPP))
    beta_DMPP <- matrix(0,nrow=length(unlist(a_DMPP)),ncol=length(mu_DMPP[[1]]))
    for(i in 1:length(unlist(a_DMPP))){
      beta_DMPP[i,] <- rmnorm(n=1,mu_DMPP[[i]],sigma2_DMPP[i]*solve(Lambda_DMPP[[i]]))
    }
    result_DMPP<-as.mcmc(as.data.frame((beta_DMPP))) 
    summary(result_DMPP)
    
    
    
    ## MH for Robust_DMPP_1
    return_MH_Robust_DMPP_1=MH_Robust_DMPP_1(niter,post_mean,X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n,post_mean,post_sd,n_hist_trials)
    # posterior of power parameters Robust_DMPP_1
    postweight_Robust_DMPP_1=return_MH_Robust_DMPP_1$s_weight 
    n_Robust_DMPP_1_post <- length(postweight_Robust_DMPP_1[,1])
    ################  Computation of model parameters Robust_DMPP_1
    Lambda_Robust_DMPP_1 = list(n_Robust_DMPP_1_post)
    mu_Robust_DMPP_1 = list(n_Robust_DMPP_1_post)
    a_Robust_DMPP_1 = list(n_Robust_DMPP_1_post)
    b_Robust_DMPP_1 = list(n_Robust_DMPP_1_post)
    return_f04 = list(n_Robust_DMPP_1_post)
    for (i in 1: n_Robust_DMPP_1_post){
      return_f04[[i]] <- f0(weight=postweight_Robust_DMPP_1[i,],X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
      Lambda_Robust_DMPP_1[[i]] =  return_f04[[i]]$Lambda_star
      mu_Robust_DMPP_1[[i]]=  return_f04[[i]]$mu_star
      a_Robust_DMPP_1[[i]] = return_f04[[i]]$a_star
      b_Robust_DMPP_1[[i]] = return_f04[[i]]$b_star
    }
    sigma2_Robust_DMPP_1 <- rinvgamma(n=length(unlist(a_Robust_DMPP_1)),shape=unlist(a_Robust_DMPP_1),unlist(b_Robust_DMPP_1))
    beta_Robust_DMPP_1 <- matrix(0,nrow=length(unlist(a_Robust_DMPP_1)),ncol=length(mu_Robust_DMPP_1[[1]]))
    for(i in 1:length(unlist(a_Robust_DMPP_1))){
      beta_Robust_DMPP_1[i,] <- rmnorm(n=1,mu_Robust_DMPP_1[[i]],sigma2_Robust_DMPP_1[i]*solve(Lambda_Robust_DMPP_1[[i]]))
    }
    result_Robust_DMPP_1<-as.mcmc(as.data.frame((beta_Robust_DMPP_1))) 
    summary(result_Robust_DMPP_1)
   
    
    #MH for Robust_DMPP_2
    return_MH_Robust_DMPP_2=MH_Robust_DMPP_2(niter,post_mean,X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n,post_mean,post_sd,n_hist_trials)
    ## posterior weights
    postweight_Robust_DMPP_2=return_MH_Robust_DMPP_2$s_weight
    n_Robust_DMPP_2_post <- length(postweight_Robust_DMPP_2[,1])
    ##Computation of other model parameters Robust_DMPP_2
    Lambda_Robust_DMPP_2 = list(n_Robust_DMPP_2_post)
    mu_Robust_DMPP_2 = list(n_Robust_DMPP_2_post)
    a_Robust_DMPP_2 = list(n_Robust_DMPP_2_post)
    b_Robust_DMPP_2 = list(n_Robust_DMPP_2_post)
    return_f04 = list(n_Robust_DMPP_2_post)
    for (i in 1: n_Robust_DMPP_2_post){
      return_f04[[i]] <- f0(weight=postweight_Robust_DMPP_2[i,],X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
      Lambda_Robust_DMPP_2[[i]] =  return_f04[[i]]$Lambda_star
      mu_Robust_DMPP_2[[i]]=  return_f04[[i]]$mu_star
      a_Robust_DMPP_2[[i]] = return_f04[[i]]$a_star
      b_Robust_DMPP_2[[i]] = return_f04[[i]]$b_star
    }
    sigma2_Robust_DMPP_2 <- rinvgamma(n=length(unlist(a_Robust_DMPP_2)),shape=unlist(a_Robust_DMPP_2),unlist(b_Robust_DMPP_2))
    beta_Robust_DMPP_2 <- matrix(0,nrow=length(unlist(a_Robust_DMPP_2)),ncol=length(mu_Robust_DMPP_2[[1]]))
    for(i in 1:length(unlist(a_Robust_DMPP_2))){
      beta_Robust_DMPP_2[i,] <- rmnorm(n=1,mu_Robust_DMPP_2[[i]],sigma2_Robust_DMPP_2[i]*solve(Lambda_Robust_DMPP_2[[i]]))
    }
    result_Robust_DMPP_2<-as.mcmc(as.data.frame((beta_Robust_DMPP_2))) 
    summary(result_Robust_DMPP_2)
    
    
    #Bayesian meta-analystic-predictive (MAP) analysis 
    MAP_data = list(response=response, Intervention=Intervention,cata=cata,cont=cont,
                    Trial=dataset$Trial,minTrial=min(dataset$Trial),maxTrial=max(dataset$Trial),N=N) 
    MAP_inits <- function() {list (beta0=0,beta1=0,beta2=0,beta3=0,isigmasq=1,tau=1,.RNG.name="base::Mersenne-Twister",.RNG.seed=80000+sim)}
    MAP_parameters.to.save<- c("beta0","beta1","beta2","beta3","sigmasq","tau")
    MAP_model<-jags.model(file="MAPmodel.bug", data=  MAP_data, inits=MAP_inits , n.chains=1, n.adapt=2*niter) 
    update(MAP_model, n.iter=2*niter) # burn in
    result_MAP<-coda.samples(MAP_model, variable.names=MAP_parameters.to.save,thin = 5, n.iter=5*niter)
    varnames(result_MAP)[1]<-"Beta0"
    varnames(result_MAP)[2]<-"Treatment_effect"
    varnames(result_MAP)[3]<-"Categorical"
    varnames(result_MAP)[4]<-"Continuous"
    varnames(result_MAP)[5]<-"sigmasq"
    varnames(result_MAP)[6]<-"tau"
    summary(result_MAP)
    
    
    #Bayesian Robust_MAP data analysis  
    Robust_MAP_data = list(response=response, Intervention=Intervention,cata=cata,cont=cont,
                           Trial=dataset$Trial,minTrial=min(dataset$Trial),maxTrial=max(dataset$Trial),N=N) 
    Robust_MAP_inits <- function() {list (beta0=0,beta1=0,beta2=0,beta3=0,isigmasq=1,tau=1,wr=2,.RNG.name="base::Mersenne-Twister",.RNG.seed=80000+sim)}
    Robust_MAP_parameters.to.save<- c("beta0","beta1","beta2","beta3","sigmasq","tau","wr_post")
    Robust_MAP_model<-jags.model(file="Robust_MAPmodel", data= Robust_MAP_data, inits=Robust_MAP_inits , n.chains=1, n.adapt=2*niter)
    update(Robust_MAP_model, n.iter=2*niter) # burn in
    result_Robust_MAP<-coda.samples(Robust_MAP_model, variable.names=Robust_MAP_parameters.to.save,thin = 5, n.iter=5*niter)
    varnames(result_Robust_MAP)[1]<-"Beta0"
    varnames(result_Robust_MAP)[2]<-"Treatment_effect"
    varnames(result_Robust_MAP)[3]<-"Categorical"
    varnames(result_Robust_MAP)[4]<-"Continuous"
    varnames(result_Robust_MAP)[5]<-"sigmasq"
    varnames(result_Robust_MAP)[6]<-"tau"
    varnames(result_Robust_MAP)[7]<-"wr_post"
    summary(result_Robust_MAP)
  
  
