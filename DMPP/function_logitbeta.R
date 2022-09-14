dlogitbeta<-function(x,alpha,beta,log){
  #Density function of logit transform of beta distribution.
  #For explanation of the density, see 
  #https://stats.stackexchange.com/questions/286042/derivation-of-the-mean-and-variance-of-the-logit-transform-of-a-beta-random-vari
  
  if(log==TRUE){
    res<-rep(0,length(x))
    for(i in 1:length(x)){
      if(x[i]< -30){
        res[i]<- -lbeta(alpha,beta)+alpha*x[i]+beta*log(expit(-x[i]))
      }else if(x[i]> 30){
        res[i]<- -lbeta(alpha,beta)+alpha*log(expit(x[i]))-beta*x[i]
      }
      else{
        res[i]<- -lbeta(alpha,beta)+alpha*log(expit(x[i]))+beta*log(expit(-x[i]))
      }
    }
  }
  res
}