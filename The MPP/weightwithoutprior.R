############  Estimate of the weight parameter without prior
st<-rep(0,5000)
w<-replicate(n_hist_trials,rep(0,5000))
for(i in 1:5000){
  w[i,] <- runif(n_hist_trials)
  return_f0 <- f0(w[i,],X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n)
  st[i]<- return_f0$ML0_SC0
  if(st[i]>=-720){
    st[i]<-st[i]}
  else{
    st[i]<--720
  }
}
post_mean<-colMeans(w*exp(st)/mean(exp(st))) ## the average of weights
post_sd<-sqrt(colMeans(exp(st)/mean(exp(st))*(w-mean(w*exp(st)/mean(exp(st))))^2)) ## sd of weights 