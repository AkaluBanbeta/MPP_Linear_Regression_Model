 
f0 <- function(weight,X0,Y0,n0,dim,a0,b0,mu0,Lambda0,X,Y,n){
     weightedsumX0TX0 <- matrix(0,nrow=dim,ncol=dim)
     weightedsumX0TY0 <- matrix(0,nrow=dim,ncol=1)
     weightedsumY0TY0 <- 0
     weightedsumn0 <- 0
   for (i in 1:n_hist_trials){
     weightedsumX0TX0<- weightedsumX0TX0+ weight[i]*t(X0[[i]])%*%X0[[i]]
     weightedsumX0TY0<- weightedsumX0TY0+ weight[i]*t(X0[[i]])%*%Y0[[i]]
     weightedsumY0TY0<- weightedsumY0TY0+ weight[i]*t(Y0[[i]])%*%Y0[[i]]
     weightedsumn0    <- weightedsumn0+ weight[i]*n0[i]
     }
    Lambda_star = t(X)%*%X + weightedsumX0TX0 + Lambda0
    mu_star=  solve(Lambda_star)%*%(t(X)%*%Y + weightedsumX0TY0+Lambda0 %*%mu0)
    a_star=a0+(weightedsumn0+n)/2
    b_star = b0 + 1/2*( t(Y)%*%Y + weightedsumY0TY0 +t(mu0)%*% Lambda0%*%mu0 - 
                        t(mu_star)%*% Lambda_star%*%mu_star)
  fact1=(n/2+weightedsumn0/2)*log(1/(2*pi))
  fact2= 1/2*logdet(Lambda0)-1/2*logdet(Lambda_star)  +a0*log(b0) - a_star*log(b_star)+lgamma(a_star)-lgamma(a0) 
  log_ML0=fact2 + fact1
    Lambda_star_sc = weightedsumX0TX0 + Lambda0
    mu_star_sc=  solve(Lambda_star_sc)%*%(weightedsumX0TY0+Lambda0 %*%mu0)
    a_star_sc=a0+(weightedsumn0)/2
    b_star_sc= b0 + 1/2*(weightedsumY0TY0 +t(mu0)%*% Lambda0%*%mu0 - 
                         t(mu_star_sc)%*% Lambda_star_sc%*%mu_star_sc)
  fact1=(n/2+weightedsumn0/2)*log(1/(2*pi))
  fact2= 1/2*logdet(Lambda0)-1/2*logdet(Lambda_star)  +a0*log(b0) - a_star*log(b_star)+lgamma(a_star)-lgamma(a0) 
  log_ML0=fact2 + fact1
  fact1_sc=(weightedsumn0/2)*log(1/(2*pi))
  fact2_sc= 1/2*logdet(Lambda0) - 1/2*logdet(Lambda_star_sc) + a0*log(b0) - a_star_sc*log(b_star_sc)+lgamma(a_star_sc)-lgamma(a0)
  log_SC0=fact2_sc + fact1_sc 
  ML0_SC0 = log_ML0 - log_SC0
  return(list(Lambda_star=Lambda_star,mu_star=mu_star,a_star=a_star,b_star=b_star,
              Lambda_star_sc=Lambda_star_sc,mu_star_sc=mu_star_sc,a_star_sc=a_star_sc,b_star_sc=b_star_sc,
              log_ML0=log_ML0,log_SC0=log_SC0,ML0_SC0=ML0_SC0))
  }
