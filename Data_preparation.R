
#packages
library(Rlab)
library(MCMCpack)
library(invgamma)
library(mnormt)
library(rjags)
library(tree)
library(msos)
library(expm)
library(msm)
library(heavy)
library(locfit)
library(mvtnorm)


#settings
n_patients<- 100                # Number of patients set per treatment arm per trial in simulation
n_hist_trials <- 3              # Number of historical trials
baseline_mean <- 10             # baseline mean of the linear regression model
sd_data <- 1                    # standard deviation of a data in a trial
sc<- 1                          # Relevant covariates with fixed coefficients, HOMOG
intervention_effect<- 0         # Treatment effect
adjustment<- 1                  # Covariates are added in the model
niter=10000                     # number of iterations


#list of scenarios in the simulation study
source("../Scenarios_list.R",local=TRUE)


set.seed(sim)  #Use the same seed for each scenario


trialnr<-rep(0,(n_hist_trials+2)*n_patients)
response<- rep(0,(n_hist_trials+2)*n_patients)
intervention<- rep(0,(n_hist_trials+2)*n_patients)
catagorical<- rep(0,(n_hist_trials+2)*n_patients)
continuous<- rep(0,(n_hist_trials+2)*n_patients)
trial_effect_indiv<- rep(0,(n_hist_trials+2)*n_patients)


### Historical trials
for (i in 1:(n_hist_trials)){
  intercept<-rnorm(n=1,mean= baseline_mean ,sd=sd_trial_effect_mean) 
  cata_effect<-rnorm(n=1,mean= cata_effect0,sd=sd_cata_effect)
  cont_effect<-rnorm(n=1,mean= cont_effect0,sd=sd_cont_effect)
  for (j in 1:n_patients){    
    trialnr[(i-1)*n_patients+j]<-i
    if(i==1){
      catagorical[(i-1)*n_patients+j] <-rbern(n=1,p_cata_01)
      continuous[(i-1)*n_patients+j] <- rnorm(1,mean_cont_01,1)}
    if(i==2){
      catagorical[(i-1)*n_patients+j] <-rbern(n=1,p_cata_02)
      continuous[(i-1)*n_patients+j] <- rnorm(1,mean_cont_02,1)}
    if(i==3){
      catagorical[(i-1)*n_patients+j] <-rbern(n=1,p_cata_03)
      continuous[(i-1)*n_patients+j] <- rnorm(1,mean_cont_03,1)}
    
    if(sc<=11){
      trial_effect_indiv[(i-1)*n_patients+j]<- intercept+cata_effect*catagorical[(i-1)*n_patients+j]+cont_effect*continuous[(i-1)*n_patients+j]}
    if(sc==12){
      trial_effect_indiv[(i-1)*n_patients+j]<- intercept+cata_effect*catagorical[(i-1)*n_patients+j]+cont_effect*continuous[(i-1)*n_patients+j]*continuous[(i-1)*n_patients+j]}
    if(sc==13){
      trial_effect_indiv[(i-1)*n_patients+j]<- intercept+cata_effect*catagorical[(i-1)*n_patients+j]+cont_effect*continuous[(i-1)*n_patients+j]+
        cont_effect*cont_effect*continuous[(i-1)*n_patients+j]*continuous[(i-1)*n_patients+j]}    
    response[(i-1)*n_patients+j]<-rnorm(n=1,trial_effect_indiv[(i-1)*n_patients+j],sd_data)
  }  
}


### Control group of the current trial 
i<- n_hist_trials+1  
intercept_current_trial<-rnorm(n=1,mean= baseline_mean ,sd=sd_trial_effect_mean)
cata_effect_current<-rnorm(n=1,mean= cata_effect0,sd=sd_cata_effect)
cont_effect_current<-rnorm(n=1,mean= cont_effect0,sd=sd_cont_effect)
for (j in 1:n_patients){  
  trialnr[(i-1)*n_patients+j]<-n_hist_trials+1  
  catagorical[(i-1)*n_patients+j] <-rbern(n=1,p_cata)
  continuous[(i-1)*n_patients+j] <- rnorm(1,mean_cont,1)
  if(sc<=11){
    trial_effect_indiv[(i-1)*n_patients+j]<- intercept_current_trial+cata_effect_current*catagorical[(i-1)*n_patients+j]+cont_effect_current*continuous[(i-1)*n_patients+j]}
  if(sc==12){
    trial_effect_indiv[(i-1)*n_patients+j]<- intercept_current_trial+cata_effect_current*catagorical[(i-1)*n_patients+j]+cont_effect_current*continuous[(i-1)*n_patients+j]*continuous[(i-1)*n_patients+j]}
  if(sc==13){
    trial_effect_indiv[(i-1)*n_patients+j]<- intercept_current_trial+cata_effect_current*catagorical[(i-1)*n_patients+j]+cont_effect_current*continuous[(i-1)*n_patients+j]+
      cont_effect_current*cont_effect_current*continuous[(i-1)*n_patients+j]*continuous[(i-1)*n_patients+j]}
  response[(i-1)*n_patients+j]<-rnorm(n=1,trial_effect_indiv[(i-1)*n_patients+j],sd_data)
} 


#Intervention group of current trial
i<-n_hist_trials+2   
for (j in 1:n_patients){  
  intervention[(i-1)*n_patients+j]<-1
  trialnr[(i-1)*n_patients+j]<- n_hist_trials+1
  catagorical[(i-1)*n_patients+j] <-rbern(n=1,p_cata)
  continuous[(i-1)*n_patients+j] <- rnorm(1,mean_cont,1)
  if(sc<=11){
    trial_effect_indiv[(i-1)*n_patients+j]<- intercept_current_trial+intervention_effect+cata_effect_current*catagorical[(i-1)*n_patients+j]+cont_effect_current*continuous[(i-1)*n_patients+j]}
  if(sc==12){
    trial_effect_indiv[(i-1)*n_patients+j]<- intercept_current_trial+intervention_effect+cata_effect_current*catagorical[(i-1)*n_patients+j]+cont_effect_current*continuous[(i-1)*n_patients+j]*continuous[(i-1)*n_patients+j]}
  if(sc==13){
    trial_effect_indiv[(i-1)*n_patients+j]<- intercept_current_trial+intervention_effect+cata_effect_current*catagorical[(i-1)*n_patients+j]+cont_effect_current*continuous[(i-1)*n_patients+j] + 
      cont_effect_current*cont_effect_current*continuous[(i-1)*n_patients+j]*continuous[(i-1)*n_patients+j]}
  response[(i-1)*n_patients+j]<-rnorm(n=1,trial_effect_indiv[(i-1)*n_patients+j],sd_data)
}  
dataset<-data.frame(Trial=trialnr,response=response,InterventionArray=intervention,catagorical=catagorical,continuous=continuous)


#Variable names   
response= dataset$response
Intervention= dataset$InterventionArray
cata= dataset$catagorical
cont= dataset$continuous  
Y0=list(n_hist_trials)
for(i in 1:n_hist_trials){
  Y0[[i]]=dataset$response[(1+(i-1)*n_patients):(i*n_patients)]}
Y=dataset$response[(n_hist_trials*n_patients+1):((n_hist_trials+2)*n_patients)]

if(adjustment==1){      # without covariates (cata and cont)
  X0=list(n_hist_trials)
  for(i in 1:n_hist_trials){
    X0[[i]]=cbind(1,dataset$InterventionArray[(1+(i-1)*n_patients):(i*n_patients)])}
  X=cbind(1,dataset$InterventionArray[(n_hist_trials*n_patients+1):((n_hist_trials+2)*n_patients)])
  X_Pooled=cbind(1,Intervention)  
}

if(adjustment==2){      # with covariates
  X0=list(n_hist_trials)
  for(i in 1:n_hist_trials){
    X0[[i]]=cbind( 1,Intervention[(1+(i-1)*n_patients):(i*n_patients)],cata[(1+(i-1)*n_patients):(i*n_patients)],cont[(1+(i-1)*n_patients):(i*n_patients)])}
  X=cbind(1,Intervention[(n_hist_trials*n_patients+1):((n_hist_trials+2)*n_patients)],cata[(n_hist_trials*n_patients+1):((n_hist_trials+2)*n_patients)],cont[(n_hist_trials*n_patients+1):((n_hist_trials+2)*n_patients)])
  X_Pooled=cbind(1,Intervention,cata,cont) 
}  

n0=rep(n_patients,n_hist_trials) # number of each historical trial
n= 2*n_patients                  # total number of patients in current trial
N=sum(n0)+n                      # total number of patients in the study
dim=ncol(X0[[1]])                # dimention of the covariate matrix



######  Initials for  Bayesian LM parametes  ######
a0  <- 0.001   ###### a0=v0/2; v0 = 0.02
b0  <- 0.001   ###### b0 = v0s0^2/2;  ssquare0 = 1
mu0 <- rep(0,dim)
Lambda0 <- diag(0.001,dim)
Lambda0_inv=solve(Lambda0)
