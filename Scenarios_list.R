#scenarios
if(sc==1){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==2){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-sd_data/8
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==3){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-sd_data/4
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==4){
  cata_effect0 <-0
  cont_effect0 <- 0
  sd_trial_effect_mean<-0
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==5){
  cata_effect0 <-0
  cont_effect0 <- 0
  sd_trial_effect_mean<-sd_data/8
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==6){
  cata_effect0 <-0
  cont_effect0 <- 0
  sd_trial_effect_mean<-sd_data/4
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==7){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0.2
  sd_cont_effect<-0
}
if(sc==8){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0.8
  sd_cont_effect<-0.8
}
if(sc==9){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==10){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==11){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==12){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0
  sd_cont_effect<-0
}
if(sc==13){
  cata_effect0 <-1
  cont_effect0 <- 0.5
  sd_trial_effect_mean<-0
  sd_cata_effect<-0
  sd_cont_effect<-0
}


if(sc<=8||sc>11){
  p_cata_01=p_cata_02=p_cata_03=p_cata=0.5
  mean_cont_01=mean_cont_02=mean_cont_03=mean_cont=0
}

if(sc==9) {
  p_cata_01=0.4
  p_cata_02=0.4
  p_cata_03=0.6
  p_cata=0.5
  mean_cont_01=mean_cont_02=mean_cont_03=mean_cont=0
}

if(sc==10) {
  p_cata_01=0.25
  p_cata_02=0.25
  p_cata_03=0.75
  p_cata=0.5
  mean_cont_01=mean_cont_02=mean_cont_03=mean_cont=0
}

if(sc==11) {
  p_cata_01=p_cata_02=p_cata_03=p_cata=0.5
  mean_cont_01=-1
  mean_cont_02=-1
  mean_cont_03=1
  mean_cont=0
}
