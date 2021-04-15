library("depmixS4")
library('lubridate')

id_all<-c(3, 7, 4,16,2121,16,20)
id<-10
infoall<-read.csv(paste0('S',id,'.csv'))
info.new = infoall[seq(1, nrow(infoall), 5), ]               
info=info.new



#write.csv(info, 'act_data.csv')

info$sqValueT<-sqrt(info$tempreture) # square root is applied to improve normality
info$sqValueA<-sqrt(info$activity)
#temp <- info$tempreture
#timet<- info$time
#temp_data <- cbind(timet, temp)
#temp_data <-as.data.frame(temp_data)
#write.csv(temp_data,'temp_data.csv')
#set.seed(1)
mod <- depmix(sqValueA ~ 1, data = info, nstates = 2, family = gaussian(), transition = ~ scale(sqValueT), instart = runif(2))
fm <- fit(mod, verbose = FALSE)
summary(fm)
summary(fm, which = "transition")



y<-depmix(sqValue~1,data=info,nstates=3,family=gaussian(),ntimes=nrow(info))
HMMs<-fit(y,verbose=F)


sf<-1/12 #sampling frequency in hour unit, i.e. here sf=5min/60min=1/12
lag<-seq(0,(nrow(info)-1))
info$sin_part<-sin(2*pi*lag*sf/24)
info$cos_part<-cos(2*pi*lag*sf/24)

set.seed(6)
y_hmc<- depmix(sqValue~1,transition=~sin_part+cos_part,data=info,nstates=3,family=gaussian(),ntimes=nrow(info))
HMM_hmc<-fit(y_hmc,verbose=F)

source('Homogeneous_HMMs.R')
HMM_results<-Homogeneous_HMMs(HMM=HMMs)

HMM_results

source('Harmonic_HMMs.R')


HMM_results_hmc<-Harmonic_HMMs(HMM=HMM_hmc,sin_part=info$sin_part,cos_part=info$cos_part)

source('figure_HMM_decoding_SP.R')
figure_HMM_decoding_SP(HMM_results,info, id)
#for fig to work colum name date on csv file must be chaned to time 

circadian_states_prob<-HMM_results_hmc$circadian_states_prob



#time<- info$time
#NROW(time)
#time <-as.character(time)

#circa<- data.frame(circadian_states_prob)
#df<- cbind(circa, time)
#ggplot()+geom_line(data=df, aes(x=time, y=state_3, group=3), color='red')+
  #geom_line(data =df, aes(x=time, y=state_2, group =2), color='green')+
  #geom_line(data =df, aes(x=time, y=state_1, group =1), color='blue')



#combs <- cbind(df, dfff) 
#ggplot()+geom_line(data=combs, aes(x=time, y=state_3, group=1), color='red')+
  #geom_line(data =combs, aes(x=time, y=state_2t, group =2), color='blue')+
  #geom_line(data =combs, aes(x=time, y=state_1, group =3), color='green')


time<-as.POSIXct(info$time, format="%Y-%m-%d %H:%M:%S")
circadian_states_prob
library('lubridate')
hour_day_start<-12
find_day_start<-function(hour_day_start,time){
  hour_time<-hour(time)
  min_time<-minute(time)
  A<-which(hour_time==hour_day_start)
  B<-which(min_time<60*sf)
  one_day_start<-A[min(which(A%in%B, arr.ind = TRUE))]
  return(one_day_start)
}

source('find_day_start.R')
one_day_start<-find_day_start(hour_day_start=hour_day_start,time=time) 
one_day_end<-one_day_start+(24/sf-1)

# one day profile
one_day_prob<-circadian_states_prob[one_day_start:one_day_end,]



p1<-one_day_prob$state_1 
p2<-one_day_prob$state_2 


p3<-one_day_prob$state_3 

rest_amount<-24*sum(p1)/nrow(one_day_prob)

index<-seq(0,24-1/12,1/12) #absolute index to compute the gravity centre of p1
center_rest<-sum(p1*index/sum(p1))

index_clocktime<-seq(hour_day_start,hour_day_start+24-sf,sf) 
clocktime<-index_clocktime 
clocktime[index_clocktime>24]<-(index_clocktime-24)[index_clocktime>24] 
center_rest_clock<-clocktime[which.min(abs(index-center_rest))]

worst_p1<-rep(mean(p1),24/sf)

source('find_perfect_p1.R')
perfect_p1<-find_perfect_p1(center_rest,rest_amount,index,sf)

RI<-(sum(p1[which(perfect_p1>0)])*sf/rest_amount-rest_amount/24)*24/(24- rest_amount)

index_position<-seq(0,24,2)
index_lable<-c('12','14','16','18','20','22',
               '0/24','2','4','6','8','10','12') #corresponding clock time

n_ML_states<- table(HMM_results_hmc$ML_states)
n_ML_states<- as.data.frame.array(n_ML_states)
n_ML_states <- t(n_ML_states)
n_ML_states<- as.data.frame(n_ML_states)
ave_state_duration <-n_states*5
trans_state_prob <- as.data.frame(HMM_results_hmc$transition_prob)
aave_trans_state_prob <-colMeans(HMM_results_hmc$transition_prob)
ave_circ_state_prob <- colMeans(HMM_results_hmc$circadian_states_prob)
ave_circ_state_prob <- as.data.frame(ave_circ_state_prob)
ave_circ_state_prob <- t(ave_circ_state_prob)

ave_state_duration <- n_ML_states * 5


ave_trans_state_prob<- t(aave_trans_state_prob)
ave_trans_state_prob<- as.data.frame(ave_trans_state_prob)
ave <- cbind(n_ML_states, ave_circ_state_prob, ave_trans_state_prob, ave_state_duration)
ave <-as.data.frame(ave)


#write_csv(ave, "ave.csv")

