#function integrating the area under charts
check_ov<-function(a,b){
  # by the above formula
c = floor((b-a)/(log(b)-log(a)))
value=ppois(c,b)  + 1 - ppois(c,a)
return(value)
}

#5%
check_ov(1,10) #~2%
check_ov(10,26) #~5%
#25%
check_ov(4,10) #~24%
check_ov(10,19)#~24%
#50%
check_ov(6,10) #~48%
check_ov(10,15) #~48%

#5%
check_ov(12,30) #~4.3%
check_ov(30,55) #~5.2%
#25%
check_ov(19,30) #26%
check_ov(30,44) # 25%
#50%
check_ov(23,30) #50%
check_ov(30,38) #49%

#5%
check_ov(57,90)
check_ov(90,131)
#25%
check_ov(70,90) #24%
check_ov(90,113)
#50%
check_ov(78,90) #51
check_ov(90,103) #50


library(ggplot2)
library(tidyverse)

#Dwell time distribution
#Based on the approximate overlap of 5%, 25% and 50%
dwell_assumed<-data.frame(mean=c(rep(10,3),rep(30,3),rep(90,3)),prop_diff=rep(c("5% overlap","25% overlap","50% overlap"),3),state1=c(1,4,6,12,19,23,57,70,78),state2=c(10,10,10,30,30,30,90,90,90),state3=c(26,19,15,55,44,38,131,113,103))

##### PLOT ####
##for the purposes of plotting the data
dwell_assumed_data1<-data.frame(scenario=c(rep("1 Mean 10; Proportion overlap: 5%",1500),rep("2 Mean 10; Proportion overlap: 25%",1500),rep( "3 Mean: 10; Proportion overlap: 50%",1500)),state=rep(c(rep(1,500),rep(2,500),rep(3,500)),3),value=NA)
dwell_assumed_data2<-data.frame(scenario=c(rep("4 Mean 30; Proportion overlap: 5%",1500),rep("5 Mean 30; Proportion overlap: 25%",1500),rep("6 Mean: 30; Proportion overlap: 50%",1500)),state=rep(c(rep(1,500),rep(2,500),rep(3,500)),3),value=NA)
dwell_assumed_data3<-data.frame(scenario=c(rep("7 Mean 90; Proportion overlap: 5%",1500),rep("8 Mean 90; Proportion overlap: 25%",1500),rep("9 Mean: 90; Proportion overlap: 50%",1500)),state=rep(c(rep(1,500),rep(2,500),rep(3,500)),3),value=NA)
dwell_assumed_data<-rbind(dwell_assumed_data1,dwell_assumed_data2,dwell_assumed_data3)
for(i in 1:9){
  dwell_assumed_data$value[(i*1500-1500+1):(i*1500)]<-c(rpois(500,lambda = dwell_assumed[i,3]),rpois(500,lambda = dwell_assumed[i,4]),rpois(500,lambda = dwell_assumed[i,5]))+1
}
#plot the data
neworder <- c("Mean 10; Proportion overlap: 5%","Mean 10; Proportion overlap: 25%","Mean: 10; Proportion overlap: 50%","Mean 30; Proportion overlap: 5%","Mean 30; Proportion overlap: 25%","Mean: 30; Proportion overlap: 50%","Mean 90; Proportion overlap: 5%","Mean 90; Proportion overlap: 25%","Mean: 90; Proportion overlap: 50%")
dwell_assumed_data$state<-as.factor(dwell_assumed_data$state)
dwell_assumed_data$scenario<-as.factor(dwell_assumed_data$scenario)
dwell_assumed_data<-dwell_assumed_data %>% arrange(mutate(dwell_assumed_data,
                                      scenario=factor(scenario,levels=neworder)),scenario)
#plot
dwell_assumed_data%>% ggplot(., aes(x=value,group=state,fill=state)) + geom_density(alpha=.3)+facet_wrap(~scenario,nrow =3 ,scales = "free_y")+theme_minimal()

