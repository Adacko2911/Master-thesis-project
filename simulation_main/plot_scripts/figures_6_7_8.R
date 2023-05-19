library(tidyverse)
library(patchwork)
library(readr)
colors=tibble(blue="#5778a4",
              orange="#e49444",
              red="#d1615d",
              teal="#85b6b2",
              green="#6a9f58",
              yellow="#e7ca60",
              purple="#a87c9f",
              pink="#f1a2a9",
              brown="#967662",
              grey="#b8b0ac")

dwell_gamma_fin <- read_csv("simulation study final/dwell_gamma_finn2.csv") %>% filter(dwell_time!="mix")
dwell_gamma_fin<-dwell_gamma_fin[complete.cases(dwell_gamma_fin),] %>% dplyr::select(-c(mean,bias_mean,bias_mean_pp))
dwell_gamma_fin<-rbind(dwell_gamma_fin,c("medhmm",4,1000,3.5,"dwell_mu_bar_1",3.48,0,0,0,0,0),c("medhmm",4,500,99.5,"dwell_mu_bar_1",99.50,0,0,0,0,0))
#plot the expected calculated from the distributional parameters distributions
dwell_gamma_fin$model<-factor(dwell_gamma_fin$model,levels = c("medhmm","mhmm"),labels = c("MEDHMM","MHMM"))
dwell_gamma_fin$state<-factor(dwell_gamma_fin$state,levels = c(3,4),labels = c("Three state scenarios", "Four state scenarios"))
dwell_gamma_fin$observations<-factor(dwell_gamma_fin$observations,levels = c(200,500,1000), labels = c("200 observations", "500 observations","1000 observations"))
dwell_gamma_fin$`Dwell time scenario`<-factor(dwell_gamma_fin$dwell_time, levels = c(1.4,3.5,19.5,99.5),labels = c("1.4","3.5","19.5","99.5"))
dwell_gamma_fin$bias_median_pp%<>%as.numeric()
## absolute bias for expected dwell time
data<-dwell_gamma_fin %>% filter(estimates %in% c("dwell_mu_bar_1","dwell_mu_bar_2","dwell_mu_bar_3","dwell_mu_bar_4")) %>%
  mutate(bias_median =abs(bias_median_pp)) %>% group_by(model,state,observations,dwell_time) %>%
  mutate(abs_mean_bias=mean(bias_median),min_abs_mean_bias=min(bias_median),max_abs_mean_bias=max(bias_median),
         abs_cov=mean(coverage_pp),min_abs_cov=min(coverage_pp),max_abs_cov=max(coverage_pp))

pdf("~/GitHub/Master-thesis-project/Simulation Study/simulation study final/outputs plots/expected_dwell_biasv2.pdf",width = 6,height =5 )
data%>%
  ggplot(aes(x=`Dwell time scenario`,weight=abs_mean_bias,ymin=min_abs_mean_bias, ymax=max_abs_mean_bias,y=abs_mean_bias, fill=model))+
  geom_bar(aes(y=abs_mean_bias),stat = "identity",position = position_dodge())+
geom_errorbar(position = position_dodge(width=0.9), width=0.3, colour="black", alpha=0.7, size=0.2)+
  facet_wrap(~state+observations, nrow = 2,scales = "free_y")+
  labs(y="Average absolute % bias", title = "")+
  scale_fill_manual(values = c(colors$blue,colors$teal),name="Model:")+theme_minimal()+theme(legend.position = "bottom")
dev.off()


#absolute bias for the transition probabilities MEDHMM and the adjusted
states<-expand.grid(s1=1:4,s2=1:4)
states<-states[-c(1,6,11,16),]
states<-states %>% apply(.,1,function(x){paste0("gamma_prob_bar_S",x[1],"toS",x[2])})

post_est_final <- read_csv("simulation study final/dwell_gamma_finn2.csv") %>% filter(dwell_time!="mix")
post_est_final<-post_est_final %>%
  dplyr::select(-c(mean,bias_mean,bias_mean_pp))%>% filter(estimates %in% states)
dwell_gamma_fin<-post_est_final
dwell_gamma_fin<-rbind(dwell_gamma_fin,c("medhmm",4,1000,3.5,"gamma_prob_bar_S1toS2",0.300,0,0,0,0,0),c("medhmm",4,500,99.5,"gamma_prob_bar_S1toS2",0.300,0,0,0,0,0))

dwell_gamma_fin$model<-factor(dwell_gamma_fin$model,levels = c("medhmm","mhmm"),labels = c("MEDHMM","MHMM"))
dwell_gamma_fin$state<-factor(dwell_gamma_fin$state,levels = c(3,4),labels = c("Three state scenarios", "Four state scenarios"))
dwell_gamma_fin$observations<-factor(dwell_gamma_fin$observations,levels = c(200,500,1000), labels = c("200 observations", "500 observations","1000 observations"))
dwell_gamma_fin$`Dwell time scenario`<-factor(dwell_gamma_fin$dwell_time, levels = c(1.4,3.5,19.5,99.5),labels = c("1.4","3.5","19.5","99.5"))
dwell_gamma_fin$bias_median_pp%<>%as.numeric()
data<-dwell_gamma_fin %>%
  mutate(bias_median =abs(bias_median_pp)) %>% group_by(model,state,observations,dwell_time) %>%
  mutate(abs_mean_bias=mean(bias_median),min_abs_mean_bias=min(bias_median),max_abs_mean_bias=max(bias_median),
         abs_cov=mean(coverage_pp),min_abs_cov=min(coverage_pp),max_abs_cov=max(coverage_pp))
pdf("~/GitHub/Master-thesis-project/Simulation Study/simulation study final/outputs plots/transition_prob_biasv2.pdf",width = 6,height =5 )
data%>%
  ggplot(aes(x=`Dwell time scenario`,weight=abs_mean_bias,ymin=min_abs_mean_bias, ymax=max_abs_mean_bias,y=abs_mean_bias, fill=model))+
  geom_bar(aes(y=abs_mean_bias),stat = "identity",position = position_dodge())+
  geom_errorbar(position = position_dodge(width=0.9), width=0.3, colour="black", alpha=0.7, size=0.15)+
  facet_wrap(~state+observations, nrow = 2,scales = "free_y")+
  labs(y="Average absolute % bias", title = "")+
  scale_fill_manual(values = c(colors$blue,colors$teal),name="Model:")+theme_minimal()+theme(legend.position = "bottom")
dev.off()


