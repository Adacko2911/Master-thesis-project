library(tidyverse)
library(ggshadow)
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
colornames=c(colors$blue,colors$orange,colors$red,colors$teal)

#PREPARATION Figure B4
dwell_final <- read_csv("simulation study final/dwell_emp_final.csv") %>% filter(dwell_time!="mix")

dwell_final$model<-factor(dwell_final$model,levels = c("medhmm","mhmm"),labels = c("MEDHMM","MHMM"))
dwell_final$state<-factor(dwell_final$state,levels = c(3,4),labels = c("Three state scenarios", "Four state scenarios"))
#dwell_final$observations<-factor(dwell_final$observations,levels = c(200,500,1000), labels = c("200 observations", "500 observations","1000 observations"))
dwell_final$`Dwell time scenario`<-factor(dwell_final$dwell_time, levels = c(1.4,3.5,19.5,99.5),labels = c("Dwell time: 1.4","Dwell time: 3.5","Dwell time: 19.5","Dwell time: 99.5"))

dwell_final$estimates<-factor(dwell_final$estimates,levels = unique(dwell_final$estimates),labels = c("1st state","2nd state","3rd state","4th state"))

dwell_final$median_emp_dwell_true<-round(dwell_final$median_emp_dwell_true)

# Figure B4 ##########################################################################################################################################################################
pdf("~/GitHub/Master-thesis-project/Simulation Study/simulation study final/outputs plots/dwell_time_plot.pdf",width =9,height = 6)
dwell_final %>%arrange(dwell_time) %>%  ggplot(.,aes(x=observations,y=median_emp_dwell,color=estimates))+
  geom_line(data=dwell_final, mapping=aes(x=observations, y=median_emp_dwell_true),size=2.1,color="black")+
  geom_line(aes(linetype=model),size=1)+geom_point(size=2)+
  scale_color_manual(values =colornames ,name = "Estimate for:")+
  scale_x_continuous(breaks = c(200,500,1000))+
  scale_y_continuous(n.breaks = 4)+
  facet_wrap(~state+`Dwell time scenario`,nrow = 2,scales = "free_y",dir = "h")+
  labs(colour="Estimate for:", linetype="Model: ",y="Empirical state duration",x="Observations")+theme_minimal()+theme(legend.position = "bottom")
dev.off()
