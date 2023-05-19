library(patchwork)
library(ggplot2)
library(tidyverse)
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
my_grey_text="#696969"


#plot 1a
cols <- c("col1" = colors$blue, "col2" = colors$orange, "col3" = colors$red, "col4"=colors$teal)
datay<-dnorm(seq(0,120,0.05), mean= 30,sd = sqrt(60))
data<-data.frame(x=seq(0,120,0.05),y1=datay,y2=dnorm(seq(0,120,0.05), mean= 60,sd = sqrt(120)))
plot1a<-ggplot(data=data,aes(x=x,y=y1))+geom_line(size=1,aes(color="col1"))+geom_line(size=1,aes(y=y2,color="col2"))+
  theme_minimal()+xlim(c(0,120))+ggtitle("Dependent Variable 1")+
  labs(x = "",
       y = "",
       color = "two-state scenario:") +
  scale_color_manual(values = cols,
                     labels=c(expression(1^"st"~"state" ),expression(2^"nd"~"state" )))+
  scale_y_continuous(limits=c(0,0.055),breaks=c(0,0.02,0.04))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour=rep("#696969",3)),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(7, 6, 6, 6),legend.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

#plot 1b
datay<-dnorm(seq(0,120,0.05), mean= 60,sd = sqrt((120+0.25*120)))
data<-data.frame(x=seq(0,120,0.05),y1=datay,y2=dnorm(seq(0,120,0.05), mean= 30,sd = sqrt((60+0.25*60))))
plot1b<-ggplot(data=data,aes(x=x,y=y1))+geom_line(size=1,aes(color="col1"))+geom_line(size=1,aes(y=y2,color="col2"))+
  theme_minimal()+xlim(c(0,120))+ggtitle("Dependent Variable 2")+#+ggtitle("Dwell time distributions")
  labs(x = "",
       y = "",
       color = "two-state scenario:") +
  scale_color_manual(values = cols,
                     labels=c(expression(1^"st"~"state" ),expression(2^"nd"~"state" )))+
  scale_y_continuous(limits=c(0,0.055),breaks=c(0,0.02,0.04))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(7, 6, 6, 6),legend.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

#plot 2a
cols <- c("col1" = colors$blue, "col2" = colors$orange, "col3" = colors$red, "col4"=colors$teal)
datay<-dnorm(seq(0,120,0.05), mean= 10,sd = sqrt(20))
data<-data.frame(x=seq(0,120,0.05),y1=datay,y2=dnorm(seq(0,120,0.05), mean= 30,sd = sqrt(60)),y3=dnorm(seq(0,120,0.05), mean= 60,sd = sqrt(120)))
plot2a<-ggplot(data=data,aes(x=x,y=y1))+geom_line(size=1,aes(color="col1"))+geom_line(size=1,aes(y=y2,color="col2"))+geom_line(size=1,aes(y=y3,color="col3"))+
  scale_y_continuous(limits=c(0,0.09),breaks=c(0,0.04,0.08))+ggtitle("Dependent Variable 1")+
  theme_minimal()+xlim(c(0,120))+
  labs(x = "",
       y = "",
       color = "three-state scenario:") +
  scale_color_manual(values = cols,
                     labels=c(expression(1^"st"~"state" ),expression(2^"nd"~"state" ),expression(3^"rd"~"state" )))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(colour=rep("#696969",3)),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(7, 6, 6, 6),legend.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

#plot 2b
datay<-dnorm(seq(0,120,0.05), mean= 60,sd = sqrt((120+0.25*120)))
data<-data.frame(x=seq(0,120,0.05),y1=datay,y2=dnorm(seq(0,120,0.05), mean= 30,sd = sqrt((60+0.25*60))),y3=dnorm(seq(0,120,0.05), mean= 10,sd = sqrt((20+0.25*20))))
plot2b<-ggplot(data=data,aes(x=x,y=y1))+geom_line(size=1,aes(color="col1"))+geom_line(size=1,aes(y=y2,color="col2"))+geom_line(size=1,aes(y=y3,color="col3"))+
  scale_y_continuous(limits=c(0,0.09),breaks=c(0,0.04,0.08))+
  theme_minimal()+xlim(c(0,120))+ggtitle("Dependent Variable 2")+#+ggtitle("Dwell time distributions")
  labs(x = "",
       y = "",
       color = "three-state scenario:") +
  scale_color_manual(values = cols,
                     labels=c(expression(1^"st"~"state" ),expression(2^"nd"~"state" ),expression(3^"rd"~"state" )))+
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(7, 6, 6, 6),legend.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

#plot 3a
cols <- c("col1" = colors$blue, "col2" = colors$orange, "col3" = colors$red, "col4"=colors$teal)
datay<-dnorm(seq(0,120,0.05), mean= 10,sd = sqrt(20))
data<-data.frame(x=seq(0,120,0.05),y1=datay,y2=dnorm(seq(0,120,0.05), mean= 30,sd = sqrt(60)),y3=dnorm(seq(0,120,0.05), mean= 60,sd = sqrt(120)),y4=dnorm(seq(0,120,0.05), mean= 98,sd = sqrt(150)))
plot3a<-ggplot(data=data,aes(x=x,y=y1))+geom_line(size=1,aes(color="col1"))+geom_line(size=1,aes(y=y2,color="col2"))+geom_line(size=1,aes(y=y3,color="col3"))+geom_line(size=1,aes(y=y4,color="col4"))+
  scale_y_continuous(limits=c(0,0.09),breaks=c(0,0.04,0.08))+
  theme_minimal()+xlim(c(0,120))+#+ggtitle("Dwell time distributions")
  labs(x = "",
       y = "",
       color = "four-state scenario:") +
  scale_color_manual(values = cols,
                     labels=c(expression(1^"st"~"state" ),expression(2^"nd"~"state" ),expression(3^"rd"~"state" ),expression(4^"th"~"state" )))+
  theme(axis.text.x = element_text(colour="#696969"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour=rep("#696969",3)),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(7, 6, 6, 6),legend.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

#plot 3b
datay<-dnorm(seq(0,120,0.05), mean= 98,sd = sqrt((150+0.25*150)))
data<-data.frame(x=seq(0,120,0.05),y1=datay,y2=dnorm(seq(0,120,0.05), mean= 60,sd = sqrt((120+0.25*120))), y3=dnorm(seq(0,120,0.05), mean= 30,sd = sqrt((60+0.25*60))),y4=dnorm(seq(0,120,0.05), mean= 10,sd = sqrt((20+0.25*20))))
plot3b<-ggplot(data=data,aes(x=x,y=y1))+geom_line(size=1,aes(color="col1"))+geom_line(size=1,aes(y=y2,color="col2"))+geom_line(size=1,aes(y=y3,color="col3"))+geom_line(size=1,aes(y=y4,color="col4"))+
  scale_y_continuous(breaks=c(0,0.04,0.08),limits=c(0,0.09))+
  theme_minimal()+xlim(c(0,120))+#+ggtitle("Dwell time distributions")
  labs(x = "",
       y = "",
       color = "four-state scenario:") +
  scale_color_manual(values = cols,
                     labels=c(expression(1^"st"~"state" ),expression(2^"nd"~"state" ),expression(3^"rd"~"state" ),expression(4^"th"~"state" )))+
  theme(axis.text.x = element_text(colour="#696969"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(7, 6, 6, 6),legend.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE))

pdf("emiss_dist_plot.pdf",width = 12,height = 8)
(plot2a+plot2b)/(plot3a+plot3b)+ plot_layout(guides = "collect") & theme(legend.title =element_text(size = 13, colour = "black") , legend.text = element_text(size = 13, colour = "black"),legend.position = 'bottom')
dev.off()
getwd()


















