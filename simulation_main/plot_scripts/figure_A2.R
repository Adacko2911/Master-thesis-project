library(tidyverse)
library(ggplot2)
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


cols <- c("col1" = colors$blue, "col2" = colors$orange, "col3" = colors$red, "col4"=colors$teal)
datay<-dlnorm(seq(0,150,0.05), meanlog = log(3),sdlog = sqrt(log(1.10093)))
data<-data.frame(x=seq(0,150,0.05),y=datay)
plot<-ggplot(data=data,aes(x=x,y=y))+geom_line(size=1,aes(color="col1"))

datay<-dlnorm(seq(0,150,0.05), meanlog = log(6),sdlog = sqrt(log(1.10093)))
data<-data.frame(x=seq(0,150,0.05),y=datay)
plot<-plot+geom_line(data=data,aes(x=x,y=y,color="col2"),size=1.2)

datay<-dlnorm(seq(0,150,0.05), meanlog = log(19.5),sdlog = sqrt(log(1.10093)))
data<-data.frame(x=seq(0,150,0.05),y=datay)
plot<-plot+geom_line(data=data,aes(x=x,y=y,color="col3"),size=1.2)
datay<-dlnorm(seq(0,150,0.05), meanlog = log(99.5),sdlog = sqrt(log(1.05902)))
data<-data.frame(x=seq(0,150,0.05),y=datay)
plot<-plot+geom_line(data=data,aes(x=x,y=y,color="col4"),size=1.2)
plot<-plot+theme_minimal()+xlim(c(0,125))+#+ggtitle("Dwell time distributions")
  labs(x = "Dwell time",
       y = "",
       color = "Mean Dwell Time") +
  scale_color_manual(values = cols, labels=c(expression("d"=="1.4"),expression("d"=="3.5"),expression("d"=="19.5"),expression("d"=="99.5"))) +
  scale_y_continuous(breaks=c(0,0.2,0.4))+
  theme(axis.text.x = element_text(colour="#696969"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour=rep("#696969",3)),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(7, 6, 6, 6),legend.title =element_text(size = 13, colour = "black") , legend.text = element_text(size = 13, colour = "black"),legend.text.align = 0
  )
pdf("lognorm_dist_plot.pdf",width = 7,height = 4)
plot
dev.off()
