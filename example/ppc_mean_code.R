### MEDHMM
library(ggplot2)
library(R.utils)
library(tidyverse)
library(patchwork)
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
#This part is suppose perform predictive posterior checks checking the mean of the values
#Let us check the variance as well (this is to be done)
dep=12
# Specify the name of a file created during simulating ppc
ppc <- readRDS("./example/outputs/ppc.rds")
group_level_ppc<-matrix(nrow=500,ncol=12)
ppc_out_list1<-list()
for(i in 1:500){
  ppc_out_list1[[i]]<-ppc[[i]][["observations"]] %>%  as.data.frame()
  names(ppc_out_list1[[i]]) <- c("patient_id", c("down","dread","worry",
                                                 "inadequate", "tired", "content",
                                                 "agitated", "irritated", "focused",
                                                 "extremely well", "full of ideas", "thought racing"))
  ppc_out_list1[[i]]<-ppc_out_list1[[i]]%>%
    group_by(patient_id) %>%
    mutate(occasion = row_number()) %>%
    gather(variable, value, -patient_id, -occasion) %>%
    mutate(rep = i) %>% ungroup()

  new<-ppc[[i]][["observations"]] %>% as.data.frame() %>% apply(.,2,mean)%>% as.matrix(.,nrow=1,by_row=T) %>% t()
  #group level pcc
  group_level_ppc[i,1:12]<-new[,-1]
}

colnames(group_level_ppc)<-c(paste0("dep_",1:dep))


#group_level_ppc
#subject_level_ppc$subject<-as.factor(subject_level_ppc$subject)
equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    if(seq[2]-seq[1] < 10^(-r)) seq else round(seq, r)
  }
}
group_statistics <- readRDS("./example/outputs/group_statistics.rds")
group_statistics <- group_statistics %>%
  mutate(variable = factor(variable, levels = c("down","dread nrest of day","worry",
                                                "inadequate", "tired", "content",
                                                "agitated", "irritated", "switch and focus",
                                                "extremely well", "full of ideas", "thoughts are racing"),
                           labels = c("down","dread","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "focused",
                                      "extremely well", "full of ideas", "thought racing")))


colnames(group_statistics)[1]<-"dep"
group_statistics<-group_statistics[,1:2]
#here the posterior predictive checks for the group level emission distributions are checked
group_level_ppc_long<-group_level_ppc %>%as.data.frame() %>%  pivot_longer(cols=dep_1:dep_12,names_to = "dep", values_to = "mean")    %>% mutate_if(is.character,as.factor)  %>% mutate(dep=fct_relevel(dep, c("dep_10","dep_11","dep_12"), after = 9))%>% arrange(dep)
#reorder and add the true values (sample values)
#extract the true values from observations
levels(group_level_ppc_long$dep)
group_level_ppc_long<-group_level_ppc_long %>% mutate(dep=factor(dep, levels = c("dep_1","dep_2", "dep_3","dep_4","dep_5","dep_6","dep_7","dep_8","dep_9","dep_10","dep_11","dep_12"), labels =  c("down","dread","worry",
                                                                                                                                                                                                   "inadequate", "tired", "content",
                                                                                                                                                                                                   "agitated", "irritated", "focused",
                                                                                                                                                                                                   "extremely well", "full of ideas", "thought racing")))
my_grey_text="#696969"
#we visualized the data
p1<-group_level_ppc_long%>%  right_join(.,group_statistics,by="dep") %>% ggplot(aes(x=mean))+geom_density(fill = my_grey_text, colour = "black", alpha = 0.1,size=0.2)+theme_minimal()+geom_vline(data=group_statistics, aes(xintercept = mean_value),linewidth = 0.5,colour=colors$red,linetype="dashed")+
  facet_wrap(~dep,nrow=2,scales = "free_y") +theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y = element_text(size = 9),axis.title.x = element_text(size = 9),axis.text.x = element_text(size = 8,colour="#696969"))+ylab("Density")+xlab("Value")
rm(ppc)

#EDHMM
ppc_data_mean_ss<-do.call(rbind,ppc_out_list1) %>%
  group_by(rep, patient_id, variable) %>%
  summarise(mean = mean(value)) %>%
  group_by(rep, variable) %>%
  arrange(desc(mean)) %>%
  mutate(id = row_number()) %>%
  group_by(variable, rep) %>%
  arrange(id) %>%mutate(subject=id) %>% ungroup()%>% dplyr::select(c(subject,variable,mean))
rm(ppc_out_list1)
gc()
patient_statistics <- readRDS("~/GitHub/Master-thesis-project/Simulation Study/empirical examples script/patient_statistics.rds")
patient_statistics <- patient_statistics %>%
  mutate(variable = factor(variable, levels = c("down","dread nrest of day","worry",
                                                "inadequate", "tired", "content",
                                                "agitated", "irritated", "switch and focus",
                                                "extremely well", "full of ideas", "thoughts are racing"),
                           labels = c("down","dread","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "focused",
                                      "extremely well", "full of ideas", "thought racing")))
colnames(patient_statistics)[1]<-"subject"
patient_statistics<-patient_statistics[,c(1:3)] %>% group_by(variable)%>%
  arrange(desc(mean_value)) %>%mutate(id = row_number()) %>% mutate(subject=id) %>% ungroup()%>% dplyr::select(c(subject,variable,mean_value))

neee<-ppc_data_mean_ss%>% right_join(.,patient_statistics,by=c("subject","variable"))
neee$subject<-as.factor(neee$subject)

p1a<-neee %>%ggplot(aes(x=mean,y=subject))+geom_boxplot(outlier.shape = "",size=0.1)+geom_point(aes(x=mean_value,y=subject),color=colors$red,shape=4,size =0.4,fill=colors$red)+facet_wrap(~variable,nrow=2,scales = "fixed")+
  theme_minimal()+ ylab("Patient no.")+xlab("Value")+scale_x_continuous(breaks=c(0,20,50,80))+xlim(c(0,80))+
  theme(panel.grid.minor = element_blank(),axis.title.y = element_text(size = 9),axis.title.x = element_text(size = 9),axis.text.x = element_text(size = 8,colour="#696969"))

pdf("./example/plots/ppc_mean_emiss_hmm.pdf",width =7 ,height =10 )
p1+p1a+plot_layout(ncol=1,heights=c(1,2))+plot_annotation(tag_levels = c('A', '1'))
dev.off()


#### MHMM
library(ggplot2)
library(R.utils)
library(tidyverse)
library(patchwork)
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
#This part is suppose perform predictive posterior checks checking the mean of the values
#Let us check the variance as well (this is to be done)
dep=12
# Specify the name of a file created during simulating ppc
ppc <- readRDS("./example/outputs/ppc_mhmm.rds" )
group_level_ppc<-matrix(nrow=500,ncol=12)
ppc_out_list1<-list()
for(i in 1:500){
  ppc_out_list1[[i]]<-ppc[[i]][["obs"]] %>%  as.data.frame()
  names(ppc_out_list1[[i]]) <- c("patient_id", c("down","dread","worry",
                                                 "inadequate", "tired", "content",
                                                 "agitated", "irritated", "focused",
                                                 "extremely well", "full of ideas", "thought racing"))
  ppc_out_list1[[i]]<-ppc_out_list1[[i]]%>%
    group_by(patient_id) %>%
    mutate(occasion = row_number()) %>%
    gather(variable, value, -patient_id, -occasion) %>%
    mutate(rep = i) %>% ungroup()

  new<-ppc[[i]][["obs"]] %>% as.data.frame() %>% apply(.,2,mean)%>% as.matrix(.,nrow=1,by_row=T) %>% t()
  #group level pcc
  group_level_ppc[i,1:12]<-new[,-1]
}

colnames(group_level_ppc)<-c(paste0("dep_",1:dep))


#group_level_ppc
#subject_level_ppc$subject<-as.factor(subject_level_ppc$subject)
equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    if(seq[2]-seq[1] < 10^(-r)) seq else round(seq, r)
  }
}

group_statistics <- readRDS("./example/outputs/group_statistics.rds")
group_statistics <- group_statistics %>%
  mutate(variable = factor(variable, levels = c("down","dread nrest of day","worry",
                                                "inadequate", "tired", "content",
                                                "agitated", "irritated", "switch and focus",
                                                "extremely well", "full of ideas", "thoughts are racing"),
                           labels = c("down","dread","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "focused",
                                      "extremely well", "full of ideas", "thought racing")))


colnames(group_statistics)[1]<-"dep"
group_statistics<-group_statistics[,1:2]
#here the posterior predictive checks for the group level emission distributions are checked
group_level_ppc_long<-group_level_ppc %>%as.data.frame() %>%  pivot_longer(cols=dep_1:dep_12,names_to = "dep", values_to = "mean")    %>% mutate_if(is.character,as.factor)  %>% mutate(dep=fct_relevel(dep, c("dep_10","dep_11","dep_12"), after = 9))%>% arrange(dep)
#reorder and add the true values (sample values)
#extract the true values from obs
levels(group_level_ppc_long$dep)
group_level_ppc_long<-group_level_ppc_long %>% mutate(dep=factor(dep, levels = c("dep_1","dep_2", "dep_3","dep_4","dep_5","dep_6","dep_7","dep_8","dep_9","dep_10","dep_11","dep_12"), labels =  c("down","dread","worry",
                                                                                                                                                                                                   "inadequate", "tired", "content",
                                                                                                                                                                                                   "agitated", "irritated", "focused",
                                                                                                                                                                                                   "extremely well", "full of ideas", "thought racing")))
my_grey_text="#696969"
#we visualized the data
p1<-group_level_ppc_long%>%  right_join(.,group_statistics,by="dep") %>% ggplot(aes(x=mean))+geom_density(fill = my_grey_text, colour = "black", alpha = 0.1,size=0.2)+theme_minimal()+geom_vline(data=group_statistics, aes(xintercept = mean_value),linewidth = 0.5,colour=colors$red,linetype="dashed")+
  facet_wrap(~dep,nrow=2,scales = "free_y") +theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y = element_text(size = 9),axis.title.x = element_text(size = 9),axis.text.x = element_text(size = 8,colour="#696969"))+ylab("Density")+xlab("Value")
rm(ppc)

#EDHMM
ppc_data_mean_ss<-do.call(rbind,ppc_out_list1) %>%
  group_by(rep, patient_id, variable) %>%
  summarise(mean = mean(value)) %>%
  group_by(rep, variable) %>%
  arrange(desc(mean)) %>%
  mutate(id = row_number()) %>%
  group_by(variable, rep) %>%
  arrange(id) %>%mutate(subject=id) %>% ungroup()%>% dplyr::select(c(subject,variable,mean))
rm(ppc_out_list1)
gc()
patient_statistics <- readRDS("./example/outputs/patient_statistics.rds")
patient_statistics <- patient_statistics %>%
  mutate(variable = factor(variable, levels = c("down","dread nrest of day","worry",
                                                "inadequate", "tired", "content",
                                                "agitated", "irritated", "switch and focus",
                                                "extremely well", "full of ideas", "thoughts are racing"),
                           labels = c("down","dread","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "focused",
                                      "extremely well", "full of ideas", "thought racing")))
colnames(patient_statistics)[1]<-"subject"
patient_statistics<-patient_statistics[,c(1:3)] %>% group_by(variable)%>%
  arrange(desc(mean_value)) %>%mutate(id = row_number()) %>% mutate(subject=id) %>% ungroup()%>% dplyr::select(c(subject,variable,mean_value))

neee<-ppc_data_mean_ss%>% right_join(.,patient_statistics,by=c("subject","variable"))
neee$subject<-as.factor(neee$subject)

p1a<-neee %>%ggplot(aes(x=mean,y=subject))+geom_boxplot(outlier.shape = "",size=0.1)+geom_point(aes(x=mean_value,y=subject),color=colors$red,shape=4,size =0.4,fill=colors$red)+facet_wrap(~variable,nrow=2,scales = "fixed")+
  theme_minimal()+ ylab("Patient no.")+xlab("Value")+scale_x_continuous(breaks=c(0,20,50,80))+xlim(c(0,80))+
  theme(panel.grid.minor = element_blank(),axis.title.y = element_text(size = 9),axis.title.x = element_text(size = 9),axis.text.x = element_text(size = 8,colour="#696969"))

pdf("./example/plots/ppc_mean_hmm.pdf",width =7 ,height =10 )
p1+p1a+plot_layout(ncol=1,heights=c(1,2))+plot_annotation(tag_levels = c('A', '1'))
dev.off()

gc()

sessionInfo()
