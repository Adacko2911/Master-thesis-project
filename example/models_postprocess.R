#libraries loading
library(medHMM)
library(mHMMbayes)
library(devtools)
library(ggplot2)
library(R.utils)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reshape)
library(ggridges)
library(ggshadow)

#source functions (please make sure all functions are loaded) in case some errors regarding
#unknown functions appear please run the script manually by opening its location
source("./example/utility_functions.R")
#load data
out_medHMM_m4_12dv_it2000_c1 <- readRDS("./example/outputs/out_medHMM_m4_12dv_it2000_c1.rds")
decoding_medhmm<-local_decoding(out_medHMM_m4_12dv_it2000_c1)
out_mHMM_m4_12dv_it2000_c1 <- readRDS("./example/outputs/out_mHMM_m4_12dv_it2000_c1.rds")
decoding_mhmm<-local_decoding(out_mHMM_m4_12dv_it2000_c1)
J=2000
B=1000
m=4
out_medHMM_m4_12dv_it2000_c2 <- readRDS("./example/outputs/out_medHMM_m4_12dv_it2000_c3.rds")
out_mHMM_m4_12dv_it2000_c2 <- readRDS("./example/outputs/out_mHMM_m4_12dv_it2000_c2.rds")

#--Inspect Gelman Rubin statistics of group-level statistics
ge_rub_two_models_emp(out_medHMM_m4_12dv_it2000_c1,out_medHMM_m4_12dv_it2000_c2,m = m,burn_in = B, J=J)
ge_rub_two_models_emp(out_mHMM_m4_12dv_it2000_c1,out_mHMM_m4_12dv_it2000_c2,m = m,burn_in = B, J=J)


### TRACE PLOTS ###########################################################################################################################
test_medhmm_chain1=out_medHMM_m4_12dv_it2000_c1
test_medhmm_chain2=out_medHMM_m4_12dv_it2000_c2
test_mhmm_chain1=out_mHMM_m4_12dv_it2000_c1
test_mhmm_chain2=out_mHMM_m4_12dv_it2000_c2

pdf("./example/plots/tarce_plots.pdf")
#============emiss_mu=====================================================================================
traceplots_emiss_mu(test_medhmm_chain1,test_medhmm_chain2,colors = c("#485296","#8F3481"), titles=rep("MEDHMM emiss mu",4),dep = 4)
traceplots_emiss_mu(test_mhmm_chain1,test_mhmm_chain2,colors = c("#00BCB4","#FF4747"), titles=rep("MHMM emiss mu",4),dep = 4)

#===========emiss_var=====================================================================================
traceplots_emiss_var(test_medhmm_chain1,test_medhmm_chain2,colors = c("#485296","#8F3481"), titles=rep("MEDHMM emiss var",4),dep = 4)
traceplots_emiss_var(test_mhmm_chain1,test_mhmm_chain2,colors = c("#00BCB4","#FF4747"), titles=rep("MHMM emiss var",4),dep = 4)

#=========gamma estimates=================================================================================
traceplots_gamma(test_medhmm_chain1,test_medhmm_chain2,colors = c("#485296","#8F3481"), titles=rep("MEDHMM gamma",16))
traceplots_gamma(test_mhmm_chain1,test_mhmm_chain2,colors = c("#00BCB4","#FF4747"), titles=rep("MHMM gamma",16))

#=========dwell mu and dwell var for mhmm=================================================================
traceplots_dwell_mu(test_medhmm_chain1,test_medhmm_chain2,c(1,4),colors = c("#485296","#8F3481"), titles=rep("MEDHMM dwell mu",4))
traceplots_dwell_var(test_medhmm_chain1,test_medhmm_chain2,ylim1=c(0,2),colors = c("#485296","#8F3481"), titles=rep("MEDHMM dwell var",4))
#=========================================================================================================

dev.off()

#### Get the group-level estimates to inspect ##########################################################################################################################
MAP_medhmm_out<-MAP_medhmm(out_medHMM_m4_12dv_it2000_c1,J=2000,B=1000,m=4)
dep_var_names<-out_medHMM_m4_12dv_it2000_c1$input$dep_labels
emiss_mu<-do.call(rbind,MAP_medhmm_out$emiss_mu_bar)%>% as.data.frame() %>% mutate(item=factor(dep_var_names, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                                                         "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                                                         "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                                                         "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                                                               labels = c("down","dreadnest of day","worry",
                                                                                                          "inadequate", "tired", "content",
                                                                                                          "agitated", "irritated", "switch and focus",
                                                                                                          "extremely well", "full of ideas", "thoughts are racing"))) %>% pivot_longer(-item,
                                                                                                                                                                                       names_to = c("Var", ".value"),
                                                                                                                                                                                       names_sep="[.]" )
emiss_var<-do.call(rbind,MAP_medhmm_out$emiss_var_bar)%>% as.data.frame() %>% mutate(item=factor(dep_var_names, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                                                           "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                                                           "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                                                           "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                                                                 labels = c("down","dreadnest of day","worry",
                                                                                                            "inadequate", "tired", "content",
                                                                                                            "agitated", "irritated", "switch and focus",
                                                                                                            "extremely well", "full of ideas", "thoughts are racing"))) %>% pivot_longer(-item,
                                                                                                                                                                                         names_to = c("Var", ".value"),
                                                                                                                                                                                         names_sep="[.]" )
emiss_varmu<-do.call(rbind,MAP_medhmm_out$emiss_varmu_bar)%>% as.data.frame() %>% mutate(item=factor(dep_var_names, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                                                               "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                                                               "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                                                               "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                                                                     labels = c("down","dreadnest of day","worry",
                                                                                                                "inadequate", "tired", "content",
                                                                                                                "agitated", "irritated", "switch and focus",
                                                                                                                "extremely well", "full of ideas", "thoughts are racing"))) %>% pivot_longer(-item,
                                                                                                                                                                                             names_to = c("Var", ".value"),
                                                                                                                                                                                             names_sep="[.]" )
data_to_display1<-rbind(emiss_mu,emiss_var,emiss_varmu)%>% mutate(state=as.integer(str_split_i(Var,pattern = "_",-1)),estimate=as.factor(str_split_i(Var,pattern = "_",1)), .after=2)  %>% dplyr::select(-Var) %>% mutate(model=rep("medhmm",nrow(.)),.before = 1)


MAP_mhmm<-MAP_mHMM(out_mHMM_m4_12dv_it2000_c1,J=2000,B=1000,m=4)
emiss_mu<-do.call(rbind,MAP_mhmm$emiss_mu_bar)%>% as.data.frame() %>% mutate(item=factor(dep_var_names, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                                                   "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                                                   "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                                                   "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                                                         labels = c("down","dreadnest of day","worry",
                                                                                                    "inadequate", "tired", "content",
                                                                                                    "agitated", "irritated", "switch and focus",
                                                                                                    "extremely well", "full of ideas", "thoughts are racing"))) %>% pivot_longer(-item,
                                                                                                                                                                                 names_to = c("Var", ".value"),
                                                                                                                                                                                 names_sep="[.]" )
emiss_var<-do.call(rbind,MAP_mhmm$emiss_var_bar)%>% as.data.frame() %>% mutate(item=factor(dep_var_names, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                                                     "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                                                     "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                                                     "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                                                           labels = c("down","dreadnest of day","worry",
                                                                                                      "inadequate", "tired", "content",
                                                                                                      "agitated", "irritated", "switch and focus",
                                                                                                      "extremely well", "full of ideas", "thoughts are racing"))) %>%pivot_longer(-item,
                                                                                                                                                                                  names_to = c("Var", ".value"),
                                                                                                                                                                                  names_sep="[.]" )
emiss_varmu<-do.call(rbind,MAP_mhmm$emiss_varmu_bar)%>% as.data.frame() %>% mutate(item=factor(dep_var_names, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                                                                                         "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                                                                                         "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                                                                                         "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                                                                                               labels = c("down","dreadnest of day","worry",
                                                                                                          "inadequate", "tired", "content",
                                                                                                          "agitated", "irritated", "switch and focus",
                                                                                                          "extremely well", "full of ideas", "thoughts are racing")))  %>% pivot_longer(-item,
                                                                                                                                                                                        names_to = c("Var", ".value"),
                                                                                                                                                                                        names_sep="[.]" )
data_to_display2<-rbind(emiss_mu,emiss_var,emiss_varmu)%>% mutate(state=as.integer(str_split_i(Var,pattern = "_",-1)),estimate=as.factor(str_split_i(Var,pattern = "_",1)), .after=2) %>% dplyr::select(-Var) %>% mutate(model=rep("mhmm",nrow(.)),.before = 1)

final_to_display_dep<-rbind(data_to_display2,data_to_display1)
final_to_display_dep$state<-as.factor(final_to_display_dep$state)

write.csv(final_to_display_dep,"./example/outputs/new2.csv")



# Figure 9 PREPARATION

#prepare estimates for the subject specific means
out_ss_mu_medhmm<-pbmcapply::pbmclapply(1:length(out_medHMM_m4_12dv_it2000_c1$PD_subj), function(s){
  data_out<-out_medHMM_m4_12dv_it2000_c1$PD_subj[[s]][(B+1):J,1:(4*12)] %>%
    as.data.frame() %>% mutate(subj=rep(s, J-B))%>%  pivot_longer(-subj ,  names_to = c("item", ".value"), names_sep="_mu_" ) %>%
    group_by(item) %>% pivot_longer(-c(item,subj),  names_to = c("state"),values_to = c("value")) %>% group_by(item,state) %>%
    summarise(mean_value = mean(value, na.rm = TRUE),median_value= median(value, na.rm = TRUE),sd_value= sd(value, na.rm = TRUE), ci_025=quantile(value, probs=c(0.025)),ci_975=quantile(value, probs=c(0.975))) %>%
    ungroup() %>% mutate(subj=rep(s, nrow(.))) %>% mutate(item=factor(item, levels = paste0("dep",1:12),labels = c("down","dreadnest of day","worry",
                                                                                                                   "inadequate", "tired", "content",
                                                                                                                   "agitated", "irritated", "switch and focus",
                                                                                                                   "extremely well", "full of ideas", "thoughts racing"))) %>%
    mutate(model=rep("medhmm",nrow(.)))

  return(data_out)
},mc.cores = 1)

out_ss_mu_medhmm2<-do.call(rbind, out_ss_mu_medhmm)
out_ss_mu_medhmm2$state<-factor(out_ss_mu_medhmm2$state, levels = c("S1","S2","S3","S4"),labels = c("euthymic","manic","mixed","depressive"))
write.csv(out_ss_mu_medhmm2,"C:/Users/ola/OneDrive - Universiteit Utrecht/Documents/GitHub/Master-thesis-project/Simulation Study/new.csv")
# bipolar_data_items has to be substituted with the data with first column subject, 2:13 the items we chose
medhmm_simulated_data_bip_sum<-data.frame(subject=out_ss_mu_medhmm2$subj,item=out_ss_mu_medhmm2$item,median=out_ss_mu_medhmm2$median_value,states=out_ss_mu_medhmm2$state)
medhmm_simulated_data_bip_sum$states<-as.factor(medhmm_simulated_data_bip_sum$states)
medhmm_simulated_data_bip_sum$subject<-as.factor(medhmm_simulated_data_bip_sum$subject)


###### FIGURE 9 ###########################################################################################################################
pdf("./example/plots/ss_emiss.pdf",width = 7,height =8)
#we could add to the data frame variable representing the model type and add shape as an extra variable so that we van compare the results
p1<-medhmm_simulated_data_bip_sum%>% filter(states=="euthymic") %>% ggplot(aes(x=factor(item, level=c("agitated","irritated", "switch and focus","extremely well","full of ideas","thoughts racing",  "down","dreadnest of day", "worry","inadequate","tired","content")),y=median,color=item))+ geom_jitter(pch="o",size=2,position = position_jitter(width = 0.2, height = 0.1),show.legend = FALSE)+theme_minimal()+geom_hline(yintercept = 50,linetype = "dashed")+scale_y_continuous(breaks = c(0,25,50,75),limits = c(0,80))+ggtitle("Euthymic state")+scale_color_manual(values = as.character(gamam_colors))+
  theme(axis.text.x = element_blank(),axis.title = element_text(color = my_grey_text,size = 8),plot.title = element_text( size = 10,face = "bold"))
p2<-medhmm_simulated_data_bip_sum%>% filter(states=="manic") %>%  ggplot(aes(x=factor(item, level=c("agitated","irritated", "switch and focus","extremely well","full of ideas","thoughts racing",  "down","dreadnest of day", "worry","inadequate","tired","content")),y=median,color=item))+ geom_jitter(pch="o",size=2,position = position_jitter(width = 0.2, height = 0.1),show.legend = FALSE)+theme_minimal()+geom_hline(yintercept = 50,linetype = "dashed")+scale_y_continuous(breaks = c(0,25,50,75),limits = c(0,80))+ggtitle("Manic state")+scale_color_manual(values = as.character(gamam_colors))+
  theme(axis.text.x = element_blank(),axis.title = element_text(color = my_grey_text,size = 8),plot.title = element_text( size = 10,face = "bold"))
p3<-medhmm_simulated_data_bip_sum %>% filter(states=="mixed") %>% ggplot(aes(x=factor(item, level=c("agitated","irritated", "switch and focus","extremely well","full of ideas","thoughts racing",  "down","dreadnest of day", "worry","inadequate","tired","content")),y=median,color=item))+ geom_jitter(pch="o",size=2,position = position_jitter(width = 0.2, height = 0.1),show.legend = FALSE)+theme_minimal()+geom_hline(yintercept = 50,linetype = "dashed")+scale_y_continuous(breaks = c(0,25,50,75),limits = c(0,80))+ggtitle("Mixed state")+scale_color_manual(values = as.character(gamam_colors))+
  theme(axis.text.x = element_blank(),axis.title = element_text(color = my_grey_text,size = 8),plot.title = element_text( size = 10,face = "bold"))

p4<-medhmm_simulated_data_bip_sum%>% filter(states=="depressive") %>% ggplot(aes(x=factor(item, level=c("agitated","irritated", "switch and focus","extremely well","full of ideas","thoughts racing",  "down","dreadnest of day", "worry","inadequate","tired","content")),y=median,color=item))+ geom_jitter(pch="o",size=2,position = position_jitter(width = 0.2, height = 0.1),show.legend = FALSE)+theme_minimal()+geom_hline(yintercept = 50,linetype = "dashed")+scale_y_continuous(breaks = c(0,25,50,75),limits = c(0,80))+ggtitle("Depressive state")+scale_color_manual(values = as.character(gamam_colors))+
  theme(axis.text.x = element_text(angle = 35, hjust=1,size=9),axis.title = element_text(color = my_grey_text,size = 9),plot.title = element_text( size = 10,face = "bold"))



ylab <- "EMA estimated mean"
xlab<- ""
p1$labels$y <- p2$labels$y<- p3$labels$y  <- p4$labels$y <- " "
p1$labels$x <- p2$labels$x<- p3$labels$x  <- p4$labels$x <- " "
(p1/p2/p3/p4)+ plot_annotation(
  title = ''
)
grid::grid.draw(grid::textGrob(ylab, x = 0.02, rot = 90,gp = gpar(col = "black", fontsize = 11)))

dev.off()

###########################################################################################################################################


## Figure C3 PREPARATION
#extract median transition matrices
#MEDHMM
data_trans_mat2<-matrix(as.numeric(MAP_medhmm_out$gamma_prob_bar[,1:(m*m)]),nrow=m,byrow = T)
row.names(data_trans_mat2)<-paste0("State",1:m)
colnames(data_trans_mat2)<-paste0("State",1:m)
g3<-melt(data_trans_mat2)
names(g3)
g3$X1 <- factor(g3$X1,
                levels = paste0("State",m:1), labels =c("depressive","mixed","manic","euthymic") )
g3$X2 <- factor(g3$X2,
                levels = paste0("State",1:m), labels =c("euthymic","manic","mixed","depressive") )
g3$value<-round(g3$value,2)
tr_plt_med<-ggplot(g3, aes(X2,X1)) + geom_tile(aes(fill= value))+coord_fixed()+geom_text(aes(label = value), color = "white", size = 5) +scale_fill_distiller(direction = 1)+guides(fill = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                                                           barheight = 5,title = "Transition probability"))+
  xlab("To state")+ylab("From state")+theme_minimal()+theme(legend.position = "none",plot.title = element_text(),plot.subtitle = element_text(size=9,hjust=0.05),panel.grid.minor =element_blank(),axis.text.x = element_text(angle = 35, hjust=1,size=9),
                                                            axis.text.y = element_text(angle = 35,size=10))+ggtitle("Group-level TPM")

#MHMM
m=4
data_trans_mat1<-matrix(as.numeric(MAP_mhmm$gamma_prob_bar[,17:(2*m*m)]),nrow=m,byrow = T)
row.names(data_trans_mat1)<-paste0("State",1:m)
colnames(data_trans_mat1)<-paste0("State",1:m)
g2<-melt(data_trans_mat1)
g2$X1 <- factor(g2$X1,
                levels = paste0("State",m:1), labels =c("depressive","mixed","manic","euthymic") )
g2$X2 <- factor(g2$X2,
                levels = paste0("State",1:m), labels =c("euthymic","manic","mixed","depressive") )
g2$value<-round(g2$value,2)


## Figure C3 ##############################################################################################################################
pdf("./example/plots/group_trandition_matricesvv2.pdf")
ggplot(g2, aes(X2,X1)) + geom_tile(aes(fill= value))+coord_fixed()+geom_text(aes(label = value), color = "white", size = 5)+
  scale_fill_distiller(direction = 0)+guides(fill = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                                               barheight = 5,title = "Transition probability"))+xlab("")+ylab("")+theme_minimal()+theme(legend.position = "none",plot.title = element_text(size=12),plot.subtitle = element_text(size=9,hjust=0.05),panel.grid.minor =element_blank(),axis.text.x = element_text(angle = 35, hjust=1,size=9))+ggtitle("Group-level TPM: MHMM")
dev.off()
###########################################################################################################################################

### Comparison plots for the two models
pdf("./example/plots/group_trandition_matrices.pdf")
tr_plt_med+tr_plt_m
dev.off()

part1<-MAP_medhmm_out$gamma_prob_bar%>%as.data.frame() %>%mutate(model="medhmm") %>%  pivot_longer(-model,
                                                                                                   names_to = c("Var", ".value"),
                                                                                                   names_sep="[.]" )
part2<-MAP_mhmm$gamma_prob_bar_no_diag%>%as.data.frame() %>%mutate(model="mhmm") %>%  pivot_longer(-model,
                                                                                                   names_to = c("Var", ".value"),
                                                                                                   names_sep="[.]" )
saveRDS(rbind(part1,part2),"./example/outputs/out_group_gamma.rds")
write.csv(rbind(part1,part2),"./example/outputs/new3.csv")



##### Figure C4 ############################################################################################################################
pdf("./example/plots/comparison_emiss_group_level2.pdf",width = 16,height = 8)
# Default bar plot
supp.labs <- c("Mean", "Variance","Var mu")
names(supp.labs) <- levels(final_to_display_dep$estimate)
dose.labs <- c("down","dreadnest of day","worry",
               "inadequate", "tired", "content",
               "agitated", "irritated", "switch and focus",
               "extremely well", "full of ideas", "thoughts racing")
names(dose.labs) <- c("down","dreadnest of day","worry",
                      "inadequate", "tired", "content",
                      "agitated", "irritated", "switch and focus",
                      "extremely well", "full of ideas", "thoughts are racing")

final_to_display_dep %>% filter(estimate!="varmu") %>% mutate(State=factor(state,levels = c("1","2","3","4"),labels = c("euthymic", "manic", "mixed","depressive")))  %>% mutate(Model=factor(model,levels = c("medhmm","mhmm"),labels = c("MEDHMM","MHMM")))%>%  ggplot(., aes(x=State, y=mean, group=Model,shape=Model,color=State)) +
  scale_shape_manual(values = c(19,23))+scale_color_manual(values =c(colors$blue,colors$orange,colors$red,colors$teal))+
  geom_errorbar(aes(ymin=ci_025, ymax=ci_975), width=0.5,linewidth=0.55,
                position=position_dodge(.6))+
  geom_point(position=position_dodge(0.6),aes(shape=Model),size=2)+scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "")) +
  guides(y = "none")+
  #scale_x_discrete(breaks=c(1,2,3,4),limits=c(0.5,4.5))+
  facet_grid(rows=vars(estimate),cols =vars(item),scales = "free_y",switch = "y",
             labeller = labeller(item=dose.labs, estimate=supp.labs))+theme_minimal()+theme(strip.text.x = element_text(size = 12),strip.text.y = element_text(size =12),axis.text.x = element_blank(),
                                                                                            axis.text.y = element_blank(),
                                                                                            axis.ticks = element_blank() ,axis.title = element_blank(),legend.text.align = 0,legend.position = 'bottom')

dev.off()

## Save output as RDS to compare
saveRDS(final_to_display_dep, "./example/outputs/out_group_emiss.rds")

#############################################################################################################################################

#### Figure 10 PREPERATION ###################################################################################################################
desired_names_gam<-c("S1toS1","S1toS2" ,"S1toS3","S1toS4", "S2toS1", "S2toS2","S2toS3","S2toS4","S3toS1","S3toS2","S3toS3", "S3toS4","S4toS1","S4toS2","S4toS3", "S4toS4")

##----------subject specific output if needed supplied here------------------------------------------------------------------------
###----MEDHMM----------------------------------------------------------------------------------------------------------------------
out_ss_gamma_medhmm<-lapply(out_medHMM_m4_12dv_it2000_c1$PD_subj, function(x){new=parent.frame()$i[]
x %>% as.data.frame()%>%  dplyr::select(desired_names_gam) %>% MAP_any_df%>% as.data.frame()%>% pivot_longer(everything() ,  names_to = c("item", ".value"), names_sep="[.]" ) %>% mutate(subj=new,model="medhmm")} )
out_ss_gamma_medhmm<-do.call(rbind, out_ss_gamma_medhmm)
###---MHMM-------------------------------------------------------------------------------------------------------------------------
out_ss_gamma_mhmm<-lapply(out_mHMM_m4_12dv_it2000_c1$PD_subj, function(x){new=parent.frame()$i[]
x %>% as.data.frame()%>%  dplyr::select(desired_names_gam) %>% MAP_any_df %>% as.data.frame() %>%  pivot_longer(-subj ,  names_to = c("item", ".value"), names_sep="[.]" )%>% mutate(subj=new,model="mhmm")%>% dplyr::select(-"S1toS4.mean")} )
out_ss_gamma_mhmm<-do.call(rbind, out_ss_gamma_mhmm)
###----
saveRDS(out_gamma1,"./example/outputs/out_ss_gamma.rds")
write.csv(out_gamma1,"./example/outputs/new4.csv")
#----------------------------------------------------------------------------------------------------------------------------------

#### Plot ss transition matrices
plot_subject_gamma<-function(subj=1){
  sub1_ed<-out_medHMM_m4_12dv_it2000_c1$PD_subj[[subj]]
  sub1_hm<-out_mHMM_m4_12dv_it2000_c1$PD_subj[[subj]] %>% as.data.frame()%>%  dplyr::select(desired_names_gam)
  map_sub1_ed<-MAP_any_df(sub1_ed)
  map_sub1_hm<-MAP_any_df(sub1_hm,mod="mhmm_gam")
  #MEDHMM
  data_trans_mat2<-matrix(as.numeric(map_sub1_ed[,1:(m*m)]),nrow=m,byrow = T)
  row.names(data_trans_mat2)<-paste0("State",1:m)
  colnames(data_trans_mat2)<-paste0("State",1:m)
  g3<-melt(data_trans_mat2)
  g3$Var1 <- factor(g3$Var1, levels = paste0("State",m:1), labels =c("depressive","mixed","manic","euthymic") )
  g3$Var2 <- factor(g3$Var2,
                    levels = paste0("State",1:m), labels =c("euthymic","manic","mixed","depressive") )
  g3$value<-round(g3$value,2)
  tr_plt_med<-ggplot(g3, aes(Var2,Var1)) + geom_tile(aes(fill= value))+coord_fixed()+geom_text(aes(label = value), color = "white", size = 5) +scale_fill_distiller(direction = 1)+guides(fill = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                                                                 barheight = 5,title = "Transition probability"))+xlab("")+ylab("")+theme_minimal()+theme(legend.position = "none",plot.title = element_text(size=12, hjust = 0.5),plot.subtitle = element_text(size=9,hjust=0.05),panel.grid.minor =element_blank(),axis.text.x = element_text(angle = 35, hjust=1,size=9),
                                                                                                                                                                                                                                                                                                          axis.text.y = element_text(angle = 35,size=9))+labs(title = paste0("Patient no. ",subj))

  #MHMM
  data_trans_mat1<-matrix(as.numeric(map_sub1_hm[,17:(2*m*m)]),nrow=m,byrow = T)
  row.names(data_trans_mat1)<-paste0("State",1:m)
  colnames(data_trans_mat1)<-paste0("State",1:m)
  g2<-melt(data_trans_mat1)
  g2$Var1 <- factor(g2$Var1,
                    levels = paste0("State",m:1), labels =c("depressive","mixed","manic","euthymic") )
  g2$Var2 <- factor(g2$Var2,
                    levels = paste0("State",1:m), labels =c("euthymic","manic","mixed","depressive") )
  g2$value<-round(g2$value,2)
  tr_plt_m<-ggplot(g2, aes(Var2,Var1)) + geom_tile(aes(fill= value))+coord_fixed()+geom_text(aes(label = value), color = "white", size = 5) +scale_fill_distiller(direction = 1)+guides(fill = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                                                               barheight = 5,title = "Transition probability"))+xlab("")+ylab("")+theme_minimal()+theme(legend.position = "none",plot.title = element_text(size=12, hjust = 0.5),plot.subtitle = element_text(size=9,hjust=0.05),panel.grid.minor =element_blank(),axis.text.x = element_text(angle = 35, hjust=1,size=9),
                                                                                                                                                                                                                                                                                                        axis.text.y = element_blank())+ggtitle("MHMM",subtitle = paste0("Subject ", subj))

  print(tr_plt_med)
}


#plot desired trace plots as exampples
#plot patient 3
data_trans_mat2<-matrix(as.numeric(out_ss_gamma_medhmm[[3]]$median),nrow=m,byrow = T)
row.names(data_trans_mat2)<-paste0("State",1:m)
colnames(data_trans_mat2)<-paste0("State",1:m)
g3<-melt(data_trans_mat2)
g3$X1 <- factor(g3$X1, levels = paste0("State",m:1), labels =c("depressive","mixed","manic","euthymic") )
g3$X2 <- factor(g3$X2,
                levels = paste0("State",1:m), labels =c("euthymic","manic","mixed","depressive") )
g3$value<-round(g3$value,2)


pp3<-ggplot(g3, aes(X2,X1)) + geom_tile(aes(fill= value))+coord_fixed()+geom_text(aes(label = value), color = "white", size = 5) +scale_fill_distiller(direction = 1)+guides(fill = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                                                    barheight = 5,title = "Transition probability"))+xlab("")+ylab("From state")+theme_minimal()+theme(legend.position = "none",plot.title = element_text(),plot.subtitle = element_text(size=11,hjust=0.05),panel.grid.minor =element_blank(),axis.text.x = element_text(angle = 35, hjust=1,size=10),
                                                                                                                                                                                                                                                                                                       axis.text.y = element_text(angle = 35,size=10),axis.title = element_text(size=12))+labs(title = paste0("patient no. ",3))

#plot patient 4
data_trans_mat2<-matrix(as.numeric(out_ss_gamma_medhmm[[4]]$median),nrow=m,byrow = T)
row.names(data_trans_mat2)<-paste0("State",1:m)
colnames(data_trans_mat2)<-paste0("State",1:m)
g3<-melt(data_trans_mat2)
g3$X1 <- factor(g3$X1, levels = paste0("State",m:1), labels =c("depressive","mixed","manic","euthymic") )
g3$X2 <- factor(g3$X2,
                levels = paste0("State",1:m), labels =c("euthymic","manic","mixed","depressive") )
g3$value<-round(g3$value,2)


pp4<-ggplot(g3, aes(X2,X1)) + geom_tile(aes(fill= value))+coord_fixed()+geom_text(aes(label = value), color = "white", size = 5) +scale_fill_distiller(direction = 1)+guides(fill = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                                                    barheight = 5,title = "Transition probability"))+xlab("To state")+ylab("")+theme_minimal()+theme(legend.position = "none",plot.title = element_text(),plot.subtitle = element_text(size=11,hjust=0.05),panel.grid.minor =element_blank(),axis.text.x = element_text(angle = 35, hjust=1,size=10),
                                                                                                                                                                                                                                                                                                     axis.text.y = element_blank(),axis.title = element_text(size=12))+labs(title = paste0("patient no. ",4))
#plot 14
data_trans_mat2<-matrix(as.numeric(out_ss_gamma_medhmm[[14]]$median),nrow=m,byrow = T)
row.names(data_trans_mat2)<-paste0("State",1:m)
colnames(data_trans_mat2)<-paste0("State",1:m)
g3<-melt(data_trans_mat2)
g3$X1 <- factor(g3$X1, levels = paste0("State",m:1), labels =c("depressive","mixed","manic","euthymic") )
g3$X2 <- factor(g3$X2,
                levels = paste0("State",1:m), labels =c("euthymic","manic","mixed","depressive") )
g3$value<-round(g3$value,2)


pp14<-ggplot(g3, aes(X2,X1)) + geom_tile(aes(fill= value))+coord_fixed()+geom_text(aes(label = value), color = "white", size = 5) +scale_fill_distiller(direction = 1)+guides(fill = guide_colourbar(barwidth = 0.5,
                                                                                                                                                                                                     barheight = 5,title = "Transition probability"))+xlab("")+ylab("")+theme_minimal()+theme(legend.position = "none",plot.title = element_text(),plot.subtitle = element_text(size=11,hjust=0.05),panel.grid.minor =element_blank(),axis.text.x = element_text(angle = 35, hjust=1,size=10),
                                                                                                                                                                                                                                                                                              axis.text.y = element_blank())+labs(title = paste0("patient no. ",14))

desired_names_dwell<-c("dur_logmu1" ,"dur_logmu2","dur_logmu3", "dur_logmu4" ,"dur_logsigma21", "dur_logsigma22", "dur_logsigma23","dur_logsigma24")
desired_names_dwell_hm<-c("S1toS1","S2toS2","S3toS3", "S4toS4")
list_final_log<-list()
list_final_exp<-list()

out_lognorm_est<-pbmcapply::pbmclapply(1:length(out_medHMM_m4_12dv_it2000_c1$PD_subj), function(subj){
  sub1_ed<-out_medHMM_m4_12dv_it2000_c1$PD_subj[[subj]] %>% as.data.frame()%>%  dplyr::select(desired_names_dwell) %>% MAP_any_df(.) %>% as.data.frame() %>% mutate(subj=subj)%>% pivot_longer(-subj ,  names_to = c("variable", ".value"), names_sep="[.]" ) %>% mutate(state=rep(1:4,2),variable=c(rep("logmu",m),rep("logsigma2",m)))
  return(sub1_ed)
},mc.cores = 1)
out_lognorm_est<-do.call(rbind,out_lognorm_est)
saveRDS(out_lognorm_est, "./example/outputs/dwell_time_medhmm.rds")

exp_g_dwell_medhmm<- MAP_medhmm(out_medHMM_m4_12dv_it2000_c1,J=2000,B=1000,m=4)[["dwell_mu_bar"]]  %>% as.data.frame()%>%mutate(model="medhmm") %>% pivot_longer(-model,  names_to = c("variable", ".value"), names_sep="[.]" ) %>% mutate(expected_median=exp(mean),.before=3)
exp_g_dwell_medhmm2<- MAP_medhmm(out_medHMM_m4_12dv_it2000_c1,J=2000,B=1000,m=4)[["dwell_var_bar"]]  %>% as.data.frame()%>%mutate(model="medhmm") %>% pivot_longer(-model,  names_to = c("variable", ".value"), names_sep="[.]" ) %>% mutate(expected_median=exp(mean),.before=3)
exp_g_dwell_mhmm<- MAP_mHMM(out_mHMM_m4_12dv_it2000_c1,J=2000,B=1000,m=4)[["dwell_mu_bar"]]%>% as.data.frame()%>%mutate(model="mhmm") %>% pivot_longer(-model,  names_to = c("variable", ".value"), names_sep="[.]" )
write.csv(exp_g_dwell_mhmm,"./example/outputs/exp_g_dwell_mhmm.csv")
write.csv(exp_g_dwell_medhmm,"./example/outputs/exp_g_dwell_medhmm.csv")
write.csv(exp_g_dwell_medhmm2,"./example/outputs/exp_g_dwell_medhmm2.csv")

#=====Plot group-level
df <- data.frame(
  state=factor(rep(c("Euthymic", "Manic","Mixed","Depressive"), each=10000)),
  distr=(c(rlnorm(10000, meanlog=exp_g_dwell_medhmm$median[1], sdlog =sqrt(exp_g_dwell_medhmm2$median[1])),
           rlnorm(10000, meanlog=exp_g_dwell_medhmm$median[2], sdlog =sqrt(exp_g_dwell_medhmm2$median[2])),
           rlnorm(10000, meanlog=exp_g_dwell_medhmm$median[3], sdlog =sqrt(exp_g_dwell_medhmm2$median[3])),
           rlnorm(10000, meanlog=exp_g_dwell_medhmm$median[4], sdlog =sqrt(exp_g_dwell_medhmm2$median[4]))))
)
dwell_time_medhmm <- readRDS("./example/outputs/dwell_time_medhmm.rds")
ss_df<-list()
for (i in 1:20) {
  dwell_time_medhmm2<-dwell_time_medhmm %>% filter(subj==i)
  ss_df[[i]]<-data.frame(subj=as.factor(rep(i,4*50000)),
                         state=factor(rep(c("Euthymic", "Manic","Mixed","Depressive"), each=50000)),
                         distr=(c(rlnorm(50000, meanlog=dwell_time_medhmm2$median[1], sdlog =sqrt(dwell_time_medhmm2$median[5])),
                                  rlnorm(50000, meanlog=dwell_time_medhmm2$median[2], sdlog =sqrt(dwell_time_medhmm2$median[6])),
                                  rlnorm(50000, meanlog=dwell_time_medhmm2$median[3], sdlog =sqrt(dwell_time_medhmm2$median[7])),
                                  rlnorm(50000, meanlog=dwell_time_medhmm2$median[4], sdlog =sqrt(dwell_time_medhmm2$median[8])))))
}

new_at<-rbind(ss_df[[4]],ss_df[[3]],ss_df[[14]])
theme_set(theme_minimal())

## FIGURE 10 ################################################################################################################################################################
pdf("./example/plots/my_dwell_time.pdf")
plot1<-df %>% ggplot(., aes(x =distr, y = factor(state,levels = c("Depressive","Mixed", "Manic","Euthymic")), fill=state) )+
  stat_density_ridges(quantile_lines = FALSE,
                      quantiles = 2,bandwidth = 0.8, alpha = 0.85,from=0, to=20)+
  scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+ xlab("State duration (in hours)")+ylab("")+scale_x_continuous(breaks = c(0,5,10,15,20),
                                                                                                                                                                                                                                                                 labels = c('0', '15', '30',"45","60"))+theme(legend.position = "none",axis.title.y = element_blank(),axis.text.y = element_text(size=11),axis.text.x = element_text(size=10))+labs(title = "Group-level dwell time distributions")
plot1
dev.off()


p2b<-ss_df[[4]]%>% ggplot(., aes(x =distr, y = factor(state,levels = c("Depressive","Mixed", "Manic","Euthymic")), fill=state) )+
  stat_density_ridges(mapping = aes(x=distr,y=factor(state,levels = c("Depressive","Mixed", "Manic","Euthymic"))),bandwidth=0.593, alpha = 0.85,from=0, to=20)+
  scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+ xlab("State duration (in hours)")+ylab("")+scale_x_continuous(breaks = c(0,5,10,15,20),
                                                                                                                                                                                                                                                                 labels = c('0 ', '15 ', '30 ',"45 ","60 "))+theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_text(size= 11),axis.text.y = element_blank(),axis.text.x = element_text(size=10))+labs(title ="patient no. 4")

p3b<-ss_df[[3]]%>% ggplot(., aes(x =distr, y = factor(state,levels = c("Depressive","Mixed", "Manic","Euthymic")), fill=state) )+
  stat_density_ridges(mapping = aes(x=distr,y=factor(state,levels = c("Depressive","Mixed", "Manic","Euthymic"))),bandwidth=0.8, alpha = 0.85,from=0, to=20)+
  scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+ xlab("")+ylab("")+scale_x_continuous(breaks = c(0,5,10,15,20),
                                                                                                                                                                                                                                        labels = c('0 ', '15 ', '30 ',"45 ","60 "))+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_text(size=11),axis.text.x = element_text(size=10))+labs(title="patient no. 3")

p4b<-ss_df[[14]]%>% ggplot(., aes(x =distr, y = factor(state,levels = c("Depressive","Mixed", "Manic","Euthymic")), fill=state) )+
  stat_density_ridges(mapping = aes(x=distr,y=factor(state,levels = c("Depressive","Mixed", "Manic","Euthymic"))), bandwidth = 0.593,alpha = 0.85,from=0, to=20)+
  scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+scale_fill_manual(name = "State", values = c(colors$blue,colors$orange,colors$red,colors$teal))+ xlab("")+ylab("")+scale_x_continuous(breaks = c(0,5,10,15,20),
                                                                                                                                                                                                                                        labels = c('0 ', '15 ', '30 ',"45 ","60 " ))+theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(size=10))+labs(title="patient no. 14")

thm <- theme(plot.title = element_text(face = 2, size = 16))



plotfirs<-wrap_elements((pp3+pp4+pp14) +plot_annotation(title = "C", theme = thm))
plotsec<-wrap_elements((p3b+p2b+p4b)+plot_annotation(title = "D", theme = thm))
#### FIGURE 10 IMPROVED########################
pdf("extra_examples1.pdf",height = 5,width = 10)
wrap_elements(wrap_elements(tr_plt_med+plot_annotation(title = "A", theme = thm))+wrap_elements(plot1+plot_annotation(title = "B", theme = thm)))
dev.off()
pdf("extra_examples2.pdf",height = 10,width = 10)
(plotfirs/ plotsec)+ plot_layout( heights = c(1.6,1))
# plot_layout(design = layout)
dev.off()


####################################################################################################################################################################

#additional matrices saved

desired_names_dwell<-c("dur_logmu1" ,"dur_logmu2","dur_logmu3", "dur_logmu4" ,"dur_logsigma21", "dur_logsigma22", "dur_logsigma23","dur_logsigma24")
desired_names_dwell_hm<-c("S1toS1","S2toS2","S3toS3", "S4toS4")
list_final_log<-list()
list_final_exp<-list()

out_lognorm_est<-pbmcapply::pbmclapply(1:length(out_medHMM_m4_12dv_it2000_c1$PD_subj), function(subj){
  sub1_ed<-out_medHMM_m4_12dv_it2000_c1$PD_subj[[subj]] %>% as.data.frame()%>%  dplyr::select(desired_names_dwell) %>% MAP_any_df(.) %>% as.data.frame() %>% mutate(subj=subj)%>% pivot_longer(-subj ,  names_to = c("variable", ".value"), names_sep="[.]" ) %>% mutate(state=rep(1:4,2),variable=c(rep("logmu",m),rep("logsigma2",m)))
  return(sub1_ed)
},mc.cores = 1)
out_lognorm_est<-do.call(rbind,out_lognorm_est)
write.csv(out_lognorm_est,"./example/outputs/new6.csv")
saveRDS(out_lognorm_est, "./example/outputs/dwell_time_medhmm.rds")

out_rate_est<-pbmcapply::pbmclapply(1:length(out_mHMM_m4_12dv_it2000_c1$PD_subj), function(subj){
  sub1_hm<-out_mHMM_m4_12dv_it2000_c1$PD_subj[[subj]] %>% as.data.frame()%>%  dplyr::select(desired_names_dwell_hm)
  map_sub1_hm<-MAP_any_df(sub1_hm) %>% as.data.frame() %>% mutate(subj=subj)%>% pivot_longer(-subj ,  names_to = c("variable", ".value"), names_sep="[.]" ) %>% mutate(state=rep(1:4))})
out_rate_est<-do.call(rbind,out_rate_est)
saveRDS(out_rate_est, "./example/outputs/dwell_time_rates_mhmm.rds")
write.csv(out_rate_est,"./example/outputs/new5.csv")



## Figure 12 PREPARATION ##################

desired_names_dwell<-c("dur_logmu1" ,"dur_logmu2","dur_logmu3", "dur_logmu4" ,"dur_logsigma21", "dur_logsigma22", "dur_logsigma23","dur_logsigma24")
desired_names_dwell_hm<-c("S1toS1","S2toS2","S3toS3", "S4toS4")
list_final_log<-list()
list_final_exp<-list()
for(subj in 1:20){
  sub1_ed<-out_medHMM_m4_12dv_it2000_c1$PD_subj[[subj]] %>% as.data.frame()%>%  dplyr::select(desired_names_dwell)
  sub1_hm<-out_mHMM_m4_12dv_it2000_c1$PD_subj[[subj]] %>% as.data.frame()%>%  dplyr::select(desired_names_dwell_hm)
  map_sub1_ed<-MAP_any_df(sub1_ed)[,9:16]
  map_sub1_hm<-MAP_any_df(sub1_hm)[,5:8]
  log_distr<-data.frame(state=c(1,2,3,4),log_mu=as.numeric(map_sub1_ed[1:4]),logsigma2=as.numeric(map_sub1_ed[5:8]))
  expo_dist<-data.frame(state=c(1,2,3,4),rate=as.numeric(map_sub1_hm))
  x=seq(0,100,0.01)
  final_df_log<-data.frame(subject=as.factor(rep(subj,length(x))),x=x)
  final_df_exp<-data.frame(subject=as.factor(rep(subj,length(x))),x=x)
  for(st in 1:4){
    final_df_log$names<-dlnorm(x,meanlog =log_distr$log_mu[st],sdlog = sqrt(log_distr$logsigma2[st]))
    colnames(final_df_log)[st+2]<-paste0("state",st)
    final_df_exp$names<-dexp(x,rate = expo_dist$rate[st])
    colnames(final_df_exp)[st+2]<-paste0("state",st)
  }

  list_final_log[[subj]]<-final_df_log
  list_final_exp[[subj]]<-final_df_exp
}
final_df_log<-do.call(rbind,list_final_log)
final_df_exp<-do.call(rbind,list_final_exp)


final_df_log2<-final_df_log %>% pivot_longer(cols = c(3:6),names_to = "states",values_to = "value") %>% filter(subject %in% c(3,4,14))
final_df_log2$subject<- as.numeric(final_df_log2$subject)
final_df_exp2<-final_df_exp %>% pivot_longer(cols = c(3:6),names_to = "states",values_to = "value")%>% filter(subject %in% c(3,4,14))
final_df_log2$subject<-factor(final_df_log2$subject,levels=unique(final_df_log2$subject),labels = c("patient no. 3","patient no. 4","patient no. 14"))
final_df_exp2$subject<-factor(final_df_exp2$subject,levels=unique(final_df_exp2$subject),labels = c("patient no. 3","patient no. 4","patient no. 14"))
final_df_log2$states<-factor(final_df_log2$states,levels=unique(final_df_log2$states),labels = c("euthymic","manic","mixed","depressive"))
final_df_exp2$states<-factor(final_df_exp2$states,levels=unique(final_df_exp2$states),labels =  c("euthymic","manic","mixed","depressive"))

##### FIGURE 12 ##############################################################################################################################################
pdf("./example/plots/expected_dwell_distr.pdf",width = 10,height = 4)
final_df_log2 %>%
  ggplot(.,aes(x=x,y=value,color=states))+
  geom_line(linetype = "solid",size=1)+
  geom_line(data = final_df_exp2 ,aes(x=x,y=value,color=states),linetype = "dashed",size=1.1)+  ylim(c(0,0.2))+
  theme_minimal()+
  scale_colour_manual(values=c(colors$orange,colors$red,colors$teal,colors$blue))+
  scale_x_continuous(breaks = c(0,4,8,12),limits = c(0,12),
                     labels = c('0 ', '12 ', '24',"36" ))+facet_wrap(~subject)+
  labs(x="State duration (in hours)",y="",color="Mood State:")+theme(legend.position = "bottom",strip.text.x = element_text(hjust = 0, margin=margin(l=0)))
dev.off()

##############################################################################################################################################################

out_switches_medhmm<-pbmcapply::pbmclapply(levels(decoding_medhmm$subject), function(sub){
  sto<-decoding_medhmm%>% filter(subject==sub)
  nobs<-nrow(sto)
  out<-data.frame(switches=length(rle(sto[,2])[[2]])-1)%>% mutate(subject=rep(as.numeric(sub),nrow(.)),observations=nobs,.before=1) %>% mutate(relative_switch=round(switches/nobs,3))
  return(out)
})
out_switches_medhmm<-do.call(rbind,out_switches_medhmm) %>% mutate(model=rep("medhmm",nrow(.)))

out_switches_mhmm<-pbmcapply::pbmclapply(levels(decoding_mhmm$subject), function(sub){
  sto<-decoding_mhmm%>% filter(subject==sub)
  nobs<-nrow(sto)
  out<-data.frame(switches=length(rle(sto[,2])[[2]])-1)%>% mutate(subject=rep(as.numeric(sub),nrow(.)),observations=nobs,.before=1) %>% mutate(relative_switch=round(switches/nobs,3))
  return(out)
})
out_switches_mhmm<-do.call(rbind,out_switches_mhmm) %>% mutate(model=rep("mhmm",nrow(.)))

out_switches<-rbind(out_switches_medhmm,out_switches_mhmm) %>% pivot_wider(names_from = model,values_from = c(-subject,-observations))

saveRDS(out_switches,"./example/outputs/out_ss_switches.rds")
write.csv(out_switches,"./example/outputs/out_ss_switches.csv")


#### PPC simulation #################################################################################################################################################################
#Note that emission distribution and more need to be taken from the new_map output
new_map<-MAP_medhmm(out_medHMM_m4_12dv_it2000_c1,B=1000,J=2000,m=4)
gamma<-matrix(c(as.numeric(new_map$gamma_prob_bar[1,(m*m+1):(2*m*m)])),byrow = T,nrow = m)
gamma_ss_var<-matrix(as.vector(new_map[["gamma_V_int_bar"]][,1:8]),nrow=4,byrow = T)

emissions<-list(        #1st dep variable
  matrix(c(new_map$emiss_mu_bar[[1]][(m+1):(2*m)],new_map$emiss_var_bar[[1]][(m+1):(2*m)]),byrow = F, nrow = m),
  #2nd dep variable
  matrix(c(new_map$emiss_mu_bar[[2]][(m+1):(2*m)],new_map$emiss_var_bar[[2]][(m+1):(2*m)]),byrow = F, nrow = m),
  #3rd dep variable
  matrix(c(new_map$emiss_mu_bar[[3]][(m+1):(2*m)],new_map$emiss_var_bar[[3]][(m+1):(2*m)]),byrow = F, nrow = m),
  #4th dep variable
  matrix(c(new_map$emiss_mu_bar[[4]][(m+1):(2*m)],new_map$emiss_var_bar[[4]][(m+1):(2*m)]),byrow = F, nrow = m),
  #5th dep variable
  matrix(c(new_map$emiss_mu_bar[[5]][(m+1):(2*m)],new_map$emiss_var_bar[[5]][(m+1):(2*m)]),byrow = F, nrow = m),
  #6th dep variable
  matrix(c(new_map$emiss_mu_bar[[6]][(m+1):(2*m)],new_map$emiss_var_bar[[6]][(m+1):(2*m)]),byrow = F, nrow = m),
  #7th dep variable
  matrix(c(new_map$emiss_mu_bar[[7]][(m+1):(2*m)],new_map$emiss_var_bar[[7]][(m+1):(2*m)]),byrow = F, nrow = m),
  #8th dep variable
  matrix(c(new_map$emiss_mu_bar[[8]][(m+1):(2*m)],new_map$emiss_var_bar[[8]][(m+1):(2*m)]),byrow = F, nrow = m),
  #9th dep variable
  matrix(c(new_map$emiss_mu_bar[[9]][(m+1):(2*m)],new_map$emiss_var_bar[[9]][(m+1):(2*m)]),byrow = F, nrow = m),
  #10th dep variable
  matrix(c(new_map$emiss_mu_bar[[10]][(m+1):(2*m)],new_map$emiss_var_bar[[10]][(m+1):(2*m)]),byrow = F, nrow = m),
  #11th dep variable
  matrix(c(new_map$emiss_mu_bar[[11]][(m+1):(2*m)],new_map$emiss_var_bar[[11]][(m+1):(2*m)]),byrow = F, nrow = m),
  #12th dep variable
  matrix(c(new_map$emiss_mu_bar[[12]][(m+1):(2*m)],new_map$emiss_var_bar[[12]][(m+1):(2*m)]),byrow = F, nrow = m))

emiss_ss_var<-lapply(new_map$emiss_varmu_bar,function(x){x[,5:8]})
dwell<-cbind(exp(new_map$dwell_mu_bar[(m+1):(2*m)]),exp(new_map$dwell_var_bar[(m+1):(2*m)]),new_map$dwell_varmu_bar[(m+1):(2*m)])

set.seed(12)
sampl<-500
posterior_checks_list<-list()
for (iter in 1:sampl) {
  posterior_checks_list[[iter]] <-mhsmm_sim(n = 20,
                                            n_t = c(386, 589, 538, 405, 431, 540, 608, 559, 399, 613, 308, 427, 434, 428, 486, 541, 197, 869, 489, 566),
                                            m = 4,
                                            n_dep = 12,
                                            emiss = emissions,
                                            gamma = gamma,
                                            gamma_ss_var = gamma_ss_var,
                                            emiss_ss_var = emiss_ss_var,
                                            dwell_distr = dwell)}

# save pcc for the MEDHMM
saveRDS(posterior_checks_list,file = "./example/outputs/ppc_medhmm.rds")

#MHMM
new_map<-MAP_mHMM(out_mHMM_m4_12dv_it2000_c1,B=1000,J=2000,m=4)
gamma<-matrix(c(as.numeric(new_map$gamma_prob_bar[1,(m*m+1):(2*m*m)])),byrow = T,nrow = m)
gamma_ss_var<-matrix(as.vector(new_map[["gamma_V_int_bar"]][,1:12]),nrow=4,byrow = T)

emissions<-list(        #1st dep variable
  matrix(c(new_map$emiss_mu_bar[[1]][(m+1):(2*m)],new_map$emiss_var_bar[[1]][(m+1):(2*m)]),byrow = F, nrow = m),
  #2nd dep variable
  matrix(c(new_map$emiss_mu_bar[[2]][(m+1):(2*m)],new_map$emiss_var_bar[[2]][(m+1):(2*m)]),byrow = F, nrow = m),
  #3rd dep variable
  matrix(c(new_map$emiss_mu_bar[[3]][(m+1):(2*m)],new_map$emiss_var_bar[[3]][(m+1):(2*m)]),byrow = F, nrow = m),
  #4th dep variable
  matrix(c(new_map$emiss_mu_bar[[4]][(m+1):(2*m)],new_map$emiss_var_bar[[4]][(m+1):(2*m)]),byrow = F, nrow = m),
  #5th dep variable
  matrix(c(new_map$emiss_mu_bar[[5]][(m+1):(2*m)],new_map$emiss_var_bar[[5]][(m+1):(2*m)]),byrow = F, nrow = m),
  #6th dep variable
  matrix(c(new_map$emiss_mu_bar[[6]][(m+1):(2*m)],new_map$emiss_var_bar[[6]][(m+1):(2*m)]),byrow = F, nrow = m),
  #7th dep variable
  matrix(c(new_map$emiss_mu_bar[[7]][(m+1):(2*m)],new_map$emiss_var_bar[[7]][(m+1):(2*m)]),byrow = F, nrow = m),
  #8th dep variable
  matrix(c(new_map$emiss_mu_bar[[8]][(m+1):(2*m)],new_map$emiss_var_bar[[8]][(m+1):(2*m)]),byrow = F, nrow = m),
  #9th dep variable
  matrix(c(new_map$emiss_mu_bar[[9]][(m+1):(2*m)],new_map$emiss_var_bar[[9]][(m+1):(2*m)]),byrow = F, nrow = m),
  #10th dep variable
  matrix(c(new_map$emiss_mu_bar[[10]][(m+1):(2*m)],new_map$emiss_var_bar[[10]][(m+1):(2*m)]),byrow = F, nrow = m),
  #11th dep variable
  matrix(c(new_map$emiss_mu_bar[[11]][(m+1):(2*m)],new_map$emiss_var_bar[[11]][(m+1):(2*m)]),byrow = F, nrow = m),
  #12th dep variable
  matrix(c(new_map$emiss_mu_bar[[12]][(m+1):(2*m)],new_map$emiss_var_bar[[12]][(m+1):(2*m)]),byrow = F, nrow = m))

emiss_ss_var<-lapply(new_map$emiss_varmu_bar,function(x){x[,5:8] %>% as.matrix()})
set.seed(12)
sampl<-500
posterior_checks_list<-pbapply::pblapply(1:500, function(s) sim_mHMM_plnorm(n_t = 600, n = 20,
                                                                            data_distr = "continuous",
                                                                            m = out_mHMM_m4_12dv_it2000_c1$input$m,
                                                                            n_dep =out_mHMM_m4_12dv_it2000_c1$input$n_dep,
                                                                            gamma = gamma,
                                                                            emiss_distr = emissions,
                                                                            var_gamma = gamma_ss_var,
                                                                            var_emiss = emiss_ss_var,
                                                                            return_ind_par = TRUE))
saveRDS(posterior_checks_list,file = "./example/outputs/ppc_mhmm.rds")


##### FIGURES C1 & C2 in separate files
