#devtools::install_github("https://github.com/emmekeaarts/medHMM",ref='dev-poisson',build = TRUE,force = TRUE)
#devtools::install_github("https://github.com/emmekeaarts/mHMMbayes",ref='dev-poisson',build = TRUE, force = TRUE)
library(mHMMbayes)
library(medHMM)
require(pbmcapply) #for parallel computing
require(tidyverse) #for data wrangling
library(RColorBrewer)
library(gridExtra)
library(grid)
library(reshape2)
library(R.utils)
source("~/GitHub/Master-thesis-project/further_research/utils_sim_medhmm_shift_poiss_map_traceplots.R")
source("~/GitHub/Master-thesis-project/further_research/script-snellius-ola/utils_mapping_scenarios.R")


crit<-expand.grid(dur_s=1:12,ind_s=3, occ_s=2, dep_s=3) %>% as.data.frame()

cores=1
J = 15
burn_in = 10
m=3
j=6 # <<<< This parameter is passed by stoppos and based on the crit correct scenario is selected


# start
mhmm_and_medhmm_all_cores_result <-pbmcapply::pbmclapply(1:nrow(cores), function(j) {


#number of set of mean duration matrices
dur_s <- crit$dur_s[j]
#number of scenarios of the number of observations per individual
ind_s <- crit$ind_s[j]
n_subj=n_distr[[ind_s]]
#number of scenarios of the number of observations per individual
occ_s <-crit$occ_s[j]
n_t=n_t_distr[[occ_s]]
#number of dependent variable scenarios
dep_s<-crit$dep_s[j]
#number of samples within the same scenario
sampl <- 1

if(dep_s<=3){
  n_dep<-1
}else{
  n_dep=2
}

#=============== MHMM + data simulation ==========================================

#dwell_distr1 is expected dwell time and I calculate the start gamma matrix out of it
dwell_distr1=matrix(unlist(dwell_distr[[dur_s]]),ncol=1)
gam_start<-data.frame(expected=dwell_distr1) %>% summarise(gamma=exp(-1/expected)-0.1) %>% summarise(gamma=gamma,rest=(1-gamma)/2)
matrigg=c(gam_start[1,1],gam_start[1,2],gam_start[1,2],
          gam_start[2,2],gam_start[2,1],gam_start[2,2],
          gam_start[3,2],gam_start[3,2],gam_start[3,1])
gamma_start=matrix(data =matrigg ,nrow=3,byrow=TRUE)


if(n_dep==1){
  emiss=matrix(unlist(dep_distr[[dep_s]]), nrow = m, byrow = FALSE)
  emiss_start=cbind(emiss[,1],emiss[,2]*1.5)
  emiss<-list(emiss)

}else{
  emiss=lapply(dep_distr[[dep_s]],function(x){matrix(unlist(x), nrow = m, byrow = FALSE)}  )
emiss_start=lapply(emiss, function(x){x[,2]<-x[,2]*1.5
cbind(x[,1],x[,2])})
}


if(dep_s<=3){
  emiss_hyp_pr <- prior_emiss_cont(
  gen = list(m = m, n_dep = n_dep),
  emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][,1]), nrow = 1)),
  emiss_K0  = list(1),
  emiss_nu  = list(1),
  emiss_V   = list(rep(10, m)),
  emiss_a0  = list(rep(0.01, m)),
  emiss_b0  = list(rep(0.01, m))
)
}else{
  emiss_hyp_pr <- prior_emiss_cont(
    gen = list(m = m, n_dep = n_dep),
    emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][[1]][,1]), nrow = 1),
                     matrix(unlist(dep_distr[[dep_s]][[2]][,1]), nrow = 1)),
    emiss_K0  = list(1, 1),
    emiss_nu  = list(1, 1),
    emiss_V   = list(rep(10, m), rep(10, m)),
    emiss_a0  = list(rep(0.01, m), rep(0.01, m)),
    emiss_b0  = list(rep(0.01, m), rep(0.01, m))
  )
}

gamma <-matrix(c(0, 0.5, 0.5,
                  0.5, 0, 0.5,
                  0.5, 0.5, 0), nrow = m, ncol = m, byrow = TRUE)

dwell_ss_var<-dwell_var_ss_s[[dur_s]]

# SIMULATE DATA

sim_data <- mhsmm_sim_shift_pois(n = n_subj,
                                        n_t = n_t,
                                        m = m,
                                        n_dep = n_dep,
                                        emiss = emiss,
                                        gamma = gamma,
                                        gamma_ss_var = gamma_ss_var,
                                        emiss_ss_var = emiss_ss_var,
                                        dwell_distr = log(dwell_distr1),
                                        dwell_ss_var=dwell_ss_var)

assign(paste0("sim_data_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s), sim_data)
#saveRDS(sim_data,paste0("~/GitHub/Master-thesis-project/further_research/results/simulation/","hmm_sim_data_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s,"_iter",J,".rds"))

#check simulated data
colnames(sim_data$states)<-c("subject","state")
summary_dwell_ss<-sim_data$states %>%as.data.frame() %>%  group_by(subject)%>%
  summarise(length=rle(state)[[1]],state=rle(state)[[2]])%>%
  group_by(subject,state) %>%summarise(mean_emp_dwell=mean(length),median_emp_dwell=median(length))
summary_dwell_ss$subject<-as.factor(summary_dwell_ss$subject)
summary_dwell_ss$state<-as.factor(summary_dwell_ss$state)
summary_dwell_ss=as.data.frame(summary_dwell_ss)

ob<-data.frame(sim_data[['observations']])
colnames(ob)<-c('subject',paste0('dep',1:n_dep))
dw<-data.frame(sim_data[['states']])
colnames(dw)<-c('subject','state')
data=cbind(ob,state=dw$state) %>% mutate(., state=as.factor(state), subject=as.factor(subject))

#check simulated emission distributions
#ggplot(data,aes(x=dep1,fill=subject))+geom_density(show.legend = FALSE,alpha=0.5)+facet_wrap(~state,ncol=1)+theme_minimal()+ggtitle("Emission distribution check")
#check simulated distribution of shifted poisson lambdas
ggplot(summary_dwell_ss,aes(x=median_emp_dwell))+geom_histogram(show.legend = FALSE,alpha=0.5)+facet_wrap(~state,ncol=1,drop = TRUE)+theme_minimal()+ggtitle("Dwell mean distribution check")+coord_cartesian(xlim=c(0,70))

# MHMM
start_time1 <- Sys.time()


mHMM_cont<- try({mHMMbayes::mHMM(s_data = sim_data$observations,
                                             data_distr = 'continuous',
                                             gen = list(m = m, n_dep = n_dep),
                                             start_val = c(list(gamma_start), emiss),
                                             emiss_hyp_prior = emiss_hyp_pr,
                                             mcmc = list(J = J, burn_in = burn_in),
                                             return_path = TRUE)
  })

end_time1 <- Sys.time()

if(!is.null(mHMM_cont)){
  stored_map_mhmm<-MAP_mHMM(case_out = mHMM_cont,iteration = j,J=J,B = burn_in,m = m)
  mhmm_traceplots_emiss_mu1<-try(traceplots_emiss_mu(model2 = mHMM_cont,dep=1))
  if(n_dep>1){
    mhmm_traceplots_emiss_mu2<-try(traceplots_emiss_mu(model2 = mHMM_cont,dep=2))
  }
  mhmm_traceplots_gamma<-try(traceplots_gamma(model1 = mHMM_cont))

  state_decoding<-try(local_decoding(out=mHMM_cont))
  true_st<-as.data.frame(sim_data$states)
  mhmm_state_decoding_table<-try(data.frame(subject=true_st[,1],true_state_decoding=true_st[,-1],vit_state_decoding=state_decoding$state))
}

mhmm_case_out<-list(MAP=stored_map_mhmm,state_decoding=mhmm_state_decoding_table,execution_time=end_time1-start_time1)

### ??? Placeholder to save the pdf with trace plots

#============ MEDHMM ==============================================================

if(dur_s<=4){
  max_dwell=200
}else if(dur_s<=8){
  max_dwell=300
}else{
  max_dwell=400
}

dwell_distr1=matrix(unlist(dwell_distr[[dur_s]]),ncol=1)

dwell_hyp_pr <- list(
  dwell_mu0 = matrix(log(as.numeric(dwell_distr[[dur_s]])), nrow = 1, ncol = 3),
  dwell_K0  = c(1),
  dwell_nu  = c(1),
  dwell_V   = rep(0.1, m)
)
if(n_dep==1){
  emiss=matrix(unlist(dep_distr[[dep_s]]), nrow = m, byrow = FALSE)
  emiss_start=cbind(emiss[,1],emiss[,2]*1.5)
  emiss<-list(emiss)

}else{
  emiss=lapply(dep_distr[[dep_s]],function(x){matrix(unlist(x), nrow = m, byrow = FALSE)}  )
  emiss_start=lapply(emiss, function(x){x[,2]<-x[,2]*1.5
  cbind(x[,1],x[,2])})
}

if(dep_s<=3){
  emiss_hyp_pr <- list(
    emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][,1]), nrow = 1)),
    emiss_K0  = list(1),
    emiss_nu  = list(1),
    emiss_V   = list(rep(10, m)),
    emiss_a0  = list(rep(0.01, m)),
    emiss_b0  = list(rep(0.01, m))
  )
}else{
  emiss_hyp_pr <- list(
    emiss_mu0 = list(matrix(unlist(dep_distr[[dep_s]][[1]][,1]), nrow = 1),
                     matrix(unlist(dep_distr[[dep_s]][[2]][,1]), nrow = 1)),
    emiss_K0  = list(1, 1),
    emiss_nu  = list(1, 1),
    emiss_V   = list(rep(10, m), rep(10, m)),
    emiss_a0  = list(rep(0.01, m), rep(0.01, m)),
    emiss_b0  = list(rep(0.01, m), rep(0.01, m))
  )
}

start_time2 <- Sys.time()
medHMM_cont_shiftpois <- try(medHMM_cont_shiftpois(s_data = sim_data$observations,
                                                   shift = 1,
                                                   gen = list(m = m, n_dep = n_dep),
                                                   start_val = c(list(gamma), emiss, list(dwell_distr1)),
                                                   emiss_hyp_prior = emiss_hyp_pr,
                                                   dwell_hyp_prior = dwell_hyp_pr,
                                                   show_progress = TRUE,
                                                   mcmc = list(J = J, burn_in = burn_in),
                                                   return_path = TRUE,
                                                   max_dwell = max_dwell))

end_time2 <- Sys.time()


if(!is.null(medHMM_cont_shiftpois)){
  stored_map_medhmm<-MAP_medhmm(case_out = medHMM_cont_shiftpois,iteration = j,J=J,B = burn_in,m = m)

  medhmm_traceplots_emiss_mu1<-try(traceplots_emiss_mu(model2 = medHMM_cont_shiftpois,dep=1))
  if(n_dep>1){
    medhmm_traceplots_emiss_mu2<-try(traceplots_emiss_mu(model2 = medHMM_cont_shiftpois,dep=2))
  }
  medhmm_traceplots_gamma<-try(traceplots_gamma(model1 = medHMM_cont_shiftpois))
  medhmm_traceplots_dwell_mu<-try(traceplots_dwell_mu(model2 = medHMM_cont_shiftpois))

  state_decoding<-try(local_decoding(out=medHMM_cont_shiftpois))
  true_st<-as.data.frame(sim_data$states)
  medhmm_state_decoding_table<-try(data.frame(subject=true_st[,1],true_state_decoding=true_st[,-1],vit_state_decoding=state_decoding$state))
}

medhmm_case_out<-list(MAP=stored_map_medhmm,state_decoding=medhmm_state_decoding_table,execution_time=end_time2-start_time2)

### ??? Placeholder to save the pdf with trace plots


mhmm_and_medhhmm_case_out<-list(MHMM=mhmm_case_out,
                                MEDHMM=medhmm_case_out)
# below if intermediate result needed
#saveRDS(mhmm_and_medhhmm_case_out, paste0("~/GitHub/Master-thesis-project/further_research/results/main_run/","mhmm_results_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s,"_iter",J,".rds"))
return(mhmm_and_medhhmm_case_out)

}, mc.cores = cores, mc.set.seed = 123L)

saveRDS(mhmm_and_medhmm_all_cores_result, paste0("~/GitHub/Master-thesis-project/further_research/results/main_run/","mhmm_results_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s,"_iter",J,".rds"))

# IGNORE

#pdf(paste0("~/GitHub/Master-thesis-project/further_research/results/main_run/","hmm_results_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s,"_iter",J,".pdf"))

#PLOT EMISS
# if(n_dep==1){
#   colnames(sim_data[["observations"]])=c("subject", "var1")
#   summary_emiss_ss<-cbind(sim_data[["observations"]] %>% as.data.frame(), sim_data$states %>% as.data.frame() %>% select(state)) %>% round(4)
#   summary_emiss_ss<-summary_emiss_ss %>%
#   pivot_longer(.,cols =2,names_to = "variable",values_to = "value")}
# else{
#   colnames(sim_data[["observations"]])=c("subject", "var1","var2")
#   summary_emiss_ss<-cbind(sim_data[["observations"]] %>% as.data.frame(), sim_data$states %>% as.data.frame() %>% select(state)) %>% round(4)
#   summary_emiss_ss<-summary_emiss_ss %>% pivot_longer(.,cols =2:3,names_to = "variable",values_to = "value")
#   }
# summary_emiss_ss$subject=as.factor(summary_emiss_ss$subject)
# summary_emiss_ss$state=as.factor(summary_emiss_ss$state)
# means_emiss<-matrix(ncol=3)
#
# for(i in 1:length(emiss)){
# means_emiss<-rbind(means_emiss, cbind(state=1:3,variable=rep(paste0("var",i),m),value=unlist(emiss[[1]])[,1]))
# }
# means_emiss<-means_emiss[-1,] %>% as.data.frame()
# means_emiss$state<-as.factor(means_emiss$state)
# means_emiss$variable<-as.factor(means_emiss$variable)
# means_emiss$value<-as.numeric(means_emiss$value)
#
# summary_emiss_ss<-summary_emiss_ss %>% group_by(subject,state,variable)%>%
# summarise(mean=mean(value),median=median(value),variance=var(value))
#
# emiss_est_matrix<-matrix(ncol=6)
# for(i in 1:length(emiss)){
# emiss_est_matrix=rbind(emiss_est_matrix, cbind(state=1:3,variable=rep(paste0("var",i),m),
# true=emiss[[i]][,1],
# est_m=stored_map[["emiss_mu_bar"]][[i]][1:3],
# est_25=stored_map[["emiss_mu_bar"]][[i]][10:12],
# est_97=stored_map[["emiss_mu_bar"]][[i]][13:15]))
# }
#
# datx= emiss_est_matrix[-1,] %>% as.data.frame()
# colnames(datx)<-c("state","variable","true","est_m","est_25","est_97")
# datx<-transform(datx, true = as.numeric(true),
# est_m = as.numeric(est_m),
# est_25=as.numeric(est_25),
# est_97=as.numeric(est_97))
# datx$state<-as.factor(datx$state)
# datx$variable<-as.factor(datx$variable)
# plot_emiss_mu_data=right_join(summary_emiss_ss,datx,by=c("state","variable"))
# plot_emiss_mu_data$state<-as.factor(plot_emiss_mu_data$state)
#
# #PLOT
# simulated_vs_estimated_emiss_mu<-plot_emiss_mu_data%>%
# ggplot(.,aes(x=median,y=state,fill=state))+geom_boxplot()+
# geom_point(aes(y=state,x=true),color="red",size=2)+
# geom_errorbar(aes(xmin=est_25, xmax=est_97),color="blue")+
# scale_fill_brewer(palette = "RdYlBu")+theme_minimal()+ facet_wrap(~variable,ncol = 1)+
# ggtitle(paste0("MHMM Simulated emiss check (MEAN):","sim_data_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s))+
# labs(x=" Median emiss distribution means",y="State")+theme(axis.text.y=element_blank())
#
#
# emiss_est_matrix<-matrix(ncol=6)
# for(i in 1:length(emiss)){
# emiss_est_matrix=rbind(emiss_est_matrix, cbind(state=1:3,variable=rep(paste0("var",i),m),
# true=emiss[[i]][,2],
# est_m=stored_map[["emiss_var_bar"]][[i]][1:3],
# est_25=stored_map[["emiss_var_bar"]][[i]][10:12],
# est_97=stored_map[["emiss_var_bar"]][[i]][13:15]))
# }
#
# datx= emiss_est_matrix[-1,] %>% as.data.frame()
# colnames(datx)<-c("state","variable","true","est_m","est_25","est_97")
# datx<-transform(datx, true = as.numeric(true),
# est_m = as.numeric(est_m),
# est_25=as.numeric(est_25),
# est_97=as.numeric(est_97))
# datx$state<-as.factor(datx$state)
# datx$variable<-as.factor(datx$variable)
# plot_emiss_var_data=right_join(summary_emiss_ss,datx,by=c("state","variable"))
# plot_emiss_var_data$state<-as.factor(plot_emiss_mu_data
# $state)
#
# simulated_vs_estimated_emiss_var<-plot_emiss_var_data%>%
# ggplot(.,aes(x=variance,y=state,fill=state))+geom_boxplot()+
# geom_point(aes(y=state,x=true),color="red",size=2)+
# geom_errorbar(aes(xmin=est_25, xmax=est_97),color="blue")+
# scale_fill_brewer(palette = "RdYlBu")+theme_minimal()+ coord_flip() + facet_wrap(~variable,nrow = 1)+
# ggtitle(paste0("Simulated emiss check (VAR):","sim_data_","dur_s_",dur_s,"ind_s_",ind_s,"occ_s_",occ_s,"dep_s_",dep_s))+
# labs(x="MHMM Median emiss distribution variances",y="State")+theme(axis.text.x=element_blank())




