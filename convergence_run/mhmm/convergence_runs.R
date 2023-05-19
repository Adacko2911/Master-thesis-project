#-------------------------------------------------------------------------------
# Convergence checks script for mhmm
#-------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
#install or load required libraries

#library(devtools)
#devtools::install_github("emmekeaarts/mHMMbayes@continuous-emiss")
library(mHMMbayes)
library(pbmcapply)
library(tidyverse)
cores<-96


# Define helper function to output only MAP estimates
MAP_mHMM<-function(case_out=NULL,iteration=1,B=NULL,J=NULL,m=NULL){
  #---emiss_mu_bar---------
  emiss_mu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_mu_bar"]])){
    median<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_025),"_ci_975")
    emiss_mu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_mu_bar)<-paste0("emiss_mu_bar",iteration)
    emiss_mu_bar_list[[i]]<-emiss_mu_bar
  }


  #---emiss_var_bar--------
  emiss_var_bar_list<-list()
  for(i in 1:length(case_out[["emiss_var_bar"]])){
    median<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
    emiss_var_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_var_bar)<-paste0("emiss_var_bar",iteration)
    emiss_var_bar_list[[i]]<-emiss_var_bar
  }
  #---emiss_varmu_bar------
  emiss_varmu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_varmu_bar"]])){
    median<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),"_median")
    mean<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),"_mean")
    sd<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),"_sd")
    ci_025<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
    emiss_varmu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_varmu_bar)<-paste0("emiss_varmu_bar",iteration)
    emiss_varmu_bar_list[[i]]<-emiss_varmu_bar
  }
  #---gamma_prob_bar------
  median<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
  gamma_prob_bar<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(gamma_prob_bar)<-paste0("gamma_prob_bar",iteration)

  #---gamma_prob_bar_no_diag------
  #preprocess the chain outputs
  self_trans_lab<-c()
  gamma_prob_bar_no_diag<-case_out[["gamma_prob_bar"]][(B+1):J,]
  for (i in 1:m) {
    self_trans_lab[i]<-paste0("S",i,"toS",i)
  }
  zeros<-matrix(rep(0,m*(J-B)),ncol=m)
  gamma_prob_bar_no_diag[,self_trans_lab] <-zeros
  n<-0
  while(n<m*m){
    gamma_prob_bar_no_diag[,(n+1):(n+m)]<-gamma_prob_bar_no_diag[,(n+1):(n+m)]%>% apply(.,1,function(x){x/sum(x)}) %>% t()
    n<-n+m
  }
  #assign the gamma without diagonal entries

  median<-gamma_prob_bar_no_diag %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-gamma_prob_bar_no_diag %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-gamma_prob_bar_no_diag %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-gamma_prob_bar_no_diag %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-gamma_prob_bar_no_diag %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
  gamma_prob_bar_no_diag<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(gamma_prob_bar_no_diag)<-paste0("gamma_prob_bar_no_diag")


  #---gamma_int_prob_bar------
  median<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
  gamma_int_bar<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(gamma_int_bar)<-paste0("gamma_int_bar",iteration)

  #--gamma_V_int_bar-----------
  #subject specific parameters of variance covariance between
  all_st<-c()
  for(st_cur in 1:m){
    all_st<-append(all_st,paste0("S",st_cur,"toS",2:m))
  }
  names_need<-c()
  for(al in 1:length(all_st)){
    names_need<-append(names_need,paste0("var_int_",all_st[al],"_with_int_",all_st[al]))
  }
  case_out_gamma_V_int_bar<-case_out[["gamma_V_int_bar"]][(B+1):J,names_need]
  median<-case_out_gamma_V_int_bar %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out_gamma_V_int_bar %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out_gamma_V_int_bar %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
  gamma_V_int_bar<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(gamma_V_int_bar)<-paste0("gamma_V_int_bar",iteration)

  #dwell_mu_bar_exp
  self_trans_lab<-c()
  gamma_prob_bar_new<-case_out[["gamma_prob_bar"]][(B+1):J,]
  for (i in 1:m) {
    self_trans_lab[i]<-paste0("S",i,"toS",i)
  }
  dwell_mu_bar_exp<-gamma_prob_bar_new[,self_trans_lab]
  dwell_mu_bar_exp<-dwell_mu_bar_exp%>% apply(., 1, function(x){new=1/(-log(x))
  return(new)}) %>% t()
  colnames(dwell_mu_bar_exp)<-paste0("dwell_mu_bar",1:m)
  median<-dwell_mu_bar_exp %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-dwell_mu_bar_exp %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-dwell_mu_bar_exp %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-dwell_mu_bar_exp %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-dwell_mu_bar_exp %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),"_ci_975")
  dwell_mu_bar_exp<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(dwell_mu_bar_exp)<-paste0("dwell_mu_bar",iteration)


  MAP<-list(emiss_mu_bar_list,emiss_var_bar_list,emiss_varmu_bar_list,gamma_prob_bar,gamma_prob_bar_no_diag,gamma_int_bar,gamma_V_int_bar,dwell_mu_bar_exp)
  names(MAP)<-c("emiss_mu_bar","emiss_var_bar","emiss_varmu_bar","gamma_prob_bar","gamma_prob_bar_no_diag","gamma_int_bar","gamma_V_int_bar","dwell_mu_bar")
  return(MAP)
}
rep_each <- function(x, times) {
  times <- rep(times, length.out = length(x))
  rep(x, times = times)
}
local_decoding<-function(out=NA){
  out<-out
  n<-out$input$n_subj
  m<-out$input$m
  J=out$input$J
  B=out$input$burn_in
  local_d=out$sample_path
  subj<-length(local_d)
  len<-out$input$n_vary
  local_d<-do.call(rbind,local_d)
  sub<-rep_each(c(1:subj),as.vector(len))
  sel_data<-local_d[,-c(1:B)] %>% as.data.frame()
  st1<-sel_data%>% apply(.,1,function(x){cl<-names(sort(summary(as.factor(x)), decreasing=T))
  most_lik<-cl[1]
  return(most_lik)})
  st<-as.numeric(st1)
  output<-data.frame(subject=sub,state=st)
  return(output)
}

#------------------------------------------------------------------------------
#Gelman rubin statistics calculation
gelman_rubin <- function(par_matrix, J=J,burn_in=burn_in){
  chains<-2
  samples<-J-burn_in
  # Coerce to matrix
  par_matrix <- as.data.frame(par_matrix)
  out<-c()
  for(param in 2:ncol(par_matrix)){
    # Mean over all samples
    all_mean <- mean(par_matrix[,param])

    # Mean of each chain
    chain_mean <- tapply(par_matrix[,param], par_matrix[,1], mean)

    # Variance of each chain
    chain_var <- tapply(par_matrix[,param], par_matrix[,1], stats::var)
    W <- (1 / chains) * sum(chain_var)
    B <- samples / (chains - 1) * sum((chain_mean - all_mean)^2)
    V <- (1 - 1 / samples) * W + (1 / samples) * B
    out[param-1]<-round(sqrt(V / W), 4)
  }
  return(round(out,1))
}
ge_rub_two_models_emp<-function(model1=NULL, model2=NULL, m=2, J=4000,burn_in=2000){

  chains_matrix_gamma<-matrix(ncol = (m^2-(m)))
  chains_matrix_gamma<-rbind(chains_matrix_gamma,model1$gamma_int_bar[(burn_in+1):J,1:(m^2-(m))],model2$gamma_int_bar[(burn_in+1):J,1:(m^2-(m))])
  chains_matrix_gamma<-na.omit(chains_matrix_gamma)
  var_data<-rbind(model1$emiss_var_bar[[1]][(burn_in+1):(J),],model2$emiss_var_bar[[1]][(burn_in+1):(J),])
  mean_data<-rbind(model1$emiss_mu_bar[[1]][(burn_in+1):J,1:m],model2$emiss_mu_bar[[1]][(burn_in+1):J,1:m])
  var_data2<-rbind(model1$emiss_var_bar[[2]][(burn_in+1):(J),],model2$emiss_var_bar[[2]][(burn_in+1):(J),])
  mean_data2<-rbind(model1$emiss_mu_bar[[2]][(burn_in+1):J,1:m],model2$emiss_mu_bar[[2]][(burn_in+1):J,1:m])
  chains_matrix_emiss<-cbind(mean_data,var_data,mean_data2,var_data2)
  chains_matrix_emiss<-na.omit(chains_matrix_emiss)
  mean_dwell<-rbind(model1$dwell_mu_bar[(burn_in+1):J,1:m],model2$dwell_mu_bar[(burn_in+1):J,1:m])
  var_dwell<-rbind(model1$dwell_var_bar[(burn_in+1):(J),],model2$dwell_var_bar[(burn_in+1):(J),])
  chains_matrix_dwell<-cbind(mean_dwell,var_dwell)
  index<-c(rep(1,J-burn_in),rep(2,J-burn_in))
  gelman_input<-cbind(index,chains_matrix_emiss,chains_matrix_dwell,chains_matrix_gamma)
  out_scenario_gelman_empirical<-gelman_rubin(gelman_input,J = J,burn_in = burn_in)
  out_scenario_gelman_empirical<-t(as.matrix(out_scenario_gelman_empirical))
  colnames(out_scenario_gelman_empirical)<-colnames(gelman_input)[-1]
  return(out_scenario_gelman_empirical)
}

#--------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#parameter specification
#-------------------------------------------------------------------------------
#----Legend---------------------------------
# There are 3 hidden states scenarios m={1,2,3} (m)
# There are 4 senarios for observation number (o)
# * 1 - 100 obs
# * 2 - 200 obs
# * 3 - 500 obs
# * 4 - 1000 obs
# There are 3 scenarios for dwell time (d)
# * 1 - 1.4 time points (self-transition in mhmm is 0.5)
# * 2 - 3.5 time points (self-transition in mhmm is 0.75)
# * 3 - 19.5 time points (self-transition in mhmm is 0.95)
# J is representation of total number of iterations in MCMC chain (J)
# B is representation of burn_in period for iterations in MCMC chain (B)
# There are 100 simulated datasets from each scenario and each of the datasets
#  will to be fitted to the mhmm model (s)

# Creating vectors
m <- 3:4
o <- 1:3
d <- 1:4
s<-1:5

# Calling expand.grid() Function to create the parameter dataframe
pars = expand.grid(m=m,o=o,d=d,s=s)
#set global variables
J  <- 4000
B  <- 2000
n_dep=2

# Define the random datasets to be inspected with convergence
set.seed(1234)
sampled_datasets<-sample(1:128,max(s))
#========= Load datasets and extract the chosen ones====================================================
#---------Three state scenario chose appropriate number of datasets-------------------------------------
#In order to run the scripts without any issues using snellius we need to adjust the below target data directory
three_state<-readRDS(paste0("./data/three_state_design_2dep_final.rds"))

three_state_chosen_datasets<-list()
for(obs in o){
  dwell_list<-list()
  for(dur in d){
    sampled_data<-list()
    for(sampl in sampled_datasets){
      sampled_data[[which(sampled_datasets==sampl)]]<-three_state[[obs]][[dur]][[sampl]]
    }
    dwell_list[[dur]]<-sampled_data
  }
  three_state_chosen_datasets[[obs]]<-dwell_list
}
rm(three_state)
#---------Four state scenario chose appropriate number of datasets---------------------------------------
#In order to run the scripts without any issues using snellius we need to adjust the below target data directory
four_state<-readRDS(paste0("./data/four_state_design_2dep_final.rds"))
four_state_chosen_datasets<-list()
for(obs in o){
  dwell_list<-list()
  for(dur in d){
    sampled_data<-list()
    for(sampl in sampled_datasets){
      sampled_data[[which(sampled_datasets==sampl)]]<-four_state[[obs]][[dur]][[sampl]]
    }
    dwell_list[[dur]]<-sampled_data
  }
  four_state_chosen_datasets[[obs]]<-dwell_list
}
rm(four_state)


three_state_extra<-readRDS(paste0("./data/three_state_design_mixed_dwell_2dep_final.rds"))
three_state_chosen_datasets2<-list()
for(obs in 1){
  dwell_list<-list()
  for(dur in 1){
    sampled_data<-list()
    for(sampl in sampled_datasets){
      sampled_data[[which(sampled_datasets==sampl)]]<-three_state_extra[[obs]][[dur]][[sampl]]
    }
    dwell_list[[dur]]<-sampled_data
  }
  three_state_chosen_datasets2[[obs]]<-dwell_list
}
rm(three_state_extra)


# Save data files
#In order to run the scripts without any issues using smellius we need to adjust the below data directory
saveRDS(list(three_state_chosen_datasets2,sampled_datasets), paste0("./data/convergence_mhmm_run_chosen_data_three_mixed_states.rds"))

# Save only chosen data sets such that they are not occupy the computer less
#In order to run the scripts without any issues using snellius we need to adjust the below target data directory
 saveRDS(list(two_state_chosen_datasets,sampled_datasets), paste0("./data/convergence_mhmm_run_chosen_data_two_states.rds"))
 saveRDS(list(three_state_chosen_datasets,sampled_datasets), paste0("./data/convergence_mhmm_run_chosen_data_three_states.rds"))
 saveRDS(list(four_state_chosen_datasets,sampled_datasets), paste0("./data/convergence_mhmm_run_chosen_data_four_states.rds"))

#==================================================================================================================================================================================================
#Outstart clean
#==================================================================================================================================================================================================
out_list1 <-pbmcapply::pbmclapply(1:nrow(pars), function(parameter_row) {

  # Gen
  case_out<-NULL
  out1<-NULL
  out2<-NULL

  # Specify simulation parameters
  m  <- pars[parameter_row,1]
  o  <- pars[parameter_row,2]
  d  <- pars[parameter_row,3]
  s  <- pars[parameter_row,4]
  #==================================================================================================================================================================================================
  #load: simulated data, hyper-priors and starting values according (m)

    if(m==3){
      #simulated data
      sim_data<-three_state_chosen_datasets
      #---Emission distribution starting values
      emiss_distr <- matrix(c(10,20*1.2,
                              30, 60*1.2,
                              60,120*1.2), nrow = m, byrow = TRUE)
      emiss_distr2 <- matrix(c(60,(120+0.25*120)*1.2,
                               30, (60+0.25*60)*1.2,
                               10,(20+20*0.25)*1.2), nrow = m, byrow = TRUE)
      #---Emission distribution starting values 2nd chain
      emiss_distr_2nd_chain <- matrix(c(10*sample(c(0.8,1.2),1),20*1.2*sample(c(0.8,1.2),1),
                                        30*sample(c(0.8,1.2),1), 60*1.2*sample(c(0.8,1.2),1),
                                        60*sample(c(0.8,1.2),1),120*1.2*sample(c(0.8,1.2),1)), nrow = m, byrow = TRUE)
      emiss_distr2_2nd_chain<- matrix(c(60*sample(c(0.8,1.2),1),(120+0.25*120)*1.2*sample(c(0.8,1.2),1),
                                        30*sample(c(0.8,1.2),1), (60+0.25*60)*1.2*sample(c(0.8,1.2),1),
                                        10*sample(c(0.8,1.2),1),(20+20*0.25)*1.2*sample(c(0.8,1.2),1)), nrow = m, byrow = TRUE)

      #---Gamma starting values
      if(d==1){
        gamma   <- matrix(c(0.5, 0.25,0.25,
                            0.25, 0.5, 0.25,
                            0.25, 0.25, 0.5), ncol = m, byrow = TRUE)
        gamma2   <- matrix(c(0.5, 0.25,0.25,
                             0.25, 0.5, 0.25,
                             0.25, 0.25, 0.5)+rnorm(m*m,mean=0,sd=0.04), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
      }else if(d==2){
        gamma   <- matrix(c(0.7, 0.15,0.15,
                            0.15, 0.7, 0.15,
                            0.15, 0.15, 0.7), ncol = m, byrow = TRUE)
        gamma2   <- matrix(c(0.7, 0.15,0.15,
                             0.15, 0.7, 0.15,
                             0.15, 0.15, 0.7)+rnorm(m*m,mean=0,sd=0.04), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
      }else if(d==3){
        gamma   <- matrix(c(0.95, 0.025,0.025,
                            0.025, 0.95, 0.025,
                            0.025, 0.025, 0.95), ncol = m, byrow = TRUE)
        gamma2   <- matrix(c(0.95, 0.025,0.025,
                             0.025, 0.95, 0.025,
                             0.025, 0.025, 0.95)+rnorm(m*m,mean=0,sd=0.03), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
      }else if(d==4){
        gamma   <- matrix(c(0.95, 0.025,0.025,
                            0.025, 0.95, 0.025,
                            0.025, 0.025, 0.95), ncol = m, byrow = TRUE)
        gamma2   <- matrix(c(0.95, 0.025,0.025,
                             0.025, 0.95, 0.025,
                             0.025, 0.025, 0.95)+rnorm(m*m,mean=0,sd=0.04), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
      }else if(d=5){
        sim_data<-three_state_chosen_datasets2
        gamma   <- matrix(c(0.75, 0.25/2,0.25/2,
                            0.025, 0.95, 0.025,
                            0.025, 0.025, 0.95), ncol = m, byrow = TRUE)
        gamma2   <- matrix(c(0.75, 0.25/2,0.25/2,
                             0.025, 0.95, 0.025,
                             0.025, 0.025, 0.95)+rnorm(m*m,mean=0,sd=0.04), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
      }
      #---Specify hyper-prior for the continuous emission distribution----------------
      hyp_pr <- list(
        emiss_mu0 = list(matrix(c(10,30,60), nrow = 1),
                         matrix(c(60,30,10), nrow = 1)),
        emiss_K0  = list(1,1),
        emiss_nu  = list(1,1),
        emiss_V   = list(rep(100, m),
                         rep(100, m)),
        emiss_a0  = list(rep(0.001,m),
                         rep(0.001,m)),
        emiss_b0  = list(rep(0.001,m),
                         rep(0.001,m)))


    }else
      if(m==4){
        #simulated data
        sim_data<-four_state_chosen_datasets
        #---Emission distribution starting values
        emiss_distr <- matrix(c(10,20*1.2,
                                30, 60*1.2,
                                60,120*1.2,
                                98,150*1.2), nrow = m, byrow = TRUE)
        emiss_distr2 <- matrix(c(98,(150+0.25*150)*1.2,
                                 60, (120+0.25*120)*1.2,
                                 30,(60+0.25*60)*1.2,
                                 10,(20+0.25*20)*1.2), nrow = m, byrow = TRUE)
        #---Emission distribution starting values 2nd chain
        emiss_distr_2nd_chain <- matrix(c(10*sample(c(0.8,1.2),1),20*1.2*sample(c(0.8,1.2),1),
                                          30*sample(c(0.8,1.2),1), 60*1.2*sample(c(0.8,1.2),1),
                                          60*sample(c(0.8,1.2),1),120*1.2*sample(c(0.8,1.2),1),
                                          98*sample(c(0.8,1.2),1),150*1.2*sample(c(0.8,1.2),1)), nrow = m, byrow = TRUE)
        emiss_distr2_2nd_chain <- matrix(c(98*sample(c(0.8,1.2),1),(150+0.25*150)*1.2*sample(c(0.8,1.2),1),
                                           60*sample(c(0.8,1.2),1), (120+0.25*120)*1.2*sample(c(0.8,1.2),1),
                                           30*sample(c(0.8,1.2),1),(60+0.25*60)*1.2*sample(c(0.8,1.2),1),
                                           10*sample(c(0.8,1.2),1),(20+0.25*20)*1.2*sample(c(0.8,1.2),1)), nrow = m, byrow = TRUE)
        #---Gamma starting values
        if(d==1){
          gamma   <- matrix(c(0.5, 0.5/3,0.5/3,0.5/3,
                              0.5/3,0.5,0.5/3,0.5/3,
                              0.5/3,0.5/3,0.5,0.5/3,
                              0.5/3,0.5/3,0.5/3,0.5), ncol = m, byrow = TRUE)
          gamma2   <- matrix(c(0.5, 0.5/3,0.5/3,0.5/3,
                               0.5/3,0.5,0.5/3,0.5/3,
                               0.5/3,0.5/3,0.5,0.5/3,
                               0.5/3,0.5/3,0.5/3,0.5)+rnorm(m*m,mean=0,sd=0.04), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
        }else if(d==2){
          gamma   <- matrix(c(0.75, 0.25/3,0.25/3,0.25/3,
                              0.25/3,0.75,0.25/3,0.25/3,
                              0.25/3,0.25/3,0.75,0.25/3,
                              0.25/3,0.25/3,0.25/3,0.75), ncol = m, byrow = TRUE)
          gamma2   <- matrix(c(0.75, 0.25/3,0.25/3,0.25/3,
                               0.25/3,0.75,0.25/3,0.25/3,
                               0.25/3,0.25/3,0.75,0.25/3,
                               0.25/3,0.25/3,0.25/3,0.75)+rnorm(m*m,mean=0,sd=0.04), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
        }else if(d==3){
          gamma   <- matrix(c(0.9, 0.1/3,0.1/3,0.1/3,
                              0.1/3,0.9,0.1/3,0.1/3,
                              0.1/3,0.1/3,0.9,0.1/3,
                              0.1/3,0.1/3,0.1/3,0.9), ncol = m, byrow = TRUE)
          gamma2   <- matrix(c(0.9, 0.1/3,0.1/3,0.1/3,
                               0.1/3,0.9,0.1/3,0.1/3,
                               0.1/3,0.1/3,0.9,0.1/3,
                               0.1/3,0.1/3,0.1/3,0.9)+rnorm(m*m,mean=0,sd=0.03), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
        }else if(d==4){
          gamma   <- matrix(c(0.9, 0.1/3,0.1/3,0.1/3,
                              0.1/3,0.9,0.1/3,0.1/3,
                              0.1/3,0.1/3,0.9,0.1/3,
                              0.1/3,0.1/3,0.1/3,0.9), ncol = m, byrow = TRUE)
          gamma2   <- matrix(c(0.9, 0.1/3,0.1/3,0.1/3,
                               0.1/3,0.9,0.1/3,0.1/3,
                               0.1/3,0.1/3,0.9,0.1/3,
                               0.1/3,0.1/3,0.1/3,0.9)+rnorm(m*m,mean=0,sd=0.04), ncol = m, byrow = TRUE) %>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
        }


        #---Specify hyper-prior for the continuous emission distribution----------------
        hyp_pr <- list(
          emiss_mu0 = list(matrix(c(10,30,60,98), nrow = 1),
                           matrix(c(98,60,30,10), nrow = 1)),
          emiss_K0  = list(1,1),
          emiss_nu  = list(1,1),
          emiss_V   = list(rep(100, m),
                           rep(100, m)),
          emiss_a0  = list(rep(0.001,m),
                           rep(0.001,m)),
          emiss_b0  = list(rep(0.001,m),
                           rep(0.001,m)))

      }


  #Fit 2nd chain MHMM
  # Try statement, to catch errors and avoid interrupting the script
  ti <- Sys.time()

  input_data_mhmm<-sim_data[[o]][[d]][[s]]
  out1<- try(mHMMbayes::mHMM_cont(s_data = input_data_mhmm$observations,
                                  gen = list(m = m, n_dep = n_dep),
                                  start_val = c(list(gamma), list(emiss_distr,emiss_distr2)),
                                  emiss_hyp_prior = hyp_pr,
                                  mcmc = list(J = J, burn_in = B),
                                  show_progress = TRUE) )
  out1[["time1"]] <- Sys.time() - ti
  if(!is.null(out1)){
  MAP1<-try(MAP_mHMM(case_out =out1,J=J,B=B,m=m,iteration=s))
}
  #Fit 2nd chain MHMM
  # Try statement, to catch errors and avoid interrupting the script
  ti <- Sys.time()

  out2<- try(mHMMbayes::mHMM_cont(s_data = input_data_mhmm$observations,
                                  gen = list(m = m, n_dep = n_dep),
                                  start_val = c(list(gamma2), list(emiss_distr_2nd_chain,emiss_distr2_2nd_chain)),
                                  emiss_hyp_prior = hyp_pr,
                                  mcmc = list(J = J, burn_in = B),
                                  show_progress = TRUE) )
  out2[["time"]] <- Sys.time() - ti
  if(!is.null(out2)){
  MAP2<-try(MAP_mHMM(case_out =out2,J=J,B=B,m=m,iteration=s))
}


  # Calculate necessary statistics to be saved
  if(!is.null(out1) & !is.null(out2)){
    gelman_rubin_out<-try(ge_rub_two_models_emp(model1 = out1,model2 = out2,m = m,J = J,burn_in = B))
  }

  # Define new output object
  case_out <- list(
    chain1=out1,
    MAP1=MAP1,
    chain2=out2,
    MAP2=MAP2,
    gelman_rubin_stats = gelman_rubin_out
  )
  # Save output
  #In order to run the scripts without any issues using snellius we need to adjust the below target directory
  try(saveRDS(case_out, paste0("./adacko/convergence_run/mhmm/outputs/convergence_run",
                               "_m",m,
                               "_nt",o,
                               "_dwell",d,
                               "dt_set",s,
                               "_it",J,
                               "_burn_in",B,
                               ".rds")))


  return(case_out)
}, mc.cores = cores, mc.set.seed = 42L)


# Calling expand.grid() Function to create the parameter dataframe
pars3 = expand.grid(m=3,o=o,d=d,s=s)
pars4 = expand.grid(m=4,o=o,d=d,s=s)


gelman_list_all3 <-pbmcapply::pbmclapply(1:nrow(pars3), function(parameter_row) {


  # Specify simulation parameters
  m  <- pars3[parameter_row,1]
  o  <- pars3[parameter_row,2]
  d  <- pars3[parameter_row,3]
  s  <- pars3[parameter_row,4]

  #In order to run the scripts without any issues using snellius we need to adjust the below data directory
  dataset<- readRDS(paste0("./convergence_run/mhmm/outputs/convergence_run",
                               "_m",m,
                               "_nt",o,
                               "_dwell",d,
                               "dt_set",s,
                               "_it",J,
                               "_burn_in",B,
                               ".rds"))
  out_gelman_row<-dataset$gelman_rubin_stats
  rownames(out_gelman_row)<- paste0("convergence_run",
                               "_m",m,
                               "_nt",o,
                               "_dwell",d,
                               "dt_set",s,
                               "_it",J,
                               "_burn_in",B)
  return(out_gelman_row)


}, mc.cores = cores, mc.set.seed = 42L)

gelman_list_all3<-do.call(rbind,gelman_list_all3)

#In order to run the scripts without any issues using snellius we need to adjust the below target directory
saveRDS(gelman_list_all3,paste0("./convergence_run/mhmm/outputs/all_gelaman_mhmm_m3.rds"))


gelman_list_all4 <-pbmcapply::pbmclapply(1:nrow(pars4), function(parameter_row) {


  # Specify simulation parameters
  m  <- pars4[parameter_row,1]
  o  <- pars4[parameter_row,2]
  d  <- pars4[parameter_row,3]
  s  <- pars4[parameter_row,4]

  #In order to run the scripts without any issues using snellius we need to adjust the below data directory
  dataset<- readRDS(paste0("./convergence_run/mhmm/outputs/convergence_run",
                               "_m",m,
                               "_nt",o,
                               "_dwell",d,
                               "dt_set",s,
                               "_it",J,
                               "_burn_in",B,
                               ".rds"))
  out_gelman_row<-dataset$gelman_rubin_stats
  rownames(out_gelman_row)<- paste0("convergence_run",
                               "_m",m,
                               "_nt",o,
                               "_dwell",d,
                               "dt_set",s,
                               "_it",J,
                               "_burn_in",B)
  return(out_gelman_row)


}, mc.cores = cores, mc.set.seed = 42L)

gelman_list_all4<-do.call(rbind,gelman_list_all4)

#In order to run the scripts without any issues using snellius we need to adjust the below target data directory
saveRDS(gelman_list_all4,paste0("./convergence_run/mhmm/outputs/all_gelaman_mhmm_m4.rds"))
#In order to run the scripts without any issues using snellius we need to adjust the below target data directory
saveRDS(sampled_datasets,paste0("./convergence_run/mhmm/outputs/randomly_chosen_datasets.rds"))


library(ggplot2)
library(gridExtra)
library(grid)

traceplots_emiss_mu<-function(model1,model2,dep=1){
  model1_em<-as.data.frame(model1$emiss_mu_bar[[dep]])
  model2_em<-as.data.frame(model2$emiss_mu_bar[[dep]])
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]] <-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("Emiss mu",dep,":dwell",d," :obs",o," :state",i," :set",s ))

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob(conv_name))))
  print(out)
}
traceplots_emiss_var<-function(model1,model2,dep=1){
  model1_em<-as.data.frame(model1$emiss_var_bar[[dep]])
  model2_em<-as.data.frame(model2$emiss_var_bar[[dep]])
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]]<-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("Emiss var",dep,":dwell",d," :obs",o," :state",i," :set",s ))

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob(conv_name))))
  print(out)
}


traceplots_gamma<-function(model1,model2){
  model1_em<-as.data.frame(model1$gamma_prob_bar)
  model2_em<-as.data.frame(model2$gamma_prob_bar)
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]]<-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+geom_line(data=model2_em, aes_string(x=1:data_points,y=model2_em[,i]),color="blue")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])+ggtitle(paste0("gamma: dwell",d," :obs",o," :state",i," :set",s))
  }
  out<-do.call(grid.arrange, c(p, list(top=textGrob(conv_name))))
  print(out)
}

# Creating vectors
m <- 3:4
o <- 1:3
d <- 1:4
s<-1:5
pars = expand.grid(m=m,o=o,d=d,s=s)
pars1<-expand.grid(m=3,o=2,d=5,s=s)
pars=rbind(pars,pars1)
#In order to run the scripts without any issues using snellius we need to adjust
#the below directory such that it is a desired convergence runs directory
pdf("./convergence_run/mhmm/outputs/trace_plots_mhmm.pdf")
for(parameter in 1:nrow(pars)){
  m  <- pars[parameter,1]
  o  <- pars[parameter,2]
  d  <- pars[parameter,3]
  s  <- pars[parameter,4]

  dataset<- readRDS(paste0("./convergence_run/mhmm/outputs/convergence_run","_m",m,
                           "_nt",o,
                           "_dwell",d,
                           "dt_set",s,
                           "_it",J,
                           "_burn_in",B,
                           ".rds"))
  conv_name<-paste0("mhmm_convergence_run","_m",m,
                    "_nt",o,
                    "_dwell",d,
                    "dt_set",s,
                    "_it",J,
                    "_burn_in",B)
  model1=try(dataset$chain1)
  model2=try(dataset$chain2)
  if(is.list(model2) && is.list(model1)){
    traceplots_emiss_mu(model1 = model1,model2 = model2,dep=1)
    traceplots_emiss_var(model1 = model1,model2 = model2,dep=1)
    traceplots_emiss_mu(model1 = model1,model2 = model2,dep=2)
    traceplots_emiss_var(model1 = model1,model2 = model2,dep=2)
    traceplots_gamma(model1 = model1,model2 = model2)
}
}
dev.off()


