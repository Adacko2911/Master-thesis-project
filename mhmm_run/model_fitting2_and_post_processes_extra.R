#----------------------------------------------------------------------------------------------------------------------------------------------
# Fit model for 3 hidden states, 2 dependent variable, 80 subjects, o length of time series and varying duration of a hidden state and save results
#----------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
#specify number of cores
cores=128
# Defined a priori parameters (not changing over the course of the study)
n_dep <- 2
n_subj<-80

#install or load required libraries
#library(devtools)
#devtools::install_github("emmekeaarts/mHMMbayes@continuous-emiss")
require(mHMMbayes)
require(pbmcapply)
require(tidyverse)

#------------------------------------------------------------------------------
ifrm <- function(obj, env = globalenv()) {
    obj <- deparse(substitute(obj))
    if(exists(obj, envir = env)) {
        rm(list = obj, envir = env)
    }
}

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


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load parameters passed by Stopos
pars = commandArgs(trailingOnly=TRUE)
#as a form of check print the line passed
print(pars)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#pass parameters from stoppos
#----Legend---------------------------------
# There are 3 hidden states scenarios m={1,2,3} (m)
# There are 4 senarios for observation number (o)
# * 1 - 200 obs
# * 2 - 500 obs
# * 3 - 1000 obs
# There are 3 scenarios for dwell time (d)
# * 1 - 1.4 time points (self-transition in mhmm is 0.5)
# * 2 - 3.5 time points (self-transition in mhmm is 0.75)
# * 3 - 19.5 time points (self-transition in mhmm is 0.95)
# * 4 - 99.5 time points (self-transition in mhmm is 0.99)
# J is representation of total number of iterations in MCMC chain (J)
# B is representation of burn_in period for iterations in MCMC chain (B)
# There are 100 simulated datasets from each scenario and each of the datasets
#  will to be fitted to the mhmm model (s)


m  <- as.numeric(pars[1])
o  <- as.numeric(pars[2])
d  <- as.numeric(pars[3])
J  <- as.numeric(pars[4])
B  <- as.numeric(pars[5])


#==================================================================================================================================================================================================
#load: simulated data, hyper-priors and starting values according that depends on m and d
#==================================================================================================================================================================================================

    #simulated data
    #In order to run the scripts without any issues using snellius we need to adjust the below data directory
    sim_data<-readRDS(paste0("./data/three_state_design_mixed_dwell_2dep_final.rds"))
    #---Emission distribution starting values
    emiss_distr <- matrix(c(10,20*1.2,
                            30, 60*1.2,
                            60,120*1.2), nrow = m, byrow = TRUE)
    emiss_distr2 <- matrix(c(60,(120+0.25*120)*1.2,
                             30, (60+0.25*60)*1.2,
                             10,(20+20*0.25)*1.2), nrow = m, byrow = TRUE)

 gamma   <- matrix(c(0.75, 0.25/2,0.25/2,
                          0.025, 0.95, 0.025,
                          0.025, 0.025, 0.95), ncol = m, byrow = TRUE)
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



#==================================================================================================================================================================================================
#Outstart clean
#==================================================================================================================================================================================================
#loop over iteration of the datasets
out_list1 <-pbmcapply::pbmclapply(1:cores, function(s) {
  # Gen
  case_out<-NULL

  # Specify simulation parameters
  m  <- as.numeric(pars[1])
  o  <- as.numeric(pars[2])
  d  <- as.numeric(pars[3])
  J  <- as.numeric(pars[4])
  B  <- as.numeric(pars[5])
  d<-d-4
  o<-o-1
  #Fit MHMM
  # Try statement, to catch errors and avoid interrupting the script
  ti <- Sys.time()

  input_data_mhmm<-sim_data[[o]][[d]][[s]]
  case_out<- try(mHMMbayes::mHMM_cont(s_data = input_data_mhmm$observations,
                                      gen = list(m = m, n_dep = n_dep),
                                      start_val = c(list(gamma), list(emiss_distr,emiss_distr2)),
                                      emiss_hyp_prior = hyp_pr,
                                      mcmc = list(J = J, burn_in = B),
                                      show_progress = FALSE,return_path = TRUE) )
  case_out[["time"]] <- Sys.time() - ti

  # Calculate necessary statistics to be saved
  if(!is.null(case_out)){
    state_decoding<-local_decoding(case_out)
    true_st<-as.data.frame(input_data_mhmm$states)
    state_decoding_table<-data.frame(subject=true_st[,1],true_state_decoding=true_st[,-1],vit_state_decoding=state_decoding$state)
    rm(state_decoding)
    #---MAP estimates---------
    MAP<-MAP_mHMM(case_out =case_out,J=J,B=B,m=m,iteration=s)


  }

  # Define new output object
  case_out <- list(
    decoding = state_decoding_table,
    MAP=MAP)
  return(case_out)
}, mc.cores = cores, mc.set.seed = 42L)

# Save output as a list
#In order to run the scripts without any issues using snellius we need to adjust the below target directory
saveRDS(out_list1, paste0("./mhmm_run/outputs/mhmm_list_output",
                          "_m",m,
                          "_nt",o+1,
                          "_dwell",d+4,
                          "_it",J,
                          "_burn_in",B,
                          ".rds"))

m  <- as.numeric(pars[1])
o  <- as.numeric(pars[2])
d  <- as.numeric(pars[3])
J  <- as.numeric(pars[4])
B  <- as.numeric(pars[5])
#drop simulated data to save memory used
ifrm(sim_data)
ifrm(case_out)

namessd<-c()
namessd<-append(namessd,paste0("st_",m,"_obs_",o,"_dwell_",d,"_mhmm"))
dwell_time_names<-c(1.4,3.5,19.5,99.5,"mixed")
observation_lenght<-c(200,500,1000)

count_switches <- function(seq.input){
  COUNT <- 0
  switch_position<-c()
  for (i in 2:length(seq.input)) {
    if (seq.input[i] != seq.input[i - 1]) {
      COUNT <- COUNT + 1
      switch_position<-c(switch_position,i)
    }
  }
  total_switches <- COUNT
  list_out<-list(switch_position=switch_position,
       total_switches=total_switches
  )
  return(list_out)
}

kappa <- function(mat) {
  N <- sum(mat)
  No <- sum(diag(mat))
  dim_rm<-as.vector(which(rowSums(acc)==colSums(acc)))
  mat<-mat[-dim_rm,-dim_rm]
  Ne <- 1 / N * sum(colSums(mat) * rowSums(mat))
  return( (No - Ne) / (N - Ne) )
}

exp_mean1=1/(-log(0.5))
exp_mean2=1/(-log(0.75))
exp_mean3=1/(-log(0.95))
exp_mean4=1/(-log(0.99))
#define standard deviations
exp_sigma1=1.10093 #such that exp_sigma1=1/3*exp_mean1
exp_sigma2=1.10093 #such that exp_sigma2=1/3*exp_mean2
exp_sigma3=1.10093 #such that exp_sigma3=1/3*exp_mean3
exp_sigma4=1.05902 #such that exp_sigma4=1/4*exp_mean4

true_emiss_varmu_bar <- rep(9,m)
true_emiss_var_bar<-c(20,60,120)
true_emiss_mu_bar<-c(10,30,60)
true_emiss_varmu_bar2 <- rep(9,m)
true_emiss_var_bar2<-c(120+0.25*120,60+0.25*60,20+20*0.25)
true_emiss_mu_bar2<-c(60,30,10)
true_dwell<-c(exp_mean2,exp_mean3,exp_mean4)
true_gamma_prob_bar<-as.vector(matrix(c(0.75, 0.25/2,0.25/2,
                                   0.05/2,0.95,0.05/2,
                                   0.01/2,0.01/2,0.99), ncol = m, byrow = TRUE))
true_gamma_prob_bar_no_diag<-as.vector(matrix(c(0, 0.5,0.5,
                                            0.5,0,0.5,
                                            0.5,0.5,0), ncol = m, byrow = TRUE))

true_gamma_V_int_bar<-rep(0.17,(m*(m-1)))


emiss_mu_bar<-list()
emiss_mu_bar2<-list()
emiss_var_bar<-list()
emiss_var_bar2<-list()
emiss_varmu_bar<-list()
emiss_varmu_bar2<-list()
gamma_prob_bar<-list()
gamma_prob_bar_no_diag<-list()
gamma_int_bar<-list()
dwell_mu_bar<-list()
#list for V gamma and names relevant
gamma_V_int_bar<-list()
#prepare list for group level decoding list
decoding_out_list<-list()
decoding_out_list_truth<-list()
out_list1[sapply(out_list1, is.null)] <- NULL
cores<-length(out_list1)
for (i in 1:cores){
  emiss_mu_bar[[i]]<-out_list1[[i]]$MAP$emiss_mu_bar[[1]]
  emiss_var_bar[[i]]<-out_list1[[i]]$MAP$emiss_var_bar[[1]]
  emiss_varmu_bar[[i]]<-out_list1[[i]]$MAP$emiss_varmu_bar[[1]]
  emiss_mu_bar2[[i]]<-out_list1[[i]]$MAP$emiss_mu_bar[[2]]
  emiss_var_bar2[[i]]<-out_list1[[i]]$MAP$emiss_var_bar[[2]]
  emiss_varmu_bar2[[i]]<-out_list1[[i]]$MAP$emiss_varmu_bar[[2]]
  gamma_prob_bar[[i]]<-out_list1[[i]]$MAP$gamma_prob_bar
  gamma_prob_bar_no_diag[[i]]<-out_list1[[i]]$MAP$gamma_prob_bar_no_diag #added that
  dwell_mu_bar[[i]]<-out_list1[[i]]$MAP$dwell_mu_bar #
  gamma_int_bar[[i]]<-out_list1[[i]]$MAP$gamma_int_bar
  gamma_V_int_bar[[i]]<-out_list1[[i]]$MAP$gamma_V_int_bar

  decoding_out_list[[i]]<-out_list1[[i]][["decoding"]] %>% group_by(subject)%>%
   summarise(length=rle(vit_state_decoding)[[1]],point=rle(vit_state_decoding)[[2]])%>%
    group_by(point) %>%summarise(mean_emp_dwell=mean(length),median_emp_dwell=median(length),
    sd_emp_dwell=sd(length), ci_025=quantile(length,probs = 0.025),ci_975=quantile(length,probs = 0.975))

  #prepare the data sets for the true and estimated empirical dwell times estimates
  decoding_out_list_truth[[i]]<-out_list1[[i]][["decoding"]]%>% group_by(subject)%>%
    summarise(length=rle(true_state_decoding)[[1]],point=rle(true_state_decoding)[[2]])%>% group_by(point) %>%
    summarise(mean_emp_dwell=mean(length),median_emp_dwell=median(length),sd_emp_dwell=sd(length),
    ci_025=quantile(length,probs = 0.025),ci_975=quantile(length,probs = 0.975))

}
emiss_mu_bar<-do.call(rbind,emiss_mu_bar)
emiss_var_bar<-do.call(rbind,emiss_var_bar)
emiss_varmu_bar<-do.call(rbind,emiss_varmu_bar)
emiss_mu_bar2<-do.call(rbind,emiss_mu_bar2)
emiss_var_bar2<-do.call(rbind,emiss_var_bar2)
emiss_varmu_bar2<-do.call(rbind,emiss_varmu_bar2)
gamma_prob_bar<-do.call(rbind,gamma_prob_bar)
gamma_int_bar<-do.call(rbind,gamma_int_bar)
gamma_V_int_bar<-do.call(rbind,gamma_V_int_bar)
gamma_prob_bar_no_diag<-do.call(rbind,gamma_prob_bar_no_diag)
dwell_mu_bar<-do.call(rbind,dwell_mu_bar)


#joint the out results form each simulation list and calculate the out results for each simulation list
empirical_dwell_time_true<-do.call(rbind,decoding_out_list_truth) %>%
  as.data.frame() %>%group_by(point) %>% group_modify(~ apply(.x,2,mean) %>%t() %>%
                                                        as.data.frame()) %>% rename("state_scenario" = "point") %>%
                                                        as.data.frame() %>% mutate(dwell_time=rep(d,m), .before = state_scenario) %>%
                                                        mutate(observations=rep(o,m), .before = dwell_time)  %>%
                                                        mutate(state=rep(m,m), .before = observations) %>%mutate(state=rep("mhmm",m), .before = state)
empirical_dwell_time_mhmm<-do.call(rbind,decoding_out_list) %>% as.data.frame()%>%
  group_by(point) %>% group_modify(~ apply(.x,2,mean) %>%t() %>% as.data.frame())%>%
  rename("state_scenario" = "point")%>% as.data.frame() %>% mutate(dwell_time=rep(d,m), .before = state_scenario) %>%
  mutate(observations=rep(o,m), .before = dwell_time)  %>%  mutate(state=rep(m,m), .before = observations) %>%
  mutate(state=rep("mhmm",m), .before = state)

all_map<-list(emiss_mu_bar=emiss_mu_bar,emiss_var_bar=emiss_var_bar,emiss_varmu_bar=emiss_varmu_bar,
              emiss_mu_bar2=emiss_mu_bar2,emiss_var_bar2=emiss_var_bar2,emiss_varmu_bar2=emiss_varmu_bar2,
              gamma_prob_bar=gamma_prob_bar,gamma_int_bar=gamma_int_bar,gamma_V_int_bar=gamma_V_int_bar,gamma_prob_bar_no_diag=gamma_prob_bar_no_diag,
              dwell_mu_bar=dwell_mu_bar)

cols_kapp<-c("Correct_decod(mean)","Miss_class_cases","Kappa_stats")
cols<-c("true_average_switches(mean)","MHMM_average_switches(mean)","MHMM_Proportion_exact_switch_point(mean)","true_average_switches(sd)","MHMM_average_switches(sd)","MHMM_Proportion_exact_switch_point(sd)","true_average_switches(2.5%)","MHMM_average_switches(2.5%)","MHMM_Proportion_exact_switch_point(2.5%)","true_average_switches(97.5%)","MHMM_average_switches(97.5%)","MHMM_Proportion_exact_switch_point(97.5%)")
decoding_perc<-data.frame(mean_occ=rep(0,cores),miss_class=rep(0,cores), kappa_est=rep(0,cores))
switches<-data.frame(switch_t=rep(0,cores),switch_est=rep(0,cores),switch_pos_prop=rep(0,cores))
for(it in 1:cores){
  miss_class<-c()
  for(sub in 1:80){
    miss_class[sub]<-sum(out_list1[[it]]$decoding[which(out_list1[[it]]$decoding[,1]==sub),2]==out_list1[[it]]$decoding[which(out_list1[[it]]$decoding[,1]==sub),3])/length(out_list1[[it]]$decoding[which(out_list1[[it]]$decoding[,1]==sub),2])<0.2
  }
  decoding_perc$miss_class[it]<-sum(miss_class)
  #obtain the correct decoding
  decoding_perc$mean_occ[it]<-sum(out_list1[[it]]$decoding[,2]==out_list1[[it]]$decoding[,3])/length(out_list1[[it]]$decoding[,2])
  #count number of swiches in the data
  switch_t<-count_switches(out_list1[[it]]$decoding[,2])
  switch_est<-count_switches(out_list1[[it]]$decoding[,3])
  switches$switch_t[it]<-((switch_t$total_switches)-79)/80
  switches$switch_est[it]<-((switch_est$total_switches)-79)/80
  #can we distingush position of the switch
  switch_pos_t<-switch_t$switch_position
  switch_pos_est<-switch_est$switch_position
  if(length(switch_pos_t)==length(switch_pos_est)){
    switches$switch_pos_prop[it]<-(sum(switch_pos_t==switch_pos_est))/(length(switch_pos_t))
  }else{
    if(length(switch_pos_t)>length(switch_pos_est)){
      switch_pos_est<-c(switch_pos_est,rep(0,length(switch_pos_t)-length(switch_pos_est)))
    }else{
      switch_pos_t<-c(switch_pos_t,rep(0,length(switch_pos_est)-length(switch_pos_t)))
    }

  }

  #table needed for the kappa
  acc<-table(pred=factor(out_list1[[it]]$decoding[,3],levels = c(1:m)),true=factor(out_list1[[it]]$decoding[,2],levels = c(1:m)))
  decoding_perc$kappa_est[it]<-kappa(acc)
}

matrix_sd_out<-t(as.matrix(apply(decoding_perc,2,mean))) %>% as.data.frame()
colnames(matrix_sd_out)<-cols_kapp
matrix_switch<-t(as.matrix(c(apply(switches,2,mean),apply(switches,2,sd))))%>% as.data.frame()
cc_int<-apply(switches, 2, function(x){quantile(x,probs=c(0.025,0.975))}) %>% as.data.frame() %>% tibble::rownames_to_column(., "estimate") %>% pivot_longer(cols = c(2:4)) %>% mutate(estimate=paste0(name,estimate)) %>% dplyr::select(estimate,value)
names_ci<-cc_int$estimate
cc_int<-cc_int$value %>% t()
matrix_switch<-cbind(matrix_switch,cc_int)
colnames(matrix_switch)<-cols

  rownames(matrix_sd_out)<-namessd
  rownames(matrix_switch)<-namessd

  #set appropriate scenario classificators
  matrix_sd_out$model<-"mhmm"
  matrix_sd_out$state<-m
  matrix_sd_out$observations<-500
  matrix_sd_out$dwell_time<-"mix"
  matrix_switch$model<-"mhmm"
  matrix_switch$state<-m
  matrix_switch$observations<-500
  matrix_switch$dwell_time<-"mix"
  matrix_sd_out<-matrix_sd_out %>% select(4:7,1:3)
  matrix_switch<-matrix_switch %>% select(13:16,1:12)
#=====================================================================================================================================
#---emiss_mu_bar-------------------------
mean<-apply(emiss_mu_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(emiss_mu_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(emiss_mu_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(emiss_mu_bar[,1:m],2,mean) %>% rbind(.,true_emiss_mu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_mu_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(emiss_mu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_mu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_mu_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(emiss_mu_bar[,1:m],2,mean) %>% rbind(.,true_emiss_mu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_mu_bar"])*100/x["true_emiss_mu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(emiss_mu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_mu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_mu_bar"])*100/x["true_emiss_mu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)

cov<- emiss_mu_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(emiss_mu_bar[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_emiss_mu_bar[true_val] & true_emiss_mu_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
emiss_mu_bar_out<-cbind(true=true_emiss_mu_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("emiss_mu_bar_",1:m)
rownames(emiss_mu_bar_out)<-rows_names

#---emiss_mu_bar2-------------------------
mean<-apply(emiss_mu_bar2[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(emiss_mu_bar2[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(emiss_mu_bar2[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(emiss_mu_bar2[,1:m],2,mean) %>% rbind(.,true_emiss_mu_bar2) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_mu_bar2"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(emiss_mu_bar2[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_mu_bar2) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_mu_bar2"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(emiss_mu_bar2[,1:m],2,mean) %>% rbind(.,true_emiss_mu_bar2) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_mu_bar2"])*100/x["true_emiss_mu_bar2"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(emiss_mu_bar2[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_mu_bar2) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_mu_bar2"])*100/x["true_emiss_mu_bar2"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)

cov<- emiss_mu_bar2[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(emiss_mu_bar2[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_emiss_mu_bar2[true_val] & true_emiss_mu_bar2[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
emiss_mu_bar2_out<-cbind(true=true_emiss_mu_bar2,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("emiss_mu_bar2_",1:m)
rownames(emiss_mu_bar2_out)<-rows_names

#---emiss_var_bar------------------------
mean<-apply(emiss_var_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(emiss_var_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(emiss_var_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(emiss_var_bar[,1:m],2,mean) %>% rbind(.,true_emiss_var_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_var_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(emiss_var_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_var_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_var_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(emiss_var_bar[,1:m],2,mean) %>% rbind(.,true_emiss_var_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_var_bar"])*100/x["true_emiss_var_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(emiss_var_bar[,(m+1):(2*m)],2,mean)%>% rbind(.,true_emiss_var_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_var_bar"])*100/x["true_emiss_var_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov_need<-rep(true_emiss_var_bar,2)
cov<- emiss_var_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(emiss_var_bar[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_emiss_var_bar[true_val] & true_emiss_var_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
emiss_var_bar_out<-cbind(true=true_emiss_var_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("emiss_var_bar_",1:m)
rownames(emiss_var_bar_out)<-rows_names

#---emiss_var_bar2------------------------
mean<-apply(emiss_var_bar2[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(emiss_var_bar2[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(emiss_var_bar2[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(emiss_var_bar2[,1:m],2,mean) %>% rbind(.,true_emiss_var_bar2) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_var_bar2"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(emiss_var_bar2[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_var_bar2) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_var_bar2"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(emiss_var_bar2[,1:m],2,mean) %>% rbind(.,true_emiss_var_bar2) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_var_bar2"])*100/x["true_emiss_var_bar2"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(emiss_var_bar2[,(m+1):(2*m)],2,mean)%>% rbind(.,true_emiss_var_bar2) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_var_bar2"])*100/x["true_emiss_var_bar2"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)

cov<- emiss_var_bar2[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(emiss_var_bar2[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_emiss_var_bar2[true_val] & true_emiss_var_bar2[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
emiss_var_bar2_out<-cbind(true=true_emiss_var_bar2,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("emiss_var_bar2_",1:m)
rownames(emiss_var_bar2_out)<-rows_names

#---emiss_varmu_bar----------------------
mean<-apply(emiss_varmu_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(emiss_varmu_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(emiss_varmu_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(emiss_varmu_bar[,1:m],2,mean) %>% rbind(.,true_emiss_varmu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_varmu_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(emiss_varmu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_varmu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_varmu_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(emiss_varmu_bar[,1:m],2,mean) %>% rbind(.,true_emiss_varmu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_varmu_bar"])*100/x["true_emiss_varmu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(emiss_varmu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_varmu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_varmu_bar"])*100/x["true_emiss_varmu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)

cov<- emiss_varmu_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(emiss_varmu_bar[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_emiss_varmu_bar[true_val] & true_emiss_varmu_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
emiss_varmu_bar_out<-cbind(true=true_emiss_varmu_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("emiss_varmu_bar_",1:m)
rownames(emiss_varmu_bar_out)<-rows_names

#---emiss_varmu_bar2----------------------
mean<-apply(emiss_varmu_bar2[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(emiss_varmu_bar2[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(emiss_varmu_bar2[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(emiss_varmu_bar2[,1:m],2,mean) %>% rbind(.,true_emiss_varmu_bar2) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_varmu_bar2"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(emiss_varmu_bar2[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_varmu_bar2) %>%apply(.,2,FUN =function(x){x["."]-x["true_emiss_varmu_bar2"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(emiss_varmu_bar2[,1:m],2,mean) %>% rbind(.,true_emiss_varmu_bar2) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_varmu_bar2"])*100/x["true_emiss_varmu_bar2"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(emiss_varmu_bar2[,(m+1):(2*m)],2,mean) %>% rbind(.,true_emiss_varmu_bar2) %>%apply(.,2,FUN = function(x){(x["."]-x["true_emiss_varmu_bar2"])*100/x["true_emiss_varmu_bar2"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)

cov<- emiss_varmu_bar2[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(emiss_varmu_bar2[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_emiss_varmu_bar2[true_val] & true_emiss_varmu_bar2[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
emiss_varmu_bar2_out<-cbind(true=true_emiss_varmu_bar2,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("emiss_varmu_bar2_",1:m)
rownames(emiss_varmu_bar2_out)<-rows_names

#---gamma_prob_bar-----------------------

mean<-apply(gamma_prob_bar[,(m*m+1):(2*m*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(gamma_prob_bar[,1:(m*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(gamma_prob_bar[,(2*m*m+1):(3*m*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(gamma_prob_bar[,1:(m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_prob_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(gamma_prob_bar[,(m*m+1):(2*m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_prob_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(gamma_prob_bar[,1:(m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_prob_bar"])*100/x["true_gamma_prob_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(gamma_prob_bar[,(m*m+1):(2*m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_prob_bar"])*100/x["true_gamma_prob_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov<- gamma_prob_bar[,(3*m*m+1):(5*m*m)] %>% as.data.frame()%>% select(order(colnames(gamma_prob_bar[,(3*m*m+1):(5*m*m)])))
new_cov<-matrix(nrow=1,ncol=m*m)
for(true_val in 1:(m*m)){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_gamma_prob_bar[true_val] & true_gamma_prob_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
gamma_prob_bar_out<-cbind(true=true_gamma_prob_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-c()
for(i in 1:m){
  for(j in 1:m){
    rows_names<-c(rows_names,paste0("gamma_prob_bar_S",i,"toS",j))
  }
}
rownames(gamma_prob_bar_out)<-rows_names


#---gamma_prob_bar_no_diag-----------------------

mean<-apply(gamma_prob_bar_no_diag[,(m*m+1):(2*m*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(gamma_prob_bar_no_diag[,1:(m*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(gamma_prob_bar_no_diag[,(2*m*m+1):(3*m*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(gamma_prob_bar_no_diag[,1:(m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar_no_diag) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_prob_bar_no_diag"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(gamma_prob_bar_no_diag[,(m*m+1):(2*m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar_no_diag) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_prob_bar_no_diag"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(gamma_prob_bar_no_diag[,1:(m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar_no_diag) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_prob_bar_no_diag"])*100/x["true_gamma_prob_bar_no_diag"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(gamma_prob_bar_no_diag[,(m*m+1):(2*m*m)],2,mean) %>% rbind(.,true_gamma_prob_bar_no_diag) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_prob_bar_no_diag"])*100/x["true_gamma_prob_bar_no_diag"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov<- gamma_prob_bar_no_diag[,(3*m*m+1):(5*m*m)] %>% as.data.frame()%>% select(order(colnames(gamma_prob_bar_no_diag[,(3*m*m+1):(5*m*m)])))
new_cov<-matrix(nrow=1,ncol=m*m)
for(true_val in 1:(m*m)){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_gamma_prob_bar_no_diag[true_val] & true_gamma_prob_bar_no_diag[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
gamma_prob_bar_no_diag_out<-cbind(true=true_gamma_prob_bar_no_diag,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-c()
for(i in 1:m){
  for(j in 1:m){
    rows_names<-c(rows_names,paste0("gamma_prob_bar_S",i,"toS",j,"_no_diag"))
  }
}
rownames(gamma_prob_bar_no_diag_out)<-rows_names

#---gamma_V_int_bar-----------------------
mean<-apply(gamma_V_int_bar[,((m*(m-1))+1):(2*((m*(m-1))))],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(gamma_V_int_bar[,1:((m*(m-1)))],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(gamma_V_int_bar[,(2*((m*(m-1)))+1):(3*((m*(m-1))))],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(gamma_V_int_bar[,1:((m*(m-1)))],2,mean) %>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_V_int_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(gamma_V_int_bar[,((m*(m-1))+1):(2*(m*(m-1)))],2,mean)%>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_V_int_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)

bias_median_pp<-apply(gamma_V_int_bar[,1:((m*(m-1)))],2,mean) %>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_V_int_bar"])*100/x["true_gamma_V_int_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(gamma_V_int_bar[,((m*(m-1))+1):(2*(m*(m-1)))],2,mean) %>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_V_int_bar"])*100/x["true_gamma_V_int_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov_need<-rep(true_gamma_V_int_bar,2)
cov<- gamma_V_int_bar[,(3*(m*(m-1))+1):(5*(m*(m-1)))] %>% as.data.frame()%>% select(order(colnames(gamma_V_int_bar[,(3*(m*(m-1))+1):(5*(m*(m-1)))])))
new_cov<-matrix(nrow=1,ncol=(m*(m-1)))
for(true_val in 1:(m*(m-1))){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_gamma_V_int_bar[true_val] & true_gamma_V_int_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
gamma_V_int_bar_out<-cbind(true=true_gamma_V_int_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)

#---dwell_mu_bar-------------------------
mean<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(dwell_mu_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(dwell_mu_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(dwell_mu_bar[,1:m],2,mean) %>% rbind(.,true_dwell) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_dwell) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(dwell_mu_bar[,1:m],2,mean) %>% rbind(.,true_dwell) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell"])*100/x["true_dwell"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_dwell) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell"])*100/x["true_dwell"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov_need<-rep(true_dwell,2)
cov<- dwell_mu_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(dwell_mu_bar[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_dwell[true_val] & true_dwell[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
dwell_mu_bar_out<-cbind(true=true_dwell,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("dwell_mu_bar_",1:m)
rownames(dwell_mu_bar_out)<-rows_names


#binding outputs
post_proc_out<-rbind(emiss_mu_bar_out,emiss_mu_bar2_out,emiss_var_bar_out,emiss_var_bar2_out,emiss_varmu_bar_out,emiss_varmu_bar2_out,gamma_prob_bar_out,gamma_prob_bar_no_diag_out, gamma_V_int_bar_out,dwell_mu_bar_out) %>% round(.,3)
post_proc_out$estimates <- rownames(post_proc_out)
post_proc_out$model<-rep("mhmm",nrow(post_proc_out))
post_proc_out$state<-rep(m,nrow(post_proc_out))
post_proc_out$observations<-rep(observation_lenght[o],nrow(post_proc_out))
post_proc_out$dwell_time<-rep("mix",nrow(post_proc_out))
rownames(post_proc_out)<-NULL
post_proc_out<-post_proc_out %>% select(11:14,10,1:9)

out_list2<-list(MAP_summary=all_map, decod_acc=matrix_sd_out, switches_est=matrix_switch,empirical_dwell_time_true=empirical_dwell_time_true,empirical_dwell_time_mhmm=empirical_dwell_time_mhmm, post_proc_est=post_proc_out)
#In order to run the scripts without any issues using snellius we need to adjust the below target directory to save outputs
saveRDS(out_list2, paste0("./mhmm_run/outputs/post_processed/mhmm_post_proc_output",
                          "_m",m,
                          "_nt",o,
                          "_dwell",d,
                          "_it",J,
                          "_burn_in",B,
                          ".rds"))

#If set of parameters is the last from parameter list hich is scenario that we consider, parse the outputs in order to obtain one table at the end
m <- 3:4
o <- 1:3
d <- 1:4
# Calling expand.grid() Function to create the parameter dataframe
pars = expand.grid(d=d,o=o,m=m)
pars<-rbind(pars,c(5,2,3))
n_dep=2
list_post_est<-list()
list_decod_acc<-list()
list_switches_est<-list()
list_dwell_emp<-list()
dwell_time_names<-c(1.4,3.5,19.5,99.5,"mix")
obs_names<-c(200,500,1000)
for(parameter in 1:nrow(pars)){
  m  <- pars[parameter,3]
  o  <- pars[parameter,2]
  d  <- pars[parameter,1]
  #In order to run the scripts without any issues using snellius we need to adjust the below data directory
  post_data<-readRDS(paste0("./mhmm_run/post_processed/mhmm_post_proc_output",
                            "_m",m,
                            "_nt",o,
                            "_dwell",d,
                            "_it",J,
                            "_burn_in",B,
                            ".rds"))
  list_post_est[[parameter]]<-post_data$post_proc_est
  list_decod_acc[[parameter]]<-post_data$decod_acc
  list_switches_est[[parameter]]<-post_data$switches_est

  state_column<-data.frame(state=rep(m,m))
  #create columns for true values
  true_cols<-data.frame(mean_emp_dwell_true=post_data$empirical_dwell_time_true$mean_emp_dwell,median_emp_dwell_true=post_data$empirical_dwell_time_true$median_emp_dwell, ci_025_emp_dwell_true=post_data$empirical_dwell_time_true$ci_025,ci_975_emp_dwell_true=post_data$empirical_dwell_time_true$ci_975)

  emp_dwell<-post_data$empirical_dwell_time_mhmm
  emp_dwell<-cbind(emp_dwell[,1],state_column,emp_dwell[,3:4],true_cols[,1:2],emp_dwell[,5:ncol(emp_dwell)],true_cols[,3:4])
  colnames(emp_dwell)<-c("model","state","observations","dwell_time", "mean_emp_dwell_true","median_emp_dwell_true","estimates","mean_emp_dwell","median_emp_dwell", "sd_emp_dwell","ci_025","ci_975","ci_025_true","ci_975_true")
  emp_dwell$estimates<-paste0("emp_dwell",1:nrow(post_data$empirical_dwell_time_true))
  emp_dwell$dwell_time<-rep(dwell_time_names[d],m)
  emp_dwell$observations<-rep(obs_names[o],m)
  emp_dwell<-emp_dwell[,c(1:4,7,5:6,8:14)]
  list_dwell_emp[[parameter]]<-emp_dwell
}
list_post_est<-do.call(rbind,list_post_est)
list_decod_acc<-do.call(rbind,list_decod_acc)
list_switches_est<-do.call(rbind,list_switches_est)
list_dwell_emp<-do.call(rbind,list_dwell_emp)
#In order to run the scripts without any issues using snellius we need to adjust the below target directory
saveRDS(list_post_est, "./mhmm_run/post_processed/list_post_est.rds")
saveRDS(list_decod_acc, "./mhmm_run/post_processed/list_decod_acc.rds")
saveRDS(list_switches_est, "./mhmm_run/post_processed/list_switches_est.rds")
saveRDS(list_dwell_emp, "./mhmm_run/post_processed/list_dwell_emp.rds")

