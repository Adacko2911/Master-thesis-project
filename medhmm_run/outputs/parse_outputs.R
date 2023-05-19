#set global variables
exp_mean1=1/(-log(0.5))
exp_mean2=1/(-log(0.75))
exp_mean3=1/(-log(0.95))
exp_mean4=1/(-log(0.99))
J=4000
B=2000
#define standard deviations
exp_sigma1=1.10093 #such that exp_sigma1=1/3*exp_mean1
exp_sigma2=1.10093 #such that exp_sigma2=1/3*exp_mean2
exp_sigma3=1.10093 #such that exp_sigma3=1/3*exp_mean3
exp_sigma4=1.05902 #such that exp_sigma4=1/4*exp_mean4
#read data
m <- 3:4
o <- 1:3
d <- 1:4
# Calling expand.grid() Function to create the parameter data frame
pars = expand.grid(m=m,o=o,d=d)
#the outputs indicators exucluding the 2 missing scenarios
pars<-rbind(pars,c(3,2,5))
pars<-pars[c(13:21,23:25),]

#start the loop that go over the outputs and parse them all once again ( this process could be omitted in case snelius code would run successfully)###################
for(s in 1:nrow(pars)){
  pars2<-pars[s,]
  m  <- as.numeric(pars2[1])
o  <- as.numeric(pars2[2])
d  <- as.numeric(pars2[3])

out_list1<-readRDS(paste0("./medhmm_run/outputs/medhmm_list_output",
                          "_m",m,
                          "_nt",o,
                          "_dwell",d,
                          "_it",J,
                          "_burn_in",B,
                          ".rds"))

namessd<-c()
namessd<-append(namessd,paste0("st_",m,"_obs_",o,"_dwell_",d,"_medhmm"))
dwell_time_names<-c(as.character(round(exp_mean1,2)),as.character(round(exp_mean2,2)),as.character(round(exp_mean3,2)),as.character(round(exp_mean4,2)),"mixed")
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

if(m==3){
  true_emiss_varmu_bar <- rep(9,m)
  true_emiss_var_bar<-c(20,60,120)
  true_emiss_mu_bar<-c(10,30,60)
  true_emiss_varmu_bar2 <- rep(9,m)
  true_emiss_var_bar2<-c(120+0.25*120,60+0.25*60,20+20*0.25)
  true_emiss_mu_bar2<-c(60,30,10)
  if(d==1){
    true_dwell_mu_bar<-rep(log(exp_mean1),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma1),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                            0.5,0,0.5,
                                            0.5,0.5,0), ncol = m, byrow = TRUE))
  }else if(d==2){
    true_dwell_mu_bar<-rep(log(exp_mean2),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma2),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                            0.5,0,0.5,
                                            0.5,0.5,0), ncol = m, byrow = TRUE))
  }else if(d==3){
    true_dwell_mu_bar<-rep(log(exp_mean3),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma3),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                            0.5,0,0.5,
                                            0.5,0.5,0), ncol = m, byrow = TRUE))
  }else if(d==4){
    true_dwell_mu_bar<-rep(log(exp_mean4),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma4),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                            0.5,0,0.5,
                                            0.5,0.5,0), ncol = m, byrow = TRUE))
  }else if(d==5){
    true_dwell_mu_bar<-c(log(exp_mean2),log(exp_mean3),log(exp_mean4))
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma2),log(exp_sigma3),log(exp_sigma4))
    true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                            0.5,0,0.5,
                                            0.5,0.5,0), ncol = m, byrow = TRUE))
  }

  true_gamma_V_int_bar<-rep(0.17,((m)*(m-2)))
}else if(m==4){
  true_emiss_varmu_bar <- rep(9,m)
  true_emiss_var_bar<-c(20,60,120,150)
  true_emiss_mu_bar<-c(10,30,60,98)
  true_emiss_varmu_bar2 <- rep(9,m)
  true_emiss_var_bar2<-c((150+0.25*150),(120+0.25*120),(60+0.25*60),(20+0.25*20))
  true_emiss_mu_bar2<-c(98,60,30,10)
  if(d==1){
    true_dwell_mu_bar<-rep(log(exp_mean1),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma1),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 1/3,1/3,1/3,
                                            1/3, 0,1/3,1/3,
                                            1/3,1/3,0,1/3,
                                            1/3,1/3,1/3,0), ncol = m, byrow = TRUE))
  }else if(d==2){
    true_dwell_mu_bar<-rep(log(exp_mean2),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma2),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 1/3,1/3,1/3,
                                            1/3, 0,1/3,1/3,
                                            1/3,1/3,0,1/3,
                                            1/3,1/3,1/3,0), ncol = m, byrow = TRUE))
  }else if(d==3){
    true_dwell_mu_bar<-rep(log(exp_mean3),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma3),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 1/3,1/3,1/3,
                                            1/3, 0,1/3,1/3,
                                            1/3,1/3,0,1/3,
                                            1/3,1/3,1/3,0), ncol = m, byrow = TRUE))
  }else if(d==4){
    true_dwell_mu_bar<-rep(log(exp_mean4),m)
    true_dwell_varmu_bar<-rep(0.2,m)
    true_dwell_var_bar<-rep(log(exp_sigma4),m)
    true_gamma_prob_bar<-as.vector(matrix(c(0, 1/3,1/3,1/3,
                                            1/3, 0,1/3,1/3,
                                            1/3,1/3,0,1/3,
                                            1/3,1/3,1/3,0), ncol = m, byrow = TRUE))
  }
  true_gamma_V_int_bar<-rep(0.17,((m)*(m-2)))

}

#post processing lists to obtaing the parameters estimates
emiss_mu_bar<-list()
emiss_mu_bar2<-list()
emiss_var_bar<-list()
emiss_var_bar2<-list()
emiss_varmu_bar<-list()
emiss_varmu_bar2<-list()
dwell_mu_bar<-list()
dwell_var_bar<-list()
dwell_varmu_bar<-list()
gamma_prob_bar<-list()
gamma_int_bar<-list()
gamma_V_int_bar<-list()
#prepare list for group level decoding list
decoding_out_list<-list()
decoding_out_list_truth<-list()
out_list1[sapply(out_list1, is.list)==FALSE] <- NULL
cores<-length(out_list1)
for (i in 1:cores){
  emiss_mu_bar[[i]]<-out_list1[[i]]$MAP$emiss_mu_bar[[1]]
  emiss_var_bar[[i]]<-out_list1[[i]]$MAP$emiss_var_bar[[1]]
  emiss_varmu_bar[[i]]<-out_list1[[i]]$MAP$emiss_varmu_bar[[1]]
  emiss_mu_bar2[[i]]<-out_list1[[i]]$MAP$emiss_mu_bar[[2]]
  emiss_var_bar2[[i]]<-out_list1[[i]]$MAP$emiss_var_bar[[2]]
  emiss_varmu_bar2[[i]]<-out_list1[[i]]$MAP$emiss_varmu_bar[[2]]
  dwell_mu_bar[[i]]<-out_list1[[i]]$MAP$dwell_mu_bar
  dwell_var_bar[[i]]<-out_list1[[i]]$MAP$dwell_var_bar
  dwell_varmu_bar[[i]]<-out_list1[[i]]$MAP$dwell_varmu_bar
  gamma_prob_bar[[i]]<-try(out_list1[[i]]$MAP$gamma_prob_bar)
  gamma_int_bar[[i]]<-try(out_list1[[i]]$MAP$gamma_int_bar)
  gamma_V_int_bar[[i]]<-try(out_list1[[i]]$MAP$gamma_V_int_bar)

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
dwell_mu_bar<-do.call(rbind,dwell_mu_bar)
dwell_var_bar<-do.call(rbind,dwell_var_bar)
dwell_varmu_bar<-do.call(rbind,dwell_varmu_bar)
if(m>2){
  gamma_prob_bar<-do.call(rbind,gamma_prob_bar)
  gamma_int_bar<-do.call(rbind,gamma_int_bar)
  gamma_V_int_bar<-do.call(rbind,gamma_V_int_bar)
  all_map<-list(emiss_mu_bar=emiss_mu_bar,emiss_var_bar=emiss_var_bar,emiss_varmu_bar=emiss_varmu_bar,
                emiss_mu_bar2=emiss_mu_bar2,emiss_var_bar2=emiss_var_bar2,emiss_varmu_bar2=emiss_varmu_bar2,
                dwell_mu_bar=dwell_mu_bar,dwell_var_bar=dwell_var_bar,dwell_varmu_bar=dwell_varmu_bar,gamma_prob_bar=gamma_prob_bar,gamma_int_bar=gamma_int_bar,gamma_V_int_bar=gamma_V_int_bar)

}else{
  all_map<-list(emiss_mu_bar=emiss_mu_bar,emiss_var_bar=emiss_var_bar,emiss_varmu_bar=emiss_varmu_bar,
                emiss_mu_bar2=emiss_mu_bar2,emiss_var_bar2=emiss_var_bar2,emiss_varmu_bar2=emiss_varmu_bar2,
                dwell_mu_bar=dwell_mu_bar,dwell_var_bar=dwell_var_bar,dwell_varmu_bar=dwell_varmu_bar)

}



#joint the out results form each simulation list and calculate the out results for each simulation list
empirical_dwell_time_true<-do.call(rbind,decoding_out_list_truth) %>%
  as.data.frame() %>%group_by(point) %>% group_modify(~ apply(.x,2,mean) %>%t() %>%
                                                        as.data.frame()) %>% rename("state_scenario" = "point") %>%
  as.data.frame() %>% mutate(dwell_time=rep(d,m), .before = state_scenario) %>%
  mutate(observations=rep(o,m), .before = dwell_time)  %>%
  mutate(state=rep(m,m), .before = observations) %>%mutate(model=rep("medhmm",m), .before = state)

empirical_dwell_time_medhmm<-do.call(rbind,decoding_out_list) %>% as.data.frame()%>%
  group_by(point) %>% group_modify(~ apply(.x,2,mean) %>%t() %>% as.data.frame())%>%
  rename("state_scenario" = "point")%>% as.data.frame() %>% mutate(dwell_time=rep(d,m), .before = state_scenario) %>%
  mutate(observations=rep(o,m), .before = dwell_time)  %>%  mutate(state=rep(m,m), .before = observations) %>%
  mutate(model=rep("medhmm",m), .before = state)



cols_kapp<-c("Correct_decod(mean)","Miss_class_cases","Kappa_stats")
cols<-c("true_average_switches(mean)","MEDHMM_average_switches(mean)","MEDHMM_Proportion_exact_switch_point(mean)","true_average_switches(sd)","MEDHMM_average_switches(sd)","MEDHMM_Proportion_exact_switch_point(sd)","true_average_switches(2.5%)","MEDHMM_average_switches(2.5%)","MEDHMM_Proportion_exact_switch_point(2.5%)","true_average_switches(97.5%)","MEDHMM_average_switches(97.5%)","MEDHMM_Proportion_exact_switch_point(97.5%)")
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
#====================================================================================================================================================================================


#set appropriate scenario classificators

matrix_sd_out$model<-"medhmm"
matrix_sd_out$state<-m
matrix_sd_out$observations<-observation_lenght[o]
matrix_sd_out$dwell_time<-dwell_time_names[d]
matrix_switch$model<-"medhmm"
matrix_switch$state<-m
matrix_switch$observations<-observation_lenght[o]
matrix_switch$dwell_time<-dwell_time_names[d]
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
cov_need<-rep(true_emiss_mu_bar,2)
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
cov_need<-rep(true_emiss_mu_bar2,2)
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
cov_need<-rep(true_emiss_varmu_bar,2)
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
cov_need<-rep(true_emiss_varmu_bar2,2)
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

#---dwell_mu_bar-------------------------
mean<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(dwell_mu_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(dwell_mu_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(dwell_mu_bar[,1:m],2,mean) %>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_mu_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_mu_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(dwell_mu_bar[,1:m],2,mean) %>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_mu_bar"])*100/x["true_dwell_mu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_mu_bar"])*100/x["true_dwell_mu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov_need<-rep(true_dwell_mu_bar,2)
cov<- dwell_mu_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(dwell_mu_bar[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_dwell_mu_bar[true_val] & true_dwell_mu_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
dwell_mu_bar_out<-cbind(true=true_dwell_mu_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("dwell_mu_bar_",1:m)
rownames(dwell_mu_bar_out)<-rows_names

#---dwell_var_bar------------------------
mean<-apply(dwell_var_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(dwell_var_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(dwell_var_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(dwell_var_bar[,1:m],2,mean) %>% rbind(.,true_dwell_var_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_var_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(dwell_var_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_dwell_var_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_var_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(dwell_var_bar[,1:m],2,mean) %>% rbind(.,true_dwell_var_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_var_bar"])*100/x["true_dwell_var_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(dwell_var_bar[,(m+1):(2*m)],2,mean)%>% rbind(.,true_dwell_var_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_var_bar"])*100/x["true_dwell_var_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov_need<-rep(true_dwell_var_bar,2)
cov<- dwell_var_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(dwell_var_bar[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_dwell_var_bar[true_val] & true_dwell_var_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
dwell_var_bar_out<-cbind(true=true_dwell_var_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("dwell_var_bar_",1:m)
rownames(dwell_var_bar_out)<-rows_names

#---dwell_varmu_bar----------------------
mean<-apply(dwell_varmu_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
median<-apply(dwell_varmu_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
se<-apply(dwell_varmu_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
bias_median<-apply(dwell_varmu_bar[,1:m],2,mean) %>% rbind(.,true_dwell_varmu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_varmu_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
bias_mean<-apply(dwell_varmu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_dwell_varmu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_varmu_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
bias_median_pp<-apply(dwell_varmu_bar[,1:m],2,mean) %>% rbind(.,true_dwell_varmu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_varmu_bar"])*100/x["true_dwell_varmu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
bias_mean_pp<-apply(dwell_varmu_bar[,(m+1):(2*m)],2,mean) %>% rbind(.,true_dwell_varmu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_varmu_bar"])*100/x["true_dwell_varmu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
cov_need<-rep(true_dwell_varmu_bar,2)
cov<- dwell_varmu_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(dwell_varmu_bar[,(3*m+1):(5*m)])))
new_cov<-matrix(nrow=1,ncol=m)
for(true_val in 1:m){
  new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
    if(x[1]<=true_dwell_varmu_bar[true_val] & true_dwell_varmu_bar[true_val]<=x[2]){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }) %>% sum(.)*100/cores
}
new_cov<-new_cov %>% round(.,2)%>% t() %>% as.data.frame() %>% rename(., coverage_pp =V1)
dwell_varmu_bar_out<-cbind(true=true_dwell_varmu_bar,mean,median,se,bias_mean,bias_mean_pp, bias_median,bias_median_pp, new_cov)
rows_names<-paste0("dwell_varmu_bar_",1:m)
rownames(dwell_varmu_bar_out)<-rows_names


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

  #---gamma_V_int_bar-----------------------
  mean<-apply(gamma_V_int_bar[,((m)*(m-2)+1):(2*((m)*(m-2)))],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.)
  median<-apply(gamma_V_int_bar[,1:((m)*(m-2))],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.)
  se<-apply(gamma_V_int_bar[,(2*((m)*(m-2))+1):(3*((m)*(m-2)))],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
  bias_median<-apply(gamma_V_int_bar[,1:((m)*(m-2))],2,mean) %>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_V_int_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
  bias_mean<-apply(gamma_V_int_bar[,((m)*(m-2)+1):(2*((m)*(m-2)))],2,mean)%>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_gamma_V_int_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
  bias_median_pp<-apply(gamma_V_int_bar[,1:((m)*(m-2))],2,mean) %>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_V_int_bar"])*100/x["true_gamma_V_int_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
  bias_mean_pp<-apply(gamma_V_int_bar[,((m)*(m-2)+1):(2*((m)*(m-2)))],2,mean) %>% rbind(.,true_gamma_V_int_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_gamma_V_int_bar"])*100/x["true_gamma_V_int_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
  cov_need<-rep(true_gamma_V_int_bar,2)
  cov<- gamma_V_int_bar[,(3*((m)*(m-2))+1):(5*((m)*(m-2)))] %>% as.data.frame()%>% select(order(colnames(gamma_V_int_bar[,(3*((m)*(m-2))+1):(5*((m)*(m-2)))])))
  new_cov<-matrix(nrow=1,ncol=((m)*(m-2)))
  for(true_val in 1:((m)*(m-2))){
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

  all_st<-c()
  for(st_cur in 1:m){
    to_st<-2:m
    to_st<-to_st[which(to_st!=st_cur)]
    all_st<-append(all_st,paste0("S",st_cur,"toS",to_st))
  }
  rows_names<-c()
  for(al in 1:length(all_st)){
    rows_names<-append(rows_names,paste0("gamma_V_int_",all_st[al],"_with_int_",all_st[al]))
  }

  rownames(gamma_V_int_bar_out)<-rows_names[-1]

  post_proc_out<-rbind(emiss_mu_bar_out,emiss_mu_bar2_out,emiss_var_bar_out,emiss_var_bar2_out,emiss_varmu_bar_out,emiss_varmu_bar2_out,dwell_mu_bar_out,dwell_var_bar_out,dwell_varmu_bar_out,gamma_prob_bar_out,gamma_V_int_bar_out) %>% round(.,1)
  post_proc_out$estimates <- rownames(post_proc_out)
  post_proc_out$model<-rep("medhmm",nrow(post_proc_out))
  post_proc_out$state<-rep(m,nrow(post_proc_out))
  post_proc_out$observations<-rep(observation_lenght[o],nrow(post_proc_out))
  post_proc_out$dwell_time<-rep(dwell_time_names[d],nrow(post_proc_out))
  rownames(post_proc_out)<-NULL
  post_proc_out<-post_proc_out %>% select(11:14,10,1:9)

  out_list2<-list(MAP_summary=all_map, post_proc_est=post_proc_out,decod_acc=matrix_sd_out, switches_est=matrix_switch,empirical_dwell_time_true=empirical_dwell_time_true,empirical_dwell_time_medhmm=empirical_dwell_time_medhmm)
  saveRDS(out_list2, paste0("./medhmm_run/outputs/post_processed/medhmm_post_proc_output",
                            "_m",m,
                            "_nt",o,
                            "_dwell",d,
                            "_it",J,
                            "_burn_in",B,
                            ".rds"))


}
J=4000
B=2000
#If set of parameters is the last from parameter list, parse the outputs in order to obtain one table at the end
#read data
m <- 3:4
o <- 1:3
d <- 1:4
# Calling expand.grid() Function to create the parameter dataframe
pars = expand.grid(d=d,o=o,m=m)
pars<-rbind(pars,c(5,2,3))
pars<-pars[c(1:19,21,23:25),] #this need to be added when the oter two runs will succesfully manage
  n_dep=2
  list_post_est<-list()
  list_decod_acc<-list()
  list_switches_est<-list()
  list_dwell_emp<-list()
  dwell_time_names<-c(1.4,3.5,19.5,99.5,"mix")
  obs_names<-c(200,500,1000)

  #parse post processed outputs and include them in one data frame
   for(parameter in 1:nrow(pars)){
    m  <- pars[parameter,3]
    o  <- pars[parameter,2]
    d  <- pars[parameter,1]

    post_data<-readRDS(paste0("./medhmm_run/outputs/post_processed/medhmm_post_proc_output",
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

    emp_dwell<-post_data$empirical_dwell_time_medhmm
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
  #Save all outputs
  saveRDS(list_post_est, "./medhmm_run/outputs/post_processed/list_post_est.rds")
  saveRDS(list_decod_acc, "./medhmm_run/outputs/post_processed/list_decod_acc.rds")
  saveRDS(list_switches_est,"./medhmm_run/outputs/post_processed/list_switches_est.rds")
  saveRDS(list_dwell_emp, "./medhmm_run/outputs/post_processed/list_dwell_emp.rds")


