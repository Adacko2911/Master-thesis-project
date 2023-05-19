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
# Calling expand.grid() Function to create the parameter dataframe
pars = expand.grid(m=m,o=o,d=d)
pars<-rbind(pars,c(3,2,5))
pars<-pars[c(1:11,13:21,23:25),]

for(s in 1:nrow(pars)){
  pars2<-pars[s,]
  m  <- as.numeric(pars2[1])
  o  <- as.numeric(pars2[2])
  d  <- as.numeric(pars2[3])
  #In order to run the scripts without any issues using snellius we need to adjust the below data directory
  out_list1<-try(readRDS(paste0("./medhmm_run/outputs/medhmm_list_output",
                                "_m",m,
                                "_nt",o,
                                "_dwell",d,
                                "_it",J,
                                "_burn_in",B,
                                ".rds")))

  namessd<-c()
  namessd<-append(namessd,paste0("st_",m,"_obs_",o,"_dwell_",d,"_medhmm"))
  dwell_time_names<-c(1.4,3.5,19.5,99.5,"mix")
  observation_lenght<-c(200,500,1000)


  if(m==3){
    true_emiss_varmu_bar <- rep(9,m)
    true_emiss_var_bar<-c(20,60,120)
    true_emiss_mu_bar<-c(10,30,60)
    true_emiss_varmu_bar2 <- rep(9,m)
    true_emiss_var_bar2<-c(120+0.25*120,60+0.25*60,20+20*0.25)
    true_emiss_mu_bar2<-c(60,30,10)
    if(d==1){
      true_dwell_mu_bar<-rep((exp_mean1),m)
      true_dwell_varmu_bar<-rep(0.2,m)
      true_dwell_var_bar<-rep(log(exp_sigma1),m)
      true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                              0.5,0,0.5,
                                              0.5,0.5,0), ncol = m, byrow = TRUE))
    }else if(d==2){
      true_dwell_mu_bar<-rep((exp_mean2),m)
      true_dwell_varmu_bar<-rep(0.2,m)
      true_dwell_var_bar<-rep(log(exp_sigma2),m)
      true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                              0.5,0,0.5,
                                              0.5,0.5,0), ncol = m, byrow = TRUE))
    }else if(d==3){
      true_dwell_mu_bar<-rep((exp_mean3),m)
      true_dwell_varmu_bar<-rep(0.2,m)
      true_dwell_var_bar<-rep(log(exp_sigma3),m)
      true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                              0.5,0,0.5,
                                              0.5,0.5,0), ncol = m, byrow = TRUE))
    }else if(d==4){
      true_dwell_mu_bar<-rep((exp_mean4),m)
      true_dwell_varmu_bar<-rep(0.2,m)
      true_dwell_var_bar<-rep(log(exp_sigma4),m)
      true_gamma_prob_bar<-as.vector(matrix(c(0, 0.5,0.5,
                                              0.5,0,0.5,
                                              0.5,0.5,0), ncol = m, byrow = TRUE))
    }else if(d==5){
      true_dwell_mu_bar<-c((exp_mean2),(exp_mean3),(exp_mean4))
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
      true_dwell_mu_bar<-rep((exp_mean1),m)
      true_dwell_varmu_bar<-rep(0.2,m)
      true_dwell_var_bar<-rep(log(exp_sigma1),m)
      true_gamma_prob_bar<-as.vector(matrix(c(0, 1/3,1/3,1/3,
                                              1/3, 0,1/3,1/3,
                                              1/3,1/3,0,1/3,
                                              1/3,1/3,1/3,0), ncol = m, byrow = TRUE))
    }else if(d==2){
      true_dwell_mu_bar<-rep((exp_mean2),m)
      true_dwell_varmu_bar<-rep(0.2,m)
      true_dwell_var_bar<-rep(log(exp_sigma2),m)
      true_gamma_prob_bar<-as.vector(matrix(c(0, 1/3,1/3,1/3,
                                              1/3, 0,1/3,1/3,
                                              1/3,1/3,0,1/3,
                                              1/3,1/3,1/3,0), ncol = m, byrow = TRUE))
    }else if(d==3){
      true_dwell_mu_bar<-rep((exp_mean3),m)
      true_dwell_varmu_bar<-rep(0.2,m)
      true_dwell_var_bar<-rep(log(exp_sigma3),m)
      true_gamma_prob_bar<-as.vector(matrix(c(0, 1/3,1/3,1/3,
                                              1/3, 0,1/3,1/3,
                                              1/3,1/3,0,1/3,
                                              1/3,1/3,1/3,0), ncol = m, byrow = TRUE))
    }else if(d==4){
      true_dwell_mu_bar<-rep((exp_mean4),m)
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
  dwell_mu_bar<-list()
  dwell_var_bar<-list()
  gamma_prob_bar<-list()
  out_list1[sapply(out_list1, is.list)==FALSE] <- NULL
  cores<-length(out_list1)
  for (i in 1:cores){

    dwell_mu_bar[[i]]<-out_list1[[i]]$MAP$dwell_mu_bar
    dwell_var_bar[[i]]<-out_list1[[i]]$MAP$dwell_var_bar
    gamma_prob_bar[[i]]<-try(out_list1[[i]]$MAP$gamma_prob_bar)


  }

  dwell_mu_bar<-do.call(rbind,dwell_mu_bar)
  dwell_var_bar<-do.call(rbind,dwell_var_bar)
  gamma_prob_bar<-do.call(rbind,gamma_prob_bar)



  #---dwell_mu_bar-------------------------
  mean<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., mean =.) %>% exp(.)
  median<-apply(dwell_mu_bar[,1:m],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., median =.) %>% exp(.)
  se<-apply(dwell_mu_bar[,(2*m+1):(3*m)],2,mean)%>% round(.,2)%>% as.data.frame() %>% rename(., sd =.)
  bias_median<-apply(dwell_mu_bar[,1:m],2,mean) %>% exp(.) %>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_mu_bar"]}) %>% round(.,2) %>% as.data.frame() %>% rename(., bias_median =.)
  bias_mean<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean) %>% exp(.) %>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN =function(x){x["."]-x["true_dwell_mu_bar"]}) %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean =.)
  bias_median_pp<-apply(dwell_mu_bar[,1:m],2,mean) %>% exp(.)%>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_mu_bar"])*100/x["true_dwell_mu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_median_pp =.)
  bias_mean_pp<-apply(dwell_mu_bar[,(m+1):(2*m)],2,mean) %>% exp(.) %>% rbind(.,true_dwell_mu_bar) %>%apply(.,2,FUN = function(x){(x["."]-x["true_dwell_mu_bar"])*100/x["true_dwell_mu_bar"]})%>% abs() %>% round(.,2)%>% as.data.frame() %>% rename(., bias_mean_pp =.)
  cov_need<-rep(true_dwell_mu_bar,2)
  cov<- dwell_mu_bar[,(3*m+1):(5*m)] %>% as.data.frame()%>% select(order(colnames(dwell_mu_bar[,(3*m+1):(5*m)])))
  new_cov<-matrix(nrow=1,ncol=m)
  for(true_val in 1:m){
    new_cov[true_val]<-cov[,(2*true_val-1):(2*true_val)] %>% apply(.,1,FUN = function(x){
      if(x[1]<=log(true_dwell_mu_bar[true_val]) & log(true_dwell_mu_bar[true_val])<=x[2]){
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



  post_proc_out<-rbind(dwell_mu_bar_out,gamma_prob_bar_out) %>% round(.,2)
  post_proc_out$estimates <- rownames(post_proc_out)
  post_proc_out$model<-rep("medhmm",nrow(post_proc_out))
  post_proc_out$state<-rep(m,nrow(post_proc_out))
  post_proc_out$observations<-rep(observation_lenght[o],nrow(post_proc_out))
  post_proc_out$dwell_time<-rep(dwell_time_names[d],nrow(post_proc_out))
  rownames(post_proc_out)<-NULL
  post_proc_out<-post_proc_out %>% select(11:14,10,1:9)
  #In order to run the scripts without any issues using snellius we need to adjust the target data directory
  write.csv(post_proc_out, paste0("./medhmm_run/outputs/post_processed/medhmm_post_proc_dwell_gamma_output",
                                  "_m",m,
                                  "_nt",o,
                                  "_dwell",d,
                                  "_it",J,
                                  "_burn_in",B,
                                  ".csv"))


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
for(parameter in 1:nrow(pars)){
  m  <- pars[parameter,3]
  o  <- pars[parameter,2]
  d  <- pars[parameter,1]
  #In order to run the scripts without any issues using snellius we need to adjust the below data directory
  post_data<-read.csv( paste0("./medhmm_run/outputs/post_processed/medhmm_post_proc_dwell_gamma_output",
                              "_m",m,
                              "_nt",o,
                              "_dwell",d,
                              "_it",J,
                              "_burn_in",B,
                              ".csv"))
  list_post_est[[parameter]]<-post_data

}
list_post_est<-do.call(rbind,list_post_est)
write.csv(list_post_est, "./medhmm_run/outputs/post_processed/list_post_est_dwell_gamma.csv")
