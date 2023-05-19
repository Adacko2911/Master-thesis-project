library(tidyverse)
library(R.utils)
library(tidyr)

##
#This script is not needed to run unless the all_gelaman_*.rds  is not created
#this parsing is just to get an information from large convergence files
##

setwd("./convergence_run/medhmm")
files <- list.files()
dbf.files <- files[grep(".rds", files, fixed=T)]

#------------------------------------------------------------------------------
#Gelman rubin
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
  if(m>2){
    chains_matrix_gamma<-matrix(ncol = (m^2-(2*m)))
    chains_matrix_gamma<-rbind(chains_matrix_gamma,model1$gamma_int_bar[(burn_in+1):J,1:(m^2-(2*m))],model2$gamma_int_bar[(burn_in+1):J,1:(m^2-(2*m))])
    chains_matrix_gamma<-na.omit(chains_matrix_gamma)
  }
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
  if(m==2){
    gelman_input<-cbind(index,chains_matrix_emiss,chains_matrix_dwell)
  }else{
    gelman_input<-cbind(index,chains_matrix_emiss,chains_matrix_dwell,chains_matrix_gamma)
  }

  out_scenario_gelman_empirical<-gelman_rubin(gelman_input,J = J,burn_in = burn_in)
  out_scenario_gelman_empirical<-t(as.matrix(out_scenario_gelman_empirical))
  colnames(out_scenario_gelman_empirical)<-colnames(gelman_input)[-1]
  return(out_scenario_gelman_empirical)
}
three_list_gelman<-list()
four_list_gelman<-list()

for(i in dbf.files){
  #In order to run the scripts without any issues using snellius we need to adjust the below target directory
  new<-try(readRDS(as.character(paste0("./convergence_run/medhmm/outputs",i))))
  if(is.list(new)){
    gel<-new[["gelman_rubin_stats"]]
    n <- 32
    length(gel)<-n

  numbers<- regmatches(i, gregexpr("[[:digit:]]+", i)) %>% unlist()
  gel<-matrix(gel,nrow=1)
  rownames(gel)<-paste0("state ",numbers[1],"obs l",numbers[2], "dwell ",numbers[3], "set ",numbers[4])
  m=numbers[1]
  if(m==1 | m==2){
    short_list_gelman<-append(short_list_gelman,list(gel))
  }else{
    long_list_gelman<-append(long_list_gelman,list(gel))
  }
  }

}

a<-do.call(rbind,three_list_gelman)
b<-do.call(rbind, four_list_gelman)
