#---------------------------------------------------
# Fit model for M hidden states and save results
#---------------------------------------------------

#Fitting the model
library(mHMMbayes)
library(medHMM)
library(tidyverse)
library(depmixS4)

#Load utility functions
source("R/utils.r")

#---read data
## Note the data is not included in the repository as it is only available on a special query
bipolar_data <- foreign::read.spss("ESM bipolar cleaned data_UU project.sav", to.data.frame = TRUE)
bipolar_data <- bipolar_data %>%
  dplyr::select(patient_id, time,
                bs_diary_5, bs_diary_13, bs_diary_22,
                bs_diary_15, bs_diary_9, bs_diary_10,
                bs_diary_7, bs_diary_11, bs_diary_17,
                bs_diary_8, bs_diary_14, bs_diary_16) %>%
  drop_na() %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id()) %>%
  ungroup() %>%
  arrange(patient_id, time) %>%
  dplyr::select(-time) %>%
  as.matrix()

#------------------------------------------------------------------------------#

out <- readRDS(paste0("./example/outputs/out_cont_gamma_prior_emiss_prior_m4_12dv_it4000_c",2,".rds"))


#------------------------------------------------------------------------------#
#Pooled data analisys
#load the data for the model
library(depmixS4)
#put together the model for the pooled hmm
mod <- depmixS4::depmix(list(dep_1~1,dep_2~1,dep_3~1,dep_4~1,dep_5~1,dep_6~1,dep_7~1,dep_8~1,dep_9~1,dep_10~1,dep_11~1,dep_12~1),data=bipolar_data,nstates=4,
                        family=list(gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian()),ntimes = c(386, 589, 538, 405, 431, 540, 608, 559, 399, 613, 308, 427, 434, 428, 486, 541, 197, 869, 489, 566))
# print the model, formula and parameter values
mod
#set seed
set.seed(123)
# fit the model by calling fit function 5 times and chose the best one
fm1 <- fit(mod,em=em.control(maxit=500))
fm2 <- fit(mod,em=em.control(maxit=500))
fm3 <- fit(mod,em=em.control(maxit=500))
fm4 <- fit(mod,em=em.control(maxit=500))
fm5 <- fit(mod,em=em.control(maxit=500))
#in this case the best model is model number 3
plot(1:5,c(BIC(fm1),BIC(fm2),BIC(fm3),BIC(fm4),BIC(fm5)),ty="b")
#obtain the transition probability parameters
summary(fm3)
summary(fm3,which="response")

#==================================================END of group level


# Def
m <- 4
n_dep <- 12
#---preprocess data
#---set seed--------------------------------------------------------------------
set.seed(123)
#--don't change
m=4
n_dep=12
J=2000
B=1000

#---Specify hyper-prior for the continuous emission distribution----------------
#gamma taken from the output
gamma_start <- int_to_prob(matrix(apply(out$gamma_int_bar[2001:4000,],2,median), nrow = 4, byrow = TRUE))

#calculate gamma_start_no_diag the gamma starting value
gamma_start_no_diag<-gamma_start
diag(gamma_start_no_diag)<-rep(0,m)
gamma_start_no_diag<-gamma_start_no_diag%>% apply(.,1,FUN=function(x){x/sum(x)}) %>% t()



# Def starting values emiss
emiss_start <- lapply(1:12, function(q) matrix(c(apply(out$emiss_mu_bar[[q]][2001:4000,],2,median),
                                               apply(out$emiss_var_bar[[q]][2001:4000,],2,median)),
                                             nrow = m, ncol = 2, byrow = FALSE))




#calculate the starting log-dwell-mean value
self_trans<-diag(gamma_start)
# expected_mean_dwell_times<-1/-log(self_trans)
expected_mean_dwell_times <- c(17.1, 9.8, 8.8, 16.3)
expected_mean_dwell_times<-1/-log(self_trans)
#set the standard log deviation to 20 to be on the save side
sd_dwell_times<-rep(100,m)

#dwell start ready-to use matrix
dwell_start<-log(cbind(expected_mean_dwell_times,sd_dwell_times))

emiss_hyp_prior <- list(
  emiss_mu0 = list(matrix(c(5, 25, 50, 50), nrow = 1), # down
                   matrix(c(5, 25, 50, 50), nrow = 1), # dread rest of day
                   matrix(c(5, 25, 50, 50), nrow = 1), # worry
                   matrix(c(5, 25, 50, 50), nrow = 1), # inadequate
                   matrix(c(25, 25, 50, 50), nrow = 1), # tired
                   matrix(c(50, 25, 5, 5), nrow = 1), # content

                   matrix(c(5, 25, 25, 10), nrow = 1), # agitated
                   matrix(c(5, 25, 25, 10), nrow = 1), # irritated
                   matrix(c(25, 25, 25, 10), nrow = 1), # switch
                   matrix(c(25, 25, 25, 10), nrow = 1), # extremely well
                   matrix(c(25, 25, 25, 10), nrow = 1), # ideas
                   matrix(c(25, 25, 25, 10), nrow = 1)  # thoughts racing
  ),
  emiss_K0  = rep(list(1),n_dep),
  emiss_nu  = rep(list(1),n_dep),
  emiss_V   = rep(list(rep(400, m)),n_dep),
  emiss_a0  = rep(list(rep(0.001, m)),n_dep),
  emiss_b0  = rep(list(rep(0.001, m)),n_dep))

dwell_hyp_prior<-list(d_mu0=log(expected_mean_dwell_times),
                      s2_0=log(rep(100, m)),
                      alpha.sigma20	= rep(0.01, m),
                      beta.sigma20	= rep(0.01, m),
                      alpha.tau20		= rep(0.01, m),
                      beta.tau20		= rep(0.01, m))

#================Run the MHMM model on the bipolar data==============================================================
out_empirical_mhmm_st4_chain1<- mHMM_cont(s_data = bipolar_data,
                                          gen = list(m = m, n_dep = n_dep),
                                          start_val = c(list(gamma_start), emiss_start),
                                          emiss_hyp_prior = emiss_hyp_prior,
                                          mcmc = list(J = J, burn_in = B),
                                          return_path = TRUE)

saveRDS(out_empirical_mhmm_st4_chain1, paste0("./example/outputs/out_mHMM_m",m,"_12dv_it2000_c1.rds"))

#================Run the SECOND CHAIN MHMM model on the bipolar data==================================================
gamma2   <- gamma_start+rnorm(m*m,mean=0,sd=0.02)
gamma2<-gamma2%>% abs() %>% t(.) %>% apply(.,2,FUN=function(x){x/sum(x)}) %>% t()
emissions2<-emiss_start %>% lapply(., function(x){ x*matrix(sample(c(0.95,1.05),m*2,replace = T),ncol = 2,byrow = T)})

out_empirical_mhmm_st4_chain3<-mHMM_cont(s_data = bipolar_data,
                                         gen = list(m = m, n_dep = n_dep),
                                         start_val = c(list(gamma2), emissions2),
                                         emiss_hyp_prior = emiss_hyp_prior,
                                         mcmc = list(J = J, burn_in = B), return_path = TRUE)

saveRDS(out_empirical_mhmm_st4_chain3, paste0("./example/outputs/out_mHMM_m",m,"_12dv_it2000_c3.rds"))

#================Run the MEDHMM model on the bipolar data=============================================================
out_empirical_medhmm_st4_chain1<- try(medHMM_cont(s_data = bipolar_data,
                                              gen = list(m = m, n_dep = n_dep),
                                              start_val = c(list(gamma_start_no_diag),  emiss_start, list(dwell_start)),
                                              emiss_hyp_prior = emiss_hyp_prior,
                                              dwell_hyp_prior = dwell_hyp_prior,
                                              show_progress = TRUE,
                                              mcmc = list(J = J, burn_in = B),
                                              # mcmc = list(J = 20, burn_in = 10),
                                              return_path = TRUE, max_dwell = NULL))

saveRDS(out_empirical_medhmm_st4_chain1, paste0("./example/outputs/out_medHMM_m",m,"_12dv_it2000_c1.rds"))

#================Run the SECOND CHAIN MEDHMM model on the bipolar data================================================
dwell_start2<-dwell_start*matrix(sample(c(0.95,1.05),m*2,replace = T),ncol = 2,byrow = T)
gamma2   <- gamma_start_no_diag
diag(gamma2)=NA
gamma2<-t(matrix(t(gamma2)[which(!is.na(gamma2))],nrow=nrow(gamma2)-1,ncol=ncol(gamma2)))+matrix(rnorm(m*(m-1),mean=0,sd=0.02),ncol=m-1) %>% abs()
gamma2<-gamma2%>% apply(.,1,FUN=function(x){x/sum(x)})  %>% as.numeric()
gamma2<-matrix(data = insert(gamma2,
                             seq(from = 1,to = length(gamma2) + 1,by = m),
                             value = rep(0, m)),
               byrow = T,
               ncol = m)

out_empirical_medhmm_st4_chain3<- medHMM_cont(s_data = bipolar_data,
                                              gen = list(m = m, n_dep = n_dep),
                                              start_val = c(list(gamma2), emissions2, list(dwell_start2)),
                                              emiss_hyp_prior = emiss_hyp_prior,
                                              dwell_hyp_prior = dwell_hyp_prior,
                                              show_progress = TRUE,
                                              mcmc = list(J = J, burn_in = B), return_path = TRUE, max_dwell = NULL)

saveRDS(out_empirical_medhmm_st4_chain3, paste0("./example/outputs/out_medHMM_m",m,"_12dv_it2000_c3.rds"))

#=====================================================================================================================
