library(tidyverse)
#dwell time distribution parsing
dwell_medhmm<-readRDS("./medhmm_run/outputs/post_processed/list_dwell_emp.rds")
dwell_mhmm<-readRDS("./mhmm_run/outputs/post_processed/list_dwell_emp.rds")
dwell_final<-rbind(dwell_medhmm,dwell_mhmm) %>% write.csv(.,"./simulation_main/post_processed_data/dwell_final.csv",row.names=FALSE)

#switches parsing
switches_medhmm<-readRDS("./medhmm_run/outputs/post_processed/list_switches_est.rds") %>% as.matrix()
switches_mhmm<-readRDS("./mhmm_run/outputs/post_processed/list_switches_est.rds") %>% as.matrix()
switches_final<-rbind(switches_medhmm,switches_mhmm) %>% write.csv(.,"./simulation_main/post_processed_data/switches_final.csv",row.names=FALSE)
#decoding parsing
decoding_medhmm<-readRDS("./medhmm_run/outputs/post_processed/list_decod_acc.rds") %>% as.matrix()
decoding_mhmm<-readRDS("./mhmm_run/outputs/post_processed/list_decod_acc.rds") %>% as.matrix()
decoding_final<-rbind(decoding_medhmm,decoding_mhmm) %>% write.csv(.,"./simulation_main/post_processed_data/decoding_final.csv",row.names=FALSE)
#parse the outputs form the maximum a  posterior estimates
post_est_medhmm<-readRDS("./medhmm_run/outputs/post_processed/list_post_est.rds")
post_est_mhmm<-readRDS("./mhmm_run/outputs/post_processed/list_post_est.rds")
post_est_final<-rbind(post_est_medhmm,post_est_mhmm) %>% write.csv(.,"./simulation_main/post_processed_data/post_est_final.csv",row.names=FALSE)

#parse together results that were separated form the post_est_final however only include comparable dwell time and gamma distribution
dwell_param_gamma_medhmm<-read.csv("./medhmm_run/outputs/post_processed/list_post_est_dwell_gamma.csv")
dwell_param_gamma_mhmm<-read.csv("./mhmm_run/outputs/post_processed/list_post_est_dwell_gamma.csv")
dwell_param_gamma_finial<-rbind(dwell_param_gamma_medhmm,dwell_param_gamma_mhmm) %>%
  dplyr::select(-c(1,2))%>% write.csv(.,"./simulation_main/post_processed_data/dwell_param_gamma_finial.csv",row.names=FALSE)
