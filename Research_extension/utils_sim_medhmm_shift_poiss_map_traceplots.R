# Data simulation: Poisson distribution
#----------Start of the data simulating file-----------------------------------------------------------------------
#utility functions
int_to_prob<-function(alfa){
  dim<-dim(alfa)
  gamma<-matrix(ncol = dim[2]+1,nrow = dim[1])
  exp_alfa<-exp(alfa)
  exp_alfa<-cbind(rep(1,dim[1]),exp_alfa)
  for(i in 1:dim[1]){
    gamma[i,]<-exp_alfa[i,]/sum(exp_alfa[i,])
  }
  return(round(gamma,4))
}
prob_to_int<-function(gamma){
  gamma<-gamma+0.00001
  if(is.matrix(gamma)){
    dim<-dim(gamma)
    alfa<-matrix(ncol = dim[2]-1,nrow = dim[1])
    for (i in 1:dim[1]) {
      for(k in 2:dim[2]){
        alfa[i,k-1]<-log(gamma[i,k]*(1/gamma[i,1]))
      }
    }
  }
  return(round(alfa,4))

}


#-------------Three state scenario------------------------------------------------------------------------------------
mhsmm_sim_shift_pois <-
  function(n,
           n_t,
           m,
           n_dep,
           emiss,
           gamma,
           gamma_ss_var,
           emiss_ss_var,
           dwell_ss_var,
           dwell_distr,
           shift=1) {
    #define output lists
    gamma_ss <- rep(list(NULL), n)
    state_sequence <-
      as.matrix(cbind(rep(1:n, each = n_t), rep(NA, n_t * n)))
    emiss_ss <- rep(list(vector("list", n_dep)), n)
    dwell_ss<-rep(list(NULL), n)
    obs_ss <-
      cbind(rep(1:n, each = n_t), matrix(rep(NA, n_t * n * n_dep), ncol = n_dep))
    gamma_no_diag <- gamma
    diag(gamma_no_diag) = NA
    gamma_no_diag <-
      t(matrix(
        t(gamma_no_diag)[which(!is.na(gamma_no_diag))],
        nrow = nrow(gamma_no_diag) - 1,
        ncol = ncol(gamma_no_diag)
      ))
    #change the form of a transition matrix
    gamma_int <- prob_to_int(gamma_no_diag)
    for (subj in 1:n) {
      #generate a set of subject specific transition matrices we introduce random subject-specific effect to the intercepts
      gamma_numeric <-
        as.numeric(t(int_to_prob(
          gamma_int + rnorm(
            n =  m * (m - 1) / 2,
            mean = 0,
            sd = sqrt(gamma_ss_var)
          )
        )))
      gamma_ss[[subj]] <-
        matrix(
          data = insert(
            gamma_numeric,
            seq(
              from = 1,
              to = length(gamma_numeric) + 1,
              by = m
            ),
            value = rep(0, m)
          ),
          byrow = T,
          ncol = m
        )
      for (dep in 1:n_dep) {
        #generate subject specific means of emission distributions of dependent variables,we introduce random subject-specific effect which is deviation from the mean of normally conditionally distributed mean
        emiss_ss[[subj]][[dep]] <- emiss[[dep]]
        emiss_ss[[subj]][[dep]][, 1] <-
          matrix(emiss[[dep]][, 1] + rnorm(
            n = m,
            mean = 0,
            sd = sqrt(emiss_ss_var)), nrow = m)
      }

      #for each subject we sample the initial state
      inst = c(1:m)
      state_sequence[subj * n_t - n_t + 1, 2] <- sample(inst, 1)
      dwell_row <- state_sequence[subj * n_t - n_t + 1, 2]
      for (dep in 1:n_dep) {
        obs_ss[subj * n_t - n_t + 1, dep + 1] <-
          rnorm(n = 1,
                mean = emiss_ss[[subj]][[dep]][dwell_row, 1],
                sd = sqrt(emiss_ss[[subj]][[dep]][dwell_row, 2]))
      }
      a = 1
      mu_dwell<-rep(NA,m)
      for(s in 1:m){
        mu_dwell[s]<-exp(dwell_distr[s,] + rnorm(1, mean = 0, sd = sqrt(dwell_ss_var[[s]])))
        mu_dwell[s]=round(max(1,mu_dwell[s]),0)
      }

      dwell_ss[[subj]] <- matrix(mu_dwell,nrow=m,ncol=1)
      while (a < n_t) {
        #draw the subject specific discrete dwell time
        mean_dwell <-
          round(rpois(1, lambda = dwell_ss[[subj]][dwell_row,1]), 0)+shift
        mean_dwell<-mean_dwell
        #condition to meet the n_t subject observations sequence length
        if (mean_dwell <= 1) {
          mean_dwell = 1
        }else{
          if (a + (mean_dwell - 1) <= n_t) {
            #repeat the most probable state
            state_given_dwell <- rep(dwell_row, mean_dwell - 1)
            #append the sequence to the observation matrix
            state_sequence[(subj * n_t - n_t + a + 1):(subj * n_t - n_t + a + mean_dwell - 1), 2] <- state_given_dwell
            #update the a(current index in the subject specific vector of observations)
            a <- a + mean_dwell - 1
          } else{
            mean_dwell <- n_t - a
            state_given_dwell <- rep(dwell_row, (mean_dwell))
            state_sequence[(subj * n_t - n_t + a + 1):(subj * n_t -n_t + a + mean_dwell), 2] <- state_given_dwell
            a <- a + mean_dwell
          }
        }
        if (a + 1 <= n_t) {
          state_sequence[subj * n_t - n_t + a + 1, 2] <-
            sample(x = 1:m,
                   size = 1,
                   prob = gamma_ss[[subj]][dwell_row, ])  #this loop is out of bound idk
          dwell_row <-
            state_sequence[subj * n_t - n_t + a + 1, 2] #the row of state specific dwell time
          a <- a + 1
        }
      }


      #for each subject and each state we sample the mean emission
      for (t_occ in 2:n_t) {
        for (dep in 1:n_dep) {
          obs_ss[subj * n_t - n_t + t_occ, dep + 1] <-
            rnorm(n = 1,
                  mean = emiss_ss[[subj]][[dep]][state_sequence[subj * n_t - n_t + t_occ, 2], 1],
                  sd = sqrt(emiss_ss[[subj]][[dep]][state_sequence[subj * n_t - n_t + t_occ, 2], 2])) %>% round(digits = 5)
        }#end of dependent variable loop
      }#end of series loop
    }#end of subject loop
    #the output is return here

    state_sequence<-matrix(as.integer(state_sequence),ncol=2)
    obs_ss[,1]<-as.integer(obs_ss[,1])
    output <- list( state_sequence,  obs_ss,gamma_ss,emiss_ss,dwell_ss)   #
    names(output) <- c( "states",  "observations","gamma_ss","emiss_ss","dwell_ss")#
    return(output)
  }



emiss_maximum_a_posteriori<-function(case_out=NULL,case_out_name="emiss_mu_bar",iteration=1,J=NULL,B=NULL,m=NULL){
  median<-case_out[[glue("{case_out_name}")]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out[[glue("{case_out_name}")]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out[[glue("{case_out_name}")]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out[[glue("{case_out_name}")]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[[glue("{case_out_name}")]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_025),"_ci_975")
  output<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(output)<-paste0(glue("{case_out_name}"),iteration)
  return(output)
}
no_emiss_maximum_a_posteriori<-function(case_out=NULL,case_out_name="dwell_mu_bar",iteration=1,J=NULL,B=NULL,m=NULL){
  median<-case_out[[glue("{case_out_name}")]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),"_median")
  mean<-case_out[[glue("{case_out_name}")]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),"_mean")
  sd<-case_out[[glue("{case_out_name}")]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),"_sd")
  ci_025<-case_out[[glue("{case_out_name}")]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[[glue("{case_out_name}")]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),"_ci_025")
  colnames(ci_975)<-paste0(colnames(ci_025),"_ci_975")
  output<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(output)<-paste0(glue("{case_out_name}"),iteration)
  return(output)
}

# Define helper function to output only MAP estimates for medhmm
MAP_medhmm<-function(case_out=NULL,iteration=1,J=NULL,B=NULL,m=NULL){
  iteration=iteration
  J=J
  B=B
  m=m
  emiss_mu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_mu_bar"]])){
    emiss_mu_bar_list[[i]]<-emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "emiss_mu_bar",iteration=iteration,J=J,B=B,m=m)
  }
  emiss_var_bar_list<-list()
  for(i in 1:length(case_out[["emiss_var_bar"]])){
    emiss_var_bar_list[[i]]<-emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "emiss_var_bar",iteration=iteration,J=J,B=B,m=m)
  }
  emiss_varmu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_varmu_bar"]])){
    emiss_varmu_bar_list[[i]]<-emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "emiss_varmu_bar",iteration=iteration,J=J,B=B,m=m)
  }
  dwell_mu_bar<-no_emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "dwell_mu_bar",iteration=iteration,J=J,B=B,m=m)
  dwell_varmu_bar<-no_emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "dwell_varmu_bar",iteration=iteration,J=J,B=B,m=m)

  gamma_prob_bar<-no_emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "gamma_prob_bar")
  gamma_int_prob_bar<-no_emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "gamma_int_bar",iteration=iteration,J=J,B=B,m=m)
  #subject specific parameters of variance covariance between
  gamma_V_names<-c()
  for(st_cur in 1:m){
    to_st<-2:m
    to_st<-to_st[which(to_st!=st_cur)]
    gamma_V_names<-append(gamma_V_names,paste0("S",st_cur,"toS",to_st))
  }
  names_need<-c()
  for(al in 1:length(gamma_V_names)){
    names_need<-append(names_need,paste0("var_int_",gamma_V_names[al],"_with_int_",gamma_V_names[al]))
  }
  names_need<-names_need[-1]
  case_out_gamma_V_int_bar<-list(gamma_V_int_bar=case_out[["gamma_V_int_bar"]][,names_need])
  gamma_V_int_bar<-no_emiss_maximum_a_posteriori(case_out=case_out_gamma_V_int_bar,case_out_name = "gamma_V_int_bar",iteration=iteration,J=J,B=B,m=m)

  MAP<-list(emiss_mu_bar_list,emiss_var_bar_list,emiss_varmu_bar_list,dwell_mu_bar,dwell_var_bar,dwell_varmu_bar,gamma_prob_bar,gamma_int_bar,gamma_V_int_bar)
  names(MAP)<-c("emiss_mu_bar","emiss_var_bar","emiss_varmu_bar","dwell_mu_bar","dwell_var_bar","dwell_varmu_bar","gamma_prob_bar","gamma_int_bar","gamma_V_int_bar")

  return(MAP)
}



MAP_mHMM<-function(case_out=NULL,iteration=1,B=NULL,J=NULL,m=NULL){
  iteration=iteration
  J=J
  B=B
  m=m
  emiss_mu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_mu_bar"]])){
    emiss_mu_bar_list[[i]]<-emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "emiss_mu_bar",iteration=iteration,J=J,B=B,m=m)
  }
  emiss_var_bar_list<-list()
  for(i in 1:length(case_out[["emiss_sd_bar"]])){
    out_temp=emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "emiss_sd_bar",iteration=iteration,J=J,B=B,m=m)^2
    rownames(out_temp)<-'emiss_var_bar'
    emiss_var_bar_list[[i]]<-out_temp
  }
  rm(out_temp)
  emiss_varmu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_varmu_bar"]])){
    emiss_varmu_bar_list[[i]]<-emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "emiss_varmu_bar",iteration=iteration,J=J,B=B,m=m)
  }

  gamma_prob_bar<-no_emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "gamma_prob_bar",iteration=iteration,J=J,B=B,m=m)

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
  rownames(gamma_prob_bar_no_diag)<-paste0("gamma_prob_bar_no_diag",iteration)

  gamma_int_bar<-no_emiss_maximum_a_posteriori(case_out=case_out,case_out_name = "gamma_int_bar",iteration=iteration,J=J,B=B,m=m)

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

  MAP<-list(emiss_mu_bar_list,emiss_var_bar_list,emiss_varmu_bar_list,gamma_prob_bar,gamma_prob_bar_no_diag,gamma_int_bar,dwell_mu_bar_exp)
  names(MAP)<-c("emiss_mu_bar","emiss_var_bar","emiss_varmu_bar","gamma_prob_bar","gamma_prob_bar_no_diag","gamma_int_bar","dwell_mu_bar")
  return(MAP)
}

#local decoding function
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

#trace plots
traceplots_emiss_mu<-function(model2,dep=1){
  model2_em<-as.data.frame(model2$emiss_mu_bar[[dep]])
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]] <-ggplot(data=model2_em,aes_string(x=1:data_points, y=model2_em[,i]))+geom_line(color="red")+
      theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model2_em)[i])

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob("Convergence: emiss mu"))))
  print(out)
}

traceplots_dwell_mu<-function(model2){
  model2_em<-exp(as.data.frame(model2$dwell_mu_bar))
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  p<-list()
  for(i in 1:ncol(model2_em)){
    p[[i]]<-ggplot(data=model2_em,aes_string(x=1:data_points, y=model2_em[,i]))+geom_line(color="red")+
      theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+
      ylab(colnames(model2_em)[i])

  }
  out<-do.call(grid.arrange,c(p, list(top=textGrob("Convergence: dwell"))))
  print(out)}


traceplots_gamma<-function(model1){
  model1_em<-as.data.frame(model1$gamma_prob_bar)
  #model2_em<-as.data.frame(model2$gamma_prob_bar)
  names<-colnames(model1_em)
  data_points<-nrow(model1_em)
  p<-list()
  for(i in 1:ncol(model1_em)){
    p[[i]]<-ggplot(data=model1_em,aes_string(x=1:data_points, y=model1_em[,i]))+geom_line(color="red")+theme_minimal()+theme(plot.title = element_text(size=8,hjust = 0.5))+xlab("iteration")+ylab(colnames(model1_em)[i])
  }
  out<-do.call(grid.arrange, c(p, list(top=textGrob("Convergence: gamma"))))
  print(out)
}
