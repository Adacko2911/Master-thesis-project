
colors=tibble(blue="#5778a4",
              orange="#e49444",
              red="#d1615d",
              teal="#85b6b2",
              green="#6a9f58",
              yellow="#e7ca60",
              purple="#a87c9f",
              pink="#f1a2a9",
              brown="#967662",
              grey="#b8b0ac")

#utility functions necessary for simulation function
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
#simulate data under the medhmm
#adjusted such that the n_t can be vector of length #subjects number
mhsmm_sim <-function(n, #individual number
           n_t, #number of occasions
           m, #number of states
           n_dep, #number of dependent variables
           emiss, #emission distribution list of matrices first column are means the second variances
           gamma, #gamma transition probability matrix for group level
           gamma_ss_var, #subject specific variance in gamma intercepts
           emiss_ss_var, #subject specific variance in emission distributions
           dwell_distr) { #dwell time distribution parameters wher the first column are means and the second one are variances before log transformation

    #define output lists
    gamma_ss <- rep(list(NULL), n)
    state_sequence <-as.matrix(cbind(rep(1:n, times = n_t), rep(NA, sum(n_t))))
    emiss_ss <- rep(list(vector("list", n_dep)), n)
    dwell_ss<-rep(list(NULL), n)
    obs_ss <-
      cbind(rep(1:n, times  = n_t), matrix(rep(NA, sum(n_t) * n_dep), ncol = n_dep))
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
            n = m*(m-2),
            mean = 0, sd = sqrt(as.numeric(gamma_ss_var))
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
            sd = sqrt(as.numeric(emiss_ss_var[[dep]]))
          ), nrow = m)
      }

      #for each subject we sample the initial state
      inst = c(1:m)
      state_sequence[sum(n_t[1:subj]) - n_t[subj] + 1, 2] <- sample(inst, 1)
      dwell_row <- state_sequence[sum(n_t[1:subj]) - n_t[subj] + 1, 2]
      for (dep in 1:n_dep) {
        obs_ss[sum(n_t[1:subj]) - n_t[subj] + 1, dep + 1] <-
          rnorm(n = 1,
                mean = emiss_ss[[subj]][[dep]][dwell_row, 1],
                sd = sqrt(emiss_ss[[subj]][[dep]][dwell_row, 2]))
      }
      a = 1
      mu_dwell<-c()
      for(s in 1:m){
        mu_dwell[s]<-dwell_distr[s, 1] + rnorm(1, mean = 0, sd = sqrt(dwell_distr[s, 3]))
        if(mu_dwell[s]<0){
          mu_dwell[s]=0
        }
      }

      dwell_ss[[subj]] <- cbind(mu_dwell,dwell_distr[, 2])
      while (a < n_t[subj]) {

        #draw the subject specific discrete dwell time
        mean_dwell <-
          round(rlnorm(1, log(dwell_ss[[subj]][dwell_row,1]), log(sqrt(dwell_ss[[subj]][dwell_row,2]))), 0)
        #condition to meet the n_t subject observations sequence length
        if (mean_dwell <= 1) {
          mean_dwell = 1
        }else{
          if (a + (mean_dwell - 1) <= n_t[subj]) {
            #replete the most probable state
            state_given_dwell <- rep(dwell_row, mean_dwell - 1)
            #append the sequence to the observation matrix
            state_sequence[(sum(n_t[1:subj]) - n_t[subj] + a + 1):(sum(n_t[1:subj]) - n_t[subj] + a + mean_dwell - 1), 2] <- state_given_dwell
            #update the a(current index in the subject specific vector of observations)
            a <- a + mean_dwell - 1
          } else{
            mean_dwell <- n_t[subj] - a
            state_given_dwell <- rep(dwell_row, (mean_dwell))
            state_sequence[(sum(n_t[1:subj]) - n_t[subj] + a + 1):(sum(n_t[1:subj]) -n_t[subj] + a + mean_dwell), 2] <- state_given_dwell
            a <- a + mean_dwell
          }
        }
        if (a + 1 <= n_t[subj]) {
          state_sequence[sum(n_t[1:subj]) - n_t[subj] + a + 1, 2] <-
            sample(x = 1:m,
                   size = 1,
                   prob = gamma_ss[[subj]][dwell_row, ])  #this loop is out of bound idk
          dwell_row <-
            state_sequence[sum(n_t[1:subj]) - n_t[subj] + a + 1, 2] #the row of state specific dwell time
          a <- a + 1
        }
      }

      #for each subject and each state we sample the mean emission
      for (t_occ in 2:n_t[subj]) {
        for (dep in 1:n_dep) {
          obs_ss[sum(n_t[1:subj]) - n_t[subj] + t_occ, dep + 1] <-
            rnorm(n = 1,
                  mean = emiss_ss[[subj]][[dep]][state_sequence[sum(n_t[1:subj]) - n_t[subj] + t_occ, 2], 1],
                  sd = sqrt(emiss_ss[[subj]][[dep]][state_sequence[sum(n_t[1:subj]) - n_t[subj] + t_occ, 2], 2]))
        }#end of dependent variable loop
      }#end of series loop
    }#end of subject loop
    #the output is return here
    output <- list( state_sequence,  obs_ss,gamma_ss,emiss_ss )   #,dwell_ss
    names(output) <- c("states", "observations","gamma_ss", "emiss_ss")  # "dwell_ss"
    return(output)
  }


#Gelman Rubin statistics calculation
#help function to obtain the convergence accuracy
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
#Gelma-rubin final
ge_rub_two_models_emp<-function(model1=NULL, model2=NULL, m=4, J=4000,burn_in=2000){
  chains_matrix_gamma<-matrix(ncol = m^2)
  chains_matrix_gamma<-rbind(chains_matrix_gamma,model1$gamma_prob_bar[(burn_in+1):J,1:(m^2)],model2$gamma_prob_bar[(burn_in+1):J,1:(m^2)])
  chains_matrix_gamma<-na.omit(chains_matrix_gamma)
  var_data<-rbind(model1$emiss_var_bar[[1]][(burn_in+1):(J),],model2$emiss_var_bar[[1]][(burn_in+1):(J),])
  mean_data<-rbind(model1$emiss_mu_bar[[1]][(burn_in+1):J,1:m],model2$emiss_mu_bar[[1]][(burn_in+1):J,1:m])
  chains_matrix_emiss<-cbind(mean_data,var_data)
  chains_matrix_emiss<-na.omit(chains_matrix_emiss)
  index<-c(rep(1,J-burn_in),rep(2,J-burn_in))
  gelman_input<-cbind(index,chains_matrix_emiss,chains_matrix_gamma)
  out_scenario_gelman_empirical<-gelman_rubin(gelman_input,J = J,burn_in = burn_in)
  return(paste("The Gelman-Rubin statistics: ",out_scenario_gelman_empirical))
}


###  #MAP functions for MHMM

MAP_mHMM<-function(case_out=NULL,iteration=1,B=NULL,J=NULL,m=NULL,ss=FALSE){
  #---emiss_mu_bar---------
  emiss_mu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_mu_bar"]])){
    median<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    emiss_mu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_mu_bar)<-paste0("emiss_mu_bar",iteration)
    emiss_mu_bar_list[[i]]<-emiss_mu_bar
  }


  #---emiss_var_bar--------
  emiss_var_bar_list<-list()
  for(i in 1:length(case_out[["emiss_var_bar"]])){
    median<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    emiss_var_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_var_bar)<-paste0("emiss_var_bar",iteration)
    emiss_var_bar_list[[i]]<-emiss_var_bar
  }
  #---emiss_varmu_bar------
  emiss_varmu_bar_list<-list()
  for(i in 1:length(case_out[["emiss_varmu_bar"]])){
    median<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    emiss_varmu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(emiss_varmu_bar)<-paste0("emiss_varmu_bar",iteration)
    emiss_varmu_bar_list[[i]]<-emiss_varmu_bar
  }
  #---gamma_prob_bar------
  median<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),".median")
  mean<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),".mean")
  sd<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),".sd")
  ci_025<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
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
  colnames(median)<-paste0(colnames(median),".median")
  mean<-gamma_prob_bar_no_diag %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),".mean")
  sd<-gamma_prob_bar_no_diag %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),".sd")
  ci_025<-gamma_prob_bar_no_diag %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-gamma_prob_bar_no_diag %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
  gamma_prob_bar_no_diag<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(gamma_prob_bar_no_diag)<-paste0("gamma_prob_bar_no_diag")


  #---gamma_int_prob_bar------
  median<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
  colnames(median)<-paste0(colnames(median),".median")
  mean<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),".mean")
  sd<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),".sd")
  ci_025<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
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
  colnames(median)<-paste0(colnames(median),".median")
  mean<-case_out_gamma_V_int_bar %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),".mean")
  sd<-case_out_gamma_V_int_bar %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),".sd")
  ci_025<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
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
  colnames(median)<-paste0(colnames(median),".median")
  mean<-dwell_mu_bar_exp %>% apply(.,2,mean) %>% t()
  colnames(mean)<-paste0(colnames(mean),".mean")
  sd<-dwell_mu_bar_exp %>% apply(.,2,sd) %>% t()
  colnames(sd)<-paste0(colnames(sd),".sd")
  ci_025<-dwell_mu_bar_exp %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
  ci_975<-dwell_mu_bar_exp %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
  colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
  colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
  dwell_mu_bar_exp<-cbind(median,mean,sd,ci_025,ci_975)
  rownames(dwell_mu_bar_exp)<-paste0("dwell_mu_bar",iteration)


  MAP<-list(emiss_mu_bar_list,emiss_var_bar_list,emiss_varmu_bar_list,gamma_prob_bar,gamma_prob_bar_no_diag,gamma_int_bar,gamma_V_int_bar,dwell_mu_bar_exp)
  names(MAP)<-c("emiss_mu_bar","emiss_var_bar","emiss_varmu_bar","gamma_prob_bar","gamma_prob_bar_no_diag","gamma_int_bar","gamma_V_int_bar","dwell_mu_bar")
  return(MAP)
}

### MAP function for MEDHMM
MAP_medhmm<-function(case_out=NULL,iteration=1,J=NULL,B=NULL,m=NULL){
    #---emiss_mu_bar---------
    emiss_mu_bar_list<-list()
    for(i in 1:length(case_out[["emiss_mu_bar"]])){
      median<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
      colnames(median)<-paste0(colnames(median),".median")
      mean<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
      colnames(mean)<-paste0(colnames(mean),".mean")
      sd<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
      colnames(sd)<-paste0(colnames(sd),".sd")
      ci_025<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
      ci_975<-case_out[["emiss_mu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
      colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
      colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
      emiss_mu_bar<-cbind(median,mean,sd,ci_025,ci_975)
      rownames(emiss_mu_bar)<-paste0("emiss_mu_bar",iteration)
      emiss_mu_bar_list[[i]]<-emiss_mu_bar
    }
    #---emiss_var_bar--------
    emiss_var_bar_list<-list()
    for(i in 1:length(case_out[["emiss_var_bar"]])){
      median<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
      colnames(median)<-paste0(colnames(median),".median")
      mean<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
      colnames(mean)<-paste0(colnames(mean),".mean")
      sd<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
      colnames(sd)<-paste0(colnames(sd),".sd")
      ci_025<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
      ci_975<-case_out[["emiss_var_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
      colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
      colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
      emiss_var_bar<-cbind(median,mean,sd,ci_025,ci_975)
      rownames(emiss_var_bar)<-paste0("emiss_var_bar",iteration)
      emiss_var_bar_list[[i]]<-emiss_var_bar
    }
    #---emiss_varmu_bar------
    emiss_varmu_bar_list<-list()
    for(i in 1:length(case_out[["emiss_varmu_bar"]])){
      median<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,median) %>% t()
      colnames(median)<-paste0(colnames(median),".median")
      mean<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,mean) %>% t()
      colnames(mean)<-paste0(colnames(mean),".mean")
      sd<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,sd) %>% t()
      colnames(sd)<-paste0(colnames(sd),".sd")
      ci_025<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
      ci_975<-case_out[["emiss_varmu_bar"]][[i]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
      colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
      colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
      emiss_varmu_bar<-cbind(median,mean,sd,ci_025,ci_975)
      rownames(emiss_varmu_bar)<-paste0("emiss_varmu_bar",iteration)
      emiss_varmu_bar_list[[i]]<-emiss_varmu_bar
    }
    #---dwell_mu_bar---------
    median<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["dwell_mu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    dwell_mu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(dwell_mu_bar)<-paste0("dwell_mu_bar",iteration)

    #---dwell_var_bar--------
    median<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["dwell_var_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    dwell_var_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(dwell_var_bar)<-paste0("dwell_var_bar",iteration)

    #---dwell_varmu_bar------
    median<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["dwell_varmu_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    dwell_varmu_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(dwell_varmu_bar)<-paste0("dwell_varmu_bar",iteration)

    #---gamma_prob_bar------
    median<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["gamma_prob_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    gamma_prob_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(gamma_prob_bar)<-paste0("gamma_prob_bar",iteration)

    #---gamma_int_prob_bar------
    median<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out[["gamma_int_bar"]][(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    gamma_int_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(gamma_int_bar)<-paste0("gamma_int_bar",iteration)

    #--gamma_V_int_bar-----------
    #subject specific parameters of variance covariance between
    all_st<-c()
    for(st_cur in 1:m){
      to_st<-2:m
      to_st<-to_st[which(to_st!=st_cur)]
      all_st<-append(all_st,paste0("S",st_cur,"toS",to_st))
    }
    names_need<-c()
    for(al in 1:length(all_st)){
      names_need<-append(names_need,paste0("var_int_",all_st[al],"_with_int_",all_st[al]))
    }
    names_need<-names_need[-1]
    case_out_gamma_V_int_bar<-case_out[["gamma_V_int_bar"]][(B+1):J,names_need]
    median<-case_out_gamma_V_int_bar %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-case_out_gamma_V_int_bar %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-case_out_gamma_V_int_bar %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.025)}) %>% t()
    ci_975<-case_out_gamma_V_int_bar %>% apply(.,2,FUN = function(x){quantile(x,0.975)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    gamma_V_int_bar<-cbind(median,mean,sd,ci_025,ci_975)
    rownames(gamma_V_int_bar)<-paste0("gamma_V_int_bar",iteration)

    MAP<-list(emiss_mu_bar_list,emiss_var_bar_list,emiss_varmu_bar_list,dwell_mu_bar,dwell_var_bar,dwell_varmu_bar,gamma_prob_bar,gamma_int_bar,gamma_V_int_bar)
    names(MAP)<-c("emiss_mu_bar","emiss_var_bar","emiss_varmu_bar","dwell_mu_bar","dwell_var_bar","dwell_varmu_bar","gamma_prob_bar","gamma_int_bar","gamma_V_int_bar")
    return(MAP)
  }
### MAP function buid for gamma with no diagonal entries
MAP_any_df<-function(df,J=2000,B=1000,mod="NA"){
    if(mod=="mhmm_gam"){
      self_trans_lab<-c()
      gamma_prob_bar_no_diag<-df
      for (i in 1:m) {
        self_trans_lab[i]<-paste0("S",i,"toS",i)
      }
      zeros<-matrix(rep(0,m*(nrow(gamma_prob_bar_no_diag))),ncol=m)
      gamma_prob_bar_no_diag[,self_trans_lab] <-zeros
      n<-0
      while(n<m*m){
        gamma_prob_bar_no_diag[,(n+1):(n+m)]<-gamma_prob_bar_no_diag[,(n+1):(n+m)]%>% apply(.,1,function(x){x/sum(x)}) %>% t()
        n<-n+m
      }
      df<-gamma_prob_bar_no_diag
    }

    median<-df[(B+1):J,] %>% apply(.,2,median) %>% t()
    colnames(median)<-paste0(colnames(median),".median")
    mean<-df[(B+1):J,] %>% apply(.,2,mean) %>% t()
    colnames(mean)<-paste0(colnames(mean),".mean")
    sd<-df[(B+1):J,] %>% apply(.,2,sd) %>% t()
    colnames(sd)<-paste0(colnames(sd),".sd")
    ci_025<-df[(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.025,na.rm = TRUE)}) %>% t()
    ci_975<-df[(B+1):J,] %>% apply(.,2,FUN = function(x){quantile(x,0.975,na.rm = TRUE)}) %>% t()
    colnames(ci_025)<-paste0(colnames(ci_025),".ci_025")
    colnames(ci_975)<-paste0(colnames(ci_975),".ci_975")
    gamma_prob_bar<-cbind(median,mean,sd,ci_025,ci_975)


    return(gamma_prob_bar)
  }


### Local decoding function
rep_each <- function(x, times) {
    times <- rep(times, length.out = length(x))
    rep(x, times = times)}


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
    st<-as.integer(st1)
    output<-data.frame(subject=as.factor(sub),state=st)
    return(output)
  }



#trace plots
traceplots_gamma<-function(model1,model2,true=NULL,titles=NULL,colors=NULL){
  model1_em<-as.data.frame(model1$gamma_prob_bar)
  model2_em<-as.data.frame(model2$gamma_prob_bar)
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  for(i in 1:ncol(model2_em)){
    plot<-ggplot(data=model1_em,aes(x=1:data_points, y=model1_em[,i]))+theme_minimal()+xlab("iteration")+ylab(colnames(model2_em)[i])
    if(!is.null(colors)){
      plot<-plot+geom_line(color=colors[1])+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color=colors[2])
    }else{
      plot<-plot+geom_line(color="red")+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color="blue")
    }

    if(!is.null(true)){
      plot<-plot+geom_hline(yintercept = true[i],size=1.5)
    }
    if(!is.null(titles)){
      plot<-plot+ggtitle(titles[i])
    }

    print(plot)
  }
}
traceplots_emiss_mu<-function(model1,model2,dep=1,true=NULL,titles=NULL,colors=NULL){
  model1_em<-as.data.frame(model1$emiss_mu_bar[[dep]])
  model2_em<-as.data.frame(model2$emiss_mu_bar[[dep]])
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  for(i in 1:ncol(model2_em)){
    plot<-ggplot(data=model1_em,aes(x=1:data_points, y=model1_em[,i]))+theme_minimal()+xlab("iteration")+ylab(colnames(model2_em)[i])
    if(!is.null(colors)){
      plot<-plot+geom_line(color=colors[1])+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color=colors[2])
    }else{
      plot<-plot+geom_line(color="red")+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color="blue")
    }
    if(!is.null(titles)){
      plot<-plot+ggtitle(titles[i])
    }
    print(plot)
  }}
traceplots_emiss_var<-function(model1,model2,dep=1,true=NULL,titles=NULL,colors=NULL){
  model1_em<-as.data.frame(model1$emiss_var_bar[[dep]])
  model2_em<-as.data.frame(model2$emiss_var_bar[[dep]])
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  for(i in 1:ncol(model2_em)){
    plot<-ggplot(data=model1_em,aes(x=1:data_points, y=model1_em[,i]))+theme_minimal()+xlab("iteration")+ylab(colnames(model2_em)[i])
    if(!is.null(colors)){
      plot<-plot+geom_line(color=colors[1])+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color=colors[2])
    }else{
      plot<-plot+geom_line(color="red")+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color="blue")
    }
    if(!is.null(titles)){
      plot<-plot+ggtitle(titles[i])
    }
    print(plot)
  }

}

traceplots_dwell_mu<-function(model1,model2,true=NULL,ylim1=NULL,titles=NULL,colors=NULL){
  model1_em<-as.data.frame(model1$dwell_mu_bar)
  model2_em<-as.data.frame(model2$dwell_mu_bar)
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  for(i in 1:ncol(model2_em)){
    plot<-ggplot(data=model1_em,aes(x=1:data_points, y=model1_em[,i]))+theme_minimal()+xlab("iteration")+ylab(colnames(model2_em)[i])
    if(!is.null(colors)){
      plot<-plot+geom_line(color=colors[1])+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color=colors[2])
    }else{
      plot<-plot+geom_line(color="red")+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color="blue")
    }
    if(!is.null(ylim1)){
      plot<-plot+ylim(ylim1)
    }
    if(!is.null(titles)){
      plot<-plot+ggtitle(titles[i])
    }

    print(plot)
  }
}
traceplots_dwell_var<-function(model1,model2,true=NULL,ylim1=NULL,titles=NULL,colors=NULL){
  model1_em<-as.data.frame(model1$dwell_var_bar)
  model2_em<-as.data.frame(model2$dwell_var_bar)
  names<-colnames(model2_em)
  data_points<-nrow(model2_em)
  for(i in 1:ncol(model2_em)){
    plot<-ggplot(data=model1_em,aes(x=1:data_points, y=model1_em[,i]))+theme_minimal()+xlab("iteration")+ylab(colnames(model2_em)[i])
    if(!is.null(colors)){
      plot<-plot+geom_line(color=colors[1])+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color=colors[2])
    }else{
      plot<-plot+geom_line(color="red")+geom_line(data=model2_em,aes(x=1:data_points,y=model2_em[,i]),color="blue")
    }
    if(!is.null(ylim1)){
      plot<-plot+ylim(ylim1)
    }
    if(!is.null(titles)){
      plot<-plot+ggtitle(titles[i])
    }
    print(plot)
  }
}

sim_mHMM_plnorm <- function(n_t, n, data_distr = 'poisson-lognormal', m, n_dep = 1,
                            start_state = NULL, q_emiss = NULL, gamma, emiss_distr, xx_vec = NULL, beta = NULL,
                            var_gamma = 0.1, var_emiss = NULL, return_ind_par = FALSE){

  #############
  # Inbuild checks for correct specification of parameters ---------------------
  #############

  if(is.null(xx_vec)){
    xx_vec <- rep(list(NULL), n_dep + 1)
    for(i in 1:(n_dep + 1)){
      xx_vec[[i]] <- rep(1,n)
    }
  } else {
    for(i in 1:(n_dep + 1)){
      if(is.null(xx_vec[[i]])) {
        xx_vec[[i]] <- rep(1,n)
      }
    }
  }
  if(is.null(beta)){
    beta <- rep(list(NULL), n_dep + 1)
    beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
    for(i in 2:(n_dep + 1)){
      if(data_distr == 'categorical'){
        beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
      } else if (data_distr == 'continuous'){
        beta[[i]] <- matrix(0, ncol = 1, nrow = m)
      } else if (data_distr == 'poisson-lognormal'){
        beta[[i]] <- matrix(0, ncol = 1, nrow = m)
      }
    }
  } else {
    if(is.null(beta[[1]])) {
      beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
    }
    for (i in 2:(n_dep + 1)){
      if (is.null(beta[[i]])) {
        if(data_distr == 'categorical'){
          beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
        } else if (data_distr == 'continuous'){
          beta[[i]] <- matrix(0, ncol = 1, nrow = m)
        } else if (data_distr == 'poisson-lognormal'){
          beta[[i]] <- matrix(0, ncol = 1, nrow = m)
        }
      }
    }
  }

  if(data_distr %in% c('poisson-lognormal','continuous')){

    # If a single value of var_gamma specified, use for all transitions
    if(length(var_gamma) == 1){
      var_gamma <- matrix(rep(var_gamma, m*(m-1)),nrow = m, byrow = TRUE)
    } else if(is.matrix(var_gamma)){
      if (dim(var_gamma)[1] != m){
        stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
      }
      if (dim(var_gamma)[2] != m-1){
        stop(paste("The betweem-subject variance matrix for the transition distribution should be a", m, "by", m-1, "matrix."))
      }
    }

    # If a single value of var_emiss specified, use for all states and n_dep
    if(is.null(var_emiss)){
      var_emiss <- rep(list(NULL), n_dep)
      for(i in 1:n_dep){
        var_emiss[[i]] <- matrix(rep(0.1, m), nrow = m, byrow = TRUE)
      }
    } else if(is.numeric(var_emiss) & length(var_emiss) == n_dep){
      arg_var_emiss <- var_emiss
      var_emiss <- rep(list(NULL), n_dep)
      for(i in 1:n_dep){
        var_emiss[[i]] <- matrix(rep(arg_var_emiss[i], m), nrow = m, byrow = TRUE)
      }
    } else if(is.list(var_emiss)){
      for(i in 1:n_dep){
        if(dim(var_emiss[[i]])[1] != m){
          stop(paste0("The number of rows of the between-subject variance for the emission distribution should be
                           equal to ",m,", the number of hidden states chosen."))
        }
        if(dim(var_emiss[[i]])[2] != 1){
          stop(paste0("The number of columns of the between-subject variance for the emission distribution should be
                           equal to one."))
        }
      }
    } else if(length(var_emiss) != n_dep){
      stop("The length of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep. Note that var_emiss can either by a list of matrices, or a numeric vector.")
    }

  }

  #############
  # Simulating the data ---------------------
  #############

  states <- matrix(ncol = 2, nrow = n_t*n)
  states[,1] <- rep(1:n, each = n_t)
  obs <- matrix(ncol = 1 + n_dep, nrow = n_t*n)
  obs[,1] <- rep(1:n, each = n_t)
  sub_gamma <- rep(list(NULL), n)
  sub_emiss <- rep(list(vector("list", n_dep)), n)
  mnl_gamma <- prob_to_int(gamma)
  if(data_distr == "categorical"){
    mnl_emiss <- rep(list(NULL), n_dep)
    for(i in 1:n_dep){
      mnl_emiss[[i]] <- prob_to_int(emiss_distr[[i]])
    }
  }
  for(j in 1:n){
    sub_gamma[[j]] <- int_to_prob(mnl_gamma + xx_vec[[1]][j] * beta[[1]] +
                                    rnorm(n = m * (m-1), mean = 0, sd = sqrt(as.numeric(var_gamma))))
    for(i in 1:n_dep){
      if(data_distr == "categorical"){
        sub_emiss[[j]][[i]] <- int_to_prob(mnl_emiss[[i]] + xx_vec[[1+i]][j] * beta[[1+i]] +
                                             rnorm(n = m * (q_emiss[i]-1), mean = 0, sd = sqrt(as.numeric(var_emiss[[i]]))))
      } else if(data_distr == "continuous"){
        sub_emiss[[j]][[i]] <- emiss_distr[[i]]
        sub_emiss[[j]][[i]][,1] <- emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
          rnorm(n = m, mean = 0, sd = sqrt(as.numeric(var_emiss[[i]])))
      } else if(data_distr == "poisson-lognormal"){
        sub_emiss[[j]][[i]] <- emiss_distr[[i]]
        sub_emiss[[j]][[i]][,1] <- exp(emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
                                         rnorm(n = m, mean = 0, sd = sqrt(as.numeric(var_emiss[[i]]))))
      }
    }

    if(n_t != 0){
      init <- solve(t(diag(m) - sub_gamma[[j]] + 1), rep(1, m))
      if (is.null(start_state)){
        states[((j-1) * n_t + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
      } else {
        states[((j-1) * n_t + 1), 2] <- start_state
      }
      if(data_distr == "categorical"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],])
        }
      } else if (data_distr == "continuous"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],2]))
        }
      } else if (data_distr == "poisson-lognormal"){
        for(i in 1:n_dep){
          obs[((j-1) * n_t + 1), (1+i)] <- rpois(1, lambda = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1])
        }
      }
      for(t in 2:n_t){
        states[((j-1) * n_t + t), 2] <- sample(x = 1:m, size = 1, prob = sub_gamma[[j]][states[((j-1) * n_t + t - 1), 2],])
        if(data_distr == "categorical"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],])
          }
        } else if (data_distr == "continuous"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],2]))
          }
        } else if (data_distr == "poisson-lognormal"){
          for(i in 1:n_dep){
            obs[((j-1) * n_t + t), (1+i)] <- rpois(1, lambda = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1])
          }
        }
      }
    }
  }

  #############
  # Returning output  ---------------------
  #############
  colnames(states) <- c("subj", "state")
  colnames(obs)    <- c("subj", paste("observation", 1:n_dep))
  if (return_ind_par == FALSE & n_t != 0){
    return(list(states = states, obs = obs))
  } else if (return_ind_par == TRUE & n_t != 0){
    return(list(states = states, obs = obs, subject_gamma = sub_gamma, subject_emiss = sub_emiss))
  } else if (n_t == 0){
    return(list(subject_gamma = sub_gamma, subject_emiss = sub_emiss))
  }
}
