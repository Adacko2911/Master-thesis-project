#----------Start of the data simulating file-----------------------------------------------------------------------
# This file is the main simulation file. Note that in order to reproduce the data that can be found in my thesis you need to run this and also
# second data_simulation_2.R R script.
library(R.utils)
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
mhsmm_sim <-
  function(n,
           n_t,
           m,
           n_dep,
           emiss,
           gamma,
           gamma_ss_var,
           emiss_ss_var,
           dwell_ss_var,
           dwell_distr) {
    #,dwell_ss_var
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
            sd = sqrt(emiss_ss_var)
          ), nrow = m)
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
      mu_dwell<-c()
      for(s in 1:m){
        mu_dwell[s]<-dwell_distr[s, 1] + rnorm(1, mean = 0, sd = sqrt(dwell_ss_var))
        if(mu_dwell[s]<0){
          mu_dwell[s]=0
        }
      }

      dwell_ss[[subj]] <- cbind(mu_dwell,dwell_distr[, 2])
      while (a < n_t) {

        #draw the subject specific discrete dwell time
        mean_dwell <-
          round(rlnorm(1, dwell_ss[[subj]][dwell_row,1], sqrt(dwell_ss[[subj]][dwell_row,2])), 0)
        #condition to meet the n_t subject observations sequence length
        if (mean_dwell <= 1) {
          mean_dwell = 1
        }else{
          if (a + (mean_dwell - 1) <= n_t) {
            #replete the most probable state
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
                  sd = sqrt(emiss_ss[[subj]][[dep]][state_sequence[subj * n_t - n_t + t_occ, 2], 2]))
        }#end of dependent variable loop
      }#end of series loop
    }#end of subject loop
    #the output is return here
    state_sequence<-matrix(as.integer(state_sequence),ncol=2)
    output <- list( state_sequence,  obs_ss,gamma_ss,emiss_ss,dwell_ss)   #
    names(output) <- c( "states",  "observations","gamma_ss","emiss_ss","dwell_ss")#
    return(output)
  }

#vary the transitions over states. Vary the duration across scenarios but keep it constant accros the states.
#12.5,37.5,62.5,87.5
#set the seed
set.seed(124)
n <- 80 #subjects
n_t <- c(200, 500, 1000) #time series length
m <- 3 #number of states

#number of set of mean duration matrices
md_s <- 4
#number of scenarios of the number of observations per individual
occ_s <- 3
#number of samples within the same scenario
sampl <- 128
#We will not vary the number of dependent variable hence n_dep is of length 1 and the emiss is 1 list
n_dep <- 2


emiss <- list(matrix(c(10,20,
                       30, 60,
                       60,120), nrow = m, byrow = TRUE),
              matrix(c(60,120+0.25*120,
                       30, 60+0.25*60,
                       10,20+20*0.25), nrow = m, byrow = TRUE)
)
emiss_ss_var <- 9
#gamma is a list of length st_s
gamma <- list(matrix(
  data = c(0, 0.5,0.5,
           0.5, 0,0.5,
           0.5,0.5,0),
  byrow = T,
  nrow = m
))
gamma_ss_var <- 0.17

#calculate according expected dwell times
exp_mean1=1/(-log(0.5))
exp_mean2=1/(-log(0.75))
exp_mean3=1/(-log(0.95))
exp_mean4=1/(-log(0.99))
#dwell_distr is a list of length md_s
dwell_distr <- list(
  matrix(data = log(c(exp_mean1, 1.10093,
                      exp_mean1,  1.10093,
                      exp_mean1,  1.10093)),
         byrow = T,

         ncol = 2,
         nrow = m
  ),
  matrix(
    data = log(c(exp_mean2, 1.10093,
                 exp_mean2,  1.10093,
                 exp_mean2,  1.10093)),
    byrow = T,
    ncol = 2,
    nrow = m
  ),

  matrix(
    data = log(c(exp_mean3, 1.10093,
                 exp_mean3,  1.10093,
                 exp_mean3,  1.10093)),
    byrow = T,
    ncol = 2,
    nrow = m
  ),
  matrix(
    data = log(c(exp_mean4, 1.05902,
                 exp_mean4,  1.05902,
                 exp_mean4, 1.05902)), #(1/-log(0.99)*1/4)^2
    byrow = T,
    ncol = 2,
    nrow = m
  )
)
dwell_ss_var=0.2
state_sc_name <- c()
nr_oc_sc_name <- c()
meandwell_sc_name <- c()
for (i in 1:m) {
  state_sc_name[i] <- paste("(states) scenario ", i, ",with m=", m)
}
for (i in 1:3) {
  nr_oc_sc_name[i] <-
    paste("(no. occasions) scenario ", i, ",with n_t=", n_t[i])
}
for (i in 1:4) {
  meandwell_sc_name[i] <- paste("(dwell time) scenario ", i)
}
three_state_design_simulated_data <- list()
md_list <- list()
sim_final <- list()
list1 <- list()
#first we got a list of state number scenarios
#----second we got number of samples within each person scenarios
#----------third we got list of state specific means scenarios
#----------------forth we got multiple samples of 1 of the scenarios
#----------------------fifth each sample consist of output of mhsmm_sim()
for (occ in 1:occ_s) {
  n_occ <- n_t[occ]
  for (dwell in 1:md_s) {
    dwel <- dwell_distr[[dwell]]
    for (iter in 1:sampl) {
      list1[[iter]] <-
        mhsmm_sim(
          n = n,
          n_t = n_occ,
          m = m,
          n_dep = 2,
          emiss = emiss,
          gamma = gamma[[1]] ,
          gamma_ss_var = gamma_ss_var,
          emiss_ss_var = emiss_ss_var,
          dwell_distr = dwel,
          dwell_ss_var=dwell_ss_var
        )
    }
    md_list[[dwell]] <- list1
    list1 <- list()
  }
  names(md_list) <- meandwell_sc_name
  three_state_design_simulated_data[[occ]] <- md_list
  md_list <- list()
}
names(three_state_design_simulated_data) <- nr_oc_sc_name

saveRDS(three_state_design_simulated_data,file = "./data/three_state_design_2dep_final.rds")
rm(three_state_design_simulated_data)
#-------------end Three state scenario------------------------------------------------------------------------------------


#-------------Four state scenario------------------------------------------------------------------------------
mhsmm_sim <-
  function(n,
           n_t,
           m,
           n_dep,
           emiss,
           gamma,
           gamma_ss_var,
           emiss_ss_var,
           dwell_distr,
           dwell_ss_var) {
    #,dwell_ss_var
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
            sd = sqrt(emiss_ss_var)
          ), nrow = m)
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
      mu_dwell<-c()
      for(s in 1:m){
        mu_dwell[s]<-dwell_distr[s, 1] + rnorm(1, mean = 0, sd = sqrt(dwell_ss_var))
        if(mu_dwell[s]<0){
          mu_dwell[s]=0
        }
      }

      dwell_ss[[subj]] <- cbind(mu_dwell,dwell_distr[, 2])
      while (a < n_t) {

        #draw the subject specific discrete dwell time
        mean_dwell <-
          round(rlnorm(1, dwell_ss[[subj]][dwell_row,1], sqrt(dwell_ss[[subj]][dwell_row,2])), 0)
        #condition to meet the n_t subject observations sequence length
        if (mean_dwell <= 1) {
          mean_dwell = 1
        }else{
          if (a + (mean_dwell - 1) <= n_t) {
            #replete the most probable state
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
                  sd = sqrt(emiss_ss[[subj]][[dep]][state_sequence[subj * n_t - n_t + t_occ, 2], 2]))
        }#end of dependent variable loop
      }#end of series loop
    }#end of subject loop
    #the output is return here
    state_sequence<-matrix(as.integer(state_sequence),ncol=2)
    output <- list( state_sequence,  obs_ss,gamma_ss,emiss_ss,dwell_ss)   #
    names(output) <- c("states", "observations","gamma_ss", "emiss_ss", "dwell_ss")  #
    return(output)
  }


#vary the transitions over states. Vary the duration across scenarios but keep it constant accros the states.
#12.5,37.5,62.5,87.5
#set the seed
set.seed(124)
n <- 80 #subjects
n_t <- c(200, 500, 1000) #time series length
m <- 4 #number of states

#number of set of mean duration matrices
md_s <- 4
#number of scenarios of the number of observations per individual
occ_s <- 3
#number of samples within the same scenario
sampl <- 128
#We will not vary the number of dependent variable hence n_dep is of length 1 and the emiss is 1 list
n_dep <- 2



emiss <- list(matrix(c(10,20,
                       30, 60,
                       60,120,
                       98,150), nrow = m, byrow = TRUE),
              matrix(c(98,(150+0.25*150),
                       60, (120+0.25*120),
                       30,(60+0.25*60),
                       10,(20+0.25*20)), nrow = m, byrow = TRUE))
emiss_ss_var <- 9
#gamma is a list of length st_s
gamma <- list(matrix(
  data = c(0, 1/3,1/3,1/3,
           1/3, 0,1/3,1/3,
           1/3,1/3,0,1/3,
           1/3,1/3,1/3,0),
  byrow = T,
  nrow = m
))
gamma_ss_var <- matrix(c(0.17,0.17,
                       0.17,0.17,
                       0.17,0.17,
                       0.17,0.17),nrow = m,ncol = m-2)
#dwell_distr is a list of length md_s
#calculate according expected dwell times
exp_mean1=1/(-log(0.5))
exp_mean2=1/(-log(0.75))
exp_mean3=1/(-log(0.95))
exp_mean4=1/(-log(0.99))
#dwell_distr is a list of length md_s
dwell_distr <- list(
  matrix(data = log(c(exp_mean1, 1.10093,
                      exp_mean1,  1.10093,
                      exp_mean1,  1.10093,
                      exp_mean1,  1.10093)),
         byrow = T,

         ncol = 2,
         nrow = m
  ),
  matrix(
    data = log(c(exp_mean2, 1.10093,
                 exp_mean2,  1.10093,
                 exp_mean2,  1.10093,
                 exp_mean2,  1.10093)),
    byrow = T,
    ncol = 2,
    nrow = m
  ),

  matrix(
    data = log(c(exp_mean3, 1.10093,
                 exp_mean3,  1.10093,
                 exp_mean3,  1.10093,
                 exp_mean3,  1.10093)),
    byrow = T,
    ncol = 2,
    nrow = m
  ),
  matrix(
    data = log(c(exp_mean4, 1.05902,
                 exp_mean4,  1.05902,
                 exp_mean4,  1.05902,
                 exp_mean4, 1.05902)), #(1/-log(0.97)*1/3)^2
    byrow = T,
    ncol = 2,
    nrow = m
  )
)
dwell_ss_var=0.2
state_sc_name <- c()
nr_oc_sc_name <- c()
meandwell_sc_name <- c()
for (i in 1:m) {
  state_sc_name[i] <- paste("(states) scenario ", i, ",with m=", m)
}
for (i in 1:3) {
  nr_oc_sc_name[i] <-
    paste("(no. occasions) scenario ", i, ",with n_t=", n_t[i])
}
for (i in 1:4) {
  meandwell_sc_name[i] <- paste("(dwell time) scenario ", i)
}
four_state_design_simulated_data <- list()
md_list <- list()
sim_final <- list()
list1 <- list()
#----second we got number of samples within each person scenarios
#----------third we got list of state specific means scenarios
#----------------forth we got multiple samples of 1 of the scenarios
#----------------------fifth each sample consist of output of mhsmm_sim()
for (occ in 1:occ_s) {
  n_occ <- n_t[occ]
  for (dwell in 1:md_s) {
    dwel <- dwell_distr[[dwell]]
    for (iter in 1:sampl) {
      list1[[iter]] <-
        mhsmm_sim(
          n = n,
          n_t = n_occ,
          m = m,
          n_dep = 2,
          emiss = emiss,
          gamma = gamma[[1]] ,
          gamma_ss_var = gamma_ss_var,
          emiss_ss_var = emiss_ss_var,
          dwell_distr = dwel,
          dwell_ss_var=dwell_ss_var
        )
    }
    md_list[[dwell]] <- list1
    list1 <- list()
  }
  names(md_list) <- meandwell_sc_name
  four_state_design_simulated_data[[occ]] <- md_list
  md_list <- list()
}
names(four_state_design_simulated_data) <- nr_oc_sc_name
saveRDS(four_state_design_simulated_data,file = "./data/four_state_design_2dep_final.rds")
rm(four_state_design_simulated_data)

#-------------end of Four state scenario------------------------------------------------------------------------------



#---------------additional example -------------------------------------------------------------------------------
#remove all variables

rm(list=ls())

mhsmm_sim <-
  function(n,
           n_t,
           m,
           n_dep,
           emiss,
           gamma,
           gamma_ss_var,
           emiss_ss_var,
           dwell_ss_var,
           dwell_distr) {
    #,dwell_ss_var
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
            sd = sqrt(emiss_ss_var)
          ), nrow = m)
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
      mu_dwell<-c()
      for(s in 1:m){
        mu_dwell[s]<-dwell_distr[s, 1] + rnorm(1, mean = 0, sd = sqrt(dwell_ss_var))
        if(mu_dwell[s]<0){
          mu_dwell[s]=0
        }
      }

      dwell_ss[[subj]] <- cbind(mu_dwell,dwell_distr[, 2])
      while (a < n_t) {

        #draw the subject specific discrete dwell time
        mean_dwell <-
          round(rlnorm(1, dwell_ss[[subj]][dwell_row,1], sqrt(dwell_ss[[subj]][dwell_row,2])), 0)
        #condition to meet the n_t subject observations sequence length
        if (mean_dwell <= 1) {
          mean_dwell = 1
        }else{
          if (a + (mean_dwell - 1) <= n_t) {
            #replete the most probable state
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
                  sd = sqrt(emiss_ss[[subj]][[dep]][state_sequence[subj * n_t - n_t + t_occ, 2], 2]))
        }#end of dependent variable loop
      }#end of series loop
    }#end of subject loop
    #the output is return here
    state_sequence<-matrix(as.integer(state_sequence),ncol=2)
    output <- list( state_sequence,  obs_ss,gamma_ss,emiss_ss,dwell_ss)   #
    names(output) <- c( "states",  "observations","gamma_ss","emiss_ss","dwell_ss")#
    return(output)
  }

#vary the transitions over states. Vary the duration across scenarios but keep it constant accros the states.
#12.5,37.5,62.5,87.5
#set the seed
set.seed(124)
n <- 80 #subjects
n_t <- c(500, 1000) #time series length
m <- 3 #number of states

#number of set of mean duration matrices
md_s <- 1
#number of scenarios of the number of observations per individual
occ_s <- 2
#number of samples within the same scenario
sampl <- 128
#We will not vary the number of dependent variable hence n_dep is of length 1 and the emiss is 1 list
n_dep <- 2


emiss <- list(matrix(c(10,20,
                       30, 60,
                       60,120), nrow = m, byrow = TRUE),
              matrix(c(60,120+0.25*120,
                       30, 60+0.25*60,
                       10,20+20*0.25), nrow = m, byrow = TRUE)
)
emiss_ss_var <- 9
#gamma is a list of length st_s
gamma <- list(matrix(
  data = c(0, 0.5,0.5,
           0.5, 0,0.5,
           0.5,0.5,0),
  byrow = T,
  nrow = m
))
gamma_ss_var <- 0.17

#calculate according expected dwell times
exp_mean1=1/(-log(0.5))
exp_mean2=1/(-log(0.75))
exp_mean3=1/(-log(0.95))
exp_mean4=1/(-log(0.99))
#dwell_distr is a list of length md_s
dwell_distr <- list(
  matrix(data = log(c(exp_mean2, 1.10093,
                      exp_mean3,  1.10093,
                      exp_mean4,  1.05902)),
         byrow = T,

         ncol = 2,
         nrow = m
  )
)
dwell_ss_var=0.2
state_sc_name <- c()
nr_oc_sc_name <- c()
meandwell_sc_name <- c()
for (i in 1:m) {
  state_sc_name[i] <- paste("(states) scenario ", i, ",with m=", m)
}
for (i in 1:2) {
  nr_oc_sc_name[i] <-
    paste("(no. occasions) scenario ", i, ",with n_t=", n_t[i])
}
for (i in 1:1) {
  meandwell_sc_name[i] <- paste("(dwell time) scenario ", i)
}
three_state_design_simulated_data <- list()
md_list <- list()
sim_final <- list()
list1 <- list()
#first we got a list of state number scenarios
#----second we got number of samples within each person scenarios
#----------third we got list of state specific means scenarios
#----------------forth we got multiple samples of 1 of the scenarios
#----------------------fifth each sample consist of output of mhsmm_sim()
for (occ in 1:occ_s) {
  n_occ <- n_t[occ]
  for (dwell in 1:md_s) {
    dwel <- dwell_distr[[dwell]]
    for (iter in 1:sampl) {
      list1[[iter]] <-
        mhsmm_sim(
          n = n,
          n_t = n_occ,
          m = m,
          n_dep = 2,
          emiss = emiss,
          gamma = gamma[[1]] ,
          gamma_ss_var = gamma_ss_var,
          emiss_ss_var = emiss_ss_var,
          dwell_distr = dwel,
          dwell_ss_var=dwell_ss_var
        )
    }
    md_list[[dwell]] <- list1
    list1 <- list()
  }
  names(md_list) <- meandwell_sc_name
  three_state_design_simulated_data[[occ]] <- md_list
  md_list <- list()
}
names(three_state_design_simulated_data) <- nr_oc_sc_name

saveRDS(three_state_design_simulated_data,file = "./data/three_state_design_mixed_dwell_2dep_final.rds")
rm(three_state_design_simulated_data)


#----------End of the data simulating file-----------------------------------------------------------------------

#remove all variables

rm(list=ls())

