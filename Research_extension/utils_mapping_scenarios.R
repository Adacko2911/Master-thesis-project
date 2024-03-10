#Utils function
#Create function to calculate between subject variance
get_logvar <- function(mu, varmu){
  logmu = log(mu)
  log(0.5*exp(-2*logmu)*(exp(2*logmu) + sqrt(4*exp(2*logmu)*varmu+exp(4*logmu))))
}


get_varmu <- function(lambda, logvar){
  logmu = log(lambda)
  abs(exp(logvar)-1)*exp(2*logmu+logvar)
}

#subjects
# scenario 1 : 25
# scenario 2 : 50
# scenario 3 : 250
# scenario 4 : 500
n_distr <- list(25, 50, 250, 500)

#time series length
# scenario 1 : 250
# scenario 2 : 500
# scenario 3 : 1000
n_t_distr <- list(250, 500, 1000)

# dwell time lambdas
# scenario 1 : 2nd state lambda=10 5% overlap
# scenario 2 : 2nd state lambda=10 25% overlap
# scenario 3 : 2nd state lambda=10 50% overlap
# scenario 4 : equal across states lambda=10
# scenario 5 : 2nd state lambda=30 5% overlap
# scenario 6 : 2nd state lambda=30 25% overlap
# scenario 7 : 2nd state lambda=30 50% overlap
# scenario 8 : equal across states lambda=30
# scenario 9 : 2nd state lambda=90 5% overlap
# scenario 10 : 2nd state lambda=90 25% overlap
# scenario 11 : 2nd state lambda=90 50% overlap
# scenario 12 : equal across states lambda=30
dwell_assumed<-data.frame(mean=c(rep(10,3),rep(30,3),rep(90,3)),overlap=rep(c("5% overlap","25% overlap","50% overlap"),3),state1=c(1,4,6,12,19,23,57,70,78),state2=c(10,10,10,30,30,30,90,90,90),state3=c(26,19,15,55,44,38,131,113,103))

dwell_distr <- list(
  matrix(dwell_assumed[1,3:5],ncol = 1),
  matrix(dwell_assumed[2,3:5],ncol = 1),
  matrix(dwell_assumed[3,3:5],ncol = 1),
  matrix(rep(10,3),ncol = 1),
  matrix(dwell_assumed[4,3:5],ncol = 1),
  matrix(dwell_assumed[5,3:5],ncol = 1),
  matrix(dwell_assumed[6,3:5],ncol = 1),
  matrix(rep(30,3),ncol = 1),
  matrix(dwell_assumed[7,3:5],ncol = 1),
  matrix(dwell_assumed[8,3:5],ncol = 1),
  matrix(dwell_assumed[9,3:5],ncol = 1),
  matrix(rep(90,3),ncol = 1)

)
names(dwell_distr)<-rep(c("5% overlap","25% overlap","50% overlap","100% overlap"),3)

# dwell time between subjects logvariance
# variance between subjects should be equal to the lambda in the second state
dwell_var_ss_s<-list(get_logvar(mu=dwell_assumed[1,3:5], varmu=10),
                     get_logvar(mu=dwell_assumed[2,3:5], varmu=10),
                     get_logvar(mu=dwell_assumed[3,3:5], varmu=10),
                     get_logvar(mu=rep(10,3), varmu=10),
                     get_logvar(mu=dwell_assumed[4,3:5], varmu=30),
                     get_logvar(mu=dwell_assumed[5,3:5], varmu=30),
                     get_logvar(mu=dwell_assumed[6,3:5], varmu=30),
                     get_logvar(mu=rep(30,3), varmu=30),
                     get_logvar(mu=dwell_assumed[7,3:5], varmu=90),
                     get_logvar(mu=dwell_assumed[8,3:5], varmu=90),
                     get_logvar(mu=dwell_assumed[8,3:5], varmu=90),
                     get_logvar(mu=rep(90,3), varmu=90))



### ??? Placeholder for gamma scenario with varying between states transitions
### I could find the code I used to calculate the gamma between subjest variance
### Probably it won't be 0.17

gamma_ss_var<-0.17


#emission distribution
# scenario 1 :  50% overlap variance=30 one variable
# scenario 2 :  25% overlap variance=30 one variable
# scenario 3 :  5% overlap variance=30 one variable
# scenario 4 : 50% overlap variance=30 two identical variables
# scenario 5 : 25% overlap variance=30 two identical variables
# scenario 6 :  5% overlap variance=30 two identical variables
# scenario 7 : 50% overlap variance=30 two different variables
# scenario 8 : 25% overlap variance=30 two different variables
# scenario 9 :  5% overlap variance=30 two different variables

emiss_ss_var<-10


fty_pct<-c(10,17.38,24.76)
fty_pct_r<-c(24.76,17.38,10)
tty_pct<-c(10,22.60,35.2)
tty_pct_r<-c(35.2,22.60,10)
f_pct<-c(10,31.47,52.94)
f_pct_r<-c(52.94,31.47,10)
overlap_names=c("50% overlap","25% overlap","5% overlap")
# dependent variable distribution
dep_univariate<-data.frame(rbind(fty_pct,tty_pct,f_pct),row.names =  overlap_names)
colnames(dep_univariate)<-c("state1","state2","state3")
dep_multivariate1<-list(dep_univariate,
                        dep_univariate)
alt_uni<-data.frame(rbind(fty_pct_r,tty_pct_r,f_pct_r),row.names = overlap_names)
colnames(alt_uni)<-c("state1","state2","state3")
dep_multivariate2<-list(dep_univariate,
                        alt_uni)

dep_distr<-list(matrix(c(dep_univariate[1,],rep(30,3)),ncol=2),
                matrix(c(dep_univariate[2,],rep(30,3)),ncol=2),
                matrix(c(dep_univariate[3,],rep(30,3)),ncol=2),
                list(matrix(c(dep_multivariate1[[1]][1,],rep(30,3)),ncol=2),matrix(c(dep_multivariate1[[2]][1,],rep(30,3)),ncol=2)),
                list(matrix(c(dep_multivariate1[[1]][2,],rep(30,3)),ncol=2),matrix(c(dep_multivariate1[[2]][2,],rep(30,3)),ncol=2)),
                list(matrix(c(dep_multivariate1[[1]][3,],rep(30,3)),ncol=2),matrix(c(dep_multivariate1[[2]][3,],rep(30,3)),ncol=2)),
                list(matrix(c(dep_multivariate2[[1]][1,],rep(30,3)),ncol=2),matrix(c(dep_multivariate2[[2]][1,],rep(30,3)),ncol=2)),
                list(matrix(c(dep_multivariate2[[1]][2,],rep(30,3)),ncol=2),matrix(c(dep_multivariate2[[2]][2,],rep(30,3)),ncol=2)),
                list(matrix(c(dep_multivariate2[[1]][3,],rep(30,3)),ncol=2),matrix(c(dep_multivariate2[[2]][3,],rep(30,3)),ncol=2)))

names(dep_distr)<-1:9

#clean up
rm(dwell_assumed)
rm(dep_univariate)
rm(dep_multivariate1)
rm(dep_multivariate2)
rm(fty_pct)
rm(fty_pct_r)
rm(tty_pct_r)
rm(tty_pct)
rm(f_pct)
rm(f_pct_r)
