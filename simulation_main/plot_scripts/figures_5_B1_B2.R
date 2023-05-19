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
library(ggshadow)

library(tidyverse)

list_post_est <- read_csv("~/GitHub/Master-thesis-project/Simulation Study/simulation study final/post_est_final.csv") %>% drop_na() %>% filter(dwell_time!="mix")

####EMISS MU BAR
data<-list_post_est

data$model<-as.factor(data$model)
data$state<-as.factor(data$state)
data$observations<-as.factor(data$observations)
data$dwell_time<-factor(data$dwell_time,levels = unique(data$dwell_time),labels=unique(data$dwell_time),ordered = TRUE)
by=c("observations","state","dwell_time")
step_plot_function<-function(data=list_post_est,parameter_name=names,labels = c("mu[1]","mu[2]","mu[3]","mu[4]"),stats="bias_median",order_conditions=c("observations","state","dwell_time"),linetype="model",stats_nam="bias", colornames=c(colors$blue,colors$orange,colors$red,colors$teal),target=0,title1="Emission distribution means: Bias",subtitle1 = " ",col_lab="Estimates:"){
  by<-order_conditions
  data<-data %>% filter(estimates %in% parameter_name)
  data$estimates<-factor(data$estimates,levels = unique(parameter_name),labels = labels)
  data$model<-as.factor(data$model)
  data$state<-as.factor(data$state)
  data$observations<-as.factor(data$observations)
  data$dwell_time<-factor(data$dwell_time,levels = unique(data$dwell_time),labels=unique(data$dwell_time),ordered = TRUE)
  stats=stats

  ### Compute internal data
  opts <- lapply(X = by, FUN = function(x) levels(data[[x]]))
  names(opts) <- by
  dgms <- do.call(expand.grid, opts)
  dgms[[".scenario"]] <- as.numeric(seq(nrow(dgms)))
  data <- merge(x = data, y = dgms)
  data <- data[order(data[[".scenario"]]), ]

  ### Compute limits and placement of nested loop plot labels
  limits <- range(data[[stats]] , na.rm = TRUE)

  delta <- diff(limits) / 10
 if(target<limits[2]){
    limits[2]<-target
  }
  placement <- vector(mode = "list", length = length(by))
  for (i in seq_along(placement)) {
    if (i == 1) {
      if(length(unique(range(data[[stats]], na.rm = TRUE)))==1){
        placement[[i]] <- c(round(limits[1], digits = 2) - 4 * delta, round(limits[1], digits = 2) - delta)
      }else{
        placement[[i]] <- c(round(limits[1], digits = 2) - 2 * delta, round(limits[1], digits = 2) - delta)
      }
    } else {
      placement[[i]] <- c(placement[[i - 1]][1] - 2 * delta, placement[[i - 1]][1] - delta)
    }
  }

  ### Rescale variables included in the nested loop plot
  for (i in seq_along(by)) {
    data[[paste0(".", by[i])]] <- scales::rescale(x = as.numeric(data[[by[i]]]), to = placement[[i]])
  }

  if(is.numeric(target)){
    data[["target"]]<-rep(target,nrow(data))
  }
  ###Plotting
  gg<-data %>%  ggplot(., mapping = ggplot2::aes(x = .scenario, y = get(stats))) +
    geom_shadowstep(aes(color=data$estimates,linetype=get(linetype)),size=1.2) +
    ggplot2::labs(x = paste0(paste(vapply(X = by, FUN = function(x) length(levels(data[[x]])), FUN.VALUE = numeric(1)), collapse = " x "), " = ", max(data[[".scenario"]]), " ordered scenarios"),y=stats_nam,linetype="Model: ")+
    ggplot2::geom_step(data = data, mapping = ggplot2::aes(x = .scenario, y = target),alpha=0.85, linetype = "dashed",size=1)+
    scale_color_manual(values =colornames ,name = col_lab)
  for (i in seq_along(by)) {
    .tmp <- rlang::sym(paste0(".", by[i]))
    gg <- gg +
      ggplot2::geom_step(mapping = ggplot2::aes(y = !!.tmp),size=1,alpha=0.5) +
      ggplot2::annotate(geom = "text", x = 1, y = placement[[i]][2] + delta /2.5, label = paste0(by[i], ": ", paste(levels(data[[by[i]]]), collapse = ", ")), hjust = 0, vjust = 0.5)+
      theme_minimal()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(hjust=0.65),
            axis.text.y = element_text(size=12),
            legend.position="top")+ggtitle(paste0(title1),subtitle = subtitle1)
  }
 return(gg)
}

library(patchwork)
####BIAS#######################################################################

#Figure 5
list_post_est2<-list_post_est %>% mutate(bias_median_pp = ifelse(bias_median<0,-bias_median_pp,bias_median_pp))
names<-c("emiss_mu_bar_1","emiss_mu_bar_2","emiss_mu_bar_3","emiss_mu_bar_4")
p1<-step_plot_function(list_post_est2 %>% filter(dwell_time=="1.4" | dwell_time =="3.5"),parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="bias_median_pp",target = 0,title1="The 1st dependent variable",col_lab="Estimates:"  )
p1a<-step_plot_function(list_post_est2 %>% filter(dwell_time=="19.5" | dwell_time =="99.5"),parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="bias_median_pp",target = 0,title1="",col_lab="Estimates:"  )
names<-c("emiss_mu_bar2_1","emiss_mu_bar2_2","emiss_mu_bar2_3","emiss_mu_bar2_4")
p2<-step_plot_function(list_post_est2 %>% filter(dwell_time=="1.4" | dwell_time =="3.5"),parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="bias_median_pp",target = 0,title1= "The 2nd dependent variable" ,col_lab="Estimates:" )
p2a<-step_plot_function(list_post_est2 %>% filter(dwell_time=="19.5" | dwell_time =="99.5"),parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="bias_median_pp",target = 0,title1="",col_lab="Estimates:"  )

pdf("~/GitHub/Master-thesis-project/Simulation Study/simulation study final/outputs plots/emiss_dep_bias2.pdf" ,height = 9, width = 11)
combined<-(p1+ p1a)/(p2+p2a)& theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
dev.off()

###Coverage######################################################################
names<-c("emiss_mu_bar_1","emiss_mu_bar_2","emiss_mu_bar_3","emiss_mu_bar_4")
step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="coverage_pp",target = 100,title1="Emission distribution means: Coverage",subtitle1 = "1st dependent variable"  )


names<-c("emiss_mu_bar2_1","emiss_mu_bar2_2","emiss_mu_bar2_3","emiss_mu_bar2_4")
step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="coverage_pp",target = 100,title1="Emission distribution means: Coverage",subtitle1 = "2nd dependent variable"  )

names<-c("emiss_var_bar_1","emiss_var_bar_2","emiss_var_bar_3","emiss_var_bar_4")
step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="coverage_pp",target = 100,title1="Emission distribution variance: Coverage",subtitle1 = "1st dependent variable"  )


names<-c("emiss_var_bar2_1","emiss_var_bar2_2","emiss_var_bar2_3","emiss_var_bar2_4")
step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="coverage_pp",target = 100,title1="Emission distribution variance: Coverage",subtitle1 = "2nd dependent variable"  )

names<-c("emiss_varmu_bar_1","emiss_varmu_bar_2","emiss_varmu_bar_3","emiss_varmu_bar_4")
step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="coverage_pp",target = 100,title1="Emission distribution ss-variance: Coverage",subtitle1 = "1st dependent variable"  )


names<-c("emiss_varmu_bar2_1","emiss_varmu_bar2_2","emiss_varmu_bar2_3","emiss_varmu_bar2_4")
step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="coverage_pp",target = 100,title1="Emission distribution ss-variance: Coverage",subtitle1 = "2nd dependent variable"  )
##SD##############################################################################
names<-c("emiss_mu_bar_1","emiss_mu_bar_2","emiss_mu_bar_3","emiss_mu_bar_4")
p1<-step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="sd",target = 0,title1="Emission distribution means: SD",subtitle1 = "1st dependent variable"  )


names<-c("emiss_mu_bar2_1","emiss_mu_bar2_2","emiss_mu_bar2_3","emiss_mu_bar2_4")
p2<-step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="sd",target = 0,title1="Emission distribution means: SD",subtitle1 = "2nd dependent variable"  )

names<-c("emiss_var_bar_1","emiss_var_bar_2","emiss_var_bar_3","emiss_var_bar_4")
p3<-step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="sd",target = 0,title1="Emission distribution variance: SD",subtitle1 = "1st dependent variable"  )


names<-c("emiss_var_bar2_1","emiss_var_bar2_2","emiss_var_bar2_3","emiss_var_bar2_4")
p4<-step_plot_function(list_post_est,parameter_name =names,labels = c("mu1","mu2","mu3","mu4"),stats ="sd",target = 0,title1="Emission distribution variance: SD",subtitle1 = "2nd dependent variable"  )

pdf("~/GitHub/Master-thesis-project/Simulation Study/simulation study final/outputs plots/emiss_sd.pdf" ,height = 15, width = 15)
combined <- (p1+p2)/(p3+p4) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
 dev.off()

####coverage results plotting
data<-list_post_est

data$model<-as.factor(data$model)
data$state<-as.factor(data$state)
data$observations<-as.factor(data$observations)
data$dwell_time<-factor(data$dwell_time,levels = unique(data$dwell_time),labels=unique(data$dwell_time),ordered = TRUE)
new_data<-data %>% filter(estimates %in% c("emiss_mu_bar_1","emiss_mu_bar_2","emiss_mu_bar_3","emiss_mu_bar_4")) %>% filter(dwell_time==1.4) %>% mutate(dwell_time=rep(100,nrow(.))) %>% mutate(coverage_pp=rep(0,nrow(.)))
p1<-data %>% filter(estimates %in% c("emiss_mu_bar_1","emiss_mu_bar_2","emiss_mu_bar_3","emiss_mu_bar_4"))%>% rbind(.,new_data)%>%
  ggplot(., aes(x=estimates, y=coverage_pp, fill=dwell_time,group=model)) +
  geom_bar (stat="identity", position =position_dodge2(width = 1, preserve = "total", padding = -0.1),alpha=0.9)+
  facet_grid(cols=vars(state),rows=vars(observations),switch = "y",scales = "free_x")+
  theme_bw()+labs(title = "Dependent variable 1", y="",x="")+
  scale_fill_manual(values=c(colors$blue,colors$orange,colors$red,colors$teal))+scale_y_continuous(position = "right")+
  theme(legend.position = "bottom",axis.text.y = element_blank(),axis.ticks.y = element_blank())

data<-list_post_est

data$model<-as.factor(data$model)
data$state<-as.factor(data$state)
data$observations<-as.factor(data$observations)
data$dwell_time<-factor(data$dwell_time,levels = unique(data$dwell_time),labels=unique(data$dwell_time),ordered = TRUE)
new_data<-data %>% filter(estimates %in% c("emiss_mu_bar2_1","emiss_mu_bar2_2","emiss_mu_bar2_3","emiss_mu_bar2_4")) %>% filter(dwell_time==1.4) %>% mutate(dwell_time=rep(100,nrow(.))) %>% mutate(coverage_pp=rep(0,nrow(.)))
p2<-data %>% filter(estimates %in% c("emiss_mu_bar2_1","emiss_mu_bar2_2","emiss_mu_bar2_3","emiss_mu_bar2_4"))%>% rbind(.,new_data)%>%
  ggplot(., aes(x=estimates, y=coverage_pp, fill=dwell_time,group=model)) +
  geom_bar (stat="identity", position =position_dodge2(width = 1, preserve = "total", padding = -0.1),alpha=0.9)+
  facet_grid(cols=vars(state),rows=vars(observations),switch = "y",scales = "free_x")+
  theme_bw()+labs(title = "Dependent variable 2")+
  scale_fill_manual(values=c(colors$blue,colors$orange,colors$red,colors$teal))+scale_y_continuous(position = "right")+
  theme(legend.position = "bottom")
pdf("~/GitHub/Master-thesis-project/Simulation Study/simulation study final/outputs plots/emiss_mu_cov.pdf" ,height = 10, width = 15)
com<-p1+p2 & theme(legend.position = "bottom")
com+plot_layout(guides = "collect")
dev.off()

