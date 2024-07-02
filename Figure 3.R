setwd("~/EpiEstim/Epi Methods/Code Final")
source("Model_synthetic structure.R")

##################Choose parameter values and simulate model
rho=0.1
k=2
N=35039
20785/35039
14254/35039
Ni=floor(35039*0.6)
Nj=floor(35039*0.4)
Ei0=20
Ej0=20
Ii0=18
Ij0=18
Ri0=0
Rj0=0



one <- data.frame(run=seq(1,6,1),dist=rep(1,6),betai=c(4*10^-5,rep(4*10^-5,5)),
                     betaj=c(4*10^-5,rep(6*10^-5,5)),
                     betaij=c(4*10^-5,1*10^-5,1*10^-6,1*10^-7,1*10^-8,0))
two <- data.frame(run=seq(1,6,1),dist=rep(2,6),betai=c(4*10^-5,rep(4*10^-5,5)),
           betaj=c(4*10^-5,rep(8*10^-5,5)),
           betaij=c(4*10^-5,1*10^-5,1*10^-6,1*10^-7,1*10^-8,0))
three <- data.frame(run=seq(1,6,1),dist=rep(3,6),betai=c(4*10^-5,rep(4*10^-5,5)),
           betaj=c(4*10^-5,rep(10*10^-5,5)),
           betaij=c(4*10^-5,1*10^-5,1*10^-6,1*10^-7,1*10^-8,0))
four <- data.frame(run=seq(1,6,1),dist=rep(4,6),betai=c(4*10^-5,rep(4*10^-5,5)),
           betaj=c(4*10^-5,rep(12*10^-5,5)),
           betaij=c(4*10^-5,1*10^-5,1*10^-6,1*10^-7,1*10^-8,0))

pars <- rbind(one, two, three, four)
pars$R0 <- rep(NA,nrow(pars))

nsim=100
Rt_threshold_list <- list()

for (i in 1:nrow(pars)){

pars$R0[i] <- R_t_func(Xi=Ni-Ei0-Ii0-Ri0, Xj=Nj-Ej0-Ij0-Rj0,
                       betai=pars$betai[i],betaj=pars$betaj[i],betaij=pars$betaij[i],
                       gamma=gamma)
make_data %>%
  simulate(
    params=c(betai=pars$betai[i],betaj=pars$betaj[i],betaij=pars$betaij[i],
             sigma=sigma,gamma=gamma,rho=rho,k=k,
             Ni=Ni,Ei0=Ei0,Ii0=Ii0,Ri0=Ri0,Nj=Nj,Ej0=Ej0,Ij0=Ij0,Rj0=Rj0),
    nsim=nsim,format="data.frame",include.data=FALSE
  ) -> sim_data

#True R_t
sim_data$R_t <- R_t_func(Xi=sim_data$Si, Xj=sim_data$Sj,
                         betai=pars$betai[i],betaj=pars$betaj[i],betaij=pars$betaij[i],gamma=gamma)
sim_data$.id <- as.numeric(sim_data$.id)
day_below_one = sim_data %>% group_by(.id)  %>% filter(R_t<1) %>% filter(row_number()==1)

##EpiEstim R_t
sim_data$reports <- sim_data$reports_i+sim_data$reports_j
T <- length(day)
t_start <- seq(14, T-6) # starting at 2 as conditional on the past observations
t_end <- t_start + 6 # adding 6 to get 7-day windows as bounds included in window

dat_list = split(sim_data, sim_data$.id)


estim_reports_df <- map_dfr(dat_list, ~estimate_R(.x$reports, 
                                                  method="parametric_si",
                                                  config = make_config(list(
                                                    t_start = t_start,
                                                    t_end = t_end,
                                                    mean_si = GI, 
                                                    std_si = sd_GI)))$R,.id=".id")
estim_reports_df$.id <- as.numeric(estim_reports_df$.id)
day_below_estimR_reports  = estim_reports_df  %>% group_by(.id)  %>% 
  filter(`Mean(R)`<1)%>% filter(row_number()==1)

thresh_error_reports <- (day_below_estimR_reports$t_start) - day_below_one$day
sim <- as.factor(rep(pars$run[i],length(thresh_error_reports)))
d_beta <- as.factor(rep(pars$dist[i],length(thresh_error_reports)))
Rt_threshold_list[[i]] <- data.frame(sim,d_beta, thresh_error_reports)
}

Rt_threshold_df <- do.call(rbind, Rt_threshold_list)

thresh.long <- melt(Rt_threshold_df, id=c("d_beta","sim"))

levels(thresh.long$d_beta) <- c("1"=TeX("$\\delta_{\\beta_i} = 2 \\times 10^{-5} $"),
                                "2"=TeX("$\\delta_{\\beta_i} = 4 \\times 10^{-5} $"),
                                "3"=TeX("$\\delta_{\\beta_i} = 6 \\times 10^{-5} $"),
                                "4"=TeX("$\\delta_{\\beta_i} = 8 \\times 10^{-5} $"))

ggplot(thresh.long, x=variable, aes(y=value, fill=sim)) +  
  guides(fill=guide_legend(title=TeX("$\\beta_{ij}$")))+geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-20, 10), breaks = seq(-20,10,5))+
  scale_fill_brewer(palette="Blues",labels=c('Single Pop.', TeX("$1 \\times 10^{-5}$"),
                                             TeX("$1 \\times 10^{-6}$"),TeX("$1 \\times 10^{-7}$"),
                                             TeX("$1 \\times 10^{-8}$"),"Separate Pops."))+
theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=12),
        panel.grid.major.x =element_blank(),
        strip.text = element_text(size = 14),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))+
  ylab(TeX("\\delta$_t$ when \\textit{\\hat{R}}$_t$<1"))+
  xlab(TeX("Population Structure $\\rightarrow$"))+theme(legend.text.align = 0)+
  facet_wrap(vars(d_beta),labeller = label_parsed, nrow = 1)+
  geom_hline(yintercept=0, linetype=3)


