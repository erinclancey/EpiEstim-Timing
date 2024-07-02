###PLOT synthetic structure
setwd("~/EpiEstim/Epi Methods/Code Final")
source("Model_synthetic structure.R")
set.seed(620)
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
betai=4*10^-5
betaj=12*10^-5
betaij=1*10^-8
nsim=1

make_data %>%
  simulate(
    params=c(betai=betai,betaj=betaj,betaij=betaij,
             sigma=sigma,gamma=gamma,rho=rho,k=k,
             Ni=Ni,Ei0=Ei0,Ii0=Ii0,Ri0=Ri0,Nj=Nj,Ej0=Ej0,Ij0=Ij0,Rj0=Rj0),
    nsim=nsim,format="data.frame",include.data=FALSE
  ) -> sim_data


##Calculate R_0 and when R_0 is less than 1
sim_data$R_0 <- rep(R_t_func(Xi=Ni-Ei0-Ii0-Ri0, Xj=Nj-Ej0-Ij0-Rj0,
                             betai=betai,betaj=betaj,betaij=betaij,gamma=gamma),nsim)
sim_data$R_0[1]
sim_data$R_t <- R_t_func(Xi=sim_data$Si, Xj=sim_data$Sj,
                         betai=betai,betaj=betaj,betaij=betaij,gamma=gamma)

sim_data$.id <- as.numeric(sim_data$.id)
sim_data %>% group_by(.id)  %>% filter(R_t<1) %>% filter(row_number()==1)


#####################Estimate R_t with EpiEstim
sim_data$ince <- sim_data$Hi+sim_data$Hj
sim_data$reports <- sim_data$reports_i+sim_data$reports_j

T <- length(day)
t_start <- seq(7, T-6) # starting at 2 as conditional on the past observations
t_end <- t_start + 6 # adding 6 to get 7-day windows as bounds included in window

dat_list = split(sim_data, sim_data$.id)

#True Incidence
estim_true_df <- map_dfr(dat_list, ~estimate_R(.x$ince, 
                                               method="parametric_si",
                                               config = make_config(list(
                                                 t_start = t_start,
                                                 t_end = t_end,
                                                 mean_si = GI, 
                                                 std_si = sd_GI)))$R, .id=".id")
estim_true_df$.id <- as.numeric(estim_true_df$.id)
estim_true_df %>% group_by(.id)  %>% 
  filter(`Mean(R)`<1)%>% filter(row_number()==1)

##Reports
estim_reports_df <- map_dfr(dat_list, ~estimate_R(.x$reports, 
                                                  method="parametric_si",
                                                  config = make_config(list(
                                                    t_start = t_start,
                                                    t_end = t_end,
                                                    mean_si = GI, 
                                                    std_si = sd_GI)))$R,.id=".id")
estim_reports_df$.id <- as.numeric(estim_reports_df$.id)
estim_reports_df  %>% group_by(.id)  %>% 
  filter(`Mean(R)`<1)%>% filter(row_number()==1)
####Agg by week
sim_data$Day <- as.Date(sim_data$day)
sim_data$week <- as.numeric(strftime(sim_data$Day , format = "%V"))
sim_week <- sim_data%>% group_by(week) %>% summarize(i=sum(reports)) %>% as.data.frame()
agg_week <- map_dfr(dat_list,~estimate_R(incid = sim_week$i,
                                         dt = 7L,
                                         dt_out = 7L,
                                         recon_opt = "naive",
                                         iter = 10L,
                                         tol = 1e-6,
                                         grid = list(precision = 0.001, min = -1, max = 1),
                                         method="parametric_si",
                                         config = make_config(list(
                                           mean_si = GI,
                                           std_si = sd_GI)))$R,.id=".id")

agg_week$.id <- as.numeric(agg_week$.id)
agg_week  %>% group_by(.id)  %>% 
  filter(`Mean(R)`<1)%>% filter(row_number()==1)
ID=nsim
subset_sim <- subset(sim_data, .id==ID)
# Ploy the simulated the process model and the case reports
A <- ggplot(subset_sim, aes(x=day)) + theme_minimal()+
  ylab("N")+xlab("t (days)")+
  geom_line(aes(y = Si+Sj, color = "Si")) + 
  geom_line(aes(y = Ei+Ej, color = "Ei")) + 
  geom_line(aes(y = Ii+Ij, color="Ii")) + 
  geom_line(aes(y = Ri+Rj, color="Ri"))+
  scale_colour_manual("", breaks = c("Si","Ei","Ii","Ri"),
                      values = c("darkred","black","steelblue","darkolivegreen4"),
                      labels=c("S","E","I","R"))+
  ggtitle("Synthetic Epidemic")+
  scale_x_continuous(limits=c(0,T),n.breaks=10)+
  scale_y_continuous(limits=c(0,35000), n.breaks=5)+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))

B <- ggplot(subset_sim, aes(x=day,y=Hi+Hj)) + theme_minimal()+ 
  geom_bar(stat="identity", color="steelblue", fill="steelblue", alpha=0.25)+
  ylab("Incidence")+xlab("t (days)")+
  ggtitle("True Incidence")+
  scale_x_continuous(limits=c(0,T),n.breaks=10)+
  scale_y_continuous(n.breaks=5)+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

C <- ggplot(subset_sim, aes(x=day,y=reports_i+reports_j)) + theme_minimal()+ 
  geom_bar(stat="identity", color="darkred", fill="darkred", alpha=0.25)+
  ylab("Reports")+xlab("t (days)")+
  ggtitle("Daily Case Reports")+
  scale_x_continuous(limits=c(0,T),n.breaks=10)+
  scale_y_continuous(n.breaks=5)+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

D <- ggplot(subset_sim, aes(x=day, y=R_t)) + theme_minimal()+
  xlab("t (days)")+geom_line(color="black")+
  ylab(TeX("\\textit{R}$_t$"))+
  ggtitle("Reproductive Number")+
  geom_hline(yintercept=1, linetype="dotted")+
  scale_x_continuous(limits=c(0,T),n.breaks=10)+
  scale_y_continuous(limits=c(0,7), n.breaks=5)+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

subset_true <- subset(estim_true_df, .id==ID)
E <- ggplot(subset_true , aes(x=t_start, y=`Mean(R)`)) + theme_minimal()+
  xlab("t (days)")+geom_line(color="steelblue")+
  ylab(TeX("\\textit{\\hat{R}}$_t$"))+
  ggtitle("Estimate from True Incidence")+
  scale_x_continuous(limits=c(0,T), n.breaks=10)+
  scale_y_continuous(limits=c(0,7), n.breaks=5)+
  geom_ribbon(aes(ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

subset_reports <- subset(estim_reports_df, .id==ID)
G <- ggplot(subset_reports, aes(x=t_start, y=`Mean(R)`)) + theme_minimal()+
  xlab("t (days)")+geom_line(color="darkred")+
  ylab(TeX("\\textit{\\hat{R}}$_t$"))+
  ggtitle("Estimate from Daily Case Reports")+
  scale_x_continuous(limits=c(0,T), n.breaks=10)+
  scale_y_continuous(limits=c(0,7), n.breaks=5)+
  geom_ribbon(aes(ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

H <- ggplot(sim_week, aes(x=week,y=i)) + theme_minimal()+
  geom_bar(stat="identity", color="darkolivegreen4", fill="darkolivegreen4", alpha=0.25)+
  ylab("Reports")+xlab("t (weeks)")+
  ggtitle("Weekly Case Reports")+
  scale_x_continuous(limits=c(0,12), n.breaks=10)+
  scale_y_continuous(n.breaks=5)+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

subset_agg <- subset(agg_week, .id==ID)
I <- ggplot(subset_agg, aes(x=t_start, y=`Mean(R)`)) + theme_minimal()+
  xlab("t (days)")+geom_line(color="darkolivegreen4")+
  ylab(TeX("\\textit{\\hat{R}}$_t$"))+
  ggtitle("Estimate from Weekly Case Reports")+
  scale_x_continuous(limits=c(0,T), n.breaks=10)+
  scale_y_continuous(limits=c(0,7), n.breaks=5)+
  geom_ribbon(aes(ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")+
  theme(plot.title = element_text(size=14),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

plot <- plot_grid(A,D,B,E,C,G,H,I, ncol = 2, nrow = 4, 
                  rel_heights=c(1,1,1,1,1,1,1), labels = c('A','', 'B','','C','','D'))
plot(plot)

