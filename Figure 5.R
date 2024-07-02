# Show Plot for a single population
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
betai=4*10^-5
betaj=4*10^-5
betaij=4*10^-5
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
sim_data$R_t <- R_t_func(Xi=sim_data$Si, Xj=sim_data$Sj,
                         betai=betai,betaj=betaj,betaij=betaij,gamma=gamma)

sim_data$.id <- as.numeric(sim_data$.id)
sim_data %>% group_by(.id)  %>% filter(R_t<1) %>% filter(row_number()==1)



###########################################
##Calculate R_0 and when R_0 is less than 1
sim_data$R_0_i <- rep(R_i_func(Xi=Ni-Ei0-Ii0-Ri0, Xj=Nj-Ej0-Ij0-Rj0,
                               betai=betai,betaij=betaij,gamma=gamma),nsim)

sim_data$R_0_j <- rep(R_i_func( Xi=Nj-Ej0-Ij0-Rj0,Xj=Ni-Ei0-Ii0-Ri0,
                                betai=betaj,betaij=betaij,gamma=gamma),nsim)

sim_data$R_i <- R_i_func(Xi=sim_data$Si, Xj=sim_data$Sj,
                         betai=betai,betaij=betaij,gamma=gamma)

sim_data$R_j <- R_i_func(Xi=sim_data$Sj, Xj=sim_data$Si,
                         betai=betaj,betaij=betaij,gamma=gamma)

#####################Estimate R_t with EpiEstim

T <- length(day)
one <- ggplot(sim_data, aes(x=day)) + theme_minimal()+
  ylab("")+xlab("")+
  geom_line(aes(y = R_t, color = "R_t", linetype="R_t"), size=1, alpha=0.5) + 
  geom_line(aes(y = R_i, color = "R_i",linetype="R_i")) + 
  geom_line(aes(y = R_j, color="R_j",linetype="R_j")) + 
  scale_linetype_manual(name="", breaks = c("R_t","R_i","R_j"),
                        values = c("solid","longdash","longdash"))+
  guides(linetype="none")+guides(color="none")+
  scale_colour_manual(name="", breaks = c("R_t","R_i","R_j"),
                      values = c("black","steelblue","darkred"),
                      labels=c(TeX("\\textit{R}$_t$"),
                               TeX("\\textit{R}$_{t_i}$"),
                               TeX("\\textit{R}$_{t_j}$")))+
  ggtitle("Single population")+
  geom_hline(yintercept=1, linetype="dotted")+
  scale_x_continuous(limits=c(0,T),n.breaks=10)+
  scale_y_continuous(limits=c(0,7), n.breaks=5)+
  theme(plot.title = element_text(size=12),
        axis.title.x =element_text(size=12),
        axis.title.y =element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        legend.position="bottom")
plot(one)
