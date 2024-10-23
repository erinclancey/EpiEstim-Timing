setwd("~/EpiEstim/Epi Methods/Code Final")
source("Model_synthetic structure.R")
setwd("~/EpiEstim/Epi Methods/Submitted Docs_revision")
source("Model_synthetic structure.R")

###Whitman County Daily Data
Whit_Co_daily <- read.csv("Whit_CO_daily.csv", header=TRUE)
Whit_Co_daily <- subset(Whit_Co_daily, J_date>229)
A <- ggplot(Whit_Co_daily, aes(x=J_date,y=Cases)) + theme_minimal()+ 
  geom_bar(stat="identity", color="darkred", fill="darkred", alpha=0.25)+
  ylab("Reports")+xlab("day of year")+
  ggtitle("Whitman Co. Daily Case Reports")+
  scale_x_continuous(limits=c(230,362),n.breaks=10)+
  scale_y_continuous(n.breaks=10)

T <- nrow(Whit_Co_daily)
t_start <- seq(7, T-6) # starting at 2 as conditional on the past observations
t_end <- t_start + 6 # adding 6 to get 7-day windows as bounds included in window
Whit_R_t <- estimate_R(Whit_Co_daily$Cases, 
                       method="parametric_si",
                       config = make_config(list(
                         t_start = t_start,
                         t_end = t_end,
                         mean_si = GI, 
                         std_si = sd_GI)))
Whit_R_t_df <- Whit_R_t$R
Whit_R_t_df$J_date <- Whit_Co_daily$J_date[-c(122:133)]
filter(Whit_R_t_df, `Mean(R)` <1)[1,]


B <- ggplot(Whit_R_t_df , aes(x=J_date, y=`Mean(R)`)) + theme_minimal()+
  ylab("R")+xlab("day of year")+geom_line(color="darkred")+
  ylab(TeX("\\textit{\\hat{R}}$_t$"))+
  ggtitle("Estimate from Daily Case Reports")+
  scale_x_continuous(limits=c(230,362), n.breaks=10)+
  scale_y_continuous(limits=c(0,4), n.breaks=5)+
  geom_ribbon(aes(ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")

#Whitman County Weekly Aggregated Data
Whit_Co_weekly <- read.csv("Whit_CO_week.csv", header=TRUE)
Whit_Co_weekly <- subset(Whit_Co_weekly, Week>33)
C <- ggplot(Whit_Co_weekly, aes(x=Week,y=Cases)) + theme_minimal()+ 
  geom_bar(stat="identity", color="darkolivegreen4", fill="darkolivegreen4", alpha=0.25)+
  ylab("Reports")+xlab("week of year")+
  ggtitle("Whitman Co. Weekly Case Reports")+
  scale_x_continuous(n.break=12)+
  scale_y_continuous(n.breaks=10)

agg_R_Whit <- estimate_R(incid = Whit_Co_weekly$Cases,
                         dt = 7L,
                         dt_out = 7L,
                         recon_opt = "naive",
                         iter = 10L,
                         tol = 1e-6,
                         grid = list(precision = 0.001, min = -1, max = 1),
                         method="parametric_si",
                         config = make_config(list(
                           mean_si = GI, 
                           std_si = sd_GI)))
agg_R_Whit_df <- agg_R_Whit$R
agg_R_Whit_df$J_date <- Whit_Co_daily$J_date[-c(121:133)]
filter(agg_R_Whit_df, `Mean(R)` <1)[1,]

D <- ggplot(agg_R_Whit_df , aes(x=J_date, y=`Mean(R)`)) + theme_minimal()+
  ylab("R")+xlab("day of year")+geom_line(color="darkolivegreen4")+
  ylab(TeX("\\textit{\\hat{R}}$_t$"))+
  ggtitle("Estimate from Weekly Case Reports")+
  scale_x_continuous(limits=c(230,362), n.breaks=10)+
  scale_y_continuous(limits=c(0,4), n.breaks=5)+
  geom_ribbon(aes(ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")

##############Whitman County Subpopulation Plot

E <- ggplot(Whit_Co_weekly, aes(x=Week,y=PUL_cases)) + theme_minimal()+ 
  geom_bar(stat="identity", color="steelblue",fill="steelblue", alpha=0.5)+
  geom_bar(aes(x=Week,y=WSU_cases),stat="identity",fill="darkred", alpha=0.35)+
  ylab("Reports")+xlab("week of year")+
  ggtitle("Whitman Co. Subpopulation Weekly Case Reports")+
  scale_x_continuous(n.breaks=12)+
  scale_y_continuous(limits=c(0,400),n.breaks=10)


WSU_R_t <- estimate_R(incid = Whit_Co_weekly$WSU_cases,
                         dt = 7L,
                         dt_out = 7L,
                         recon_opt = "naive",
                         iter = 10L,
                         tol = 1e-6,
                         grid = list(precision = 0.001, min = -1, max = 1),
                         method="parametric_si",
                         config = make_config(list(
                           mean_si = GI, 
                           std_si = sd_GI)))
WSU_R_t_df <- WSU_R_t$R
WSU_R_t_df$J_date <- Whit_Co_daily$J_date[-c(121:133)]
filter(WSU_R_t_df, `Mean(R)` <1)[1,]

PUL_R_t <- estimate_R(incid = Whit_Co_weekly$PUL_cases,
                      dt = 7L,
                      dt_out = 7L,
                      recon_opt = "naive",
                      iter = 10L,
                      tol = 1e-6,
                      grid = list(precision = 0.001, min = -1, max = 1),
                      method="parametric_si",
                      config = make_config(list(
                        mean_si = GI, 
                        std_si = sd_GI)))
#PUL_R_t_df <- PUL_R_t$R[-c(1:45),]
#PUL_R_t_df$J_date <- Whit_Co_daily$J_date[-c(1:58)]
PUL_R_t_df <- PUL_R_t$R
PUL_R_t_df$J_date <- Whit_Co_daily$J_date[-c(121:133)]
filter(PUL_R_t_df,J_date>280 & `Mean(R)` <1)[1,]

G <- ggplot()+
  geom_line(data = WSU_R_t_df, aes(x=J_date, y=`Mean(R)`), color = "darkred") +
  geom_line(data = PUL_R_t_df, aes(x=J_date, y=`Mean(R)`), color = "steelblue") +
  ylab("R")+xlab("day of year")+ theme_minimal()+
  ylab(TeX("\\textit{\\hat{R}}$_t$"))+
  ggtitle("Whitman Co. Subpopuation Estimates")+
  scale_x_continuous(limits=c(230,362), n.breaks=10)+
  scale_y_continuous(limits=c(0,8), n.breaks=5)+
  geom_ribbon(data = WSU_R_t_df,aes(x=J_date, ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_ribbon(data = PUL_R_t_df,aes(x=J_date, ymin=`Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`),alpha=0.1)+
  geom_hline(yintercept=1, linetype="dotted")

plot3 <- plot_grid(A,B,C,D,E,G, ncol = 2, nrow =3 , rel_heights=c(1,1,1,1,1,1), labels = c('A','', 'B','','C'))
plot(plot3)
