###########Load Packages
library(pomp)
library(tidyverse)
library(cowplot)
library(EpiEstim)
library(reshape2)
library(latex2exp)
library(moments)
library(mousetrap)
library(foreach)
library(iterators)
library(parallel)
library(rngtools)
library(doParallel)
library(doRNG)
registerDoParallel()
registerDoRNG(2488820)
library(ggpubr)
library(bayestestR) 
library(coda)
library(MCMCglmm) 
library(reshape2)
library(patchwork)
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)
library(tidyr)
library(tidyselect)
library(tibble)
library(readr)
library(stringr)
library(forcats)
library(lhs)

set.seed(620)

# dir.create("tmp")
#options(pomp_cdir="./tmp")



#############################Code Model in Pomp
SEIR.compart <- Csnippet("
  double dN_SiEi = rbinom(Si,1-exp((-betai*Ii-betaij*Ij)*dt));
  double dN_SjEj = rbinom(Sj,1-exp((-betaj*Ij-betaij*Ii)*dt));

  double dN_EiIi = rbinom(Ei,1-exp(-sigma*dt));
  double dN_EjIj = rbinom(Ej,1-exp(-sigma*dt));
  
  double dN_IiRi = rbinom(Ii,1-exp(-gamma*dt));
  double dN_IjRj = rbinom(Ij,1-exp(-gamma*dt));
  
  Si -= dN_SiEi;
  Sj -= dN_SjEj;
  
  Ei += dN_SiEi - dN_EiIi;
  Ej += dN_SjEj - dN_EjIj;
  
  Ii += dN_EiIi - dN_IiRi;
  Ij += dN_EjIj - dN_IjRj;
  
  Ri += dN_IiRi;
  Rj += dN_IjRj;
  
  Hi += dN_SiEi;
  Hj += dN_SjEj;
  
  Ci += dN_EiIi;
  Cj += dN_EjIj;
")

SEIR_rinit <- Csnippet("
  Si = Ni-Ei0-Ii0-Ri0;
  Sj = Nj-Ej0-Ij0-Rj0;
  Ei = Ei0;
  Ej = Ej0;
  Ii = Ii0;
  Ij = Ij0;
  Ri = Ri0;
  Rj = Rj0;
  Hi = Ei;
  Hj = Ej;
  Ci = Ii;
  Cj = Ij;
")

rmeas <- Csnippet("
  reports_i = rnbinom_mu(k,rho*Ci);
  reports_j = rnbinom_mu(k,rho*Cj);
  ")

dmeas <- Csnippet("if (betai<0||betaj<0||betaij<0||rho<0||rho>0.5||k<0) {
  lik = (give_log) ? R_NegInf : 0.0;
} else {
  lik=dnbinom_mu(reports_i,k,(rho*Ci),give_log)+dnbinom_mu(reports_j,k,(rho*Cj),give_log);
}")


accumvars=c("Hi","Hj","Ci","Cj")
statenames=c("Si","Sj","Ei","Ej","Ii","Ij","Ri","Rj","Hi","Hj","Ci","Cj")
paramnames=c("betai","betaj","betaij","sigma","gamma","rho","k",
             "Ni","Ei0","Ii0","Ri0","Nj","Ej0","Ij0","Rj0")


day <- seq(0, 80, 1)
Data <- as.data.frame(day)
Data$reports_i <- rep(NA, length(day))
Data$reports_j <- rep(NA, length(day))

make_data <- pomp(data=Data,
                  times="day",t0=0,
                  rprocess=euler(SEIR.compart,delta.t=1/7),
                  rinit=SEIR_rinit,
                  rmeasure=rmeas,
                  accumvars=accumvars,
                  statenames=statenames,
                  paramnames=paramnames)


################################################# Functions

###Calculate the generation interval
sigma=1/3.59
gamma=1/3.56
GI=1/gamma + 1/sigma
sd_GI=sqrt(1/gamma^2 + 1/sigma^2)


####Function to calculate R_t
R_t_func <- function(Xi,Xj,betai,betaj,betaij,gamma){
  (Xi*betai + Xj*betaj + 
     sqrt(Xi^2*betai^2 + 4*Xi*Xj*betaij^2 - 2*Xi*Xj*betai*betaj + Xj^2*betaj^2))/
    (2*gamma)
}

####Function to calculate R_t for subpopulations
R_i_func <- function(Xi,Xj,betai,betaij,gamma){
  (Xi*betai + Xj*betaij)/gamma
}


