####################
### ~   Setup   ~ ##
####################

library(maptools)
library(spdep)
library(dplyr)
library(MASS)
library(tidyr)

source("sim_functions.R")

logit <- function(p) {
  log(p/(1-p))
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}

## Read in geography of Kenya
adm1 = readRDS("Data/KEN_adm1.rds")

nb.RK = poly2nb(adm1, queen=FALSE, row.names=adm1@data$NAME_1)
graph <- nb2mat(nb.RK, style="B")
graph <- graph * -1
diag(graph) <- card(nb.RK)
# Scale the ICAR model
gvar <- exp(mean(log(diag(ginv(graph)))))
Ks <- gvar*graph # equivalent to scale.model=T


########################
### ~ Simulate Data ~ ##
########################

## Spatial random effects
tau <- 45
eg <- eigen(tau * Ks)
set.seed(12305044)
y <- rnorm(46, 0, sqrt(1/eg$values[-47]))
s <- eg$vectors[,-47] %*% y

## Unstructured random effects
set.seed(545)
eps <- rnorm(47, 0, sd=sqrt(1/90))

## Temporal random effects
x1 <- c(0.7, 0.4, 0.2, 0, -0.1, -0.3, -0.6)
x1 <- x1 - mean(x1)
x2 <- c(0.6, 0.5, 0.15, 0, -0.15, -0.4, -0.6)
x2 <- x2 - mean(x2)
x3 <- c(0.65, 0.5, 0.1, 0, -0.2, -0.4, -0.6)
x3 <- x3 - mean(x3)

## Fixed effects
ba <- c(-1.9, -2.94, -5.29) # intercepts

## Vector of mother's ages (1:5 FBH:SBH)
mother.ageFBH <- c(rep(15:19, each=4*20*2),
                   rep(20:24, each=4*20*2),
                   rep(25:29, each=4*20*2),
                   rep(30:34, each=4*15*2),
                   rep(35:39, each=4*10*2),
                   rep(40:44, each=4*8*2),
                   rep(45:49, each=4*7*2))
mother.ageSBH <- c(rep(15:19, each=20*20*2),
                   rep(20:24, each=20*20*2),
                   rep(25:29, each=20*20*2),
                   rep(30:34, each=20*15*2),
                   rep(35:39, each=20*10*2),
                   rep(40:44, each=20*8*2),
                   rep(45:49, each=20*7*2))
true.womenSBH <- table(c(mother.ageSBH))

hazard.dataFBH.summarized <- data.frame(period=NA,
                                        agegp=NA,
                                        N=NA,
                                        Y=NA,
                                        region=NA)
hazard.dataSBH.summarized <- data.frame(period=NA,
                                        agegp=NA,
                                        N=NA,
                                        Y=NA,
                                        region=NA)


## Fertility probabilities f(m)
prob.birth <- c(rep(0,14), rep(152, 5), rep(269, 5), rep(238, 5), rep(206, 5),
                rep((33+82+162)/3,15))/1000

## Years, periods
y <- 1975:2009
p <- rep(1:7, each=5)

## Initialize some lists and data frames to store date
sbh.small <- list()
cs <- list()
c_CEB <- list()
CEBa <- list()
fertility.dataFBHsum <- data.frame(m.s=NA, m.b=NA, CEB=NA, FP=NA, region=NA)

set.seed(435810)
start_time <- Sys.time()
for(r in 1:47) {
  print(r)
  # determine prob.death
  # columns correspond to time periods, rows correspond to age of mother (15:49)
  prob.death.logit <- matrix(0, nrow=35, ncol=7)
  for(i in 1:nrow(prob.death.logit)) {
    if(i ==1) {
      prob.death.logit[i, ] <- ba[1] + x1 + s[r] + eps[r]
    } else if(i < 6) {
      prob.death.logit[i, ] <- ba[2] + x2 + s[r] + eps[r]
    } else {
      prob.death.logit[i, ] <- ba[3] + x3 + s[r] + eps[r]
    }
  }
  prob.death <- expit(prob.death.logit)
  
  # FBH data
  fbh <- create_sim(m.s = mother.ageFBH, 
                    prob.birth, prob.death, y.s=2010, y, p)
  fbh$count <- 1
  
  # get children ever born for women of a specific age at survey time and birth time
  suppressMessages(fertility.dataFBH <- fbh %>% group_by(m.s, m.b) %>% summarise(CEB=n()))
  
  # get number of women of each age
  fpFBH <- data.frame(m.s = 15:49, FP=as.vector(table(mother.ageFBH)))
  fertility.dataFBH <- merge(fertility.dataFBH, fpFBH)
  
  # add region, then rbind to the total FBH dataframe across all regions
  fertility.dataFBH$region <- r
  fertility.dataFBHsum <- rbind(fertility.dataFBHsum, fertility.dataFBH)
  
  # create hazard data
  hazard.dataFBH <- overallHaz.shorter(fbh)
  
  # add in time period to hazard data
  hazard.dataFBH$period <- cut(hazard.dataFBH$years, seq(1975, 2010, 5),
                               include.lowest = T, right=F)
  
  # create binomial counts dataframe by time period and age group (from hazard data)
  suppressMessages(tmp <- hazard.dataFBH %>% group_by(period, agegp) %>% summarise(N=sum(N),Y=sum(Y)))
  tmp$agegp <- as.numeric(tmp$agegp)
  tmp$period <- as.numeric(tmp$period)
  tmp$region <- r
  tmp <- data.frame(tmp)
  
  # rbind to the total summarized FBH dataframe across all regions
  hazard.dataFBH.summarized <- rbind(hazard.dataFBH.summarized, tmp)
  
  # SBH data
  sbh <- create_sim(m.s = mother.ageSBH, 
                    prob.birth, prob.death, y.s=2010, y, p)
  
  sbh$count <- 1
  
  # get children ever born for women of a specific age at survey time and birth time
  suppressMessages(CEBa.tmp <- sbh %>% group_by(m.s, m.b) %>% summarise(CEB=n()))
  
  # create hazard data
  hazard.dataSBH <- overallHaz.shorter(sbh)
  
  # add in time period to hazard data
  hazard.dataSBH$period <- cut(hazard.dataSBH$years, seq(1975, 2010, 5),
                               include.lowest = T, right=F)
  
  # create binomial counts dataframe by time period and age group (from hazard data)
  suppressMessages(tmp <- hazard.dataSBH %>% group_by(period, agegp) %>% summarise(N=sum(N),Y=sum(Y)))
  tmp$agegp <- as.numeric(tmp$agegp)
  tmp$period <- as.numeric(tmp$period)
  tmp$region <- r
  tmp <- data.frame(tmp)
  
  # rbind to the total summarized FBH dataframe across all regions
  hazard.dataSBH.summarized <- rbind(hazard.dataSBH.summarized, tmp)
  
  # FBH --> SBH for SBH women
  sbh.true <- data.frame(id=unique(sbh$id))
  
  # get number of births and number of deaths for each mother in sbh dataframe
  sbh.true$no.births <- tapply(rep(1,nrow(sbh)), sbh$id, sum)
  sbh.true$no.deaths <- tapply(!is.na(sbh$death), sbh$id, sum)
  
  # add in mother's age at survey and year of survey
  sbh.true$m.s <- tapply(sbh$m.s, sbh$id, unique)
  sbh.true$y.s <- 2010
  
  # group by mother's age at survey, number of births, and number of deaths. total number of women
  # in the given categoriy is given in the 'count' column
  suppressMessages(sbh.true <- as.data.frame(sbh.true %>% 
                                               group_by(m.s, y.s, no.births, no.deaths) %>% 
                                               summarise(count =n())))
  
  # get children ever born and children died for mother's age at survey, for this region,
  # and add to sbh.small
  suppressMessages(sbh.small[[r]] <- sbh.true %>% group_by(m.s) %>% summarise(CEB = sum(no.births*count),
                                                                              CD = sum(no.deaths*count)))
  
  # add in the true number of women of that given age at survey
  sbh.small[[r]]$nwomen <- true.womenSBH[-1]
  
  # relative true birth probabilities c
  cs[[r]] <- lapply(16:49, function(ms) {
    prob.birth[15:(ms-1)]/sum(prob.birth[1:(ms-1)])
  })
  
  # multiplying c by CEB in the SBH data 
  #  (i.e. second type of birth history information known)
  c_CEB[[r]] <- matrix(0, nrow=nrow(sbh.small[[r]]), ncol=nrow(sbh.small[[r]]))
  for(ms in 16:49) {
    for(a in 1:(ms+1-16)) {
      c_CEB[[r]][ms-15, a] <- sbh.small[[r]]$CEB[ms-15] * cs[[r]][[ms-15]][ms - 16 + 2 - a]
    }
  }
  
  # true CEB 
  #  (i.e. first type of birth information known)
  CEBa[[r]] <- spread(CEBa.tmp, m.b, CEB, fill=0)
  CEBa[[r]] <- as.matrix(CEBa[[r]][,-1])
}
end_time <- Sys.time()

# Remove placeholder rows
hazard.dataFBH.summarized <- hazard.dataFBH.summarized[-1, ]
hazard.dataSBH.summarized <- hazard.dataSBH.summarized[-1, ]
fertility.dataFBHsum <- fertility.dataFBHsum[-1, ]
sbh.small.all <- do.call("rbind", sbh.small)


## Create precision matrix for RW2 over time
nperiods <- 7
K <- matrix(0, nrow=nperiods, ncol=nperiods)
for(per in 1:nperiods) {
  if(per == 1) {
    K[1, ] <- c(1, -2, 1, rep(0, nperiods-3))
  } else if(per == 2) {
    K[2, ] <- c(-2, 5, -4, 1, rep(0, nperiods-4))
  } else if(per == nperiods - 1) {
    K[nperiods-1, ] <- c(rep(0, nperiods-4), 1, -4, 5, -2)
  } else if(per == nperiods) {
    K[nperiods, ] <- c(rep(0, nperiods-3), 1, -2, 1)
  } else {
    K[per, ] <- c(rep(0, per-3), 1, -4, 6, -4, 1,
                  rep(0, nperiods - 2 - per))
  }
}
gvar <- exp(mean(log(diag(ginv(K)))))
K <- gvar*K # scale.model=T

## Randomizing starting values
set.seed(343)
eps.init <- eps + rnorm(47, 0, (1/100))
Params = list( "z_int"=ba[1]+0.1, 
               "o_int" = ba[2]-0.1,
               "f_int" = ba[3] + 0.2,
               "lkappat"=log(12),
               "lkappas"=log(4),
               "lkappaeps"=log(100),
               "phiz"=x1 + 0.05,
               "phio"=x2 + 0.05,
               "phif"=x3 - 0.05,
               "Us"=as.vector(round(s,3)),
               "eps"=eps.init)

## Data for TMB
Data = list("cd"=as.integer(sbh.small.all$CD),
            "p"=c(rep(1:7,each=5))[-1], 
            "y"=1:34, 
            "cd_fbh"=hazard.dataFBH.summarized$Y,  
            "ceb_fbh"=hazard.dataFBH.summarized$N,
            "Kt"=K,
            "Ks"=Ks,
            "At"=matrix(1, nrow=1,ncol=7),
            "As"=matrix(1, nrow=1, ncol=47),
            "ut"=1,
            "us"=1,
            "diagnum"=INLA:::inla.set.f.default()$diagonal)

## Save data, params and other helpful quantities
save(Data, Params, graph, x1, x2, x3, ba, tau, ba, s, c_CEB, cs, eps,
     hazard.dataFBH.summarized, hazard.dataSBH.summarized, CEBa,
     fertility.dataFBHsum,
     file="Data/sim_data-FINAL.rda")


load("Data/sim_data-FINAL.rda")


## Use TMB to fit the model
U_spat <- "
#include <TMB.hpp>
#include <vector>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
return R_IsNA(asDouble(x));
}

template<class Type>
Type expit(Type x) {
  return 1./(1.+exp(-x));
}

// Negative log posterior
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data
  DATA_VECTOR(cd); // cd in SBH
  DATA_IVECTOR(p);  // mapping of periods
  DATA_IVECTOR(y);  // mapping of years
  DATA_VECTOR(cd_fbh); // cd in FBH 
  DATA_VECTOR(ceb_fbh); // ceb in FBH
  DATA_MATRIX(CEBa); // CEB(a) matrix (sometimes c*CEB)
  DATA_MATRIX(Kt); // RW2
  DATA_MATRIX(Ks); // ICAR structure matrix
  DATA_MATRIX(At); // [1, ...., 1] n time pts
  DATA_MATRIX(As); // [1, ...., 1] n regions
  DATA_SCALAR(ut); // for use in pc prior (for temporal)
  DATA_SCALAR(us); // for use in pc prior (for spatial)
  DATA_SCALAR(diagnum); // INLA number to add to diagonal of structure matrix
  DATA_INTEGER(sbhI); // SBH indicator (1 if use SBH, 0 otherwise)
  DATA_INTEGER(fbhI); // SBH indicator (1 if use SBH, 0 otherwise)
  
  // Parameters
  PARAMETER(z_int);
  PARAMETER(o_int);
  PARAMETER(f_int);
  PARAMETER(lkappat); // log precision temporal
  PARAMETER(lkappas); // log precision spatial
  PARAMETER(lkappaeps); // log precision unstructure spatial error
  
  
  // Random effects
  PARAMETER_VECTOR(phiz); // unconstrained
  PARAMETER_VECTOR(phio); // unconstrained
  PARAMETER_VECTOR(phif); // unconstrained
  PARAMETER_VECTOR(Us); // spatial unconstrained
  PARAMETER_VECTOR(eps); // unstructured
  
  // Transform parameters
  Type kappaeps = exp(lkappaeps);
  Type kappat = exp(lkappat);
  Type kappas = exp(lkappas);
  matrix<Type> Qtmpt = kappat * Kt;
  for(int i=0;i<7;i++) {
    Qtmpt(i,i) += diagnum;
  }
  Eigen::SparseMatrix<Type> Qt = asSparseMatrix(Qtmpt);
  
  matrix<Type> Qtmps = kappas * Ks;
  for(int i=0;i<47;i++) {
    Qtmps(i,i) += diagnum;
  }
  Eigen::SparseMatrix<Type> Qs = asSparseMatrix(Qtmps);

  // x --> x* (using 2.30 in Rue & Held)
  // temporal
  matrix<Type> Sigt = invertSparseMatrix(Qt);
  vector<Type> A_times_invSigt = At.transpose() * (1/Sigt.sum());
  Type phiz_sum = phiz.sum();
  vector<Type> phizstar = phiz - Sigt * (A_times_invSigt) * phiz_sum;
  Type phio_sum = phio.sum();
  vector<Type> phiostar = phio - Sigt * (A_times_invSigt) * phio_sum;
  Type phif_sum = phif.sum();
  vector<Type> phifstar = phif - Sigt * (A_times_invSigt) * phif_sum;
  REPORT(phizstar);
  REPORT(phiostar);
  REPORT(phifstar);

  Type Us_sum = Us.sum();
  matrix<Type> Sigs = invertSparseMatrix(Qs);
  vector<Type> A_times_invSigs = As.transpose() * (1/Sigs.sum());
  vector<Type> Usstar = Us - Sigs * (A_times_invSigs) * Us_sum;
  REPORT(Usstar);
  

  // Objective funcction
  
  int M = 34; //c_CEB.rows();
  int T = Kt.rows();
  int R = Ks.rows();
  
  vector<Type> nll_comp(3);
  nll_comp.setZero();
  
  // Derived quantities
  vector<Type> oneqz(T);
  oneqz=z_int + phizstar;
  vector<Type> oneqo(T);
  oneqo=o_int + phiostar;
  vector<Type> oneqf(T);
  oneqf=f_int + phifstar;
  REPORT(oneqz);
  REPORT(oneqo);
  REPORT(oneqf);
  
  // Prior for lkappat and lkappas: kappa ~ pcprior(u, alpha = 0.01)
  Type lambdat = -log(0.01)/ut;
  nll_comp(2) += lambdat * exp(-lkappat/Type(2.0));
  nll_comp(2) += lkappat/Type(2.0);
  Type lambdas = -log(0.01)/us;
  nll_comp(2) += lambdas * exp(-lkappas/Type(2.0));
  nll_comp(2) += lkappas/Type(2.0);
  // Prior for lkappaeps: kappa ~ Gamma(shape=1, inv-scale=0.005)
  nll_comp(2) += kappaeps*Type(0.005);
  nll_comp(2) -= lkappaeps;
  
  // Prior on intercepts
  nll_comp(2) += Type(0.5) *  0.01 * (z_int) * (z_int); // N(0,10^2)
  nll_comp(2) += Type(0.5) *  0.01 * (o_int) * (o_int); // N(0,10^2)
  nll_comp(2) += Type(0.5) *  0.01 * (f_int) * (f_int); // N(0,10^2)
  

  // RW prior
  nll_comp(1) += GMRF(Qt) (phiz);
  nll_comp(1) += GMRF(Qt) (phio);
  nll_comp(1) += GMRF(Qt) (phif);
  
  // ICAR prior
  nll_comp(1) += GMRF(Qs) (Us);

  // iid on region prior
  for (int r=1; r<=R; r++) {
    nll_comp(1) += Type(0.5) * eps(r-1) * eps(r-1) * kappaeps;
  }
  nll_comp(1) -= Type(47.0)*Type(0.5) * lkappaeps;

  // Data likelihood
  
  // SBH
  if(sbhI==1) {
    int idxsbh=1;
    for(int r=1; r<=R; r++) {
      Type regeff = Usstar(r-1) + eps(r-1);
      vector<Type> qstara(M);
      qstara.setZero();
      qstara(0) = expit(oneqz(p(y(M-1)-1)-1) + regeff);
      for(int a=2;a<=M;a++) {
        qstara(a-1) = 1;
        for(int l=1;l<=a;l++) {
          if(l==1) {
            qstara(a-1) *= (1-expit(regeff + oneqz(p(y(M-(a-l)-1)-1)-1)));
          } else if(l < 6) {
            qstara(a-1) *= (1-expit(regeff + oneqo(p(y(M-(a-l)-1)-1)-1)));
          } else {
            qstara(a-1) *= (1-expit(regeff + oneqf(p(y(M-(a-l)-1)-1)-1)));
          }
        }
        qstara(a-1) = 1 - qstara(a-1);
      }
      for (int i=1;i<=M;i++) {
        vector<Type> CEBarow = CEBa.row(M*(r-1) + i - 1);
        vector<Type> CEBa_star = CEBarow  * qstara;
        nll_comp(0) -= dpois(cd(idxsbh-1), CEBa_star.sum() , true );
        idxsbh += 1;
      }
    }
  }
  
  
  // FBH
  if(fbhI==1) {
    int idx=1;
    for (int r=1;r<=R;r++) {
      Type regeff = Usstar(r-1) + eps(r-1);
      for (int t=1;t<=T;t++) {
        int ub;
        if(t==1) {
          ub=2;
        } else {
          ub=3;
        }
        for (int agegp=1;agegp<=ub;agegp++) {// since it only goes to 2 for 1st period
          if (agegp==1) {
            nll_comp(0) -= dbinom(cd_fbh(idx-1), ceb_fbh(idx-1), expit(oneqz(t-1) + regeff), true);
          } else if (agegp==2) {
            nll_comp(0) -= dbinom(cd_fbh(idx-1), ceb_fbh(idx-1), expit(oneqo(t-1) + regeff), true);
          } else {
            nll_comp(0) -= dbinom(cd_fbh(idx-1), ceb_fbh(idx-1), expit(oneqf(t-1) + regeff), true);
          }
          idx += 1;
        }
      } 
    }
  }
   
  
  ADREPORT(kappat);
  ADREPORT(kappas);
  ADREPORT(kappaeps);
  ADREPORT(phizstar);
  ADREPORT(phiostar);
  ADREPORT(phifstar);
  ADREPORT(Usstar);
  
  Type nll = nll_comp.sum();
  return nll;
}
"

library(TMB)
write(U_spat,file="TMB/U_spat.cpp")
compile( "TMB/U_spat.cpp")
dyn.load( dynlib("TMB/U_spat") )

# Indicators for including SBH & FBH
Data$sbhI <- 1
Data$fbhI <- 1

## 1. With true CEB(a)
CEBa.corrected <- lapply(CEBa, function(x) {
  output <- x
  for(i in 1:nrow(x)) {
    output[i,1:i ] <- rev(x[i, 1:i])
  }
  output
})

Data$CEBa <- do.call("rbind", CEBa.corrected)

ObjSBH.true = MakeADFun( data=Data, parameters=Params, 
                    random=c("phiz", "phio", "phif", "Us", "eps"), DLL="U_spat" )
t1 <- proc.time()
OptSBH.true <- nlminb(start = ObjSBH.true$par, objective = ObjSBH.true$fn,
                      gradient = ObjSBH.true$gr)
proc.time() - t1

## 2. With true fertility probabilities CEB * c(a)
Data$CEBa <- do.call("rbind", c_CEB)
ObjSBH.truec = MakeADFun( data=Data, parameters=Params,
                         random=c("phiz", "phio", "phif", "Us", "eps"), DLL="U_spat" )
t1 <- proc.time()
OptSBH.truec <- nlminb(start = ObjSBH.truec$par, objective = ObjSBH.truec$fn,
                      gradient = ObjSBH.truec$gr)
proc.time() - t1


## 3. With estimated fertility probabilities CEB * hat{c}(a)
# Step 1:  estimate c(a)
fertility.dataFBHsum$mgroup <- cut(fertility.dataFBHsum$m.b, c(seq(15,35,5), 50), include.lowest = T, right=F)
fertility.dataFBHsumALL <- fertility.dataFBHsum %>% group_by(mgroup) %>% 
  summarise(CEBtotal=sum(CEB), FPtotal=sum(FP)) # we "know" that it is constant over age groups and same in regions
fertility.fit <- glm(cbind(fertility.dataFBHsumALL$CEBtotal, 
                           fertility.dataFBHsumALL$FPtotal- fertility.dataFBHsumALL$CEBtotal) ~ -1 + fertility.dataFBHsumALL$mgroup, 
                     family=binomial)
f.est <- expit(coef(fertility.fit))
f.est <- c(rep(0,14), rep(f.est[1:4], each=5), rep(f.est[5], 15))

# Translate fhat to chat and multiply by CEB in SBH data
cs.est <- list()
c_CEB.est <- list()
for(r in 1:47) {
  cs.est[[r]] <- lapply(16:49, function(ms) {
    f.est[15:(ms-1)]/sum(f.est[1:(ms-1)])
  })
  
  c_CEB.est[[r]] <- matrix(0, nrow=34, ncol=34)
  for(ms in 16:49) {
    for(a in 1:(ms+1-16)) {
      c_CEB.est[[r]][ms-15, a] <- sum(CEBa[[r]][ms-15,]) * cs.est[[r]][[ms-15]][ms - 16 + 2 - a]
    }
  }
}

# Step 2: fit model
Data$CEBa <- do.call("rbind", c_CEB.est)
ObjSBH.estc = MakeADFun( data=Data, parameters=Params,
                          random=c("phiz", "phio", "phif", "Us", "eps"), DLL="U_spat" )
t1 <- proc.time()
OptSBH.estc <- nlminb(start = ObjSBH.estc$par, objective = ObjSBH.estc$fn,
                       gradient = ObjSBH.estc$gr)

proc.time() - t1


## Without SBH
Data$sbhI <- 0
ObjFBH = MakeADFun( data=Data, parameters=Params, 
                    random=c("phiz", "phio", "phif", "Us", "eps"), DLL="U_spat" )
t1 <- proc.time()
OptFBH <- nlminb(start = ObjFBH$par, objective = ObjFBH$fn,
                      gradient = ObjFBH$gr)
proc.time() - t1


## Results
repSBHtrue <- sdreport(ObjSBH.true, getJointPrecision = T, bias.correct = T)
summary(repSBHtrue)
repSBHtruec <- sdreport(ObjSBH.truec, getJointPrecision = T, bias.correct = T)
summary(repSBHtruec)
repSBHestc <- sdreport(ObjSBH.estc, getJointPrecision = T, bias.correct = T)
summary(repSBHestc)
repFBH <- sdreport(ObjFBH, getJointPrecision = T, bias.correct = T)
summary(repFBH)
