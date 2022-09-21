
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

