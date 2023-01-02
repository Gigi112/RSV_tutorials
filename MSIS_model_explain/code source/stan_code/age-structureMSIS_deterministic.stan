functions {
  real[] msis (real t, real[] y, real[] theta,
             real[] x_r, int[] x_i) {
               
// state, the volumes in each compartment, y;
// theta, variables used to compute f, which depend on the model parameters;
// x_r, real variables used to evaluate f, which only depend on fixed data;
// x_i, integer values used to evaluate f, which only depend on fixed data.
      int agegroups = x_i[1];
      int q=x_i[2];
      
      real birthrate[agegroups] = x_r[1:agegroups];
      real um = x_r[agegroups+1];
      real rho1 = x_r[agegroups+2];
      real rho2 = x_r[agegroups+3];
      real gamma1 = x_r[agegroups+4];
      real gamma2 = x_r[agegroups+5];
      real gamma3 = x_r[agegroups+6];
      real sigma1 = x_r[agegroups+7];
      real sigma2 = x_r[agegroups+8];
      real sigma3 = x_r[agegroups+9];
      real u[agegroups]= x_r[(agegroups+10):(9+2*agegroups)];
      real c2[agegroups*agegroups]=x_r[(10+2*agegroups):(9+2*agegroups+agegroups*agegroups)];
      //real omega = x_r[2*agegroups+agegroups*agegroups+10];
      
      real beta = theta[1];
      real b1 = theta[2];
      real phi = theta[3];
      real omega = theta[4];
      
      real dydt[9*agegroups];
      real lambda[agegroups];
      real InfectN[agegroups];
       
      real season_txn = (1+b1*cos(2*pi()*(t-phi*12)/12));
      
      for (k in 1:agegroups) {
      InfectN[k] = (y[k+2*agegroups]+rho1*y[k+4*agegroups]+rho2*y[k+6*agegroups]+rho2*y[k+8*agegroups]);
      } 
      
      for (a in 1:agegroups) {
      
       lambda[a] = season_txn*beta*gamma1/sum(y)*sum(to_vector(c2[((a-1)*agegroups+1):(a*agegroups)]).*to_vector(InfectN));
      
       dydt[a] =  log(birthrate[a]+1)/12*sum(y) - (omega+u[a]+um)*y[a]; //M
       dydt[a+agegroups] =  omega*y[a] - lambda[a]*y[a+agegroups] - (um+u[a])*y[a+agegroups]; //S0
       dydt[a+2*agegroups] =  lambda[a]*y[a+agegroups] - gamma1*y[a+2*agegroups] - (um+u[a])*y[a+2*agegroups];//I1
       dydt[a+3*agegroups] =  gamma1*y[a+2*agegroups] - sigma1*lambda[a]*y[a+3*agegroups] - (um+u[a])*y[a+3*agegroups];//S1
       dydt[a+4*agegroups] =  sigma1*lambda[a]*y[a+3*agegroups] - gamma2*y[a+4*agegroups] - (um+u[a])*y[a+4*agegroups];//I2
       dydt[a+5*agegroups] =  gamma2*y[a+4*agegroups] - sigma2*lambda[a]*y[a+5*agegroups] - (um+u[a])*y[a+5*agegroups];//S2
       dydt[a+6*agegroups] =  sigma2*lambda[a]*y[a+5*agegroups] - gamma3*y[a+6*agegroups] - (um+u[a])*y[a+6*agegroups];//I3
       dydt[a+7*agegroups] =  gamma3*y[a+6*agegroups] + gamma3*y[a+8*agegroups] - sigma3*lambda[a]*y[a+7*agegroups] - (um+u[a])*y[a+7*agegroups];//S3
       dydt[a+8*agegroups] =  sigma3*lambda[a]*y[a+7*agegroups] - gamma3*y[a+8*agegroups] - (um+u[a])*y[a+8*agegroups];
       
       if (a >1 ){
		    dydt[a] = dydt[a]+ u[a-1]*y[a-1];
		    dydt[a+agegroups] = dydt[a+agegroups] + u[a-1]*y[a+agegroups-1];
		    dydt[a+2*agegroups] = dydt[a+2*agegroups] + u[a-1]*y[a+2*agegroups-1];
		    dydt[a+3*agegroups] =dydt[a+3*agegroups]+u[a-1]*y[a+3*agegroups-1];
		    dydt[a+4*agegroups] = dydt[a+4*agegroups]+u[a-1]*y[a+4*agegroups-1];
		    dydt[a+5*agegroups] = dydt[a+5*agegroups]+u[a-1]*y[a+5*agegroups-1];
		    dydt[a+6*agegroups] = dydt[a+6*agegroups]+u[a-1]*y[a+6*agegroups-1];
		    dydt[a+7*agegroups] = dydt[a+7*agegroups]+u[a-1]*y[a+7*agegroups-1];
		    dydt[a+8*agegroups] = dydt[a+8*agegroups]+u[a-1]*y[a+8*agegroups-1];
		}
       
       }
       
       
      
      return dydt;
  }
}

// Fixed data is declared in the data block:
data {
  int<lower=1> n_months;
  int<lower=1> agegroups;
  int q;
  real y0[9*agegroups];
  real t0;
  real ts[n_months];
  int N;
  real birthrate[agegroups];
  real um;
  real rho1;
  real rho2;
  real gamma1;
  real gamma2;
  real gamma3;
  real sigma1;
  real sigma2;
  real sigma3;
  int<lower=1> hosp_cases[n_months];
  int hosp_age[8];
  real hosp1[agegroups];
  real hosp2[agegroups];
  real hosp3[agegroups];
  real c2[agegroups*agegroups];
  real u[agegroups];
}


transformed data {
  real x_r[204];
  int x_i[2];
  x_i[1] = agegroups;
  x_i[2] = q;
  
  x_r[1:agegroups]= birthrate;
  x_r[agegroups+1]=um;
  x_r[agegroups+2]=rho1;
  x_r[agegroups+3]=rho2;
  x_r[agegroups+4]=gamma1;
  x_r[agegroups+5]=gamma2;
  x_r[agegroups+6]=gamma3;
  x_r[agegroups+7]=sigma1;
  x_r[agegroups+8]=sigma2;
  x_r[agegroups+9]=sigma3;
  x_r[(agegroups+10):(9+2*agegroups)]=u;
  x_r[(10+2*agegroups):(9+2*agegroups+agegroups*agegroups)]=c2;
}

parameters {
  real<lower=0> beta;
  real<lower=0,upper=1> b1;
  real<lower=0,upper=2*pi()> phi;
  real<lower=0> omega;
  real<lower=0,upper=1> report_ratio;}

transformed parameters{
  real y[n_months, 9*agegroups]; // y is the states matrix that has row length=n_months and columns=9.
  real rel_inf[n_months,agegroups];
  real lambda[n_months, agegroups];
    // outcomes
  matrix[n_months,agegroups] HOSP_AGE;
  vector[8] output_hosp_age;
  vector[n_months] total_children;
  vector<lower = 0>[n_months] output_hosp; // overall case incidence by month
  real theta[4];
  
  {
    theta[1] = beta;
    theta[2] = b1;
    theta[3] = phi;
    theta[4] = omega;

    
    y = integrate_ode_bdf(msis, y0, t0, ts, theta, x_r, x_i,1.0E-10, 1.0E-10, 1.0E3);

  }
  
  
  for(i in 1:n_months){
    
        for(j in 1:9*agegroups){
    if (y[i,j] <= 0.0) y[i,j] = 1e-12;}
    
      for (k in 1:agegroups) {
      rel_inf[i,k] = (y[i,k+2*agegroups]+rho1*y[i,k+4*agegroups]+rho2*y[i,k+6*agegroups]+rho2*y[i,k+8*agegroups]);
      }
      
      for (a in 1:agegroups) {
      
       lambda[i,a] = (1+b1*cos(2*pi()*(i-phi*12)/12))*beta*gamma1/sum(y[i,])*sum(to_vector(c2[((a-1)*agegroups+1):(a*agegroups)]).*to_vector(rel_inf[i,]));
              HOSP_AGE[i,a] = 
			hosp1[a]*lambda[i,a]*y[i,a+agegroups]
			+hosp2[a]*lambda[i,a]*sigma1*y[i,a+3*agegroups]	
			+hosp3[a]*lambda[i,a]*sigma2*y[i,a+5*agegroups]
			+hosp3[a]*lambda[i,a]*sigma3*y[i,a+7*agegroups];
       }

    output_hosp[i] = report_ratio*sum(HOSP_AGE[i,]);
    total_children[i]=HOSP_AGE[i,1]+HOSP_AGE[i,2]+HOSP_AGE[i,3]+
    HOSP_AGE[i,4]+HOSP_AGE[i,5]+HOSP_AGE[i,6]+HOSP_AGE[i,7]+HOSP_AGE[i,8];
    }

    for (b in 1:8) {
    output_hosp_age[b] = sum(HOSP_AGE[,b])/sum(total_children);}
    
}

model {
  //priors
  beta ~ normal(3,1) T[1,5]; //truncated at 0
  b1 ~ normal(0.2,0.05) T[0.05,1]; //truncated at 0,1
  phi ~ normal(1.4,0.5) T[0,2*pi()];
  omega ~ lognormal(-1,0.6) T[0,5];
  report_ratio ~ beta(2,2); //truncated at 0,1
  
  //sampling distribution
  for(i in 1:n_months){
  target += poisson_lpmf(hosp_cases[i]| output_hosp[i]);}
  
  target += multinomial_lpmf (hosp_age|output_hosp_age);
}
