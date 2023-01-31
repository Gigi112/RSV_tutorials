model_string<-"
model {
for(k in 1:n.group){  
  for (j in 1:n.age){
      for (i in 1:n.date) { 
    
   lambda[i,j,k] <-   rd0[j,k]+
               exp(epi[epi.year[i]]) +
               exp(delta[month[i]]) +
               rsv[i,j,k]*rd2[j,k] + 
               flu[i,j,k]*rd1[epi.year[i],j,k]
  y[i,j,k] ~ dnegbin(prob[i,j,k],r)
  prob[i,j,k]<- r/(r+lambda[i,j,k])  ## likelihood 
}
  # baseline hospitalizations (must be greater than or equal to 0)
  rd0[j,k] <- exp(beta0[j,k])
  beta0[j,k]~ dnorm(mu0,tau0)
  
  for (p in 1:n.year) { 
    # epi-year effects
    # random slope
    # coefficient of influenza-associated respiratory hospitalization
    # this coefficient varies annually
    rd1[p,j,k] <- exp(beta1[p,j,k]) # ensure positive coefs
    beta1[p,j,k] ~ dnorm(beta1_mean[p,j,k],tau.flu)
    beta1_mean[p,j,k] <- gamma_flu[k]+ omega_flu[j]+ xi_flu[p]
    }
  
    # coefficient of RSV-associated respiratory hospitalization
    # this coefficient depends on SES and age 
    rd2[j,k] <- exp(beta2[j,k]) # ensure positive
    beta2[j,k] ~ dnorm(beta2_mean[j,k],tau.rsv)
    beta2_mean[j,k]<- gamma[k]+omega[j]
    
  
  # Month dummy variable, account for other pathogens and general seasonality 
}}
## hyperpriors
for (p in 1:n.year) { 
epi[p] ~ dnorm(0, tau.epi) }
  for (m in 1:12){
    delta[m] ~ dnorm(0,disp.m)
}
for(k in 1:n.group){
gamma_flu[k] ~ dnorm(0,tau3)
gamma[k] ~ dnorm(0, tau4)
}
for(j in 1:n.age){
omega[j] ~ dnorm(0, tau5)
omega_flu[j] ~ dnorm(0,tau6)
}
for (p in 1:n.year) { 
xi_flu[p] ~ dnorm(0,tau7)
}
  r ~ dunif(0,250)
  mu0 ~ dnorm(0,0.0001)
  tau0 ~ dgamma(0.01, 0.01)
  tau1 ~ dgamma(0.01, 0.01)
  tau2 ~ dgamma(0.01, 0.01)
  tau3 ~ dgamma(0.01, 0.01)
  tau4 ~ dgamma(0.01, 0.01)
  tau5 ~ dgamma(0.01, 0.01)
  tau6 ~ dgamma(0.01, 0.01)
  tau7 ~ dgamma(0.01, 0.01)
  tau8 ~ dgamma(0.01, 0.01)
  tau.epi ~ dgamma(0.01, 0.01)
  tau.flu ~ dgamma(0.01, 0.01)
  tau.rsv ~ dgamma(0.01, 0.01)
  disp.m ~ dgamma(0.01, 0.01)
}
"