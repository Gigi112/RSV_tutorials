fitmodel <-  function(parameters,dat) {   # takes the parameter values and dataset as inputs
 
  protrans <- parameters[1] # estimate parameter related to R0 (baseline transmission rate)
  b1 <- parameters[2] # estimate parameter related to seasonal amplitude
  trans <- parameters[3] # estimate parameter related to seasonal peak timing
  DMD <- parameters[4] # estimate parameter related to the duration of maternal immunity
  Amp <- exp(b1) #ensure positive 
  baseline.txn.rate <- exp(protrans) #ensure positive
  phi <-  (2*pi*(exp(trans))) / (1+exp(trans)) # transform to its scale in the model
  durx <- exp(DMD) #ensure positive
  
  #we need this for the function and model to recognize these parameters
  durx <<- durx 
  baseline.txn.rate <<- baseline.txn.rate
  Amp <<- Amp
  phi <<- phi
  
  
  # Run transmission model with initial conditions 
  # and timesteps defined above, and parameter values from function call
  results <- ode(y=yinit.vector, t=my_times,  
                 func=simple_model, 
                 parms=c(parmset,
                         baseline.txn.rate=baseline.txn.rate,
                         Amp=Amp,
                         phi=phi,
                         DurationMatImmunityDays=durx))
  
  t0 <- parmset$t0 
  # make sure the evaluation period are the same
  results<- tail(results,t0) # get rid of burn-in period
  St <- results[,-1]
  
  I1 <- St[,grep('I1', colnames(St))]
  I2 <- St[,grep('I2', colnames(St))]
  I3 <- St[,grep('I3', colnames(St))]
  I4 <- St[,grep('I4', colnames(St))]
  S0 <- St[,grep('S0', colnames(St))]
  S1 <- St[,grep('S1', colnames(St))]
  S2 <- St[,grep('S2', colnames(St))]
  S3 <- St[,grep('S3', colnames(St))]
  
  
  b=baseline.txn.rate/(parmset$dur.days1/30.44)
  beta=(b/100)/(sum(yinit.matrix)^(1-parmset$q))*parmset$contact
  
  #Force of infection
  lambda1=matrix(0,nrow=t0,ncol=parmset$N_ages)
  for (t in 1:t0)
  {lambda1[t,]<-as.vector((1+Amp*cos(2*pi*(t-phi*12)/12))*
                            ((I1[t,]+rho1*I2[t,]+rho2*I3[t,]+rho2*I4[t,])
                             %*%beta)/sum(St[t,]))}
  
  
  H1=matrix(0,nrow=t0,ncol=parmset$N_ages)#Number of hospitalizations by age
  for (i in 1:parmset$N_ages){
    H1[,i]=hosp1[i]*S0[,i]*lambda1[,i]+
      hosp2[i]*sigma1*S1[,i]*lambda1[,i]+
      hosp3[i]*sigma2*S2[,i]*lambda1[,i]+
      hosp3[i]*sigma3*S3[,i]*lambda1[,i]}
  
  H <- rowSums(H1)
  
  agedist <- c(sum(colSums(H1)[1:3]),sum(colSums(H1)[4:6]),
               sum(colSums(H1)[7:9]),
               sum(colSums(H1)[10:12]),
               colSums(H1)[13:21])/sum(H1)
  
  LLall <- sum(dpois(x = sim$Hosp_sim, lambda =H, log = TRUE)) 
  # fit to timeseries (number of cases follows poisson distribution)
  LLmulti <- dmultinom(x= sim$agedist_Sim,prob = agedist_Sim,log = T) 
  # fit to age distribution 
  # age distribution (proportion) follows mutinomial distribution
  
  #prior
  durprior <- dgamma(x=durx,22,5,log=T)
  # give a prior for the duration of maternal immunity  
  
  #total Loglikelihood (because of log, we sum up)
  LL <- LLall+LLmulti+durprior
  
  return(LL)
}