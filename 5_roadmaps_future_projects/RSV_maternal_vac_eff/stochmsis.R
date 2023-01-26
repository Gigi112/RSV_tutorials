### stochastic MSIS model to record flows of individuals
### start from no age structure

stochmsis = function(fixed,parms){ ## set =1 for fixed infectious period, =0 for exponentially distributed
  ### input parameter
  lambda=parms$fixed_lambda 
  # transmission rate in general population at each time point
  ## time-varying transmission rate derived from deterministic model in quasi-equilibrium
  
  birth=parms$birth
  # instead of birth[t] = rbinom(1,N[t-1],birthrate[t-1])
  # we treat birth number as an input over time to simulate a birth cohort
  # for the clinical trial.
  # check whether this is correct.
  
  #########################################
  # duration of infectiousness 
  #Duration in days
  dur.days1 <- 10 #days; first time infection
  dur.days2 <- 7 #days; second time infection
  dur.days3 <- 5 #days; third time and subsequent infection
  ###########################################  
  
  ###########################################
  # 1/duration of maternal immunity (DAYS)
  DurationNaturalProtection= 30
  ###########################################
  
  #############################################
  #Relaive risk of infection following 1st, 2nd, 3rd+ infections
  sigma1=0.76
  sigma2=0.6
  sigma3=0.4
  #############################################
  
  # um= -0.0002227 ## death rate
  # assume no death for the first try. Do we need death in the model?
  # for the whole clinical trial period, deaths should be negligible
  
  tmax = 365 ### max number of days to simulate
  dt = 1 ## time step for approximate method 1 day instead of 1/20 day
  
  N = M = S0 = I1 = S1 = I2 = S2 = I3 = S3 = I4 = rep(0,tmax)
  # set up a vector to store the number of individuals in each compartment 
  
  Waning = Infect1 = Infect2 = Infect3 = Infect4 = 
    Recover1 =Recover2 = Recover3 = Recover4 = rep(0,tmax)
  # create vectors to store the transition flows at each time point
  
  N[1] = 1000 ### initial number of people
  M[1] = 1000
  S0[1] = I1[1] = S1[1] = I2[1] = S2[1] = I3[1] = S3[1] = I4[1] = 0
  Infect1[1] = Infect2[1] =Infect3[1] =Infect4[1] =0
  Waning[1] = Recover1[1] =Recover2[1] = Recover3[1]= Recover4[1]= 0

  birth[1] = 1000 # in the first time point, we assumed 1000 infants into the birth cohort
  # we can change this value later.
  
  ### Start loop
  
  time = rep(0,tmax)
  t = 2
  
  while (t<=tmax){

    ### Rates
    
    if (fixed==1){ ### assuming fixed protected/infectious period
      Infect1[t] = rbinom(1,S0[t-1],lambda[t-1])
      Infect2[t] = rbinom(1,S1[t-1],sigma1*lambda[t-1])
      Infect3[t] = rbinom(1,S2[t-1],sigma2*lambda[t-1])
      Infect4[t] = rbinom(1,S3[t-1],sigma3*lambda[t-1])
        
      if (t>30){
        Waning[t] = birth[t-30]
      }  
      if (t>10){
        Recover1[t] = Infect1[t-10]}
        if (t>7){
          Recover2[t] = Infect2[t-7]
        }
      if (t>5){
        Recover3[t] = Infect3[t-5]
        Recover4[t] = Infect4[t-5]
      }
    } else{ ### assuming exponentially distributed protected/infectious periods
      ### number of events is Poisson distributed if we assume occurrence at a constant rate over a short time step
      Infect1[t] = rpois(1,S0[t-1]*lambda[t-1]*dt) 
      Infect2[t] = rpois(1,S1[t-1]*sigma1*lambda[t-1]*dt)
      Infect3[t] = rpois(1,S2[t-1]*sigma2*lambda[t-1]*dt)
      Infect4[t] = rpois(1,S3[t-1]*sigma3*lambda[t-1]*dt)
      
      Recover1[t] = rpois(1,(1/dur.days1)*I1[t-1]*dt)
      Recover2[t] = rpois(1,(1/dur.days2)*I2[t-1]*dt)
      Recover3[t] = rpois(1,(1/dur.days3)*I3[t-1]*dt)
      Recover4[t] = rpois(1,(1/dur.days3)*I4[t-1]*dt)
      
      Waning[t] = rpois(1,(1/DurationNaturalProtection)*M[t-1]*dt)
      
    }
    
    ### Transmission dynamics
    
    M[t] = max(M[t-1]+birth[t]-Waning[t],0)
    S0[t] = max(S0[t-1]+Waning[t]-Infect1[t],0)
    I1[t] = max(I1[t-1]+Infect1[t]-Recover1[t],0)
    S1[t] = max(S1[t-1]+Recover1[t]-Infect2[t],0)
    I2[t] = max(I2[t-1]+Infect2[t]-Recover2[t],0)
    S2[t] = max(S2[t-1]+Recover2[t]-Infect3[t],0)
    I3[t] = max(I3[t-1]+Infect3[t]-Recover3[t],0)
    S3[t] = max(S3[t-1]+Recover3[t]+Recover4[t]-Infect4[t],0)
    I4[t] = max(I4[t-1]+Infect4[t]-Recover4[t],0)
    
    N[t] = (M+S0+I1+S1+I2+S2+I3+S3+I4)[t]
    
    tf = t + 1
    t = t+1
  }
  
  if (fixed==1){
    NewCase = rep(0,(tf-1))
    for (i in 1:(tf-1)) {
      NewCase[i] = rbinom(1,Infect1[i],parms$hosp1)
      # hosp1 is the probability of hospitalization for
      # first time RSV infection
    }
  } else{
    NewCase = rep(0,max(1,floor(tf*dt)-1)) # create a vector that has the same length of the observed days
    for (i in 1:max(1,floor(tf*dt)-1)){
      NewCase[i] = rbinom(1,sum(Infect1[((i-1)*(1/dt)+1):(i/dt)]),parms$hosp1)
      # is this correct? I am afraid that if I use poisson to
      # get the hospitalizations in each time step (1/20 day)
      # the the number of RSV hospitalizations will likely be 0
      # therefore, I use binomial distribution on a daily basis
      }
  }
  
  return(NewCase)
}

### To test model, uncomment the following lines

### some random parameters to test the model
# input=list(fixed_lambda=rep(1.093694,365),
# birth=c(seq(1000,0,-5),rep(0,360)),
# hosp1=0.08)

### when fixed set =1 for fixed latent/infectious period, 
### when fixed set =0 for exponentially distributed
# stochmsis(1,parms = input)
# stochmsis(0,parms = input)