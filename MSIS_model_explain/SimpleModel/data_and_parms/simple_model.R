simple_model <- function(t,y,parms,time.step='month'){
  
  # read in initial states and their names
  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  # unify the time unit of parameter inputs
  if(parms$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }
  
  # waning rate of maternal immunity (by time step)
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  
  # aging rate (by time step)
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)}
  
  # rate of recovery of first infection
  gamma1= 1/(parms$dur.days1/length.step) 
  # rate of recovery of second infection
  gamma2= 1/(parms$dur.days2/length.step)
  # rate of recovery of third infection
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  
  #gamma3 stands for rate of recovery from subsequent infection
  
  # Relative risk of infection (2nd)
  sigma1=parms$sigma1 
  # Relative risk of infection (3rd)
  sigma2=parms$sigma2
  # Relative risk of infection (4th+)
  sigma3=parms$sigma3
  
  # Relative infectiousness (2nd)
  rho1=parms$rho1 
  # Relative infectiousness (3rd+)
  rho2=parms$rho2
  
  #Pull out the states for the model as vectors
  M <-  States[,'M'] # protected by maternal immunity
  S0 <-  States[,'S0'] # purely susceptible population
  I1 <-  States[,'I1'] # first time infection (infectious)
  
  S1 <-  States[,'S1'] 
  # susceptible population with build-up immunity
  I2 <-  States[,'I2'] # second time infection
  
  S2 <-  States[,'S2']
  # susceptible population with lower risk of re-infection
  I3 <-  States[,'I3'] # third time infection
  
  S3 <-  States[,'S3']
  # susceptible population with lowest risk of re-infection
  I4 <-  States[,'I4'] # subsequent time infection
  
  N_ages <- length(M) # the number of age groups
  
  ## parameter related to force of infection ################
  # per capita transmission probability
  baseline.txn.rate=parms$baseline.txn.rate
  # transmission probability per unit time
  b <- baseline.txn.rate/ (parms$dur.days1/length.step) 
  q=parms$q # q depends on transmission type 
  # (whether depends on population density or not)
  contact=parms$contact # c2 is the contact matrix
  # transmission probability per unit time in each age group
  beta <- (b/100)/(sum(yinit.matrix)^(1-q))*contact 
  # 100 is a scaling factor for the contact matrix we choose
  # (see Ginny's paper and Matlab code for details)
  # this does not matter because most likely
  # you will need to estimate baseline.txn.rate
  
  Amp=parms$Amp # seasonal amplitude
  phi=parms$phi # seasonal phase shift
  #seasonality
  seasonal.txn <- (1+Amp*cos(2*pi*(t-phi*period)/period))
  
  # seasonal transmission probability
  beta_a_i <- seasonal.txn * beta 
  infectiousN <- (I1+rho1*I2+rho2*I3+rho2*I4)/sum(States)
  # for frequency dependent transmission
  
  lambda <- infectiousN %*% beta_a_i # force of transmission
  lambda <- as.vector(lambda) # vectorize force of transmission
  
  # create a matrix to record the changing variables
  dy <- matrix(NA, nrow=N_ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- 
    log(parms$PerCapitaBirthsYear[t,]+1)/period
  # get period birth rate from annual birth rate
  # see the following page for birth rate calculation
  
  #um is death rate
  um=parms$um
  #mu represents aging to the next class
  Aging.Prop <- c(0,mu[1:(N_ages-1)])
  
  dy[,'M'] <- period.birth.rate*sum(States) - 
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N_ages-1)]) 
  
  dy[,'S0'] <- omega*M -
    lambda*S0 -
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N_ages-1)]) 
  
  dy[,'I1'] <-   lambda*S0 - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N_ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - 
    sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N_ages-1)]) 
  
  dy[,'I2'] <- sigma1*lambda*S1 - 
    gamma2*I2-(mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N_ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - 
    sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N_ages-1)]) 
  
  dy[,'I3'] <- sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N_ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 -
    sigma3*lambda*S3 -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N_ages-1)]) 
  
  dy[,'I4'] <- sigma3*lambda*S3 - 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N_ages-1)]) 
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}