maternalvac_cohort <- function(t,y,parms, time.step='week'){
  
  yinit.matrix <- parms$yinit.matrix
  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  if(parms$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }
  
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  RR=parms$RR # relative risk of infection after vaccination in infants
  waningV=1/(parms$DurationVE/length.step)
  rrM=1 # relative risk of infection after vaccination in pregnant women
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(parms$WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/lenth.step
  gamma2= 1/(parms$dur.days2/length.step)  
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  #gamma4=gamma3 stands for rate of recovery from subsequent infection
  
  rho1=parms$rho1#Relative infectiousness for 2nd infections
  rho2=parms$rho2 #Relative infectiousness for subsequent infections
  
  sigma1=parms$sigma1#Relaive risk of infection following 1st infections
  sigma2=parms$sigma2#Relaive risk of infection following 2nd infections
  sigma3=parms$sigma3#Relaive risk of infection following subsequent infections
  
  #Pull out the states  for the model as vectors
  Mv <-  States[,'Mv']
  
  Sv <-  States[,'Sv']
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  
  N.ages <- length(Mv)
  
  Vm <-  States[,'Vm'] # vaccination of pregnant women
  
  ##fixed lambda based on lambda from burn-in period##################
  lambda <- parms$lambda[t,]
  ##########ODE equations##########################  
  
  dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #B=annual birth rate
  
  #https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1004591.s016&type=supplementary
  #mu represents aging to the next class
  #um is death rate
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  um <- parms$um
  
  # newborns with vaccination of mothers is 
  # protected by natureal maternal immunity
  dy[,'Mv'] <- period.birth.rate*parms$pop - 
    (omega+(mu+um))*Mv +
    Aging.Prop*c(0,Mv[1:(N.ages-1)])
  
  # after the maternal immunity wanes, 
  # newborns with vaccination of mothers has lower risk of infections 
  dy[,'Sv'] <- omega*Mv - 
    waningV*Sv-
    RR*lambda*Sv -
    (mu + um)*Sv + 
    Aging.Prop*c(0,Sv[1:(N.ages-1)]) 
  
  dy[,'S0'] <- waningV*Sv-
    lambda*S0 -
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N.ages-1)]) 
  
  dy[,'I1'] <-   lambda*S0+ RR*lambda*Sv - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - 
    sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- sigma1*lambda*S1 - 
    gamma2*I2-
    (mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - 
    sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 +
    waningV*Vm-
    sigma3*lambda*S3 -
    c(rep(0,18),period.birth.rate[1]*sum(States),rep(0,2)) -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)]) 
  
  # pregnant women (age 20-39) will enter the Vm compartment; the number is equal to newborn (approximately because most of births are single live birth)
  # I did not add any reduced susceptibility or infectiousness at this point and therefore Vm is no different than S3
  # change rrM will allow us to evaluate the effect of reduce risk in pregnant
  # women on vaccine efficacy
  dy[,'Vm'] <- c(rep(0,18),period.birth.rate[1]*sum(States),rep(0,2)) -
    rrM*sigma3*lambda*Vm -
    (mu + um)*Vm -
    waningV*Vm+ 
    Aging.Prop*c(0,Vm[1:(N.ages-1)]) 
  
  dy[,'I4'] <- sigma3*lambda*S3+rrM*sigma3*lambda*Vm - 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}