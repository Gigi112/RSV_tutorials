waning_immunity_model <- function(t,y,parms, time.step='month'){
  
  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  if(parms$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }
  
  # if(t>=322 & t<= 335){
  #   I_ex_t= sum(States)*parms$ex_proportion*parms$ex_increase_new[t-321] # ex_increase_new is a vector indicating how 
  #   # travel restritions and virus activities in other regions affects virus importation
  # }else{
  #   I_ex_t <- sum(States)*parms$ex_proportion 
  #   ## external infections
  #   ##  I_ex_t=1/100000*St ; I_ex_t=5/100000*St ; I_ex_t=10/100000*St and etc
  # }
  # 
  
  I_ex_t=parms$I_ex*parms$airtravel[t]
  
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/lenth.step
  gamma2= 1/(parms$dur.days2/length.step)  
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  #??Is this right?? Yes, gamma3 stands for rate of recovery from subsequent infection
  
  if(parms$AllowWaning=='Yes'){
  theta1= 1/(parms$dur.immunity1/length.step) # waning immunity from S3 to W3;from S2 to W2
  theta2= 1/(parms$dur.immunity2/length.step) }
  else{theta1= 0 # waning immunity from S3 to W3;from S2 to W2
  theta2= 0}# waning immunity from W3 to S2;from W2 to S1
  
  #Pull out the states  for the model as vectors
  M <-  States[,'M']
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  
  W1 <-  States[,'W1']
  W2 <-  States[,'W2']
  W3 <-  States[,'W3']
  
  N.ages <- length(M)
  
  baseline.txn.rate[t] <- parms$R0*parms$NPIs[t]
  
  ###Check the standardization of beta and overall structure of lambda here
  #how does'baseline txn rate' figure in here?
  ##???##################
  seasonal.txn <- (1+parms$b1*cos(2*pi*(t-parms$phi*period)/period))# seasonality waves
  b <- baseline.txn.rate[t]/ (parms$dur.days1/length.step) # transmission probability per unit time
  beta <-  (b/100)/(sum(yinit.matrix)^(1-parms$q))*parms$c# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix
  
  beta_a_i <- seasonal.txn * beta/sum(States)
  infectiousN <- I1 + rho1*I2 + rho2*I3 + rho2*I4+ rho2*I_ex_t
  
  lambda <- infectiousN %*% beta_a_i 
  lambda <- as.vector(lambda)
  ##########?????????????????##########################  
  
  dy <- matrix(NA, nrow=N_ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #B=annual birth rate
  
  #https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1004591.s016&type=supplementary
  #mu represents aging to the next class
  #um is death rate
  um=parms$um
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  
  dy[,'M'] <- period.birth.rate*sum(States) - 
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N.ages-1)]) 
  
  dy[,'S0'] <- omega*M -
    lambda*S0 -
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N.ages-1)]) +
    theta2*W1
  
  dy[,'I1'] <-   lambda*S0 - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - 
    sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) +
    theta2*W2-
    theta1*S1 
  
  dy[,'W1'] <-  theta1*S1-
    sigma1*lambda*W1 -
    (mu + um)*W1 + 
    Aging.Prop*c(0,W1[1:(N.ages-1)])-
    theta2*W1
  
  dy[,'I2'] <- sigma1*lambda*S1 - 
    gamma2*I2-
    (mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) +
    sigma1*lambda*W1
  
  dy[,'S2'] <- gamma2*I2 - 
    sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N.ages-1)])+
    theta2*W3 -
    theta1*S2 
  
  dy[,'W2'] <-  theta1*S2-
    sigma2*lambda*W2 -
    (mu + um)*W2 + 
    Aging.Prop*c(0,W2[1:(N.ages-1)])-
    theta2*W2
  
  dy[,'I3'] <- sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) +
    sigma2*lambda*W2
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 -
    sigma3*lambda*S3 -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)]) -
    theta1*S3 
  
  dy[,'W3'] <-  theta1*S3-
    sigma3*lambda*W3 -
    (mu + um)*W3 + 
    Aging.Prop*c(0,W3[1:(N.ages-1)])-
    theta2*W3
    
  
  dy[,'I4'] <- sigma3*lambda*S3 +
    sigma3*lambda*W3- 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}
