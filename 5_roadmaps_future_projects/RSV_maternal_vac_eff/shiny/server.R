library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")
library("tibble")
# get all prespecified parameters. This is only for 
# simulation purpose (mostly from US, NL and Kenya info). 
# Need to take a futher look at it and
# decide whether all parameters are appropriate.
parms_sim <- readRDS("./data_and_parms/parms_sim.rds")

source("./burnin.R")  # R script for burn-in period
source("./placebo.R") # R script for placebo cohort
source("./vaccine.R") # R script for vaccine cohort
al = nrow(parms_sim$yinit.matrix)
states <- ncol(parms_sim$yinit.matrix)
hosp1=c(.18*rep(.40,3),0.08*rep(.39,3),0.07*rep(.21,3),0.06*rep(.20,3),0.06*0.16,0.05*rep(.14,3),0.02*rep(0.05,al-16))#proportion of first LRIs that are hospitalized
hosp2=.4*hosp1#proportion of second LRIs that are hospitalized
hosp3=c(rep(0,al-2),0.00001,0.00004)#proportion of third LRIs that are hospitalized

# Define server 
#
server <- function(input, output) {
  
  sim_mat_placebo_vaccine <- reactive({ 
    
  annual_birth_rate<- input$annual_birth_rate
  yinit_matrix=as.array(parms_sim$yinit.matrix)
  yinit_vector <- as.vector(parms_sim$yinit.matrix) #Vectorize the ynit matrix
  # Create array that has the labels by age, State and use this to name the yinit.vector
  name.array <- array(NA, dim=dim(yinit_matrix))
  for(i in 1:dim(name.array)[1]){
    for(j in 1:dim(name.array)[2]){
      name.array[i,j] <- paste(dimnames(yinit_matrix)[[1]][i],dimnames(yinit_matrix)[[2]][j]  )
    }
  }
  
  name.vector <- as.vector(name.array)
  names(yinit_vector) <- name.vector
  
  tmax=input$t_obs+24*52 # burn-in period as 24 years
  PerCapitaBirthsYear <-  matrix(c(annual_birth_rate,rep(0,nrow(yinit_matrix)-1)),
                                 byrow = T,nrow=tmax,ncol=nrow(yinit_matrix))
  
  parameters_burnin <- list(
    PerCapitaBirthsYear=PerCapitaBirthsYear,
    DurationMatImmunityDays=input$DurationNaturalProtection,
    b1=input$b1,
    phi=input$phase,
    baseline.txn.rate = input$baseline.txn.rate,
    time.step='week',
    yinit.matrix=yinit_matrix)
  parms_sim["yinit.matrix"] <- NULL
  parms_input <- parms_sim
  parameters_burnin_append <- append(parms_input,parameters_burnin)
  
  ## Time frame
  times_burnin <- seq(1, 24*52, by = 1)
  
  ## Solve using ode (General Solver for Ordinary Differential Equations)
  out_burnin <- ode(
    y = yinit_vector,
    t = times_burnin,
    func = burnin_quasi_equilibrium,
    parms = parameters_burnin_append) 

    ## get hospitalizations after burn-in to evaluate whether it reaches quasi-equilibium
  #get the relavent variables to calculate RSV hospitalizations
  I1_burnin <- out_burnin[,grep('I1', colnames(out_burnin))]
  I2_burnin <- out_burnin[,grep('I2', colnames(out_burnin))]
  I3_burnin <- out_burnin[,grep('I3', colnames(out_burnin))]
  I4_burnin <- out_burnin[,grep('I4', colnames(out_burnin))]
  S0_burnin <- out_burnin[,grep('S0', colnames(out_burnin))]
  S1_burnin <- out_burnin[,grep('S1', colnames(out_burnin))]
  S2_burnin <- out_burnin[,grep('S2', colnames(out_burnin))]
  S3_burnin <- out_burnin[,grep('S3', colnames(out_burnin))]
  
  # transmission per contact per unit time
  b= input$baseline.txn.rate/(parms_input$dur.days1/30.44) 
  # transmission * contact pattern
  beta=(b/100)/(sum(parms_input$yinit.matrix)^(1-parms_input$q))*parms_input$c 
   
  if(parameters_burnin$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parameters_burnin$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }
  
  lambda1_burnin=matrix(0,nrow=max(times_burnin),ncol=nrow(yinit_matrix))#Force of infection
  total_pop <- rowSums(out_burnin[,-1])
  seasonal.txn <- c()
  period.birth.rate <- c()
  # number of primary infectious individuals
  infectiousN_burnin <- matrix(0,nrow=max(times_burnin),ncol=nrow(yinit_matrix)) 
  beta_a_i=array(0,dim=c(max(times_burnin),nrow(yinit_matrix),nrow(yinit_matrix)))
  for (t in 1:max(times_burnin)){
    seasonal.txn[t] <- (1+input$b1*cos(2*pi*(t-input$phase*period)/period))# seasonality waves
    beta_a_i[t,,] <- seasonal.txn[t]*beta/total_pop[t]
    infectiousN_burnin[t,] <- I1_burnin[t,] + parms_input$rho1*I2_burnin[t,] + parms_input$rho2*I3_burnin[t,] + parms_input$rho2*I4_burnin[t,]
    lambda1_burnin[t,] <-  infectiousN_burnin[t,] %*% beta_a_i[t,,]}
    
    H_burnin=matrix(0,nrow=max(times_burnin),ncol=nrow(yinit_matrix))
    #Number of hospitalizations by age
      for (i in 1:nrow(yinit_matrix)){
      H_burnin[,i]=round(hosp1[i]*S0_burnin[,i]*lambda1_burnin[,i] +
      hosp2[i]*parms_input$sigma1*S1_burnin[,i]*lambda1_burnin[,i]+
      hosp3[i]*parms_input$sigma2*S2_burnin[,i]*lambda1_burnin[,i]+
      hosp3[i]*parms_input$sigma3*S3_burnin[,i]*lambda1_burnin[,i]
      )}
    
  fixed_lambda <- tail(lambda1_burnin,input$t_obs)
  #parameters_whole <- append(parms_sim,parameters_input)
  times_evaluation <- seq(1, input$t_obs, by = 1)
  
  y.matrix.placebo <- array(0, dim=c(al, states ))
  
  dimnames(y.matrix.placebo)[[1]] <- dimnames(yinit_matrix)[[1]]
  dimnames(y.matrix.placebo)[[2]] <- dimnames(yinit_matrix)[[2]]
  y.matrix.placebo[1,1] <- 2025
  
  parameters_placebo <- list(
    PerCapitaBirthsYear=PerCapitaBirthsYear,
    DurationMatImmunityDays=input$DurationNaturalProtection,
    lambda=fixed_lambda,
   time.step='week',
    yinit.matrix=y.matrix.placebo,
   pop=total_pop[length(times_burnin)])
  
  parameters_placebo_append <- append(parms_input,parameters_placebo)
  
  y.vector.placebo <- as.vector(y.matrix.placebo)
  names(y.vector.placebo)<- name.vector
  
  out_placebo <- ode(
    y = y.vector.placebo,
    t = times_evaluation,
    func = placebo_cohort,
    parms = parameters_placebo_append)
  
  y.matrix.vaccine <- array(0, dim=c(al,states +2))
  
  dimnames(y.matrix.vaccine)[[1]] <- dimnames(yinit_matrix)[[1]]
  dimnames(y.matrix.vaccine)[[2]] <- c("Mv","Sv",
    dimnames(yinit_matrix)[[2]][2:9],"Vm")
  
  y.matrix.vaccine [1,1] <- 2025
  
  y.vector.vaccine<- as.vector(y.matrix.vaccine)
  
  # Create array that has the labels by age, State and use this to name the yinit.vector
  name.array.vaccine <- array(NA, dim=dim(y.matrix.vaccine))
  for(i in 1:dim(name.array.vaccine)[1]){
    for(j in 1:dim(name.array.vaccine)[2]){
      name.array.vaccine[i,j] <- paste(dimnames(y.matrix.vaccine)[[1]][i],dimnames(y.matrix.vaccine)[[2]][j]  )
    }
  }
  
  name.vector.vaccine <- as.vector(name.array.vaccine)
  names(y.vector.vaccine)<- name.vector.vaccine
  
  parameters_vaccine <- list(
    PerCapitaBirthsYear=PerCapitaBirthsYear,
    DurationMatImmunityDays=input$DurationNaturalProtection,
    RR=0.3,
    DurationVE=input$DurationMatVac,
    lambda=fixed_lambda,
    time.step='week',
    yinit.matrix=y.matrix.vaccine,
    pop=total_pop[length(times_burnin)])
  
  parameters_vaccine_append <- append(parms_input,parameters_vaccine)
  
  out_vaccine <- ode(
    y = y.vector.vaccine,
    t = times_evaluation,
    func = maternalvac_cohort,
    parms = parameters_vaccine_append)
  
  #get the relavent variables to calculate RSV hospitalizations
  I1_placebo <- out_placebo[,grep('I1', colnames(out_placebo))]
  I2_placebo <- out_placebo[,grep('I2', colnames(out_placebo))]
  I3_placebo <- out_placebo[,grep('I3', colnames(out_placebo))]
  I4_placebo <- out_placebo[,grep('I4', colnames(out_placebo))]
  S0_placebo <- out_placebo[,grep('S0', colnames(out_placebo))]
  S1_placebo <- out_placebo[,grep('S1', colnames(out_placebo))]
  S2_placebo <- out_placebo[,grep('S2', colnames(out_placebo))]
  S3_placebo <- out_placebo[,grep('S3', colnames(out_placebo))]
  
  I1_vaccine <- out_vaccine[,grep('I1', colnames(out_vaccine))]
  I2_vaccine <- out_vaccine[,grep('I2', colnames(out_vaccine))]
  I3_vaccine <- out_vaccine[,grep('I3', colnames(out_vaccine))]
  I4_vaccine <- out_vaccine[,grep('I4', colnames(out_vaccine))]
  Sv_vaccine <- out_vaccine[,grep('Sv', colnames(out_vaccine))]
  S0_vaccine <- out_vaccine[,grep('S0', colnames(out_vaccine))]
  S1_vaccine <- out_vaccine[,grep('S1', colnames(out_vaccine))]
  S2_vaccine <- out_vaccine[,grep('S2', colnames(out_vaccine))]
  S3_vaccine <- out_vaccine[,grep('S3', colnames(out_vaccine))]
  
  total_Pop <- rowSums(out_vaccine[,-1])
  
  num_births <- c()
  period.birth.rate <- c()
   for (t in 1:max(times_evaluation)){
  period.birth.rate[t] <- parameters_vaccine_append$PerCapitaBirthsYear[t,1]/period #B=annual birth rate
  num_births[t] <- round(period.birth.rate[t]*total_Pop[t])
  }
  
  H_placebo=matrix(0,nrow=max(times_evaluation),ncol=nrow(yinit_matrix))#Number of hospitalizations by age
  H_vaccine=matrix(0,nrow=max(times_evaluation),ncol=nrow(yinit_matrix))#Number of hospitalizations by age
  for (i in 1:nrow(yinit_matrix)){
    H_placebo[,i]=hosp1[i]*S0_placebo[,i]*fixed_lambda[,i]
      #                   +
      # hosp2[i]*parms_input$sigma1*S1_placebo[,i]*lambda1_placebo[,i]+
      # hosp3[i]*parms_input$sigma2*S2_placebo[,i]*lambda1_placebo[,i]+
      # hosp3[i]*parms_input$sigma3*S3_placebo[,i]*lambda1_placebo[,i]
    H_vaccine[,i]=hosp1[i]*Sv_vaccine[,i]*fixed_lambda[,i]*parameters_vaccine$RR+
                          hosp1[i]*S0_vaccine[,i]*fixed_lambda[,i]
      #                   +
      # hosp2[i]*parms_input$sigma1*S1_vaccine[,i]*lambda1_vaccine[,i]+
      # hosp3[i]*parms_input$sigma2*S2_vaccine[,i]*lambda1_vaccine[,i]+
      # hosp3[i]*parms_input$sigma3*S3_vaccine[,i]*lambda1_vaccine[,i]
      
  }
  return(list("H_placebo"=as.data.frame(H_placebo),
              "H_vaccine"=as.data.frame(H_vaccine),
              "num_births"=2025*4.33,
              "H_burnin"=H_burnin[781:max(times_burnin),]))
  })
  
  output$Plot_placebo_vaccine <- renderPlotly({
      library(xts)
      timemax=input$t_obs
      dates <- seq(as.Date("2024-02-27"), length=timemax, by="weeks")
      
      H_placebo <- sim_mat_placebo_vaccine()$H_placebo 
      H_vaccine <- sim_mat_placebo_vaccine()$H_vaccine
      num_births <- sim_mat_placebo_vaccine()$num_births
      
      H_placebo_survival <- c() # choose an observation of 180 days
      H_vaccine_survival <- c() # choose an observation of 180 days
      
      # by week
      Cummulative_H_placebo <- c()
      Cummulative_H_placebo[1] <- H_placebo[1,1]
      H_week_placebo <- c()
      H_obs_t_p <- c() # time of observation in hospitalized patients H*(t-1)
      H_obs_t_p[1] <- 0
      for (i in 2:(input$t_obs-1)) {
        H_week_placebo[i] <- H_placebo[i,floor(i/4.33)+1]
        H_obs_t_p[i] <- H_obs_t_p[i-1]+H_week_placebo[i]*(i-1)
        Cummulative_H_placebo[i] <- Cummulative_H_placebo[i-1]+H_week_placebo[i]
      }
      
      for (i in 1:(input$t_obs-1)) {
      H_placebo_survival[i] <- (num_births[input$time_enroll]-Cummulative_H_placebo[i])/num_births[input$time_enroll] # choose to observe at day 1
        # should change the value later
      }
      
      Cummulative_H_vaccine <- c()
      Cummulative_H_vaccine[1] <- H_vaccine[1,1]
      H_week_vaccine <- c()
      H_obs_t_v <- c() # time of observation in hospitalized patients H*(t-1)
      H_obs_t_v[1] <- 0
      for (i in 2:(input$t_obs-1)) {
        H_week_vaccine[i] <- sum(H_vaccine[i,floor(i/4.33)+1])
        H_obs_t_v[i] <- H_obs_t_v[i-1]+H_week_vaccine[i]*(i-1)
        Cummulative_H_vaccine[i] <- Cummulative_H_vaccine[i-1]+H_week_vaccine[i]
      }
      
      for (i in 1:(input$t_obs-1)) {
        H_vaccine_survival[i] <- (num_births[input$time_enroll]-Cummulative_H_vaccine[i])/num_births[input$time_enroll] # choose to observe at day 1
        # should change the value later
      }
         
      # vaccine efficacy
      # VE by number of hospitalizations = 1-(H_vac/H_placebo)
      # 90 days approximately equals to 13 weeks
      VE_90_hosp= (1- (Cummulative_H_vaccine[13]/Cummulative_H_placebo[13]))*100
      # VE by person-time of observations = 1-(H_vac/person-time-obs_v)/(H_placebo/person-time-obs_p)  
      # person-time-obs_v = (tot_pop - cum_H)*t + H*(t-1)
      survive_obs_t_v <- (num_births[input$time_enroll]-Cummulative_H_vaccine[13])*13 #(tot_pop - cum_H)*t
      survive_obs_t_p <- (num_births[input$time_enroll]-Cummulative_H_placebo[13])*13
      
      VE_90_persontime=(1-(Cummulative_H_vaccine[13]/(H_obs_t_v[13]+survive_obs_t_v))/
        (Cummulative_H_placebo[13]/(H_obs_t_p[13]+survive_obs_t_p)))*100
      
      # plot survival curves between H_placebo and H_vaccine
      HospPlot_placebo_vaccine <- data.frame("H_placebo_survival"=c(1,H_placebo_survival),
                                  "H_vaccine_survival"=c(1,H_vaccine_survival),
                                  Date=dates)
      
      #adjust the margin
      margin <- list(pad=4,
                     autoexpand = T,
                     l = 100,
                     r = 105,
                     t = 105)
      
      
      fig <- plot_ly(HospPlot_placebo_vaccine, type = 'scatter', mode = 'lines')%>%
        add_trace(x = ~Date, y = ~H_placebo_survival, name = 'Placebo', line = list(color = '#66c2a5'))%>%
        add_trace(x = ~Date, y = ~H_vaccine_survival, name = 'Vaccine', line = list(color = '#fc8d62'))%>% 
        add_annotations(
          x= HospPlot_placebo_vaccine$Date[14],
          y= 0.98,
          xref = "x",
          yref = "y",
          text = paste0("90 days VE by case ratio=",round(VE_90_hosp,1),"%"),
          showarrow = F
        )%>% 
        add_annotations(
          x= HospPlot_placebo_vaccine$Date[14],
          y= 0.96,
          xref = "x",
          yref = "y",
          text = paste0("90 days VE by person-time at risk=",round(VE_90_persontime,1),'%'),
          showarrow = F
        ) %>% add_annotations(
          x= HospPlot_placebo_vaccine$Date[6],
          y= 0.98,
          xref = "x",
          yref = "y",
          text = paste0("vaccine\ngroup=",Cummulative_H_vaccine[13]),
          showarrow = F
        )%>% 
        add_annotations(
          x= HospPlot_placebo_vaccine$Date[6],
          y= 0.97,
          xref = "x",
          yref = "y",
          text = paste0("placebo\ngroup=",Cummulative_H_placebo[13]),
          showarrow = F
        )%>%
        add_segments(x = HospPlot_placebo_vaccine$Date[14], xend = HospPlot_placebo_vaccine$Date[14], y = 0.96, yend = 1)
      
         config(fig, toImageButtonOptions = list(format= 'png', # one of png, svg, jpeg, webp
                                              filename= 'custom_image',
                                              height= 450,
                                              width= 900,
                                              scale= 2.5 ))%>%
        layout(margin = margin,
               legend=list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = 0.5,y=1.2),
               yaxis = list(title = 'Survival probability against RSV hospitalizations',showgrid = FALSE,showline= T),
               xaxis = list(dtick = "M3",
                            ticklabelmode="period",showgrid = FALSE,showline= T))
    })
  
  output$Plot_hosp <- renderPlotly({
    library(xts)
    
    dates <- seq((as.Date("2024-02-27")-9*52*7+input$t_obs*7), length=9*52, by="weeks")
    
    H_burnin <- sim_mat_placebo_vaccine()$H_burnin 
    
    # total hospitalizations in children under 2 years of age
    Hosp_Plot <- data.frame("Children"=round(rowSums(H_burnin[,1:6])),
                                "Total"=round(rowSums(H_burnin)),
                                Date=dates)
    #adjust the margin
    margin <- list(pad=4,
                   autoexpand = T,
                   l = 100,
                   r = 105,
                   t = 105)
    
    
    fig <- plot_ly(Hosp_Plot, type = 'scatter', mode = 'lines')%>%
      add_trace(x = ~Date, y = ~Children, name = 'Hospitalizations in children', line = list(color = '#66c2a5'))%>%
      add_trace(x = ~Date, y = ~Total, name = 'Hospitalizations in total populations', line = list(color = '#fc8d62'))%>%
      layout(shapes = 
        list(type = "rect",
             fillcolor = "blue", line = list(color = "blue"), opacity = 0.1,
             x0 = "2024-02-27", x1 = (as.Date("2024-02-27")+90), xref = "x",
              y0 = 0, y1 = 1000, yref = "y"),
    annotations =list(
      x = "2024-02-27",
      y = 1000,
      text = "90 days of\nVE observational\nperiod",
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      ax = -150,
      ay = 30))
    
    # measuring at peak has the lowest VE
    # because this model potential assume leaky
    
      config(fig, toImageButtonOptions = list(format= 'png', # one of png, svg, jpeg, webp
                                            filename= 'custom_image',
                                            height= 450,
                                            width= 900,
                                            scale= 2.5 ))
    
  })
}


