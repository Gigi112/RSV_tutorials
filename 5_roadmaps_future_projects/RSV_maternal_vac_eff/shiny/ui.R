library("shiny")
library("deSolve")
library("cowplot")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")
library("shinyWidgets")
library("plotly")

ui <- fluidPage(
  titlePanel("Modeling maternal vaccination against RSV hospitalizations in infants"),
  hr(),
  p(div(HTML("Disclaimer: This simulation is for research and educational purposes only and is not intended to be a tool for epidemic prediction. 
  There are many uncertainties and debates about the details of RSV infection and transmission. This model has limitations on capturing the constantly changing non-pharmaceutical interventions against COVID-19.
             This work is licensed under a <a href=https://creativecommons.org/licenses/by-sa/4.0/> Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License </a>"))),
  sidebarLayout(
    sidebarPanel(
    #  p(div(HTML("<em>Evaluate the impact of various parameters...</em>"))),
      column(width=6,
                 sliderInput("baseline.txn.rate",
                         "Approximate R0:",
                         min = 9, max = 10.5, value = 9, step=0.1
             ),
             sliderInput("annual_birth_rate",
                         "Annual birth rate:",
                         min =  0.0080, max = 0.0160,value = 0.0094,  step=0.0002
             ),
             sliderInput("phase",
                         "Peak timing of epidemics:",
                         0, 2*pi, 3.62, step=0.01
             ),
             sliderInput("b1",
                         "Amplitude of transmission seasonality:",
                         0.1, 0.3, 0.2, step=0.01
             ),
             sliderInput("time_enroll",
                           "Time of enrollment:",
                           1, 12,1, step=1
             )),  
      column(width=6,
             sliderInput("DurationMatVac",
                         "Duration of protection of maternal vaccination:",
                         90, 180, 90, step=30
             ), 
             sliderInput("DurationNaturalProtection",
                         "Duration of protection of naturally acquired maternal immunity:",
                          20, 150, 30,step=1
             ), 
             sliderInput("t_obs",
                         "Period of observation (weeks):",
                         26, 52*2, 27, step=1
             )),
    width=4
  ),
    
    mainPanel(navbarPage("Output:",
                         
       tabPanel("Survival curves",
       fluidPage(
       fluidRow(
       h3("Simulated survival curves between vaccinated and placebo groups"),
       p("The plot shows the percentage of infants who have not yet been hospitalized"),
       plotlyOutput("Plot_placebo_vaccine"),
       br(),
       br() ))) ,
       tabPanel("Transmission dynamics",
                fluidPage(
                  fluidRow(
                    plotlyOutput("Plot_hosp"))))))
))
      