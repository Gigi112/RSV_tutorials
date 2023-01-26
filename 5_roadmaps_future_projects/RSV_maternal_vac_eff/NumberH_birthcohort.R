# function to calculate hospitalizations in a birth cohort to 90 days
# number of hospitalizations depends on the season of observation
Hosp_90 <- function(t,H1){ 
  Hosp <- round(H1[t,1]+H1[t+1,2]+H1[t+2,3])
  return(Hosp)
  }

# function to caculate hospitalizations in a birth cohort to 180 days
Hosp_180 <- function(t,H1){ 
  Hosp <- round(H1[t,1]+H1[t+1,2]+H1[t+2,3]+H1[t+3,4]+H1[t+4,5]+H1[t+5,6])
  return(Hosp)
}