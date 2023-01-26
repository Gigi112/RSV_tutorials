### function to calculate the average age of hospitalization in children under 5

average_age_H <- function(age_midpoint,Hosp_under5,time.step){
  
  if(parms$time.step=='month'){
    period=12
  }else if(parms$time.step=='week'){
    period=52
  }
  
  #multiply age midpoint with hospitalization matrix
  age_H_matrix <- sweep(Hosp_under5, MARGIN=2, age_midpoint, '*')
  
  # average age of hospitalization by year 
  average_H <- sum(age_H_matrix[1:period,])/sum(Hosp_under5[1:period,])

  
  return(average_H)
}