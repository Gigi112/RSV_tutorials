library(rstan)
library(gridExtra)
rstan_options (auto_write = TRUE)

setwd("~/project/Rstan/WA")

data_msis<- readRDS("data_msis_wa_annual.rds")

print(paste("number of cores:", future::availableCores()))

model <- stan_model("age-structureMSIS_deterministic.stan")

init_fun <- function() list(beta=3.17,b1=0.198,phi=1.46,omega=0.61,hosp_ratio=0.76)

fit_msis_poisson <- sampling(model,
                             data = data_msis,
                             iter = 2000,
			     warmup = 1000,
                             chains = 4, 
                             seed = 112,sample_file = file.path(getwd(), 'WA_MSIS_12_deterministic.csv'),
                             cores=4,
                             open_progress=T,
			     init= init_fun,
                	     control = list(stepsize=0.1))