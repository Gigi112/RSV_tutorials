# RSV tutorials
This Github repository contains tutorial materials for Pitzer-Weinberger's lab transitioning, mostly on work related to RSV.
\
\
We introduced several topics of interest. In each directory, you will find:

 1. an R project. I separated each topic as an independent R project.

 2. a PDF file, which has the same name as the directory. Read this file first to get an idea of the outline and prerequisite readings/materials of the topics. I deposited the prerequisite readings/materials either by links or PDFs in the same directory. Once you click the link in the tutorial PDF, you will be redirected to the prerequisite readings.

 3. a RMD file. This file is the source file for the PDF, you can read and run the codes for each chunk. **Note: you may need to combine several code chunks to have the correct version of code** because I separated them for the sake of creating a PDF. 

 4. Other codes and data files that are necessary to run the RMD file. 

 5. Dan and Ginny's lecture notes, which present the theoretical knowledge of the tutorial materials. 

The order of the tutorials was listed as follows:  

- Directory 0_HCUP_data_clean. This directory introduces how to convert the State Inpatient/Emergency Department Databases from ACS files to parquet files. With parquet data files, we can easily conduct data analysis using dplyr type of syntax in r. 

- Directory 1_MSIS_model_explain. This directory describes the transmission dynamic model of RSV. It includes background knowledge of transmission dynamic model of RSV (a kind of mathematical model), simulation of the transmission dynamic model, maximum a posterior/maximum likelihood estimation, STAN programming language for Bayesian analysis of ordinary differential equations. It also introduces how to create contact matrices based on the user-defined age groups and how to run R script on HPC (clusters).

- Directory 2_RSV_timing_explain. This directory demonstrates how to utilize harmonic regression to identify the peak timing of seasonal epidemics and to use the second-derivative method to pinpoint the onset timing of irregular RSV epidemics following the relaxation of mitigation measures against COVID-19 pandemic.

- Directory 3_JAGS_Bayesian_model. This directory explains how to estimate parameters in Bayesian framework using JAGS (Just Another Gibbs Sampler) model. It first gives an simple example then expands to censored data and spatial analysis. It also talks about how to create a neighboring matrix in this directory.

- Directory 4_RSV_disease_burden. I have not finished working on this but it will focus on using Hierarchical Bayesian regression to estimate the "true" burden of RSV hospitalizations in older adults.

- Directory 5_roadmaps_future_projects. Road maps for two projects that I have just started: (1) cost-effectiveness of  vaccines against RSV hospitalizations in older adults in the United States and (2) factors affecting maternal RSV vaccine efficacy in clinical trials.

