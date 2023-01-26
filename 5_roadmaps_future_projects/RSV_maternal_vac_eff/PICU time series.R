library(readxl)
PICU_UK <- read_excel("~/Library/CloudStorage/Box-Box/immuno-epi/PICU data/PICU_UK.xlsx")
library(lubridate)

PICU_UK$month <- floor_date(PICU_UK$Date_PICU_admission, "month")
PICU_UK$RSV <- 1

library(dplyr)
RSV_series <- PICU_UK %>% 
  group_by(month) %>% 
  summarise(case=sum(RSV))
RSV_series <- RSV_series[order(RSV_series$month),]

plot(RSV_series$month,RSV_series$case,type="l")

PICU_NL <- read_excel("~/Library/CloudStorage/Box-Box/immuno-epi/PICU data/PICU2.xlsx")
PICU_NL$month <- floor_date(PICU_NL$Date_PICU_admission_PICU, "month")
PICU_NL$RSV <- 1

RSV_series_NL <- PICU_NL %>% 
  group_by(month) %>% 
  summarise(case=sum(RSV))
RSV_series_NL <- RSV_series_NL[order(RSV_series_NL$month),]

plot(RSV_series_NL$month,RSV_series_NL$case,type="l")

## PICU cases are a bit sparse

## PICU age distribution UK
PICU_UK$age_at_ICU <- floor(PICU_UK$Age_PICUadmission_days/30)
hist(PICU_UK$age_at_ICU,breaks = 13,probability = T)

## PICU age distribution NL
PICU_NL$age_at_ICU <- floor(PICU_NL$Age_PICUadmission_days_PICU/30)
hist(PICU_NL$age_at_ICU,breaks = 13,probability = T)

