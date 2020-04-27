#####################################################################################################
########## SECTION 1: EMPIRICAL R(t) ESTIMATES BY COUNTY ############################################
#####################################################################################################

#Authors: LTC Nick Clark, PhD, USMA D/MATH COL Matt Dabkowski, PhD, USMA D/SE, Dusy Turner
#Date: Created 23 MAR 2020, Updated 1 APR 2020
#Purpose: This script estimates the daily effective reproduction number (R(t)) of COVID-19 by US county. 

#Description: As described in Nishiura and Chowell (2009), R(t) is defined as the actual average # 
#of secondary cases per primary case at calendar time t. R(t) shows time-dependent variation due
#to the decline in susceptible individuals (intrinsic factors) and the implementation of control 
#measures (extrinsic factors). If R(t) < 1, it suggests that the epidemic is in decline and 
#may be regarded as being under control at time t (vice versa, if R(t) > 1).

#The key equation in Nishiura and Chowell (2009) is (33), which uses daily incidence data 
#(i.e. using j_i incident cases infected between time t_i and time t_(i+1) and discretized a 
#generation time distribution w_i. As per Nishiura and Chowell (2009), generation time is
#the time from infection of a primary case to the infection of a secondary case by the primary 
#case. Given the novelty of COVID-19, the generation time distribution is currently unknown. 

#With this in mind, Zhao et al.'s (2020) "Preliminary estimation of the basic reproduction number 
#of novel coronavirus..." assumed the probability distribution for the serial interval (SI) (the time 
#interval between infection and subsequent transmission) followed the MERS' SI distribution, which 
#was Gamma with mean=7.6 days and sd=3.4 days [as per (Assiri et al., 2013)]. Alternatively, Zhao et al.
#(2020) also posits a mean=8.4 days and sd=3.8 days, which corresponded to SARS [as per (Lipsitch et al., 
#2003)]. 

#Additionally, herd immunity (1-1/R0) is a function of estimates of R0. Thw papers were used to
#establish this range: "Estimation of the reproductive number of novel coronavirus..." (Zhang et al.),
#and "The reproductive number of COVID-19..." (Liu at al.). From these, we settled on a COVID-19 R0
#range of 2 to 3, implying herd immunity may be achieved at 50% to 67% of the population being exposed. 
#Furthermore, we used the ML estimate of 2.28 for R0 presented in (Zhang et al.) for our "best guess,"
#implying a "best guess" herd immunity level of 56.1%.

##############################################################################################################
#UPDATES: 
#1 APR 2020: Nishiura, Lintona, and Akhmetzhanov's (2020) "Serial interval of novel coronavirus (COVID-19) infections"
#provides a more extensive analysis of the generation time distribution. In their analysis, "the median
#serial interval of the best-fit Weibull distribution model was estimated at with a mean and SD of 4.8
#days and 2.3 days." This is important, as (1) "The difference between these distributions suggests that
#using serial intervals estimates from SARS data will result in overestimation of the COVID-19 basic reproduction 
#number," and (2) "the serial interval of COVID-19 is close to or shorter than its median incubation period. 
#This suggests that a substantial proportion of secondary transmission may occur prior to illness onset."


#Outputs: County-level R(t) values by day

library(tidyverse)
library(extraDistr)
library(mixdist) #used to recoved the Weibull parameters from the mean and sd

jhu_all<-read_csv(paste0("https://raw.githubusercontent.com/nick3703/Parametric-Modeling-for-Time-Varying-Reproducibility-Number/master/",lubridate::today(),"_JHU.csv"))




#This function accepts the cumulative confirmed COVID-19 cases by county by day.
#It then estimates R(t) for that county by day.
rt.func.v3_dusty<-function(dat,mean.Weibull=4.8,sd.Weibull=2.3){
  r.vals<-c()
  confirmed = dat$confirmed
  #get the Weibull parameters from mixdist's weibullpar function
  mGT.params<-weibullpar(mean.Weibull, sd.Weibull, loc = 0)
  alpha<-mGT.params[2] # called shape in weibullpar, alpha in a discrete Weilbull
  beta<-mGT.params[1] # called scale in weibullpar, beta in a discrete Weibull
  #the extraDistr package uses an altrnative parameterization of the Weibull (q, beta) from
  #Nakagawa and Osaki (1975) where q = exp(-alpha^-beta), so...
  q<-exp(-as.numeric(alpha)^(-as.numeric(beta)))
  #Discretize Weibull via the extraDistr package's ddweibull function
  w<- ddweibull(0:1000, as.numeric(q), as.numeric(beta), log = FALSE)
  df<-data.frame(counts=t(c(0,confirmed)))
  # df<-data.frame(counts=t(dat[-(1:4)]))
  # first.case.index<-min(which(df>0))
  # total.cases<-max(df)
  growth<-diff(unlist(df))
  # growth<-diff(df$counts)
  growth<-pmax(growth, 0) # eliminate any erroneous downward shifts in the cumulative counts
  #Estimate R(t) from equation (33) of Nishiura and Chowell (2009)
  for(k in 2:length(growth)){
    r.vals[k-1]<-growth[k]/(sum(growth[1:k]*rev(w[1:k])))
  }
  r_val_df <- tibble(
    day = 1:length(r.vals), 
    r.vals = r.vals, 
    daily_new_cases = growth[-1])
    #population = dat$population[1], 
    #locked_down = dat$locked_down[-1])
  #Output the results
  #return(r.vals)
  return(r_val_df)
}

############
#Try it out# 
############

daily <-
  jhu_all %>% 
  #load_data_files()$jhu_all %>%
  #left_join(county_pop) %>% 
  #left_join(stay_at_home) %>% 
  #mutate(locked_down = if_else(date>effective_date,1,0)) %>% 
  filter(!is.na(country_region), county_name!="Unassigned") %>%
  group_by(fips) %>% 
  arrange(fips,date) %>% 
  # mutate(daily_confirmed = confirmed-lag(confirmed))  %>% 
  # mutate(daily_confirmed = replace_na(daily_confirmed,confirmed[1])) %>%
  # mutate(daily_confirmed = if_else(daily_confirmed<0,0,daily_confirmed))  %>% 
  filter(!(confirmed==0 & row_number()==1)) %>% 
  filter(confirmed!=0) %>%
  # filter(daily_confirmed!=0) %>% ## dirty data
  group_by(date,fips) %>% filter(confirmed==max(confirmed)) %>% 
  group_by(fips) %>%
  filter(n()>1) %>% 
  filter(country_region == "US") 
  # group_split()

rt_est <- 
  daily %>% ## cases confirmed by fips
  group_by(fips) %>% 
  nest() %>%  ## one row per fips with nested df of confirmed cases
  mutate(ans = map(data, rt.func.v3_dusty)) %>% ## applies the rt function to each dataframe
  unnest(ans) ## expands back to one dataframe

rt_wide<-rt_est%>%pivot_wider(id_cols=fips,names_from=day,values_from=r.vals)

write.csv(rt_wide,paste0(lubridate::today(),"ByCounty.csv"))
