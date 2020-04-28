# Parametric-Modeling-for-Time-Varying-Reproducibility-Number

Purpose: This script estimates the daily effective reproduction number (R(t)) of COVID-19 by US county. 

Description: As described in Nishiura and Chowell (2009), R(t) is defined as the actual average # 
of secondary cases per primary case at calendar time t. R(t) shows time-dependent variation due
to the decline in susceptible individuals (intrinsic factors) and the implementation of control 
measures (extrinsic factors). If R(t) < 1, it suggests that the epidemic is in decline and 
may be regarded as being under control at time t (vice versa, if R(t) > 1).

The key equation in Nishiura and Chowell (2009) is (33), which uses daily incidence data (i.e. using j_i incident cases infected between time t_i and time t_(i+1) and discretized a 
generation time distribution w_i. As per Nishiura and Chowell (2009), generation time is
the time from infection of a primary case to the infection of a secondary case by the primary 
case. Given the novelty of COVID-19, the generation time distribution is currently unknown. 

With this in mind, Zhao et al.'s (2020) "Preliminary estimation of the basic reproduction number 
of novel coronavirus..." assumed the probability distribution for the serial interval (SI) (the time 
interval between infection and subsequent transmission) followed the MERS' SI distribution, which 
was Gamma with mean=7.6 days and sd=3.4 days [as per (Assiri et al., 2013)]. Alternatively, Zhao et al.(2020) also posits a mean=8.4 days and sd=3.8 days, which corresponded to SARS [as per (Lipsitch et al., 
2003)]. 

Nishiura, Lintona, and Akhmetzhanov's (2020) "Serial interval of novel coronavirus (COVID-19) infections"
provides a more extensive analysis of the generation time distribution. In their analysis, "the median
serial interval of the best-fit Weibull distribution model was estimated at with a mean and SD of 4.8
days and 2.3 days." This is important, as (1) "The difference between these distributions suggests that
using serial intervals estimates from SARS data will result in overestimation of the COVID-19 basic reproduction 
number," and (2) "the serial interval of COVID-19 is close to or shorter than its median incubation period. 
This suggests that a substantial proportion of secondary transmission may occur prior to illness onset."

## Update 23 April

Added Rt by County which calculates the empirical Rt value for each county in America.  Please note that county level estimates for rural counties are unlikely to be valid due to small sample sizes.

## Update 27 April

Naming convention has stabilized.  In order to pull daily run:

### JHU Data
read_csv(paste0("https://raw.githubusercontent.com/nick3703/Parametric-Modeling-for-Time-Varying-Reproducibility-Number/master/",lubridate::today(),"_JHU.csv"))

### Daily State R(t)
read_csv(paste0("https://raw.githubusercontent.com/nick3703/Parametric-Modeling-for-Time-Varying-Reproducibility-Number/master/",lubridate::today(),"ByState.csv"))

### Daily County R(t)

read_csv(paste0("https://raw.githubusercontent.com/nick3703/Parametric-Modeling-for-Time-Varying-Reproducibility-Number/master/",lubridate::today(),"ByCounty.csv"))

