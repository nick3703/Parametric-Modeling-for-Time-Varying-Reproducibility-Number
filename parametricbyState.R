#Authors: LTC Nick Clark, PhD, USMA D/MATH COL Matt Dabkowski, PhD, USMA D/SE, Rodney X. Sturdivant, MAJ Dusty Turner
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


##############################################################################################################
#UPDATES: 
#1 APR 2020: Nishiura, Lintona, and Akhmetzhanov's (2020) "Serial interval of novel coronavirus (COVID-19) infections"
#provides a more extensive analysis of the generation time distribution. In their analysis, "the median
#serial interval of the best-fit Weibull distribution model was estimated at with a mean and SD of 4.8
#days and 2.3 days." This is important, as (1) "The difference between these distributions suggests that
#using serial intervals estimates from SARS data will result in overestimation of the COVID-19 basic reproduction 
#number," and (2) "the serial interval of COVID-19 is close to or shorter than its median incubation period. 
#This suggests that a substantial proportion of secondary transmission may occur prior to illness onset."

#Outputs: State-level R(t) values by day
install.packages('tidyverse')
install.packages('extraDistr')
install.packages('mixdist')
intall.packages('ciTools')
library(tidyverse)
library(extraDistr)
library(mixdist) #used to recoved the Weibull parameters from the mean and sd
library(ciTools)

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

#Read in the daily incidence data from usafacts.org and collect the dates
#us.data<-read_csv("https://usafactsstatic.blob.core.windows.net/public/data/covid-19/covid_confirmed_usafacts.csv")
jhu_all<-read_csv(paste0("https://raw.githubusercontent.com/nick3703/Parametric-Modeling-for-Time-Varying-Reproducibility-Number/master/",lubridate::today(),"_JHU.csv"))
#test.data<-us.data%>%filter(stateFIPS > 0)
#county_pop<-read_csv("county_pop.csv")

jhu_all<-jhu_all%>%mutate(date=as.Date(date,format="%m/%d/%Y"))

jhu_state<-jhu_all %>% group_by(province_state,date)%>%
  summarize(confirmed=sum(confirmed))


rt_est <- 
  jhu_state %>% ## cases confirmed by fips
  group_by(province_state) %>% 
  nest() %>%  ## one row per fips with nested df of confirmed cases
  mutate(ans = map(data, rt.func.v3_dusty)) %>% ## applies the rt function to each dataframe
  unnest(ans) ## expands back to one dataframe


state.names<-read_csv("https://raw.githubusercontent.com/nick3703/Parametric-Modeling-for-Time-Varying-Reproducibility-Number/master/StateNames.csv")

rt_est<-rt_est%>%
  filter(province_state %in%state.names$State)

list_states<-unique(rt_est$province_state)

output<-data.frame(state=NA,b0=NA,b1=NA,ucb=NA,curr.r=NA,
                   curr.emp.r=NA,latest.cases=NA)


p<-list()
for(j in 1:length(list_states)){
  state<-rt_est%>%filter(province_state==list_states[j])
  output[j,]$state=list_states[j]
  # Fit models with Gamma errors
  gam.glm <- glm(r.vals+.1~day, family = Gamma(link = "log"),data=state)
  output[j,]$latest.cases=state[nrow(state),]$daily_new_cases
  #summary(gam.glm)    # summary, estimates
  coefs<-gam.glm$coefficients
#gam.glm$fitted.values
  output[j,]$b0=coefs[1]
  output[j,]$b1=coefs[2]
  df = data.frame(day = state$day, r.vals = state$r.vals)   # Note can add new data here if desired
  df1 <- df %>%
    add_ci(gam.glm, names = c("lwr","upr")) %>%
    add_pi(gam.glm, names = c("plwr","pupr"))%>%
    mutate(type = "parametric")
  output[j,]$ucb=df1[nrow(df1),]$upr
  output[j,]$curr.r<-df1[nrow(df1),]$pred
  output[j,]$curr.emp.r<-df1[nrow(df1),]$r.vals


output<-output%>%mutate(current.rt=exp(b0+b1*max(rt_est$day)))


# Plot residuals
df.res = data.frame(x = state$day, y = gam.glm$residuals)
ggplot(df.res, aes(x = x, y = y) )+
  geom_point(aes(x=x,y=y))+
  geom_hline(yintercept = 0)+
  labs(x="Time",y="Residual")

# Plot predicted, confidence/prediction intervals
df = data.frame(day = state$day, r.vals = state$r.vals)   # Note can add new data here if desired
df1 <- df %>%
  add_ci(gam.glm, names = c("lwr","upr")) %>%
  add_pi(gam.glm, names = c("plwr","pupr"))%>%
  mutate(type = "parametric")

df1%>%ggplot(aes(x=day,y=r.vals))+
  geom_point()+geom_line(aes(x=day,y=pred),lwd=2,color="red")+
  theme_bw()

write.csv(output,paste0(lubridate::today(),"ByState.csv"))

p[[j]]<-ggplot(df1, aes(x = day, y = pred))+
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .4) +
  #geom_ribbon(aes(ymin = plwr, ymax = pupr), alpha = 0.2) +
  #geom_line(aes(x=day, y=pred)) +
  #facet_wrap(~province_state) +
  labs(title = paste0("Model fit- ",list_states[j],": JHU New Cases - ",round(output[j,]$latest.cases,4)))+
  ylim(c(0,5))+
  geom_hline(yintercept=1,color="red")+
  theme_bw()+theme(
    plot.title = element_text(color="black", size=14, face="bold"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")) +
      xlab("Days Since March 22nd") + ylab("R(t) Estimate")

}

pdf(paste0(lubridate::today(),"plots.pdf"))
for (i in 1:53) {
  print(p[[i]])
}
dev.off()
