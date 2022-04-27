#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Group 2: RSV
# Members: Krista Stenberg, Ally Dalby, Julia Spychalski, Seibi Kobara, and Aaron Holton

# Readme
# Codes for the primary objective
# step 1: run functions
# step 2: results (Rn, epi curve, averted)

# Codes for the secondary objective
# (need to run step 1 and 2 above)
# step 3: run functions
# step 4: simulation (takes long)
# step 5: visualize
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## -------------------------------------------------------------------------------
# packages
options(scipen = 999)
library(EpiModel)
library(tidyverse)
library(magrittr)
library(gplots)
library(egg)


#---------------------------------------------------------------------------------
# Primary objective
# step1: functions
## -------------------------------------------------------------------------------
# Hospitalization rate calculation
# in 2012 rate/100,000 from goldstein2018
age=c("<1","1","2","3","4","5","6","7","8","9","10","11","12-17","18-49","50-64","65-")
hosp_count= c(28692,6757,2801,1364,661,320,200,88,76,34,48,40,149,465,564,1174)
rate_=c(1334,311,128.4,61.4,29.2,14.1,8.9,3.9,3.4,1.5,2.1,1.7,1.1,0.6,1.7,5.3)
data=data.frame(age,hosp_count,rate=rate_/100000)

# yearly
data2 = data.frame()
for (i in 1:16){
    if (i <= 12){
        data2[i,1]= data[i,1]
        data2[i,2]= data[i,2]
        data2[i,3]= data[i,3]
        }
    else if (i ==13){
        for (j in 13:18){
            data2[j,1]= j-1
            data2[j,2]= data[i,2]/6
            data2[j,3]= data[i,3]/6
        }
    }
    else if (i == 14){
        for (j in 19:50){
            data2[j,1]= j-1
            data2[j,2]= data[i,2]/32
            data2[j,3]= data[i,3]/32
        }
    }
    else if (i == 15){
        for (j in 51:65){
            data2[j,1]= j-1
            data2[j,2]= data[i,2]/15
            data2[j,3]= data[i,3]/15
        }
    }
    else if (i == 16){
        for (j in 66){
            data2[j,1]= "65-"
            data2[j,2]= data[i,2]
            data2[j,3]= data[i,3]
        }
    }
}
names(data2)=c("age","hosp_count","rate")

data2 %<>% mutate(case_count=hosp_count/rate)

# collapse to fit our age categories
data3 = data.frame()
for (i in 1:7){
    if (i==1){
        data3[i,1]= "0-5m"
        data3[i,2]= data2[1,2]/2
        data3[i,3]= data2[1,4]/2
    } 
    if (i==2){
        data3[i,1]= "6-12m"
        data3[i,2]= data2[1,2]/2
        data3[i,3]= data2[1,4]/2
    } 
    if (i==3) {
        data3[i,1]= "1-4"
        data3[i,2]= sum(data2[2:5,2])
        data3[i,3]= sum(data2[2:5,4])
    }
    if (i==4) {
        data3[i,1]= "5-19"
        data3[i,2]= sum(data2[6:20,2])
        data3[i,3]= sum(data2[6:20,4])
    }
    if (i==5) {
        data3[i,1]= "20-39"
        data3[i,2]= sum(data2[21:40,2])
        data3[i,3]= sum(data2[21:40,4])
    }
    if (i==6) {
        data3[i,1]= "40-59"
        data3[i,2]= sum(data2[41:60,2])
        data3[i,3]= sum(data2[41:60,4])
    }
    if (i==7) {
        data3[i,1]= "60-"
        data3[i,2]= sum(data2[61:nrow(data2),2])
        data3[i,3]= sum(data2[61:nrow(data2),4])
    }
}
names(data3) =c("age","hosp_count","case_count")

data3 %<>% mutate(rate=hosp_count/case_count)
data3 %<>% mutate(case_non_hospitalized = case_count - hosp_count)
demog = data3

# function for a dynamic model
estimates=function(omega_,psi_){
   
# parameter settings
#omega_ = 0  # vaccine coverage
#psi_   = 0   # vaccine efficacy (reduction in transmission and hospitalization)
R0_    = 3
gamma_ = 0.1  # infectious length
rho_   = 1/(365)    # immunity length
US_pop_= 332578200

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# age-contact matrix M

contact=read.csv("contact_matrix.csv")
M=contact[,-c(1)]
names(M)=c("0-5 months","6-12 months", 
           "1-4 yrs","5-19 yrs","20-39 yrs","40-59 yrs","60+ yrs")
rownames(M)=c("0-5 months","6-12 months", 
              "1-4 yrs","5-19 yrs","20-39 yrs","40-59 yrs","60+ yrs")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mortality rate calculation

# annual mortality in US 2019
# 1-4 : 23.3 per 100,000
# 5-14: 13.4 per 100,000
# 25-34: 128.8 
# 35-44: 199.2
# 45-54: 392.4
# 55-64: 883.3
# 65-74: 1764.6
# 75-84:4308.3
# 85 +: 13228.6

# for 0-5month: 
death06_= 23.3/(4*2*100000*365) # per 1 person daily assuming that the mortality is equal to 1 year old
# for 6-12months
death1_ = 23.3/(4*100000*365)
# for 1 to 4 years
death5_ = 23.3/(100000*365)
# for 5 to 19 years
death10_ = 69.7*(2/3)/(100000*365)
# for 20 to 39
death30_ = (69.7*(1/3) + 128.8 + 199.2*(1/2))/(100000*365)
# for 40 to 59
death40_ = (199.2*(1/2) + 392.4 + 883.3*(1/2))/(100000*365)
# for 60 over
death60_ = (888.3*(1/2) + 1764.6 + 4308.3 + 13228.6)/(100000*365)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# birthrate 11.0 per 1000 annually
birthrate_ = 11/(1000*365)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# beta (coefficient of effec contact)

eiganvalue=5.822557704881415
# beta for non vaccinated
beta_ = R0_ * gamma_/eiganvalue
# beta for vaccinated
beta.v_ = R0_ *(1-psi_)* gamma_/eiganvalue

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model
# description
# compartment
# v06 : vaccinated 0-6months
# h06: infection from vaccinated group
# s06 : susceptible (unvaccinated)
# i06 : infection from unvaccinated group
# r06 : recovered from infection

# age stages
# 0-5 months : 06
# 6-12 months: 1
# 1-4 years  : 5
# 5-19 years : 10
# 20-39 years: 30
# 40-59 years: 40
# 60 years + : 60

RSV = function(t,t0,params){
    with(as.list(c(t0,params)),{
        # total population
        num=v06.num + s06.num + h06.num + nh06.num + r06.num +
            v1.num  + s1.num  +  h1.num +  nh1.num +  r1.num  + 
                      s5.num  +  h5.num +  nh5.num +  r5.num  +
                      s10.num +  h10.num +  nh10.num +  r10.num  +
                      s30.num +  h30.num +  nh30.num +  r30.num  +
                      s40.num +  h40.num +  nh40.num +  r40.num  +
                      s60.num +  h60.num +  nh60.num +  r60.num
        
        # each population
        num06 = v06.num + s06.num +  h06.num +  nh06.num +  r06.num
        num1  = v1.num  + s1.num  +  h1.num  +  nh1.num  +  r1.num
        num5  =           s5.num  +  h5.num  +  nh5.num  +  r5.num
        num10 =           s10.num +  h10.num +  nh10.num +  r10.num 
        num30 =           s30.num +  h30.num +  nh30.num +  r30.num  
        num40 =           s40.num +  h40.num +  nh40.num +  r40.num
        num60 =           s60.num +  h60.num +  nh60.num +  r60.num
       
         # define death rate
         death06 = death06_
         death1  = death1_
         death5  = death5_
         death10 = death10_                   
         death30 = death30_
         death40 = death40_
         death60 = death60_
         
         # hospitalization risk
         # potential mechanism is vaccine reduce the rate of hospitalization. Although vaccine               effect ends at 1 year old, the hospitalization rate among those who are older than                1 year is smaller than younger. This expects that vaccine prevent a the overall                   hospitalization rate/hospitalization averted.
         # hosp 0 -12 months,1-4 years were derived from Shi2019, and 5-60+years were from                   Goldstein2018 (see supplemental, the bottom of this file)
        
         hosp06  = 26.3/1000 # shi2017
         hosp1  =  11.3/1000 # shi2017
         hosp5  =  1.4/1000  # shi2017
         hosp06v = hosp06*(1-psi_)
         hosp1v =  hosp1*(1-psi_)

         hosp10  = demog[4,4] # goldstein2018
         hosp30  = demog[5,4] # goldstein2018
         hosp40  = demog[6,4] # goldstein2018
         hosp60  = demog[7,4] # goldstein2018
        
        # birth rate
         mu=birthrate_
        
        # differential equations
        # age 0-5months
        force_inf_v06 =      beta.v/num *(M[1,1]*(nh06.num) + 
                                          M[1,2]*(nh1.num)  +  
                                          M[1,3]*(nh5.num) + 
                                          M[1,4]*(nh10.num)+
                                          M[1,5]*(nh30.num)+
                                          M[1,6]*(nh40.num)+
                                          M[1,7]*(nh60.num))
        force_inf_06  =      beta  /num *(M[1,1]*(nh06.num) + 
                                          M[1,2]*(nh1.num)  +  
                                          M[1,3]*(nh5.num) + 
                                          M[1,4]*(nh10.num)+
                                          M[1,5]*(nh30.num)+
                                          M[1,6]*(nh40.num)+
                                          M[1,7]*(nh60.num))
        
        dV06  = num*mu*omega - force_inf_v06*v06.num - age06*v06.num - death06*v06.num
        dS06  = num*mu*(1-omega) - force_inf_06*s06.num - age06*s06.num - 
          death06*s06.num + rho*r06.num
        dH06  =  force_inf_v06*v06.num*hosp06v + force_inf_06*s06.num*hosp06 - age06*h06.num -
          death06*h06.num -gamma*h06.num
        dNH06 = force_inf_v06*v06.num*(1-hosp06v) + force_inf_06*s06.num*(1-hosp06) -
           age06*nh06.num - death06*nh06.num - gamma*nh06.num
        dR06  = gamma*h06.num + gamma*nh06.num - age06*r06.num - death06*r06.num - rho*r06.num  
        
        # age 6-12 months
        force_inf_v1 =      beta.v/num * (M[2,1]*(nh06.num) + 
                                          M[2,2]*(nh1.num)  +  
                                          M[2,3]*(nh5.num) + 
                                          M[2,4]*(nh10.num)+
                                          M[2,5]*(nh30.num)+
                                          M[2,6]*(nh40.num)+
                                          M[2,7]*(nh60.num))
        force_inf_1  =      beta  /num  *(M[2,1]*(nh06.num) + 
                                          M[2,2]*(nh1.num)  +  
                                          M[2,3]*(nh5.num) + 
                                          M[2,4]*(nh10.num)+
                                          M[2,5]*(nh30.num)+
                                          M[2,6]*(nh40.num)+
                                          M[2,7]*(nh60.num))
        
        dV1  = age06*v06.num - force_inf_v1*v1.num - age1*v1.num - death1*v1.num
        dS1  = age06*s06.num - force_inf_1*s1.num - age1*s1.num - 
          death1*s1.num + rho*r1.num
        dH1  = age06*h06.num + force_inf_v1*v1.num*hosp1v + force_inf_1*s1.num*hosp1 - 
          age1*h1.num -death1*h1.num -gamma*h1.num
        dNH1 = age06*nh06.num + force_inf_v1*v1.num*(1-hosp1v) + force_inf_1*s1.num*(1-hosp1) -
           age1*nh1.num - death1*nh1.num - gamma*nh1.num
        dR1  = age06*r06.num + gamma*h1.num + gamma*nh1.num - age1*r1.num - death1*r1.num - rho*r1.num 
        
        # age 1-4 years
        force_inf_5  =      beta  /num  *(M[3,1]*(nh06.num)+ 
                                          M[3,2]*(nh1.num)+  
                                          M[3,3]*(nh5.num) + 
                                          M[3,4]*(nh10.num)+
                                          M[3,5]*(nh30.num)+
                                          M[3,6]*(nh40.num)+
                                          M[3,7]*(nh60.num))
        
        dS5  = age1*v1.num + age1*s1.num + rho*r5.num - force_inf_5*s5.num - age5*s5.num - 
          death5*s5.num 
        dH5  = age1*h1.num + force_inf_5*s5.num*hosp5 - age5*h5.num - death5*h5.num -gamma*h5.num
        dNH5 = age1*nh1.num + force_inf_5*s5.num*(1-hosp5) - age5*nh5.num - 
          death5*nh5.num - gamma*nh5.num
        dR5  = age1*r1.num + gamma*h5.num + gamma*nh5.num - age5*r5.num - 
          death5*r5.num - rho*r5.num
        
        # age 5-19 years
        force_inf_10  =      beta  /num *(M[4,1]*(nh06.num) + 
                                          M[4,2]*(nh1.num)  +  
                                          M[4,3]*(nh5.num) + 
                                          M[4,4]*(nh10.num)+
                                          M[4,5]*(nh30.num)+
                                          M[4,6]*(nh40.num)+
                                          M[4,7]*(nh60.num))
        
        dS10  = age5*s5.num + rho*r10.num - force_inf_10*s10.num - 
          age10*s10.num - death10*s10.num 
        dH10  = age5*h5.num + force_inf_10*s10.num*hosp10 - 
          age10*h10.num - death10*h10.num -gamma*h10.num
        dNH10 = age5*nh5.num + force_inf_10*s10.num*(1-hosp10) - 
          age10*nh10.num - death10*nh10.num - gamma*nh10.num
        dR10  = age5*r5.num + gamma*h10.num + gamma*nh10.num - 
          age10*r10.num - death10*r10.num - rho*r10.num
        
        # age 20-39 years
        force_inf_30  =      beta  /num *(M[5,1]*(nh06.num) + 
                                          M[5,2]*(nh1.num)  +  
                                          M[5,3]*(nh5.num) + 
                                          M[5,4]*(nh10.num)+
                                          M[5,5]*(nh30.num)+
                                          M[5,6]*(nh40.num)+
                                          M[5,7]*(nh60.num))
        
        dS30  = age10*s10.num + rho*r30.num - force_inf_30*s30.num - 
          age30*s30.num - death30*s30.num 
        dH30  = age10*h10.num + force_inf_30*s30.num*hosp30 - 
          age30*h30.num - death30*h30.num -gamma*h30.num
        dNH30 = age10*nh10.num + force_inf_30*s30.num*(1-hosp30) - 
          age30*nh30.num - death30*nh30.num - gamma*nh30.num
        dR30  = age10*r10.num + gamma*h30.num + gamma*nh30.num - 
          age30*r30.num - death30*r30.num - rho*r30.num
        
        
        # age 40-59 years
        force_inf_40  =      beta  /num *(M[6,1]*(nh06.num) + 
                                          M[6,2]*(nh1.num)  +  
                                          M[6,3]*(nh5.num) + 
                                          M[6,4]*(nh10.num)+
                                          M[6,5]*(nh30.num)+
                                          M[6,6]*(nh40.num)+
                                          M[6,7]*(nh60.num))
        
        dS40  = age30*s30.num + rho*r40.num - force_inf_40*s40.num - 
          age40*s40.num - death40*s40.num 
        dH40  = age30*h30.num + force_inf_40*s40.num*hosp40 - 
          age40*h40.num - death40*h40.num -gamma*h40.num
        dNH40 = age30*nh30.num + force_inf_40*s40.num*(1-hosp40) - 
          age40*nh40.num - death40*nh40.num - gamma*nh40.num
        dR40  = age30*r30.num + gamma*h40.num + gamma*nh40.num - 
          age40*r40.num - death40*r40.num - rho*r40.num
        
        # age 60 + years
        force_inf_60  =      beta  /num *(M[7,1]*(nh06.num) + 
                                          M[7,2]*(nh1.num)  +  
                                          M[7,3]*(nh5.num) + 
                                          M[7,4]*(nh10.num)+
                                          M[7,5]*(nh30.num)+
                                          M[7,6]*(nh40.num)+
                                          M[7,7]*(nh60.num))
        
        dS60  = age40*s40.num + rho*r60.num - force_inf_60*s60.num - death60*s60.num 
        dH60  = age40*h40.num + force_inf_60*s60.num*hosp60 - death60*h60.num -gamma*h60.num
        dNH60 = age40*nh40.num + force_inf_60*s60.num*(1-hosp60) - 
          death60*nh60.num - gamma*nh60.num
        dR60  = age40*r40.num + gamma*h60.num + gamma*nh60.num - 
          death60*r60.num - rho*r60.num
        
        # output
        list(c(
            # compartment
            dV06,
            dS06,
            dH06,
            dNH06,
            dR06,
            
            dV1,
            dS1,
            dH1,
            dNH1,
            dR1,
            
            dS5,
            dH5,
            dNH5,
            dR5,
            
            dS10,
            dH10,
            dNH10,
            dR10,
            
            dS30,
            dH30,
            dNH30,
            dR30,
            
            dS40,
            dH40,
            dNH40,
            dR40,
            
            dS60,
            dH60,
            dNH60,
            dR60,
            
            # incidence
      si.flow.06_hos     = force_inf_v06*v06.num*hosp06v + force_inf_06*s06.num*hosp06,
      si.flow.06_non_hos = force_inf_v06*v06.num*(1-hosp06v) + force_inf_06*s06.num*(1-hosp06),
            
      si.flow.1_hos      = force_inf_v1*v1.num*hosp1v + force_inf_1*s1.num*hosp1,
      si.flow.1_non_hos  = force_inf_v1*v1.num*(1-hosp1v) + force_inf_1*s1.num*(1-hosp1)
            ))
    })
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
params=param.dcm(beta=beta_, 
                 beta.v=beta.v_,
                 age06 = 1/(365*0.5),
                 age1  = 1/(365*0.5),
                 age5  = 1/(365*4),
                 age10 = 1/(365*25),
                 age30 = 1/(365*20),
                 age40 = 1/(365*20),
                 
                 omega=omega_,
                 psi  =psi_,    
                 gamma=gamma_,
                 rho= rho_
)

# RSV cases, hospitalization cases were adjusted to the US population. The original data captured RSV cases in the State Inpatients Databases of the Healthcare Cost and Utilization Project (HCUP), which is accounting for 54% of the US population. It is unfair to calculate the RSV cases in the US by 54% weight, but it is justified because there is no other resources that have report RSV cases across age groups in the US population.

# In addition, the data from HCUP is annual cases and hospitalization. It is also unfair to divide 365 to obtain case at time 0. But, since there is no other sources that provide trend data telling the reasonable initial cases, it is justified to use these values.

demog %<>% mutate(hosp_US_time0=hosp_count/0.54 /365)
demog %<>% mutate(non_hosp_cases_US_time0=case_non_hospitalized/0.54/365)

init=init.dcm(v06.num =0, s06.num =US_pop_*0.6/100, 
                          h06.num =demog[1,6], nh06.num =demog[1,7], r06.num =0,
              v1.num  =0, s1.num  =US_pop_*0.6/100, 
                          h1.num  =demog[2,6], nh1.num  =demog[2,7], r1.num  =0,
                          s5.num  =US_pop_*4.8/100, 
                          h5.num  =demog[3,6], nh5.num  =demog[3,7], r5.num  =0,
                          s10.num =US_pop_*18.9/100,
                          h10.num =demog[4,6],  nh10.num =demog[4,7], r10.num =0,
                          s30.num =US_pop_*27.2/100,
                          h30.num =demog[5,6], nh30.num =demog[5,7], r30.num =0,
                          s40.num =US_pop_*25.2/100,
                          h40.num =demog[6,6], nh40.num =demog[6,7], r40.num =0,
                          s60.num =US_pop_*22.7/100,
                          h60.num =demog[7,6], nh60.num =demog[7,7], r60.num =0,
              
              si.flow.06_hos=0,
              si.flow.06_non_hos=0,
              si.flow.1_hos=0,
              si.flow.1_non_hos=0
)
control=control.dcm(nstep=365*3,new.mod=RSV)
sim=dcm(params,init,control)

sim=mutate_epi(sim,num=v06.num + s06.num + h06.num + nh06.num + r06.num +
                        v1.num  + s1.num  +  h1.num +  nh1.num +  r1.num  + 
                        s5.num  +  h5.num +  nh5.num +  r5.num  +
                       s10.num +  h10.num +  nh10.num +  r10.num  +
                       s30.num +  h30.num +  nh30.num +  r30.num  +
                       s40.num +  h40.num +  nh40.num +  r40.num  +
                       s60.num +  h60.num +  nh60.num +  r60.num)

sim=mutate_epi(sim,incidence06 =si.flow.06_hos + si.flow.06_non_hos)
sim=mutate_epi(sim,incidence1  =si.flow.1_hos + si.flow.1_non_hos)
sim=mutate_epi(sim,incidence0_12months=incidence06 + incidence1)

sim=mutate_epi(sim,
              hosp_incidence06 = si.flow.06_hos, 
              hosp_incidence1  = si.flow.1_hos)
sim=mutate_epi(sim, hosp_incidence0_12months = hosp_incidence06 + hosp_incidence1)

df=data.frame(sim)

# Rn
df %<>% mutate(Rn=R0_ * (s06.num+s1.num+s5.num+s10.num+s30.num+s40.num+s60.num)/num )

# incidence
c.incidence06=c()
c.incidence1=c()
c.hosp06=c()
c.hosp1=c()

timestep=365*3
for (i in 1:timestep){
        c.incidence06[i]  = sum(df$incidence06[seq(1:i)])
        c.incidence1[i]   = sum(df$incidence1[seq(1:i)])
    
        c.hosp06[i]    =sum(df$hosp_incidence06[seq(1:i)])
        c.hosp1[i]     =sum(df$hosp_incidence1[seq(1:i)])
}

df %<>% mutate(cum_incidence06=c.incidence06)
df %<>% mutate(cum_incidence1=c.incidence1)
df %<>% mutate(cum_incidence0_12months= cum_incidence06+cum_incidence1)
df %<>% mutate(cum_hosp06    =c.hosp06)
df %<>% mutate(cum_hosp1     =c.hosp1)
df %<>% mutate(cum_hosp0_12months =cum_hosp06 +cum_hosp1)
return(df)
}

avert_summary=function(days){
    # case
    # 0-12 months
    # non vaccine
    scenario1= df_com %>% filter(scenario==1) %>% filter(time==days) %>%
    dplyr::select(cum_incidence0_12months)%>% pull()
    # vaccine
    scenario2= df_com %>% filter(scenario==2) %>% filter(time==days) %>%
    dplyr::select(cum_incidence0_12months)%>% pull()
    # averted
    value0_12=(scenario1-scenario2)/scenario1
    
    # 0-5 months
    # non vaccine
    scenario1= df_com %>% filter(scenario==1) %>% filter(time==days) %>%
    dplyr::select(cum_incidence06)%>% pull()
    # vaccine
    scenario2= df_com %>% filter(scenario==2) %>% filter(time==days) %>%
    dplyr::select(cum_incidence06)%>% pull()
    # averted
    value0_5=(scenario1-scenario2)/scenario1
    
    # 6-12 months
    # non vaccine
    scenario1= df_com %>% filter(scenario==1) %>% filter(time==days) %>%
    dplyr::select(cum_incidence1)%>% pull()
    # vaccine
    scenario2= df_com %>% filter(scenario==2) %>% filter(time==days) %>%
    dplyr::select(cum_incidence1)%>% pull()
    # averted
    value1=(scenario1-scenario2)/scenario1
    
    print("Case averted")
    print(paste("0-12 months:",value0_12,"0-5 months:",value0_5,"6-12 months",value1))
    
    # hospitalization
    # 0-12 months
    # non vaccine
    scenario1= df_com %>% filter(scenario==1) %>% filter(time==days) %>%
    dplyr::select(cum_hosp0_12months)%>% pull()
    # vaccine
    scenario2= df_com %>% filter(scenario==2) %>% filter(time==days) %>%
    dplyr::select(cum_hosp0_12months)%>% pull()
    # averted
    hvalue0_12=(scenario1-scenario2)/scenario1
    
    # 0-5 months
    # non vaccine
    scenario1= df_com %>% filter(scenario==1) %>% filter(time==days) %>%
    dplyr::select(cum_hosp06)%>% pull()
    # vaccine
    scenario2= df_com %>% filter(scenario==2) %>% filter(time==days) %>%
    dplyr::select(cum_hosp06)%>% pull()
    # averted
    hvalue0_5=(scenario1-scenario2)/scenario1
    
    # 6-12 months
    # non vaccine
    scenario1= df_com %>% filter(scenario==1) %>% filter(time==days) %>%
    dplyr::select(cum_hosp1)%>% pull()
    # vaccine
    scenario2= df_com %>% filter(scenario==2) %>% filter(time==days) %>%
    dplyr::select(cum_hosp1)%>% pull()
    # averted
    hvalue1=(scenario1-scenario2)/scenario1
    print("Hospitalization averted")
    print(paste("0-12 months:",hvalue0_12,"0-5 months:",hvalue0_5,"6-12 months",hvalue1))
}

My_Theme = ggplot2::theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  panel.background = element_rect(fill = "white", colour = "grey50"))


#---------------------------------------------------------------------------------
# Primary objective
# step2: functions
## -------------------------------------------------------------------------------
# prepare dataframe
df1=estimates(omega_=0,psi_=0) 
df1 %<>% mutate(scenario=1)
df2=estimates(omega=0.5,psi=0.8) 
df2 %<>% mutate(scenario=2)
df_com=rbind(df1,df2)


## -------------------------------------------------------------------------------
png(filename="fig_Rn.png",
    units = "in",
    width = 7.5,
    height= 5,
    res=300)

df_com %>% filter(scenario ==1) %>%
    ggplot(aes(x=time, y=Rn),col="gray20")+
    geom_line()+
    geom_hline(yintercept=1.0, color="red",linetype=2)+
    xlab("Days") + 
    ylab("Net reproducthe number")+
    My_Theme

dev.off()


## -------------------------------------------------------------------------------
# range of rn
df_com %>% filter(scenario ==1) %>% summarize(range(Rn))

# convergence
df_com %>% filter(scenario ==1) %>% filter(time==365*3) %>% pull(Rn)


## -------------------------------------------------------------------------------
df_com %>% filter(scenario==1) %>% 
  ggplot(aes(x=time)) +
  geom_line(aes(y=s06.num,color="s06")) + 
  geom_line(aes(y=s1.num,color="s1"))


## -------------------------------------------------------------------------------
png(filename="fig_cases_total.png",
    units = "in",
    width = 7,
    height= 7,
    res=300)

df_com %>% filter(scenario==1) %>%
    ggplot(aes(x=time))+
    geom_line(aes(y=incidence0_12months),col="gray20")+
    xlab("Days")+
    ylab("RSV incidence")+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
# total peak 1st wave
first=df_com %>% filter(scenario==1) %>% filter(time %in% seq(0,280)) %>% summarize(max=max(incidence0_12months))

# total peak 2nd wave
second=df_com %>% filter(scenario==1) %>% filter(time %in% seq(280,450)) %>% summarize(max=max(incidence0_12months))

# total peal 3rd wave
third=df_com %>% filter(scenario==1) %>% filter(time %in% seq(450,650)) %>% summarize(max=max(incidence0_12months))

print(paste("1st peak:",first,"2nd peak:",second,"3rd peak:",third))


## -------------------------------------------------------------------------------
# total peak 1st wave
first=df_com %>% filter(scenario==1) %>% filter(time %in% seq(0,280)) %>% summarize(max=max(incidence06))

# total peak 2nd wave
second=df_com %>% filter(scenario==1) %>% filter(time %in% seq(280,450)) %>% summarize(max=max(incidence06))

# total peal 3rd wave
third=df_com %>% filter(scenario==1) %>% filter(time %in% seq(450,650)) %>% summarize(max=max(incidence06))

print(paste("1st peak:",first,"2nd peak:",second,"3rd peak:",third))


## -------------------------------------------------------------------------------
# total peak 1st wave
first=df_com %>% filter(scenario==1) %>% filter(time %in% seq(0,280)) %>% summarize(max=max(incidence1))

# total peak 2nd wave
second=df_com %>% filter(scenario==1) %>% filter(time %in% seq(280,450)) %>% summarize(max=max(incidence1))

# total peal 3rd wave
third=df_com %>% filter(scenario==1) %>% filter(time %in% seq(450,650)) %>% summarize(max=max(incidence1))

print(paste("1st peak:",first,"2nd peak:",second,"3rd peak:",third))


## -------------------------------------------------------------------------------
png(filename="fig_cases_subcategory.png",
    units = "in",
    width = 8,
    height= 7,
    res=300)

df_com %>% filter(scenario==1) %>%
    ggplot(aes(x=time))+
    geom_line(aes(y=incidence06, color="b",linetype="b"))+
    geom_line(aes(y=incidence1, color="c",linetype="c"))+
    scale_color_discrete(name="Age group",labels=c("0-5 months","6-12 months"))+
    scale_linetype_discrete(name="Age group",labels=c("0-5 months","6-12 months"))+
    xlab("Days")+
    ylab("RSV incidence")+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_cases_final_report.png",
    units = "in",
    width = 10,
    height= 7,
    res=300)
df_com %>% filter(scenario==1) %>%
    ggplot(aes(x=time))+
    geom_line(aes(y=incidence0_12months, col="a",linetype="a"))+
    geom_line(aes(y=incidence06, color="b",linetype="b"))+
    geom_line(aes(y=incidence1, color="c",linetype="c"))+
    scale_color_discrete(name="Age group",labels=c("0-12 months","0-5 months","6-12 months"))+
    scale_linetype_discrete(name="Age group",labels=c("0-12 months","0-5 months","6-12 months"))+
    xlab("Days")+
    ylab("RSV incidence")+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_hosp_total.png",
    units = "in",
    width = 7,
    height= 7,
    res=300)

df_com %>% filter(scenario==1) %>%
    ggplot(aes(x=time))+
    geom_line(aes(y=hosp_incidence0_12months),col="gray20")+
    xlab("Days")+
    ylab("RSV hospitalization")+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_hosp_subcategory.png",
    units = "in",
    width = 8,
    height= 7,
    res=300)

df_com %>% filter(scenario==1) %>%
    ggplot(aes(x=time))+
    geom_line(aes(y=hosp_incidence06, color="b",linetype="b"))+
    geom_line(aes(y=hosp_incidence1, color="c",linetype="c"))+
    scale_color_discrete(name="Age group",labels=c("0-5 months","6-12 months"))+
    scale_linetype_discrete(name="Age group",labels=c("0-5 months","6-12 months"))+
    xlab("Days")+
    ylab("RSV hospitalization")+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_hosp_final_report.png",
    units = "in",
    width = 10,
    height= 7,
    res=300)

df_com %>% filter(scenario==1) %>%
    ggplot(aes(x=time))+
    geom_line(aes(y=hosp_incidence0_12months, col="a",linetype="a"))+  
   geom_line(aes(y=hosp_incidence06, color="b",linetype="b"))+
    geom_line(aes(y=hosp_incidence1, color="c",linetype="c"))+
    scale_color_discrete(name="Age group",labels=c("0-12 months","0-5 months","6-12 months"))+
    scale_linetype_discrete(name="Age group",labels=c("0-12 months","0-5 months","6-12 months"))+
    xlab("Days")+
    ylab("RSV hospitalization")+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
# total peak 1st wave
first=df_com %>% filter(scenario==1) %>% filter(time %in% seq(0,280)) %>% summarize(max=max(hosp_incidence0_12months))

# total peak 2nd wave
second=df_com %>% filter(scenario==1) %>% filter(time %in% seq(280,450)) %>% summarize(max=max(hosp_incidence0_12months))

# total peal 3rd wave
third=df_com %>% filter(scenario==1) %>% filter(time %in% seq(450,650)) %>% summarize(max=max(hosp_incidence0_12months))

print(paste("1st peak:",first,"2nd peak:",second,"3rd peak:",third))


## -------------------------------------------------------------------------------
# total peak 1st wave
first=df_com %>% filter(scenario==1) %>% filter(time %in% seq(0,280)) %>% summarize(max=max(hosp_incidence06))

# total peak 2nd wave
second=df_com %>% filter(scenario==1) %>% filter(time %in% seq(280,450)) %>% summarize(max=max(hosp_incidence06))

# total peal 3rd wave
third=df_com %>% filter(scenario==1) %>% filter(time %in% seq(450,650)) %>% summarize(max=max(hosp_incidence06))

print(paste("1st peak:",first,"2nd peak:",second,"3rd peak:",third))


## -------------------------------------------------------------------------------
# total peak 1st wave
first=df_com %>% filter(scenario==1) %>% filter(time %in% seq(0,280)) %>% summarize(max=max(hosp_incidence1))

# total peak 2nd wave
second=df_com %>% filter(scenario==1) %>% filter(time %in% seq(280,450)) %>% summarize(max=max(hosp_incidence1))

# total peal 3rd wave
third=df_com %>% filter(scenario==1) %>% filter(time %in% seq(450,650)) %>% summarize(max=max(hosp_incidence1))

print(paste("1st peak:",first,"2nd peak:",second,"3rd peak:",third))


## -------------------------------------------------------------------------------
png(filename="fig_case_total_vac.png",
    units = "in",
    width = 7,
    height= 7,
    res=300)

df_com %>% ggplot(aes(x=time))+
    geom_line(aes(y=incidence0_12months, color=as.factor(scenario)))+
    xlab("Days")+
    ylab("RSV incidence")+
    scale_color_discrete(name="Vaccine scenario",labels=c("No vaccine"," 50% coverage+ 80% efficacy"))+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_case_0_5_vac.png",
    units = "in",
    width = 7,
    height= 7,
    res=300)

df_com %>% ggplot(aes(x=time))+
    geom_line(aes(y=incidence06, color=as.factor(scenario),linetype=as.factor(scenario)))+
    xlab("Days")+
    ylab("RSV incidence")+
    scale_color_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    scale_linetype_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_case_6-12_vac.png",
    units = "in",
    width = 7,
    height= 7,
    res=300)

df_com %>% ggplot(aes(x=time))+
    geom_line(aes(y=incidence1, color=as.factor(scenario),linetype=as.factor(scenario)))+
    xlab("Days")+
    ylab("RSV incidence")+
    scale_color_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    scale_linetype_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_hosp_total_vac.png",
    units = "in",
    width = 8,
    height= 3,
    res=300)

df_com %>% ggplot(aes(x=time))+
    geom_line(aes(y=hosp_incidence0_12months, color=as.factor(scenario)))+
    xlab("Days")+
    ylab("RSV hospitalization")+
    scale_color_discrete(name="Vaccine scenario",labels=c("No vaccine"," 50% coverage+ 80% efficacy"))+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_hosp_0_5_vac.png",
    units = "in",
    width = 8,
    height= 3,
    res=300)

df_com %>% ggplot(aes(x=time))+
    geom_line(aes(y=hosp_incidence06,
                  color=as.factor(scenario),linetype=as.factor(scenario)))+
    xlab("Days")+
    ylab("RSV hospitalization")+
    scale_color_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    scale_linetype_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
png(filename="fig_hosp_6-12_vac.png",
    units = "in",
    width = 8,
    height= 3,
    res=300)

df_com %>% ggplot(aes(x=time))+
    geom_line(aes(y=hosp_incidence1,
                  color=as.factor(scenario),linetype=as.factor(scenario)))+
    xlab("Days")+
    ylab("RSV incidence")+
    scale_color_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    scale_linetype_discrete(name="Vaccine scenario",labels=c("No vaccine","50% coverage+ 80% efficacy"))+
    My_Theme
dev.off()


## -------------------------------------------------------------------------------
# averted(time)  input must be days
avert_summary(365)
avert_summary(365*2)
avert_summary(365*3)


#---------------------------------------------------------------------------------
# Secondary objective
# step3: functions
## -------------------------------------------------------------------------------
# only interested in hospitalization 0-12 months
avert_hosp=function(days,hosp_in_scenario){
    # non vaccine
    df=estimates(omega_=0,psi_=0) 
    hosp_ref= df %>% filter(time==days) %>%
    dplyr::select(cum_hosp0_12months)%>% pull()
    # vaccine scenario
    hosp_scenario=hosp_in_scenario
    # averted
    value=(hosp_ref - hosp_scenario)/hosp_ref
    return(value)
}


#---------------------------------------------------------------------------------
# Secondary objective
# step 4: simulation (takes long)

## -------------------------------------------------------------------------------
# simulation 
omega=c(seq(0,1,0.1))
psi=c(seq(0.4,1,0.1))
setting=data.frame(expand.grid(omega,psi))

avert_summary=data.frame()
for (i in 1:nrow(setting)){
  # parameter sampling
  current_set=setting[i,]
  current_omega=current_set[,1]
  current_psi=current_set[,2]
  
  # simulation
  df=estimates(omega_=current_omega,psi=current_psi)
  hosp_cases= df %>% filter(time==365*3) %>%
    dplyr::select(cum_hosp0_12months)%>% pull()
  value=avert_hosp(365*3,hosp_cases)
  
  # summary
  avert_summary[i,1]=current_omega
  avert_summary[i,2]=current_psi
  avert_summary[i,3]=value
}
names(avert_summary) = c("omega","psi","hospital_averted")


#---------------------------------------------------------------------------------
# Secondary objective
# step 5: visualize
## -------------------------------------------------------------------------------
png(filename="fig_sensit.png",
    units = "in",
    width = 9,
    height= 7,
    res=300)

# contour plot
ggplot(avert_summary,aes(x=omega,y=psi)) +
  geom_raster(aes(fill=hospital_averted)) + 
  scale_fill_gradientn(colours=c("skyblue2","plum2","khaki2"),name="Hospitalization averted")+
  geom_contour(aes(z=hospital_averted),colour="black",size=1,alpha=0.5,breaks = c(0.1,0.2,0.3,0.4,0.5,0.6))+
  geom_text(aes(x=0.15,y=1.03),label="0.1",colour="black")+
  geom_text(aes(x=0.31,y=1.03),label="0.2",colour="black")+
  geom_text(aes(x=0.46,y=1.03),label="0.3",colour="black")+
  geom_text(aes(x=0.62,y=1.03),label="0.4",colour="black")+
  geom_text(aes(x=0.78,y=1.03),label="0.5",colour="black")+
  geom_text(aes(x=0.93,y=1.03),label="0.6",colour="black")+
  xlab("Vaccine coverage")+
  ylab("Vaccine efficacy")
dev.off()

