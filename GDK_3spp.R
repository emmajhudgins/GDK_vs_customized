## Fitting of GDKu, GDKic, GDKs, GDKs+ic, and GDKs+ic+c models for 3 United States Invasive Forest Pests##
## Written as part of Hudgins et al. "Comparing customized to generalized models for United States Forest Pests". in prep. J Ecol. ##

## Code written by Emma J. Hudgins
## PhD Candidate, McGill University
## Montreal, QC, Canada
## emma.hudgins@mail.mcgill.ca


rm(list=ls()) 
library(pdist)
plot_mode="validate" # keep this set to validate (will update to forecast later on if fit_mode='fit_to_forecast')
setwd('~/Documents/GitHub/GDK_vs_customized/') # update to your own working directory


#Modelling switches - testing different GDK layers #
spp=49 # species to fit (49=hwa, 51=gm, 54=bbd)
write_out=T # write presences and fitted parameters to csv
mk_random=F #make null model draws for r2om (output already included so not necessary)
ic=T #intercept correction
temp_t=T # temperature threshold for HWA
source="real" # start from 'centroid' or best guess of initial invasion ('real')
nm_spp<-rep(0,64)
nm_spp[c(49,51,54)]<-c("HWA", "GM", "BBD")
fit_mode="validate" ##use "fit_to_forecast" to fit to 2005 instead of 2000 (most accurate forecasts)
setwd("~/Desktop/OneDrive - McGill University/Grad/scripts/")

#Read in Data
data<-read.csv('countydatanorm_march.csv', stringsAsFactors = FALSE) # spatial data
data2<-read.csv('spdat_clean_gdk.csv', stringsAsFactors = FALSE) # species data, see Hudgins et al. corrigendum for information on why ALB (spp=3) cannot be accurately fit
host.density2<-read.csv('hostden_clean_gdk.csv') # tree host density by pest species
prez<-read.csv('prez_clean_gdk.csv') # invasible host range (from FIA)
prez2<-read.csv('prez2_clean_gdk.csv') # pest presences in 2000 for the 3 species examined here, in 2009 for all other species
L<-rep(0,64)
for (sppp in 1:64)
{
  L[sppp]<-length(which(prez[,sppp]!=0))
}

if (source=="centroid")
{sources<-as.list(read.csv('Psources_notypos.csv')[,1])}
if (source=="real")
{load('Psources_closest_march.Rdata')}

#Climatic threshold#
current_temp<-current_temp2<-0
if (temp_t==T)
{
  temp2<-temp<-read.csv('bio6_5yr.csv') 
  current_temp<-current_temp2<-temp2[,1]
}
host.density2<-read.csv("hostden_clean_gdk.csv", stringsAsFactors = FALSE)
currpopden<-as.matrix(read.csv("currpopden_clean_gdk.csv", stringsAsFactors = FALSE)) # human population density with projections to 2030




## Distance matrix pre-calculation
Tr1<-function(x)
{
  sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
}
T1<-exp(-1*sapply(1:3372, Tr1)/50000)

GDK_cor=function(par)
{
  pars<-rep(0,31)
  if (ic==F)
  {
    pars[c(1,21,22,4,18,20,8)]<-par
  }
  if (ic==T)
  {
    pars[1]<-par[1] # intercept
    pars[c(21,22,4,18,20,8)]<-c(0.000538410692229749, 0.299034706404549, -0.525670755351726, 15.6132848183217,-0.163552592351765, 0.323831382884772) # fitted values from GDKu for other important pars as published
    if (temp_t==T)
    {
      pars[31]<-par[2]
    }
  }
  par<-pars
  par[21]<-abs(par[21]) # threshold must be positive
  par[22]<-abs(par[22])+1 # growth rate must be greater than 1

  if (spp==49) #HWA
  {
    YEAR<<-50 # time of spread in fitting phase
    total_YEAR<<-55 # time including validation
    pfulz<-read.csv('hwa_3yrs.csv')[,3] # observed presences in 2005 (validation)
  }
  if (spp==51) #GM
  {
    YEAR<<-135
    total_YEAR<<-140
    pfulz<-read.csv('gm_3yrs.csv')[,3]
  }        
  if (spp==54) #BBD
  {
    YEAR<<-110
    total_YEAR<<-115
    pfulz<-read.csv('bbd_3yrs.csv')[,3]
  }
  if (fit_mode=="fit_to_forecast") # fit on 2005
  {  
    total_time<-(total_YEAR/5)+5
    fityr=total_YEAR/5
  }
  if (fit_mode=="validate") # fit on 2000, predict on 2005
  {
    fityr<-YEAR/5
    total_time<-total_YEAR/5
  }
  
  #Pest Parameters
  Pfull<<-matrix(0, 3372, total_time) # matrix will fill with indices of predicted pest presences within invasible range at each timestep (appended to zeros to maintain constant matrix size)
  pfull2<<-matrix(0,3372, (total_YEAR/5)) #  will fill with indices of observed pest presences(appended to zeros to maintain constant matrix size)
  pfull2[,(total_YEAR/5)]<<-as.matrix(pfulz) # obs distribution in 2005
  pfull2[,YEAR/5]<<-prez2[,spp] # obs distribution in 2000
  for (sppp in (total_YEAR/5))
  {pfull2[,sppp]<<-c(intersect(prez[which(prez[,spp]!=0),spp], pfull2[which(pfull2[,sppp]!=0),sppp]), rep(0,(3372-length(intersect(prez[which(prez[,spp]!=0),spp], pfull2[which(pfull2[,sppp]!=0),sppp])))))} # restrict to invasible range based on FIA
  Pfull_time<<-Pfull
  
  #time-invariant dispersal predictors
  constpD=rep(0,64)
  constpD=matrix(rep(constpD),3372,64, byrow=TRUE)
  constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]),3372,64)+par[19]*host.density2
  constpD<-as.numeric(constpD)+constpD2
  constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]),3372,64)+par[5]*host.density2
  
  #Pest Parameters
  Psource=sources[[spp]] # initial invsion location
  Discovery=2000-YEAR
  T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]] # crop distance matrix to invasible range
  vecP<-rep(0,L[spp]) # vector of pest propagule pressure
  for (rrr in 1:length(Psource))
  {vecP[which(prez[,spp]==Psource[rrr])]=1} # maintain maximum propagule pressure at source
  for (time in 1:total_time)
  {
    if (temp_t==T)
    {
      if (time<7)
      {
        current_temp<-temp2[,1]
        current_temp2<-temp2[,1]
      }
      if (time>=7)
      {
        current_temp2<-temp2[,time-5]
      }
    }
    vecP[which(prez[,spp]==Psource)]=1 # maintain max propagule pressure in source
    if (temp_t==T)
    {vecP[which(current_temp[prez[,spp]]<par[31])]<-0}
    vecP[which(prez[,spp]==Psource)]=1
    Pnext<-rep(0,L[spp]) # vecP for next timestep
    qq<-0 # dispersal kernel
    column<-(((Discovery+5*(time-1))-1790)/5)+1 # which year of human population density to consider
    qq<-matrix(rep(constpD[prez[which(vecP>=par[21]),spp],spp]+par[8]*currpopden[prez[which(vecP>=par[21]),spp],column], L[spp]), nrow=length(which(vecP>=par[21])), ncol=L[spp])
    zzz<-matrix(rep(constpD3[prez[1:L[spp],spp],spp]+par[20]*currpopden[prez[1:L[spp],spp],column], L[spp]), nrow=L[spp], ncol=L[spp], byrow=TRUE) # add in current human populations
    qq<-(2*par[1]*exp((zzz[which(vecP>=par[21]),]+qq)))/(1+exp(zzz[which(vecP>=par[21]),]+qq)) 
    qq<-T2[which(vecP>=par[21]),]^qq # add in distance
    #scale dispersal kernel
    if (length(which(vecP>=par[21]))>1){qq<-qq/rowSums(qq)}
    if (length(which(vecP>=par[21]))==1){qq<-qq/sum(qq)}
    qq[which(qq<0.001)]=0
    Pnext=(vecP[which(vecP>=par[21])])%*%(qq) # dispersal into and out of all sites
    if (temp_t==T)
    {
      Pnext[prez[which(current_temp[prez[,spp]]<par[31]),spp]]<-0
    }
   Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
    if (time>1) # do not allow extirpations
      {
        dddd<-which(!(Pfull_time[1:length(which(Pfull_time[,time-1]!=0)),time-1]%in%Pfull_time[1:length(which(Pfull_time[,time]!=0)),time]))
        ffff<-which(prez[1:length(which(prez[,spp]!=0)), spp]%in%Pfull_time[dddd,time-1])
        Pnext[ffff]<-par[21]
      Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
    }
    current_temp<-current_temp2  
    
    if (plot_mode=="forecast")# set to observed distribution before forecasting (set false absences to threshold and remove false presences)
    {
      if (time==(total_YEAR/5))
      {
        dddd<-which(prez[1:length(which(prez[,spp]!=0)),spp]%in%pfull2[1:length(which(pfull2[,total_YEAR/5]!=0)),total_YEAR/5])
        cccc<-which(!(prez[1:length(which(prez[,spp]!=0)),spp]%in%pfull2[1:length(which(pfull2[,total_YEAR/5]!=0)),total_YEAR/5]))
        eeee<-which(Pnext[dddd]<par[21])
        Pnext[dddd[eeee]]<-par[21]
        Pnext[cccc]<-0
      }
    }
    Pnext[which(Pnext>=par[21])]=Pnext[which(Pnext>=par[21])]*par[22] # growth
    Pnext[which(Pnext>=1)]<-1 # maximum propagule pressure at 1
    vecP=Pnext
    vecP[which(prez[,spp]==Psource)]=1
    Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
  
  if (time>1)
  {
  if (length(all.equal(Pfull_time[,time], Pfull_time[,time-1]))==1)# penalize parameter sets that do not spread over time
  {
    if (all.equal(Pfull_time[,time], Pfull_time[,time-1])==TRUE)
    {
      return (1000000000)
    }
  }
  }
  }
  Pfull<<-Pfull_time
  MET=function(spp) #MET calculation
  {
    dii<-2*sum(dist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), upper=F))
    dij<-sum(pdist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), cbind(data$X_coord[pfull2[1:length(which(pfull2[,spp]!=0)),spp]], data$Y_coord[pfull2[1:length(which(pfull2[,spp]!=0)),spp]]))@dist)
    djj<-2*sum(dist(cbind(data$X_coord[pfull2[1:length(which(pfull2[,spp]!=0)),spp]], data$Y_coord[pfull2[1:length(which(pfull2[,spp]!=0)),spp]])))
    return(((1/(length(which(Pfull[,spp]!=0))*length(which(pfull2[,spp]!=0))))*dij)-((1/(2*length(which(Pfull[,spp]!=0))^2))*dii)-((1/(2*length(which(pfull2[,spp]!=0))^2))*djj))
  }
  return (sum(tapply(fityr, fityr, MET)))
}

##Optimization based on MET (very sensitive to initial conditions - beware)##

if (ic==T)
{
  if (temp_t==F)
  {
    model<-optim(par=c(2), fn=GDK_cor, control=list(trace=100, maxit=1000))
  }
  if (temp_t==T)
  {
    model<-optim(par=c(1.44, -74.9), fn=GDK_cor, control=list(trace=100, maxit=1000, parscale=c(1,100)))
  }
  yy<-model$value # continue refitting until output is the same as input to avoid local minima
  xx<-1000000
  while (yy!=xx)
  {
    yy<-model$value
    if (temp_t==T)
    {
      model<-optim(par=model$par, fn=GDK_cor, control=list(trace=100, maxit=1000, parscale=c(1,100)))
    }
    if (temp_t==F)
    {
      model<-optim(par=model$par, fn=GDK_cor, control=list(trace=100, maxit=1000))  
    }
    xx<-model$value
  }
  met_fit<-GDK_cor(model$par)
}

if (ic==F) #get MET for GDKu pars
{
  pars=c(1.74213662335515,0.000538410692229749, 0.299034706404549, -0.525670755351726, 15.6132848183217,-0.163552592351765, 0.323831382884772)
  met_fit<-GDK_cor(c(pars)) # fitting MET score
  if (fit_mode=="fit_to_forecast")
  {
    plot_mode="forecast"
    GDK_cor(par=pars)
  }
}
if (fit_mode=="fit_to_forecast")
{
  plot_mode="forecast"
  GDK_cor(par=c(model$par)) # Pfull now contains forecast)
}

if (write_out==T) # write out best fitting parameters and distribution over time 
{
  if (ic==T)
  {
  write.csv(model$par, file=paste("par_GDK", nm_spp,ifelse(source=="real", "s", ""),  ifelse(ic, "ic", "u"),  ifelse(temp_t==T, "c", ""),ifelse(fit_mode=='fit_to_forecast', "forecast", "fit"),".csv", sep='_'))
  }
    write.csv(Pfull, file=paste("presences_GDK", nm_spp,ifelse(source=="real", "s", ""),  ifelse(ic, "ic", "u"),  ifelse(temp_t==T, "c", ""),ifelse(fit_mode=='fit_to_forecast', "forecast", "fit"),".csv", sep='_'))
}

## R2om calculation ##
if (fit_mode=="validate")
{
  pred<-Pfull[1:length(which(Pfull[,total_YEAR/5]!=0)),total_YEAR/5]# predicted presences
  obs<-pfull2[1:length(which(pfull2[,total_YEAR/5]!=0)),total_YEAR/5] # observed presences
  if (mk_random==T)
  {
    source('random_nn_spp_2019.R')
    random_nn<-random_nn_spp(spp)
    write.csv(random_nn, file=paste0("random_nn_",spp,".csv"), row.names=F)
    r2<-1-(om_mse(pred,obs)/mean(random_nn))
  }
  if (mk_random==F)
  {
    source('r2_om_mse.R')
    r2<-r2_om_mse(pred,obs,spp)
  }
  r2
}