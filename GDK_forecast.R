## Code to fit relevant SDK layers for all species in the Hudgins et al. (2017) dataset##
## Written as part of Hudgins et al. "Comparing customized to generalized models for United States Forest Pests". in prep. J Ecol. ##

## Code written by Emma J. Hudgins
## PhD Candidate, McGill University
## Montreal, QC, Canada
## emma.hudgins@mail.mcgill.ca

## We used GDK_3spp.R, not this script, to fit focal species from publication (49,51,54)

rm(list=ls()) 
library(pdist)
temp_t=F # minimum temperature threshold, should be set to true for spp=48
hum_t=F # maximum humidity threshold, should be set to true for spp=44
start="centroid" # spread initiation point, either 'real' or 'centroid'
spp=1 # species of interest
#setwd('~/Documents/GitHub/GDK_vs_customized/') # update to your own working directory
ic=F # fit intercept correction?

#Read in Data
data<-read.csv('countydatanorm_march.csv', stringsAsFactors = FALSE) # spatial data, including pest distributions in 2009 (fyi: other species are fit on 2005 to make timesteps equal)
data2<-read.csv('spdat_clean_gdk.csv', stringsAsFactors = FALSE) # species data, see Hudgins et al. corrigendum for information on why ALB (spp=3) cannot be accurately fit
prez<-read.csv('prez_clean_gdk.csv') # invasible range by pest species (based on FIA)
prez2<-read.csv('prez2_clean_2009.csv') # pest presences in 2009
host.density2<-read.csv('hostden_clean_gdk.csv') # tree host density by pest species
L<-rep(0,64)
for (sppp in 1:64)
{
  L[sppp]<-length(which(prez[,sppp]!=0))
}
if (spp==48)
{
  temp_t=T
}
if (spp==44)
{
  hum_t=T
}
#Climatic thresholds#
current_temp<-current_temp2<-0
if (temp_t==T)
{
  temp2<-temp<-read.csv('bio6_5yr.csv') 
  current_temp<-current_temp2<-temp2[,1]
}
hum=0
if (hum_t==T)
{
  hum<-read.csv('interp_rel_humidity.csv')[,1]
}
currpopden<-as.matrix(read.csv("currpopden_clean_gdk.csv", stringsAsFactors = FALSE)) # human population density with projections to 2030

if (start=="centroid")
{
  sources<-as.list(read.csv('Psources_notypos.csv')[,1])
}
if (start=="real")
{
  load("Psources_closest_march.Rdata")
}

Tr1<-function(x) # distance matrix pre-calculation
{
  sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
}
T1<-exp(-1*sapply(1:3372, Tr1)/50000)

GDKic=function(par)
{
  pars<-rep(0,32)
  pars[c(1)]<-par[1] # intercept
  if (temp_t==T)
  {
    pars[31]<-par[2]
  }
  if (hum_t==T)
  {
    pars[32]<-par[2]
  }
  pars[c(21,22,4,18,20,8)]<-c(0.000538410692229749, 0.299034706404549, -0.525670755351726, 15.6132848183217,-0.163552592351765, 0.323831382884772) #best fit for GDKu on 5-year timescale
  par<-pars
  par[21]<-abs(par[21]) # threshold must be positive
  par[22]<-abs(par[22])+1 # growth rate must be greater than 1
  
  #time-invariant components of dispersal kernel
  constpD=rep(0,64)
  constpD=matrix(rep(constpD),3372,64, byrow=TRUE)
  constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]),3372,64)+par[19]*host.density2
  constpD<-as.numeric(constpD)+constpD2
  constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]),3372,64)+par[5]*host.density2
  
  Pfull<<-Pfull_time<<-matrix(0, 3372, 64) # matrix will fill with indices of predicted pest presences within invasible range at each timestep (appended to zeros to maintain constant matrix size)
  
  #Pest Parameters
  YEAR=data2$YEAR[spp] # number of years of spread
  Psource=sources[[spp]] # source of initial invasion
  Discovery<-2009-YEAR # discovery year
  
  T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]] # distance matrix cropped to invasible area
  
  vecP<-rep(0,L[spp]) # vector of pest propagule pressure
  for (rrr in 1:length(Psource))
  {vecP[which(prez[,spp]==Psource[rrr])]=1} # maintain maximum propagule pressure at source
  for (time in 1:((YEAR/5)+5)) # time loop including 25 years of forecasting
  {
    if (temp_t==T)
    {
      if (time<(floor(YEAR/5)-4))
      {
        current_temp<-temp2[,1]
        current_temp2<-temp2[,1]
      }
      if (time>=(floor(YEAR/5)-4))
      {
        current_temp2<-temp2[,time-floor(YEAR/5)+6]
      }
    }
    vecP[which(prez[,spp]==Psource)]=1 # maintain max propagule pressure in source
    if (temp_t==T)
    {vecP[which(current_temp[prez[,spp]]<par[31])]<-0}
    if (hum_t==T)
    {vecP[which(hum[prez[,spp]]>par[32])]<-0}
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
      Pnext[which(current_temp[prez[,spp]]<par[31])]<-0
    }
    if (hum_t==T)
    {
      Pnext[which(hum[prez[,spp]]>par[32])]<-0
    }
    Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) # record presences
    if (time>1) # do not allow extirpations
    {
      dddd<-which(!(Pfull_time[1:length(which(Pfull_time[,time-1]!=0)),time-1]%in%Pfull_time[1:length(which(Pfull_time[,time]!=0)),time]))
      ffff<-which(prez[1:length(which(prez[,spp]!=0)), spp]%in%Pfull_time[dddd,time-1])
      Pnext[ffff]<-par[21]
      if (time==floor(YEAR/5)) # set to observed distribution before forecasting (set false absences to threshold and remove false presences)
      {
        Pfull_test<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21]))))
        dddd<-which(prez[1:length(which(prez[,spp]!=0)),spp]%in%prez2[1:length(which(prez2[,spp]!=0)),spp])
        cccc<-which(!(prez[1:length(which(prez[,spp]!=0)),spp]%in%prez2[1:length(which(prez2[,spp]!=0)),spp]))
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
    current_temp<-current_temp2
    if (ic==T)
    {
    if ((spp %in% c(7,12,59,28))==F)
    {
      if (length(which(Pfull_time[,time]!=0))<(L[spp]*0.9)) # if pest has not spread throughout 90% of its range
      {
        if (length(all.equal(Pfull_time[,time], Pfull_time[,time-1]))==1)
        {
          if (all.equal(Pfull_time[,time], Pfull_time[,time-1])==TRUE) # penalize parameter sets where pest does not spread
          {
            return (1000000000)
          }
        }
      }
    }
    if ((spp %in% c(7,59,28))==T) # species with very small invasible ranges require a more relaxed penalty
    {
      if (length(which(Pfull_time[,time]!=0))<(L[spp]-12)) # if pest has not spread throughout all but a small number of sites
      {
        if (length(all.equal(Pfull_time[,time], Pfull_time[,time-1]))==1)
        {
          if (all.equal(Pfull_time[,time], Pfull_time[,time-1])==TRUE) # penalize parameter sets where pest does not spread
          {
            return (1000000000)
          }
        }
      }
    }
    }
  }
  Pfull[,spp]<<-Pfull_test
  MET=function(spp) # MET calculation
  {
    dii<-2*sum(dist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), upper=F))
    dij<-sum(pdist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]]))@dist)
    djj<-2*sum(dist(cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]])))
    return(((1/(length(which(Pfull[,spp]!=0))*length(which(prez2[,spp]!=0))))*dij)-((1/(2*length(which(Pfull[,spp]!=0))^2))*dii)-((1/(2*length(which(prez2[,spp]!=0))^2))*djj))
  }
  return (sum(tapply(spp, spp, MET)))
}

##Optimization based on MET (very sensitive to initial conditions - beware)##
if (ic==T)
{
if(temp_t==T)
{
  model<-optim(par=c(1.76, -124.3), fn=GDKic, control=list(trace=100, maxit=1000, parscale=c(1,100))) # parscale helps with differential magnitude of temp threshold and intercept

}
if(hum_t==T)
{
  model<-optim(par=c(2.1,70), fn=GDKic, control=list(trace=100, maxit=1000, parscale=c(1,100)))
}
if(temp_t==F & hum_t==F)
{
  model<-optim(par=c(2), fn=GDKic, control=list(trace=100, maxit=1000))
  if (model$par==2) # if stuck at initial conditions, try another set, and so on
  {
    model<-optim(par=c(1), fn=GDKic, control=list(trace=100, maxit=1000))
    if (model$par==1)
    {
      model<-optim(par=c(4), fn=GDKic, control=list(trace=100, maxit=1000))
    }
    if (model$par==4)
    {
      model<-optim(par=c(0.1), fn=GDKic, control=list(trace=100, maxit=1000))
    }
    if (model$par==0.1)
    {
      model<-optim(par=c(0.25), fn=GDKic, control=list(trace=100, maxit=1000))
    }
  }
}
yy<-model$value # continue refitting until output is the same as input to avoid local minima
xx<-1000000
while (yy!=xx)
{
  yy<-model$value
  model<-optim(par=model$par, fn=GDKic, control=list(trace=100, maxit=1000))
  xx<-model$value
}
GDKic(model$par)
}

if (ic==F)
{
GDKic(1.74213662335515) # overall optimal intercept for 5-year GDK
}


##save output##
write.csv(Pfull_time, file=paste("presences_GDK", ifelse(ic, "ic", "u"), ifelse(start=="real", "real", "centroid"),spp,"csv",sep='.'), row.names=F) # vector of presences appended to string of zeros to maintain matrix size, including forecast
write.csv(model$par, file=paste("par_GDK", ifelse(ic, "ic", "u"), ifelse(start=="real", "real", "centroid"),spp,"csv",sep='.'), row.names=F) # optimal parameter set 


