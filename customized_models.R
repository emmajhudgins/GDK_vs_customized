## Code to fit customized models via forward selection for 3 United States Invasive Forest Pests##
## Written as part of Hudgins et al. "Comparing customized to generalized models for United States Forest Pests". in prep. J Ecol. ##

## Code written by Emma J. Hudgins
## PhD Candidate, McGill University
## Montreal, QC, Canada
## emma.hudgins@mail.mcgill.ca

rm(list=ls()) 
library(pdist)
i=1 # see below for meaning of i (species/model combination)
#setwd('~/Documents/GitHub/GDK_vs_customized/') # update to your own working directory

#Read in Data
data<-read.csv('countydatanorm_march.csv', stringsAsFactors = FALSE) # spatial data
data2<-read.csv('spdat_clean_gdk.csv', stringsAsFactors = FALSE) # species data, see Hudgins et al. corrigendum for information on why ALB (spp=3) cannot be accurately fit
prez<-read.csv('prez_clean_gdk.csv') # invasible range by pest species (based on FIA)
prez2<-read.csv('prez2_clean_2009.csv') # pest distributions in 2009 (fyi: these 3 species are fit on 2005 to make timesteps equal)
host.density2<-read.csv('hostden_clean_gdk.csv') # tree host density by pest species
currpopden<-as.matrix(read.csv("currpopden_clean_gdk.csv", stringsAsFactors = FALSE)) # human population density with projections to 2030
L<-rep(0,64)
for (sppp in 1:64)
{
  L[sppp]<-length(which(prez[,sppp]!=0))
}
AHS<-read.csv('AHS_grid_scaled.csv') #firewood and campground variables from the American Housing Survey

#Climatic threshold#
current_temp<-current_temp2<-temp<-temp2<-0

## Distance matrix pre-calculation
Tr1<-function(x)
{
  sqrt((data$X_coord-data$X_coord[x])^2+(data$Y_coord-data$Y_coord[x])^2)
}
T1<-exp(-1*sapply(1:3372, Tr1)/50000)


while(i<=12 && i!=-99)
{
  temp_t=F
  cc<-i%%3 #iterate 4 model types * 3 species
  if (cc==0)
  {
    cc=3
  }
  if (cc==1)
  {
  #temperature threshold switch, turned off for GM and BBD
  temp_t=T
  temp2<-temp<-read.csv('bio6_5yr.csv') 
  current_temp<-current_temp2<-temp2[,1]
  }
  spp<-c(49,51,54)[cc]
  if (i %in% c(1,2,3))
  {
    start_mode="centroid"
    dispersal_mode="normal"
  }
  if (i %in% c(4,5,6))
  {
    start_mode="centroid"
    dispersal_mode="lepto"
  }
  if (i %in% c(7,8,9))
  {
    start_mode="real"
    dispersal_mode="normal"
  }
  if (i %in% c(10,11,12))
  {
    start_mode="real"
    dispersal_mode="lepto"
  }
  if (start_mode=="real")
  {
    load('Psources_closest_march.Rdata') #'best guess' spread initiation
  }
  if (start_mode=="centroid")
  {
    sources<-as.list(read.csv('Psources_notypos.csv')[,1]) # host centroid spread initiation
  }
  
  # Spread  simulation fitting funtion
  customized=function(par)
  {
    pars<-rep(0,31)
    pars[others]<-par
    par<-pars
    par[22]<-abs(par[22])+1 #prevent negative growth
    par[21]<-abs(par[21]) #prevent inverse distance dispersal
    
    #time-invariant components of dispersal kernel
    constpD=matrix(rep(0),3372,64, byrow=TRUE) 
    constpD2<-matrix(rep(par[9]*data[,18]+par[10]*data[,16]+par[16]*data[,19]+par[17]*data[,20]+par[18]*data[,21]+par[25]*AHS$num_woodfuel+par[26]*AHS$num_seasonal+par[29]*AHS$campground),3372,64)+par[19]*host.density2
    constpD<-as.numeric(constpD)+constpD2
    constpD3<-matrix(rep(par[4]*data[,19]+par[12]*data[,20]+par[3]*data[,21]+par[13]*data[,18]+par[15]*data[,16]+par[27]*AHS$num_woodfuel+par[28]*AHS$num_seasonal+par[30]*AHS$campground),3372,64)+par[5]*host.density2
    
    Pfull<<-Pfull_time<<-matrix(0, 3372, 64) # matrix will fill with indices of predicted pest presences within invasible range at each timestep (appended to zeros to maintain constant matrix size)
    
    #HWA is spp #49, GM #51, BBD #54
    if (spp==49)
    {
      YEAR=50 # number of years to fit on
      prez2<-read.csv('hwa_5yr_fitting_july.csv')[,3:12] # observed distribution in 5-year timesteps to 2000
    }
    if (spp==51)
    {
      YEAR=135
      prez2<-read.csv('gm_5yr_fitting_july.csv')[,3:29]
    }
    if (spp==54)
    {
      YEAR=110
      prez2<-read.csv('bbd_5yr_fitting_july.csv')[,3:24]
    }
    fitting_yrs<-which(colSums(prez2)!=0)
    for (sppp in 1:ncol(prez2))
    {prez2[,sppp]<-c(intersect(prez[which(prez[,spp]!=0),spp], prez2[which(prez2[,sppp]!=0),sppp]), rep(0,(3372-length(intersect(prez[which(prez[,spp]!=0),spp], prez2[which(prez2[,sppp]!=0),sppp])))))} # restrict to known invasible distribution via FIA
    
    #Pest Parameters
    Psource=sources[[spp]] # source of initial invasion
    Discovery<-2000-YEAR # discovery year
    
    T2<-T1[prez[1:L[spp],spp],prez[1:L[spp],spp]] # distance matrix cropped to invasible area
    
    vecP<-rep(0,L[spp]) # vector of pest propagule pressure
    for (rrr in 1:length(Psource))
    {vecP[which(prez[,spp]==Psource[rrr])]=1} # maintain maximum propagule pressure at source
    for (time in 1:(floor(YEAR/5)+1))
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
      {
        vecP[which(current_temp[prez[,spp]]<par[31])]<-0
      }
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
      
      Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) # record presences
      if (time>1) # do not allow extirpations
      {
        dddd<-which(!(Pfull_time[1:length(which(Pfull_time[,time-1]!=0)),time-1]%in%Pfull_time[1:length(which(Pfull_time[,time]!=0)),time]))
        ffff<-which(prez[1:length(which(prez[,spp]!=0)), spp]%in%Pfull_time[dddd,time-1])
        Pnext[ffff]<-par[21]
      }
    
      Pnext[which(Pnext>=par[21])]=Pnext[which(Pnext>=par[21])]*par[22] # growth
      Pnext[which(Pnext>=1)]<-1 # maximum propagule pressure at 1
      vecP=Pnext
      vecP[which(prez[,spp]==Psource)]=1
      Pfull[,time]<<-Pfull_time[,time]<<-c(prez[which(Pnext>=par[21]),spp], rep(0, 3372-length(which(Pnext>=par[21])))) 
    }
    if (length(all.equal(Pfull_time[,time], Pfull_time[,time-1]))==1)# penalize parameter sets that do not spread over time
    {
      if (all.equal(Pfull_time[,time], Pfull_time[,time-1])==TRUE)
      {
        return (1000000000)
      }
    }
    
    MET=function(spp) # MET calculation in meters
    {
      dii<-2*sum(dist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), upper=F))
      dij<-sum(pdist(cbind(data$X_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]], data$Y_coord[Pfull[1:length(which(Pfull[,spp]!=0)),spp]]), cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]]))@dist)
      djj<-2*sum(dist(cbind(data$X_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]], data$Y_coord[prez2[1:length(which(prez2[,spp]!=0)),spp]])))
      return(((1/(length(which(Pfull[,spp]!=0))*length(which(prez2[,spp]!=0))))*dij)-((1/(2*length(which(Pfull[,spp]!=0))^2))*dii)-((1/(2*length(which(prez2[,spp]!=0))^2))*djj))
    }
    return (sum(tapply(fitting_yrs, fitting_yrs, MET)))
  }
  
  otherscale<-NULL
  others<-c(1,21,22)
  
  if (temp_t==T)
  {
    others<-c(1,21,22,31)
  }
  step_num=0
  
  #Number of independent fitting timesteps per species
  if (cc==1)
  {
    fityrs=5
  }
  if (cc==2)
  {
    fityrs=15
  }
  if (cc==3)
  {
    fityrs=6
  }
  #Base model (intercept only)
  if (temp_t==T)
  {
    model<-optim(par=c(1.74213662335515,0.000538410692229749, 0.299034706404549,-100), fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1,100))) }
  if (temp_t==F)
  {
    model<-optim(par=c(1.74213662335515,0.000538410692229749, 0.299034706404549), fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1))) }
  ## restart optim to assure convergence to same parameter set
  yy<-model$value
  xx<-1000000
  while (yy!=xx)
  {
    yy<-model$value
    if (temp_t==T)
    {
      model<-optim(par=model$par, fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1,100)))   
      xx<-model$value
    }
    if (temp_t==F)
    {
      model<-optim(par=model$par, fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1)))   
      xx<-model$value
    }
  }
  
  ## Write output to files
  write.csv(Pfull, file=paste("presences_customized", step_num,i,cc,"csv",sep='.'))
  write.table(model$par, file=paste("par_customized", step_num,i,cc,"csv",sep='.'))
  
  #Forward stepwise variable selection  
  if (model$value>1)
  {
    qqq<-c(4,5,8,9,10,12,13,15,16,17,18,19,20,23,24,25,26,27,28,29,30)
    improvement<-1000000
    m<-replicate(30,model)
    METS<-rep(NA,30)
    while(improvement>(5000*fityrs)) #keep variables that improve MET by at least 5km on average across fitting years
    {
      step_num<-step_num+1
      otherscale<-c(otherscale,1)
      win<-0
      for (step in qqq)
      {
        others<-c(others,step)
        if (temp_t==T)
        {
          m[[step]]<-optim(par=c(model$par, 0), fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1,100,otherscale)))
        }
        if (temp_t==F)
        {
          m[[step]]<-optim(par=c(model$par, 0), fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1,otherscale)))
        }
        yy<- m[[step]]$value
        xx<-1000000
        while (yy!=xx)
        {
          yy<- m[[step]]$value
          if (temp_t==T)
          {
            m[[step]]<-optim(par= m[[step]]$par, fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1,100,otherscale)))
          }
          if (temp_t==F)
          {
            m[[step]]<-optim(par= m[[step]]$par, fn=customized, control=list(trace=100, maxit=1000, parscale=c(0.1,0.0001,1,otherscale)))
          }
          xx<- m[[step]]$value
        }
        METS[step]<-m[[step]]$value
        others<-others[1:(length(others)-1)]
      }
      # variable causing largest decrease in MET 'wins' and is added to new model
      win[step_num]<-which(METS==min(METS, na.rm=TRUE))[1]
      if ((win[step_num])!=0)
      {
        improvement=model$value-m[[win[step_num]]]$value
        if (m[[win[step_num]]]$value<(model$value-(5000*fityrs))) # only if above threshold for inclusion
        {
          model<-m[[win[step_num]]]
          others=c(others,win[[step_num]])
          write.csv(Pfull, file=paste("presences_customized", step_num,i,cc,"csv",sep='.'))
          write.csv(cbind(model$par, others), file=paste("par_customized", step_num,i,cc,"csv",sep='.'))
        }
      }
      if (win[step_num]==0)
      {
        improvement=0
      }
    }
  }
  i=i+1
}

