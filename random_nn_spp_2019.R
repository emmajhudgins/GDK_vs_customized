random_nn_spp<-function(spp)
{
  data<-read.csv('data_minus5_july.csv', stringsAsFactors = FALSE)
  data2<-read.csv('spdat_clean_gdk.csv', stringsAsFactors = FALSE)
  prez<-read.csv('prez_clean_gdk.csv')
  if (spp==49)
  {
    start=1950
    pred=2005
    cutoff=2000
    YEAR=50
    total_YEAR=55
    prez2<-read.csv('hwa_5yr_fitting_july.csv')[,3:12]
    fitting_yrs<<-which(colSums(prez2)!=0)
    pfulz<-read.csv('hwa_3yrs.csv')[,3]
  }
  if (spp==51)
  {
    start=1865
    pred=2005
    cutoff=2000
    YEAR=135
    total_YEAR=140
    prez2<-read.csv('gm_5yr_fitting_july.csv')[,3:29]
    fitting_yrs<<-which(colSums(prez2)!=0)
    pfulz<-read.csv('gm_3yrs.csv')[,3]
  }
  if (spp==54)
  {
    start= 1890
    pred=2005
    cutoff=2000
    YEAR=110
    total_YEAR=115
    prez2<-read.csv('bbd_5yr_fitting_july.csv')[,3:24]
    fitting_yrs<<-which(colSums(prez2)!=0)
    pfulz<-read.csv('bbd_3yrs.csv')[,3]
  }
  prez2<-cbind(prez2, pfulz)
  for (sppp in 1:ncol(prez2))
  {prez2[,sppp]<-c(intersect(prez[which(prez[,spp]!=0),spp], prez2[which(prez2[,sppp]!=0),sppp]), rep(0,(3372-length(intersect(prez[which(prez[,spp]!=0),spp], prez2[which(prez2[,sppp]!=0),sppp])))))}
  total_time<-total_YEAR/5
  ffff<-prez2[1:length(which(prez2[,total_YEAR/5]!=0)),total_YEAR/5]
  hhhh<-prez[1:length(which(prez[,spp]!=0)),spp]
  source('om_mse.R')
  random_nn<-rep(0,10000)
  create_rand<-function(i)
  {
    rr<-sample(x=hhhh, size=length(ffff), replace=FALSE)
    if (i%%1000==0)
    {
      cat(i/100,"% complete\n")
    }
    return(om_mse(rr,ffff))
  }
  random_nn<-unlist(lapply(1:10000, create_rand))
  return(random_nn)
}
