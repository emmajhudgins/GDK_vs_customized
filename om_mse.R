library(pdist)
om_mse<-function(pred, obs)
{
  data<-read.csv('data_minus5_july.csv')
  library(optmatch)
  all_dist<-(as.matrix(pdist(cbind(data$X_coord[pred], data$Y_coord[pred]), cbind(data$X_coord[obs], data$Y_coord[obs]))))^2
  if (nrow(all_dist)>ncol(all_dist))
  {
    all_dist<-t(all_dist)
    row.names(all_dist)<-1:length(obs)
    colnames(all_dist)<-(length(obs)+1):(length(obs)+length(pred))
    match<-pairmatch(all_dist)
    matched_obs<-as.character(match[1:length(obs)])
    matched_pred<-as.character(match[(length(obs)+1):(length(obs)+length(pred))])
    all_dist<-t(all_dist)
  }
  if (nrow(all_dist)<=ncol(all_dist))
  {
    row.names(all_dist)<-1:length(pred)
    colnames(all_dist)<-(length(pred)+1):(length(pred)+length(obs))
    match<-pairmatch(all_dist)
    matched_pred<-as.character(match[1:length(pred)])
    matched_obs<-as.character(match[(length(pred)+1):(length(pred)+length(obs))])
  }
  min_dist<-function(i)
  {
    if (is.na(matched_pred[i])==FALSE)
    {
      return (all_dist[i, which(matched_obs==matched_pred[i])])
    }
    if (is.na(matched_pred[i])==TRUE)
    {
      return (NA)
    }
  }
  mindist<-(tapply(1:length(matched_pred), 1:length(matched_pred), min_dist))
  for (i in 1:length(matched_pred))
  {
    if (is.na(matched_pred[i])==TRUE)
    {
      mindist=c(mindist,mean(all_dist[i,]))
    }
  }
  for (i in 1:length(matched_obs))
  {
    if (is.na(matched_obs[i])==TRUE)
    {
      mindist=c(mindist,mean(all_dist[,i]))
    }
  }
 return(mean(mindist, na.rm=T))
}