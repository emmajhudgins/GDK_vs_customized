r2_om_mse<-function(pred, obs, spp)
{
  source('om_mse.R')
  dist_obs<-om_mse(pred,obs)
  random_nn<-read.csv(paste0('random_nn_', spp,".csv"))[,1]
  return(1-(dist_obs/mean(random_nn)))
}