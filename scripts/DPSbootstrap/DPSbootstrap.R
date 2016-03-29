DPSbootstrap <- function(numSamples=3,phi=exp_drop2,phi_init=c(28,0.0004),data=paq_data,tps=paq_tps,lambda=0,traj_range=c(-50000,50000),cohort=cohort){

  num_ts <- dim(data)[1]
  
  if (length(cohort)==1){cohort <- rep(cohort,num_ts)}
  
  pars <- mat.or.vec(numSamples,length(phi_init))
  
  losses <- c()
  
  cohorts <- names(table(cohort))
  
  size_cohort <- as.numeric(table(cohort))

  for (i in 1:numSamples){
    
    print(i)
    
    bootstrapSample <- c()
    
    for (j in 1:length(cohorts)){
    
      ind <- which(cohort == cohorts[j])
      
      bootstrapSample <- c(bootstrapSample,sample(ind,size_cohort[j],replace = T))
      
      #print(length(bootstrapSample))
        
    }
    
    tmp <- optim(phi_init,loss_traj_arg,data=data[bootstrapSample,],tps=tps[bootstrapSample,],lambda=lambda,traj_fun=phi,traj_range=traj_range)

    pars[i,] <- tmp$par
    
    losses[i] <- tmp$value
    
  }
  
  results <- list()
  
  results$pars <- pars
  
  results$losses <- losses
  
  return(results)

}

paq_dps_boot <- DPSbootstrap(numSamples = 200,data=mmse_paq@mat_data,tps=mmse_paq@mat_tps,cohort=mmse_paq@cohort)

aac_dps_boot <- DPSbootstrap(numSamples = 200,data=mmse_aac@mat_data,tps=mmse_aac@mat_tps,cohort=mmse_aac@cohort)

print("PAQUID DPS bootstrap CI for theta 1")
print(quantile(paq_dps_boot$par[,1],0.025))
print(quantile(paq_dps_boot$par[,1],0.975))

print("PAQUID DPS bootstrap CI for theta 2")
print(quantile(paq_dps_boot$par[,2],0.025))
print(quantile(paq_dps_boot$par[,2],0.975))


print("AAC DPS bootstrap CI for theta 1")
print(quantile(aac_dps_boot$par[,1],0.025))
print(quantile(aac_dps_boot$par[,1],0.975))

print("AAC DPS bootstrap CI for theta ")
print(quantile(aac_dps_boot$par[,2],0.025))
print(quantile(aac_dps_boot$par[,2],0.975))