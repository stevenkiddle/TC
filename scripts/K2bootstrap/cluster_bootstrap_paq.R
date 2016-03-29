
library(methods)

library(gtools)

num=as.numeric(system("echo $SGE_TASK_ID",intern=T))


setClass("multiTS", representation(numTS = "numeric", tps = "list", data = "list",mat_tps = "matrix",mat_data = "matrix" , mat_data_2 = "matrix", clusters = "numeric" , shifts = "numeric" , dilations = "numeric",subj_id = "vector",cohort = "character"))

mergeMultiTS <- function(multiTS1,multiTS2){
# requires raw data, i.e. not interpolated or filtered
  
    cohort1full <- (length(multiTS1@cohort) == multiTS1@numTS)

    cohort2full <- (length(multiTS2@cohort) == multiTS2@numTS)
  
    tmp_cohort <- c()
  
    if (cohort1full & cohort2full){
      
      tmp_cohort <- c(multiTS1@cohort,multiTS2@cohort)
    
    } else if (!cohort1full & cohort2full){
    
      tmp_cohort <- c(rep(multiTS1@cohort,multiTS1@numTS),multiTS2@cohort)
      
    } else if (cohort1full & !cohort2full){
    
      tmp_cohort <- c(multiTS1@cohort,rep(multiTS2@cohort,multiTS2@numTS))
      
    } else if (!cohort1full & !cohort2full){
    
      tmp_cohort <- c(rep(multiTS1@cohort,multiTS1@numTS),rep(multiTS2@cohort,multiTS2@numTS))
      
    }
  
    multiTS3 <- new("multiTS", numTS = (multiTS1@numTS + multiTS2@numTS) , tps = c(multiTS1@tps,multiTS2@tps) , data = c(multiTS1@data,multiTS2@data), clusters = c(multiTS1@clusters,multiTS2@clusters), shifts = c(multiTS1@shifts,multiTS2@shifts), dilations = c(multiTS1@dilations,multiTS2@dilations), cohort = tmp_cohort , subj_id = c(multiTS1@subj_id,multiTS2@subj_id)) 
 
}

plot.multiTS <- function(sim_data,cluster_num=0,yl=c(0,1),xl=c(0,2000),first_x=0,shift=0,xlabel="Time",ylabel="Quantitative measure",recolor=F){
  # Function to plot multiTS data object, only plots raw data, i.e. 'tps' and 'data'
  # 
  # Inputs : 
  #
  # sim_data = multiTS object
  # cluster_num = optional argument, provide an integer and plot will only show the data for that cluster number
  # yl = y-axis limits
  # shift = optional parameter, a vector of shifts for each time series, to look at alignment
  # xlabel = label of the x-axis
  # ylabel = label of the y-xis
  
  last_time_point <- max(unlist(sim_data@tps))
  
  num_ts <- length(sim_data@tps)
  
  if (sum(shift==0)==1){
    
    shift <- rep(0,num_ts)
    
    xl <- c(first_x,last_time_point)
    
  } 
      
  
  if (cluster_num==0){
    
    cols <- sim_data@clusters
    
    if (recolor == T){
    
      cols <- 1:num_ts
        
    }
    
    plot(sim_data@tps[[1]]+shift[1],sim_data@data[[1]],xlim=xl,ylim=yl,type="l",col=cols[1],xlab=xlabel,ylab=ylabel)
    
    for (i in 2:num_ts){
      
      lines(sim_data@tps[[i]]+shift[i],sim_data@data[[i]],xlim=xl,ylim=yl,type="l",col=cols[i])
      
    }
    
  } else {
    
    ind <- which(sim_data@clusters==cluster_num)
    
    if (recolor){
     
      plot(sim_data@tps[[ind[1]]]+shift[ind[1]],sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="l",col=1,xlab=xlabel,ylab=ylabel)
    
      for (i in 2:length(ind)){
      
        lines(sim_data@tps[[ind[i]]]+shift[ind[i]],sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="l",col=i)
      
      }      
       
    } else {
      
      plot(sim_data@tps[[ind[1]]]+shift[ind[1]],sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="l",col=sim_data@clusters[ind[1]],xlab=xlabel,ylab=ylabel)
    
      for (i in 2:length(ind)){
      
        lines(sim_data@tps[[ind[i]]]+shift[ind[i]],sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="l",col=sim_data@clusters[ind[i]])
      
      }
      
    }
    
  }
  
}

loss_traj_arg <- function(traj_args,data,tps,lambda,traj_fun,traj_range){
  # calculate the loss for a given set of parameters
  
  delta <- optimal_offset(data,tps,lambda,traj_fun,traj_args,traj_range)
  
  result <- loss(data,tps,delta,lambda,traj_fun,traj_args)
      
  return(result)
  
}

optimal_offset <- function(data,tps,lambda,traj_fun,traj_args,traj_range,alpha=2){
# Estimate delta for given fun_args, bare in mind, may be local rather than global minimum if not unimodal
  
  if (length(dim(tps)[1])){
    
    num_ts <- dim(tps)[1]
    
    delta <- mat.or.vec(num_ts,1)
  
    for (i in 1:num_ts){
  
      tmp <- optimise(f= localLoss,interval=traj_range,x=data[i,],t=tps[i,],lambda=lambda,FUN = traj_fun,fun_args = traj_args)
    
      delta[i] <- tmp$minimum
    
    }
  
  } else {
    
    tmp <- optimise(f= localLoss,interval=traj_range,x=data,t=tps,lambda=lambda,FUN = traj_fun,fun_args = traj_args)
    
    delta <- tmp$minimum
      
  }
  
  return(delta)
  
}

loss <- function(data,tps,delta,lambda,traj_fun,traj_args){
# cost of time shifted fit to given trajectory for given fun_args and delta
  
  if (length(dim(tps)[1])){
    
    num_ts <- dim(tps)[1]
  
    tmp <- 0
  
    for (i in 1:num_ts){
    
      tmp <- tmp + localLoss(delta[i],data[i,],tps[i,],lambda,traj_fun,traj_args)
  
    }
  
  } else {
    
    tmp <- localLoss(delta,data,tps,lambda,traj_fun,traj_args)
    
  }
  
  return(tmp)
  
}

localLoss <- function(delta,x,t,lambda,FUN,fun_args){
# cost of time shifted fit to FUN
  
  x <- x[!is.na(x)]
  t <- t[!is.na(t)]
  
  cost <- sum( (x - FUN(t + delta, fun_args))^2 ) + lambda*delta^2
  
  return(cost)
  
}

line <- function(t,fun_args){fun_args[1]*t + fun_args[2]}

quadratic <- function(t,fun_args){fun_args[1]*t^2 + fun_args[2]*t + fun_args[3]}

cubic <- function(t,fun_args){fun_args[1]*t^3 + fun_args[2]*t^2 + fun_args[3]*t + fun_args[4]}

exp_drop2 <- function(t,fun_args){fun_args[1] - exp(fun_args[2]*t )}

exp_drop1 <- function(t,fun_args,theta1){theta1 - exp(fun_args[1]*t )}

logistic3 <- function(t,fun_args){fun_args[3]+fun_args[1]/(1 + exp(-fun_args[2]*(t)))}

multiTS.matrix <- function(multiTS){

  num_ts <- multiTS@numTS
  
  num_tp <- mat.or.vec(num_ts,1)
  
  for (i in 1:num_ts){
    
    num_tp[i] <- length(multiTS@tps[[i]])
      
  }
  
  max_tps <- max(num_tp)
  
  matrix_tps <- array(dim=c(num_ts,max_tps))
  matrix_data <- array(dim=c(num_ts,max_tps))
  
  for (i in 1:num_ts){
    
    matrix_tps[i,1:num_tp[i]] <- multiTS@tps[[i]]
    matrix_data[i,1:num_tp[i]] <- multiTS@data[[i]]
      
  }
  
  matrices <- list(matrix_tps,matrix_data)
  
  return(matrices)
  
}

temporalClustering2 <- function(data,tps,K=2,lambda=0,phi=exp_drop2,phi_init=c(28,0.0004),min_iter=10,max_iter=30,traj_range=c(-50000,50000),junk=T){
  
  cluster <- sample(K,dim(tps)[1],replace = T) # randomly assign to clusters
  
  last_tp <- apply(tps,1,max,na.rm=T)
  
  while(!all(table(cluster)>0)) cluster <- sample(K,dim(tps)[1],replace = T) # makes sure at least one time series is assigned to each cluster

  for (iter in 1:max_iter){
  
    print(iter)
  
    cluster_old <- cluster
  
    fits <- vector("list",K)
    deltas <- vector("list",K)
    
    # Calculate means
    for (k in 1:K){
    
      fits[[k]] <- optim(phi_init,loss_traj_arg,data=data[(cluster == k) | (cluster == 0),],tps=tps[(cluster == k) | (cluster == 0),],lambda=lambda,traj_fun=phi,traj_range=traj_range)
      deltas[[k]] <- optimal_offset(data=data,tps=tps,lambda=lambda,traj_fun=phi,traj_arg=fits[[k]]$par,traj_range=traj_range)
  
    }
   
    shift_last <- numeric(dim(tps)[1])
    
    # Re-assign time series to clusters
    for (i in 1:dim(tps)[1]){
      
      losses <- array(dim = c(K,1))
      
      for (k in 1:K){
    
        losses[k] <- localLoss(deltas[[k]][i],data[i,],tps[i,],lambda=lambda,FUN =phi,fun_args=fits[[k]]$par)
      
      }
  
      cluster[i] <- which.min(losses) # reassign to the best fitting cluster model
      
      if (junk){
      
        shift_last <- last_tp[i] + deltas[[cluster[i]]][i]
        
        if (shift_last < 0){
        
          cluster[i] <- 0
            
        }
          
      }
    
    }
  
    #plot(traj_range[1]:traj_range[2],phi(traj_range[1]:traj_range[2],fits[[1]]$par),ylab="MMSE",xlab="Time in days",ylim=c(0,30))
    #for (k in 1:K){
      
      #points(traj_range[1]:traj_range[2],phi(traj_range[1]:traj_range[2],fits[[k]]$par),ylab="MMSE",xlab="Time",col=k)
    
    #}
  
    #print(table(cluster))
    
    if (all(cluster[cluster != 0] == cluster_old[cluster != 0]) & iter >= min_iter){break}
    
    if (iter < max_iter) {cluster[cluster == 0] <- sample(K,length(which(cluster==0)),replace = T)}

  }
  
  result <- list()
  
  result[[1]] <- cluster
  
  result[[2]] <- fits
  
  final_delta <- deltas[[1]]
  
  for (k in 2:K){
    
    final_delta[cluster==k] <- deltas[[k]][cluster==k]

  }
  
  result[[3]] <- final_delta
  
  return(result)
  
}

temporalClustering <- function(data,tps,K=2,lambda=0,phi=exp_drop2,phi_init=c(28,0.0004),max_iter=30,traj_range=c(-50000,50000),fix_first=T,phi_range=c(0,0.005)){
  
  if (fix_first){
    
    k1_fit <- optim(phi_init,loss_traj_arg,data=data,tps=tps,lambda=lambda,traj_fun=phi,traj_range=traj_range)
    
    k1_theta1 <- k1_fit$par[1]
    
    phi_new <- function(t,args){phi(t,c(k1_theta1,args))}
      
  }
    
  cluster <- sample(K,dim(tps)[1],replace = T) # randomly assign to clusters
  
  while(!all(table(cluster)>0)) cluster <- sample(K,dim(tps)[1],replace = T) # makes sure at least one time series is assigned to each cluster

  for (iter in 1:max_iter){
  
    print(iter)
  
    cluster_old <- cluster
  
    fits <- vector("list",K)
    deltas <- vector("list",K)
    args <- vector("list",K)
    
    # Calculate means
    for (k in 1:K){
    
      if (fix_first){
      
        fits[[k]] <- optimise(loss_traj_arg,interval=phi_range,data=data[cluster == k,],tps=tps[cluster == k,],lambda=lambda,traj_fun=phi_new,traj_range=traj_range,tol = sqrt(.Machine$double.eps))
        args[[k]] <- c(k1_theta1,fits[[k]]$minimum)
        
      } else {
        
        fits[[k]] <- optim(phi_init,loss_traj_arg,data=data[cluster == k,],tps=tps[cluster == k,],lambda=lambda,traj_fun=phi,traj_range=traj_range)
        args[[k]] <- fits[[k]]$par
  
      }
      
      deltas[[k]] <- optimal_offset(data=data,tps=tps,lambda=lambda,traj_fun=phi,traj_arg=args[[k]],traj_range=traj_range)
      
  
    }
   

    
    # Re-assign time series to clusters
    for (i in 1:dim(tps)[1]){
      
      losses <- array(dim = c(K,1))
      
      for (k in 1:K){
    
        losses[k] <- localLoss(deltas[[k]][i],data[i,],tps[i,],lambda=lambda,FUN =phi,fun_args=args[[k]])
      
      }
  
      cluster[i] <- which.min(losses) # reassign to the best fitting cluster model
    
    }
  
    #plot(traj_range[1]:traj_range[2],phi(traj_range[1]:traj_range[2],fits[[1]]$par),ylab="MMSE",xlab="Time in days",ylim=c(0,30))
    #for (k in 1:K){
      
      #points(traj_range[1]:traj_range[2],phi(traj_range[1]:traj_range[2],fits[[k]]$par),ylab="MMSE",xlab="Time",col=k)
    
    #}
  
    #print(table(cluster))
    
    if (all(cluster == cluster_old)){break}

  }
  
  result <- list()
  
  result <- vector("list",4)
  
  result[[1]] <- cluster
  
  result[[2]] <- args
  
  result[[4]] <- fits
  
  final_delta <- deltas[[1]]
  
  for (k in 2:K){
    
    final_delta[cluster==k] <- deltas[[k]][cluster==k]

  }
  
  result[[3]] <- final_delta
  
  return(result)
  
}








tempClusterBoot <- function(data=paq_data,tps=paq_tps,numSamples=3,phi=exp_drop2,phi_init=c(28,0.0004),K=2,lambda=0,traj_range=c(-50000,50000),max_iter=30,cohort,fix_first=F){

  num_ts <- dim(tps)[1]
  
  pars <- array(dim=c(numSamples*K,length(phi_init)+1))
  
  pars <- data.frame(pars)
  
  names(pars) <- c("k","theta1","theta2")
  
  cohorts <- names(table(cohort))
  
  size_cohort <- as.numeric(table(cohort))
  
  if (length(cohorts) == 1){size_cohort<-num_ts}

  for (i in 1:numSamples){
    
    print(i)
    
    bootstrapSample <- c()
    
    for (j in 1:length(cohorts)){
    
      ind <- which(cohort == cohorts[j])
      
      bootstrapSample <- c(bootstrapSample,sample(ind,size_cohort[j],replace = T))
      
      #print(length(bootstrapSample))
        
    }
    
    bootstrap_res <- temporalClustering(data[bootstrapSample,],tps[bootstrapSample,],K=2,fix_first=fix_first)
    
    for (k in 1:K){
      
      pars[(i-1)*K+k,] <- c(k,bootstrap_res[[2]][[k]])
      
    }
    
    
    
  }
  
  return(pars)

}


load("mmse.RData")

paq_cluster_boot <- tempClusterBoot(data = mmse_paq@mat_data,tps=mmse_paq@mat_tps,numSamples = 5,cohort=rep(1,421),fix_first=F)

save(paq_cluster_boot,file=paste("paq_boot_1_",num,".RData",sep=""))

paq_cluster_boot <- tempClusterBoot(data = mmse_paq@mat_data,tps=mmse_paq@mat_tps,numSamples = 5,cohort=rep(1,421),fix_first=F)

save(paq_cluster_boot,file=paste("paq_boot_2_",num,".RData",sep=""))
