library(methods)

args <- commandArgs(trailingOnly = TRUE)

sim_num <- as.numeric(args[1])

iter=as.numeric(system("echo $SGE_TASK_ID",intern=T))

setClass("multiTS", representation(numTS = "numeric", tps = "list", data = "list",mat_tps = "matrix",mat_data = "matrix" , mat_data_2 = "matrix", mat_data_3 = "matrix", clusters = "numeric" , shifts = "numeric" , dilations = "numeric",subj_id = "vector",cohort = "character"))

setClass("TCsimRes", slots = c(params = "vector",shifts = "data.frame",clusters = "data.frame",loss="data.frame"))

loss_traj_arg <- function(traj_args,data,tps,lambda,traj_fun,traj_range){
  
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
  
    store <- mat.or.vec(num_ts,1)
  
    for (i in 1:num_ts){
  
      tmp_tmp <- localLoss(delta[i],data[i,],tps[i,],lambda,traj_fun,traj_args)
    
      store[i] <- tmp_tmp
    
      tmp <- tmp + tmp_tmp
  
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
  
  ss <- sum( (x - FUN(t + delta, fun_args))^2 )
  
  cost <- ss + lambda*delta^2
  
  return(cost)
  
}

psi <- function(t,fun_args){fun_args[1]*t + fun_args[2]}

phi <- function(t,fun_args){fun_args[1]*t^2 + fun_args[2]*t + fun_args[3]}

omega <- function(t,fun_args){fun_args[1]*t^3 + fun_args[2]*t^2 + fun_args[3]*t + fun_args[4]}

logistic_shift <- function(t,fun_args){fun_args[1]/(1 + exp(-fun_args[2]*(t - fun_args[3])))}

exp_drop3 <- function(t,fun_args){fun_args[1] - exp(fun_args[2]*t + fun_args[3])}

exp_drop2 <- function(t,fun_args){fun_args[1] - exp(fun_args[2]*t )}

exp_drop1 <- function(t,fun_args){30 - exp(fun_args[1]*t )}

logistic4 <- function(t,fun_args){fun_args[4]+fun_args[1]/(1 + exp(-fun_args[2]*(t - fun_args[3])))}

logistic3 <- function(t,fun_args){fun_args[3]+fun_args[1]/(1 + exp(-fun_args[2]*(t)))}

plot.multiTS <- function(sim_data,cluster_num=0,yl=c(0,1),first_x=0,shift=0,xlabel="Time",ylabel="Quantitative measure",recolor=F){
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
    
  } else {

    start_end <- mat.or.vec(num_ts,2)
  
    for (i in 1:num_ts){
     
      start_end[i,1] <- sim_data@tps[[i]][1]+shift[i]
      
      start_end[i,2] <- sim_data@tps[[i]][length(sim_data@tps[[i]])]+shift[i]
      
    }
    
    if (first_x == 0){
      
      xl <- c(min(start_end[,1]),max(start_end[,2]))
  
    } else {
      
      xl <- c(first_x,max(start_end[,2]))
      
    }
      
  }
  
  if (cluster_num==0){
    
    plot(sim_data@tps[[1]]+shift[1],sim_data@data[[1]],xlim=xl,ylim=yl,type="b",pch=c(rep(1,length(sim_data@tps[[1]])-1),19),col=sim_data@clusters[[1]],xlab=xlabel,ylab=ylabel)
    
    for (i in 2:num_ts){
      
      lines(sim_data@tps[[i]]+shift[i],sim_data@data[[i]],xlim=xl,ylim=yl,type="b",pch=c(rep(1,length(sim_data@tps[[i]])-1),19),col=sim_data@clusters[[i]])
      
    }
    
  } else {
    
    ind <- which(sim_data@clusters==cluster_num)
    
    if (recolor){
     
      plot(sim_data@tps[[ind[1]]]+shift[ind[1]],sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="b",pch=c(rep(1,length(sim_data@tps[[ind[1]]])-1),19),col=1,xlab=xlabel,ylab=ylabel)
    
      for (i in 2:length(ind)){
      
        lines(sim_data@tps[[ind[i]]]+shift[ind[i]],sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="b",pch=c(rep(1,length(sim_data@tps[[ind[i]]])-1),19),col=i)
      
      }      
       
    } else {
      
      plot(sim_data@tps[[ind[1]]]+shift[ind[1]],sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="b",pch=c(rep(1,length(sim_data@tps[[ind[1]]])-1),19),col=sim_data@clusters[ind[1]],xlab=xlabel,ylab=ylabel)
    
      for (i in 2:length(ind)){
      
        lines(sim_data@tps[[ind[i]]]+shift[ind[i]],sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="b",pch=c(rep(1,length(sim_data@tps[[ind[i]]])-1),19),col=sim_data@clusters[ind[i]])
      
      }
      
    }
    
  }
  
}

simulateExp <- function(tps,num_time_series=421,shift_range=c(-3000,1000),exp_drop_par=array(c(30,0.0005),dim=c(1,2)),num_clusters=2){
  
  sim_data <- simulateData2(num_time_series = num_time_series,time_of_visits =tps,exp_drop_par = exp_drop_par,num_clusters = num_clusters,shift_range=shift_range,noise = 1.5)
  
  #plot.multiTS(sim_data,yl=c(0,30))
  
  sim <- multiTS.matrix(sim_data)
  
  ind <- which(sim@mat_data<0,arr.ind=T)
  
  sim@mat_data[ind] <- 0
  sim@mat_data_2[ind] <- 0
  sim@mat_data_3[ind] <- 0
  
  sim@mat_data <- round(sim@mat_data)
  
  sim@mat_data_2 <- sim@mat_data
  
  ind2 <- which(sim@mat_data>30,arr.ind=T)
  
  sim@mat_data[ind2] <- 30
  
  return(sim)
  
}

simulateTP <- function(num_time_series = 10,poisson_gamma = 5,geom_p = 0.005,extra_tp=FALSE){
  # Function to simulate time points. Based on Poisson distribution for number of time points, and geometric distribution for spacing between time points
  
  # Simulate the number of time points in each time series
  num_tps <- 2 + rpois(num_time_series,poisson_gamma)
  
  if (extra_tp){
    
    num_tps <- num_tps + 1
    
  }
  
  time_of_visits <- vector("list", num_time_series)
  
  # Simulate the visit 'dates'
  for (i in 1:num_time_series){
    
    tmp_times <- mat.or.vec(num_tps[i],1)
    
    if (num_tps[i]>1){
      
      for (j in 1:(num_tps[i]-1)){
        
        tmp_times[j + 1] <- tmp_times[j] + rgeom(1,geom_p) + 1 
        
      }	
      
    } 
    
    time_of_visits[[i]] <- tmp_times
    
  }
  
  return(time_of_visits)
  
}

simulateData2 <- function(num_time_series,time_of_visits,exp_drop_par=array(c(30,0.0005),dim=c(1,2)),num_clusters = 2,noise=2,shift_range=c(-3000,1000),dilation_sd=1,trajectory_types = "sigmoid",curve_signs = "positive"){
  # Function to simulate data from linear, cubic or sigmoid trajectories for a given set of time points
  
  new_tps <- time_of_visits
  
  # assign time series to clusters
  
  cluster_id <- mat.or.vec(num_time_series,1)
  
  cluster_id <- sample(num_clusters,num_time_series,replace=TRUE)

  
  # generate data
  #
  
  data <- vector("list", num_time_series)
  
  time_shifts <- numeric(num_time_series)
  
  #dilations <- 2^(rnorm(num_time_series,0,dilation_sd))
  
  #last_time_point <- max(unlist(time_of_visits))
  
  #exp_drop_par <- mat.or.vec(num_clusters,2)
  
  #if (num_clusters==1){
    
    #exp_drop_par[1,] <- c(30,0.0005)
    
  #} else {
    
  #  exp_drop_par[1,] <- c(30,0.0007)
  #  exp_drop_par[2,] <- c(29,0.0004)
      
  #}
  
  length_ts <- unlist(lapply(time_of_visits,max))
  
  available <- logical(num_time_series)
  
  available[]<- T
  
  for (i in 1:num_clusters){
      
    cluster_member <- which(cluster_id==i)
    
    time_shifts[cluster_member] <- sample(shift_range[2*i - 1]:shift_range[2*i],length(cluster_member),replace=TRUE)
      
    for (j in cluster_member){
      
      #print(paste(i,j,time_shifts[j]))
      
      TOD <- (log(exp_drop_par[i,1]+1)/exp_drop_par[i,2]) - time_shifts[j]
      
      shorter_tod <- length_ts < TOD
      
      available_shorts <- shorter_tod & available
      
      #print(length(which(shorter_tod > 0)))
      
      #ßprint(length(which(available_shorts > 0)))
      
      if (length(which(available_shorts > 0))) {
      
        ts <- sample((1:num_time_series)[available_shorts],1)
        
      } else {
      
        ts <- sample((1:num_time_series)[shorter_tod],1)
          
      }
      
      #print(ts)
      
      error <- rnorm(length(time_of_visits[[ts]]),sd=noise)

      available[ts] <- F
      
      
      data[[j]] <- exp_drop2(time_of_visits[[ts]]+time_shifts[j],exp_drop_par[i,])+error
      
      new_tps[[j]] <- time_of_visits[[ts]]
        
    }
      
  }
  
  sim_data <- new("multiTS", numTS = num_time_series , clusters = cluster_id , shifts = time_shifts , tps = new_tps , data = data)
  
  return(sim_data)
  
}

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
  
  matrices <- list(matrix_tps,matrix_data,matrix_data,matrix_data)
  
  multiTS@mat_tps <- matrix_tps
  multiTS@mat_data <- matrix_data
  multiTS@mat_data_2 <- matrix_data
  multiTS@mat_data_3 <- matrix_data
  
  return(multiTS)
  
}


multiTS.matrixPlot <- function(sim_matrix){

  plot(sim_matrix@mat_tps[1,],sim_matrix@mat_data[1,],type="l",ylim=c(0,30),xlim=c(0,8000))

  for (i in 2:dim(sim_matrix@mat_tps)[1]){
    lines(sim_matrix@mat_tps[i,],sim_matrix@mat_data[i,],col=i)
  }
  
}

temporalClustering <- function(data,tps,K=2,lambda=0,phi=exp_drop2,phi_init=c(28,0.0004),max_iter=30,traj_range=c(-50000,50000),fix_first=T,phi_range=c(0,0.005),max_theta1 = 30){
  
  if (fix_first){
    
    k1_fit <- optim(phi_init,loss_traj_arg,data=data,tps=tps,lambda=lambda,traj_fun=phi,traj_range=traj_range)
    
    print(k1_fit$par[1])
    
    k1_theta1 <- min(k1_fit$par[1],max_theta1)
    
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
      
        fits[[k]] <- optimise(loss_traj_arg,interval=phi_range,data=data[cluster == k,],tps=tps[cluster == k,],lambda=lambda,traj_fun=phi_new,traj_range=traj_range,tol=sqrt(.Machine$double.eps))
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

discrimK2 <- function(data,tps,result,traj_range){
  # idea: plot difference in local loss vs length of time series (or number of TP)
  
  num_ts <- dim(tps)[1]
  
  par1 <- result[[2]][[1]]
  par2 <- result[[2]][[2]]
  
  delta1 <- optimal_offset(data,tps,lambda=0,traj_fun=exp_drop2,traj_args=par1,traj_range=traj_range)
  delta2 <- optimal_offset(data,tps,lambda=0,traj_fun=exp_drop2,traj_args=par2,traj_range=traj_range)
  
  local_loss <- data.frame(array(dim=c(num_ts,3)))
  colnames(local_loss) <- c("assigned","unassigned","abs_diff")
  
  num_tp <- c()
  
  for (i in 1:num_ts){

      num_tp[i] <- length(which(!is.na(tps[i,])))
    
      local_loss[i,1] <- localLoss(delta1[i],data[i,],tps[i,],lambda=0,FUN= exp_drop2,fun_args=par1)
      local_loss[i,2] <- localLoss(delta2[i],data[i,],tps[i,],lambda=0,FUN= exp_drop2,fun_args=par2)
      local_loss[i,3] <- abs(local_loss[i,result[[1]][i]] - local_loss[i,3-result[[1]][i]])
    
  }
  

  
  #length(which(local_loss[,3]>25))
  
  #aac_discrim25_ind <- which(local_loss[,3]>25)

  
  #load("mmse_data.RData")
  
  #aac_discrim50_subj <- data.frame(mmse_cpath_adni_aibl@subj_id[aac_discrim50_ind],mmse_cpath_adni_aibl@cohort[aac_discrim50_ind])
  #aac_discrim25_subj <- data.frame(mmse_cpath_adni_aibl@subj_id[aac_discrim25_ind],mmse_cpath_adni_aibl@cohort[aac_discrim25_ind])
  
  #save(aac_discrim25_subj,aac_discrim50_subj,file="aac_discrim.RData")
  
  #pdf("discrim_plot.pdf",width=5,height=5)
  #par(mar=c(4,4,0.5,0.5))
  
  #plot(aac_tps[ind[1],]+res_aac[[3]][ind[1]]-(2-res_aac[[1]][ind[1]])*5650,aac_data[ind[1],],col=1,type="l",xlim=c(0,12000),ylim=c(0,30),xlab="Aligned time in days",ylab="MMSE")
  
  #for (i in 2:length(aac_discrim25_ind)){
    
   # print(i)
    
   # lines(aac_tps[aac_discrim25_ind[i],]+res_aac[[3]][aac_discrim25_ind[i]]-(2-res_aac[[1]][aac_discrim25_ind[i]])*5650,aac_data[aac_discrim25_ind[i],],col=i)
  
  #}
  #  lines(0:5000,exp_drop2(0:5000,par2),lwd=5)
  #lines((5650:18000)-5650,exp_drop2(5650:18000,par1),col="red",lwd=5)
  
  #legend("topright",c("AAC k=1","AAC k=2"),lty=1,lwd=5,col=1:2)
  #dev.off()
  
  return(local_loss)
  
}


load("mmse.RData")



thetas <- array(c(29,29,29,29,29,0.0005,0.0003,0.0005,0.0006,0.0008,29,29,29,29,29,0.0005,0.0001,0.0001,0.001,0.001),dim=c(5,4))

# PAQ
#shift_ranges <- data.frame(min1=c(-4000,-2000,-4000,-4000,-4000),max1=c(4000,8000,4000,4000,3000),min2=c(-4000,0,0,-6000,-6000),max2=c(4000,28000,28000,2000,2000))
#mmse <- mmse_paq

# AAC
shift_ranges <- ranges <- data.frame(min1=c(0,0,0,0,0),max1=c(6000,10000,6000,5000,3500),min2=c(0,10000,10000,-1500,-1500),max2=c(6000,30000,30000,3000,3000))
mmse <- mmse_aac

ranges <- 10^4 * c(1,2,5)

results <- vector("list",9)

sim <- simulateExp(tps = mmse@tps,num_time_series=2412,exp_drop_par = t(array(thetas[sim_num,],dim=c(2,2))),shift_range=as.numeric(shift_ranges[sim_num,]),num_clusters=2)


range <- 3

  # rounded, ceiling
  tmp <-  temporalClustering(sim@mat_data,sim@mat_tps,K=2,lambda=0,phi=exp_drop2,phi_init=c(28,0.0004),max_iter=30,traj_range=c(-ranges[range],ranges[range]))
  
  simRes <- new("TCsimRes", params = c(tmp[[2]][[1]],tmp[[2]][[2]]), shifts = data.frame(real_shift=sim@shifts,est_shift=tmp[[3]]), clusters = data.frame(real_clusters=sim@clusters,est_clusters=tmp[[1]]),loss = discrimK2(data = sim@mat_data,tps = sim@mat_tps,result = tmp,traj_range=c(-ranges[range],ranges[range])))



if (!file.exists("aac_sim")){dir.create("aac_sim")}

setwd("aac_sim")

folder_name <- paste("theta",sim_num,sep="")

if (!file.exists(folder_name)){dir.create(folder_name)}

setwd(folder_name)

file_name <- paste("sim",iter,".RData",sep="")
  
save(simRes,file=file_name)

