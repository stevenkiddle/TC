iter=as.numeric(system("echo $SGE_TASK_ID",intern=T))

library(methods)

setClass("multiTS", representation(numTS = "numeric", tps = "list", data = "list",mat_tps = "matrix",mat_data = "matrix" , mat_data_2 = "matrix", mat_data_3 = "matrix", clusters = "numeric" , shifts = "numeric" , dilations = "numeric",subj_id = "vector",cohort = "character"))

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
      
      #print(length(which(available_shorts > 0)))
      
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


load("mmse.RData")
#load("paquid.RData")

# dim 1: theta, dim 2 : bias sources, dim 3: range 10000,20000,50000, dim 4 : parameters
sim_thetas <- array(dim=c(5,3,3,2))

#multiTS1 <- mmse_paq
multiTS1 <- mmse_aac

print(paste(iter,"starting"))

ranges <- 10^4 * c(1,2,5)

# PAQUID
#theta_ranges <- data.frame(thetas=seq(1,10,2)  * 10^-4,min=c(0,-2000,-4000,-4000,-5000),max=c(28000,8000,4000,3000,3000))
# AAC
theta_ranges <- data.frame(thetas=seq(1,10,2)  * 10^-4,min=c(10000,0,0,-0,-1000),max=c(30000,10000,6000,4000,3500))


for (i1 in 1:5){
  
  sim <- simulateExp(tps = multiTS1@tps,num_time_series=2412,shift_range=as.numeric(theta_ranges[i1,2:3]),exp_drop_par =  array(c(29,theta_ranges[i1,1]),dim=c(1,2)),num_clusters = 1)
  
  for (i2 in 1:3){
    
    print(paste(i1,i2))
  
    tmp <- optim(c(28,0.0004),loss_traj_arg,data=sim@mat_data_3,tps=sim@mat_tps,lambda=0,traj_fun=exp_drop2,traj_range=c(-ranges[i2],ranges[i2]))
  
    sim_thetas[i1,1,i2,] <- tmp$par
  

    tmp <- optim(c(28,0.0004),loss_traj_arg,data=sim@mat_data_2,tps=sim@mat_tps,lambda=0,traj_fun=exp_drop2,traj_range=c(-ranges[i2],ranges[i2]))
  
    sim_thetas[i1,2,i2,] <- tmp$par


    tmp <- optim(c(28,0.0004),loss_traj_arg,data=sim@mat_data,tps=sim@mat_tps,lambda=0,traj_fun=exp_drop2,traj_range=c(-ranges[i2],ranges[i2]))
  
    sim_thetas[i1,3,i2,] <- tmp$par
  
  }
  
}

file_name <- paste("aac_sim/sim",iter,".RData",sep="")
  
save(sim_thetas,file=file_name)

