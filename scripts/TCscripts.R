# By Steven Kiddle with help from Chris Wallace
# steven.kiddle@kcl.ac.uk
# http://core.brc.iop.kcl.ac.uk/


# S4 class 'multiTS' (for storing multiple time series, with relevant information)
setClass("multiTS", representation(numTS = "numeric", tps = "list", data = "list",mat_tps = "matrix",mat_data = "matrix" , mat_data_2 = "matrix", clusters = "numeric" , shifts = "numeric" , dilations = "numeric",subj_id = "vector",cohort = "character"))

# numTS = number of time series
# tps = list containing numeric vectors, 
#       each of which contains the time points of an individual (i.e. tau_i)
# data = list containing numeric vectors, each of which contains the longitudinal
#        univariate data for an individual (i.e. x_i(tau_i) )
# mat_tps = matrix representation of tps, columns represent time points, 
#           NAs are present if individual has fewer time ponts than the
#           number of columns in the matrix
# clusters = numeric vector containing cluster labels (integers), i.e. c(i)
# offsets = numeric vector containing individual offsets, i.e. delta_i
# dilations = unused option to store dilations (individual level rates)
# subj_id = participant identifiers, may not be unique without comparing to cohort labels
# cohort = cohort labels


mergeMultiTS <- function(multiTS1,multiTS2){
# Function for merging two 'multiTS' objects (i.e. multiTS1 and multiTS2)
  
    # boolean variables, is cohort vector the length of the number of time series?
    # can be shorter if all individuals are from the same cohort, i.e. just "cohort1"
    cohort1full <- (length(multiTS1@cohort) == multiTS1@numTS)

    cohort2full <- (length(multiTS2@cohort) == multiTS2@numTS)
  
    tmp_cohort <- c()
  
    # code to make correct combined cohort vector in light of the above
    if (cohort1full & cohort2full){
      
      tmp_cohort <- c(multiTS1@cohort,multiTS2@cohort)
    
    } else if (!cohort1full & cohort2full){
    
      tmp_cohort <- c(rep(multiTS1@cohort,multiTS1@numTS),multiTS2@cohort)
      
    } else if (cohort1full & !cohort2full){
    
      tmp_cohort <- c(multiTS1@cohort,rep(multiTS2@cohort,multiTS2@numTS))
      
    } else if (!cohort1full & !cohort2full){
    
      tmp_cohort <- c(rep(multiTS1@cohort,multiTS1@numTS),rep(multiTS2@cohort,multiTS2@numTS))
      
    }
    
    # Create and return combined 'multiTS'
    multiTS3 <- new("multiTS", numTS = (multiTS1@numTS + multiTS2@numTS) , tps = c(multiTS1@tps,multiTS2@tps) , data = c(multiTS1@data,multiTS2@data), clusters = c(multiTS1@clusters,multiTS2@clusters), shifts = c(multiTS1@shifts,multiTS2@shifts), dilations = c(multiTS1@dilations,multiTS2@dilations), cohort = tmp_cohort , subj_id = c(multiTS1@subj_id,multiTS2@subj_id)) 
 
}

plot.multiTS <- function(sim_data,cluster_num=0,yl=c(0,1),xl=c(0,2000),first_x=0,shift=0,
                         xlabel="Time",ylabel="Quantitative measure",recolor=F,ltype=1,days2years=F){
  # Function to plot multiTS data object, i.e. 'tps' and 'data'. Used to generate Kiddle et al., (2016) Figures 2c-d.
  # 
  # Inputs : 
  #
  # sim_data = 'multiTS' object, names is legacy, doesn't need to be simulated.
  # cluster_num = optional argument, provide an integer and plot will only show the data for that cluster number
  # yl = y-axis limits, otherwise will show between 0 and 1
  # xl = x-axis limits, only becomes active when shifts are given
  # first_x = if shifts aren't given, what should be the first x? (defualts to zero)
  # shift = optional parameter, a vector of shifts for each time series, can use offsets (delta_i)
  # xlabel = label of the x-axis
  # ylabel = label of the y-xis
  # recolor = set to TRUE if you don't want individuals to be coloured by their cluster label
  # ltype = type of line, 1 = solid, >1 dashed/dotted etc.
  # days2years = convert days to years
  
  last_time_point <- max(unlist(sim_data@tps))
  
  num_ts <- length(sim_data@tps)
  
  if (sum(shift==0)==1){
    
    shift <- rep(0,num_ts)
    
    xl <- c(first_x,last_time_point)
    
  } 
  
  if (ltype[1] == 1){ltype = rep(1,num_ts)}
      
  if (days2years) {xl <- xl/365.25}
  
  if (cluster_num==0){
    
    cols <- sim_data@clusters
    
    if (recolor == T){
    
      cols <- 1:num_ts
        
    }
  
    if (days2years) {
      
      plot((sim_data@tps[[1]]+shift[1])/365.25,sim_data@data[[1]],xlim=xl,ylim=yl,type="l",col=cols[1],xlab=xlabel,ylab=ylabel,lty=ltype[1])
    
    } else {
      
      plot(sim_data@tps[[1]]+shift[1],sim_data@data[[1]],xlim=xl,ylim=yl,type="l",col=cols[1],xlab=xlabel,ylab=ylabel,lty=ltype[1])
    
    }   
      
    for (i in 2:num_ts){
      
      if (days2years) {
      
        lines((sim_data@tps[[i]]+shift[i])/365.25,sim_data@data[[i]],xlim=xl,ylim=yl,type="l",col=cols[i],lty=ltype[i])
        
      } else {
        
        lines(sim_data@tps[[i]]+shift[i],sim_data@data[[i]],xlim=xl,ylim=yl,type="l",col=cols[i],lty=ltype[i])
      
      }
        
    }
    
  } else {
    
    ind <- which(sim_data@clusters==cluster_num)
    
    if (recolor){
      
      if (days2years) {
     
        plot((sim_data@tps[[ind[1]]]+shift[ind[1]])/365.25,sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="l",col=1,xlab=xlabel,ylab=ylabel,lty=ltype[ind[1]])
    
      } else {
        
        plot(sim_data@tps[[ind[1]]]+shift[ind[1]],sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="l",col=1,xlab=xlabel,ylab=ylabel,lty=ltype[ind[1]])
        
      }
        
        
      for (i in 2:length(ind)){
      
        if (days2years) {
          
          lines((sim_data@tps[[ind[i]]]+shift[ind[i]])/365.25,sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="l",col=i,lty=ltype[ind[i]])
      
        } else {
          
          lines(sim_data@tps[[ind[i]]]+shift[ind[i]],sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="l",col=i,lty=ltype[ind[i]])
      
        }
          
      }      
       
    } else {
      
      if (days2years) {
      
        plot((sim_data@tps[[ind[1]]]+shift[ind[1]])/365,25,sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="l",col=sim_data@clusters[ind[1]],xlab=xlabel,ylab=ylabel,lty=ltype[ind[1]])
    
      } else {
        
        plot(sim_data@tps[[ind[1]]]+shift[ind[1]],sim_data@data[[ind[1]]],xlim=xl,ylim=yl,type="l",col=sim_data@clusters[ind[1]],xlab=xlabel,ylab=ylabel,lty=ltype[ind[1]])
    
      }
        
      for (i in 2:length(ind)){
        
        if (days2years) {
      
          lines((sim_data@tps[[ind[i]]]+shift[ind[i]])/365.25,sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="l",col=sim_data@clusters[ind[i]],lty=ltype[ind[i]])
      
        } else {
          
          lines(sim_data@tps[[ind[i]]]+shift[ind[i]],sim_data@data[[ind[i]]],xlim=xl,ylim=yl,type="l",col=sim_data@clusters[ind[i]],lty=ltype[ind[i]])
      
        }    
          
      }
      
    }
    
  }
  
}

localLoss <- function(delta,x,t,lambda,FUN,fun_args){
# cost of time shifted (offset) fit to FUN (Eqn 1)
# 
# delta = individual offset
# x = numeric vector containing univariate observations (i.e. x_i(tau_i))
# t = numeric vector containing the time of observations (i.e. tau_i)
# lambda = legacy option, essentially a ridge regression penalty
# FUN = trajectory model (i.e. phi(t;theta))
# fun_args = trajectory parameters (i.e. theta)
  
  x <- x[!is.na(x)]
  t <- t[!is.na(t)]
  
  cost <- sum( (x - FUN(t + delta, fun_args))^2 ) + lambda*delta^2
  
  return(cost)
  
}

optimal_offset <- function(data,tps,lambda,traj_fun,traj_args,
                           traj_range,alpha=2){
# Estimate optimal delta for all individuals given fun_args (Eqn 2 applied in parellel), 
# bare in mind, may be local rather than global minimums if not unimodal
# 
# data = matrix format data, see 'multiTS' definition
# tps = matrix format timepoints, see 'multiTS' definition
# lambda = legacy option, essentially a ridge regression penalty
# traj_fun = trajectory model (i.e. phi(t;theta))
# traj_args = trajectory parameters (i.e. theta)
# traj_range = delta range, i.e. interval of allowed deltas
  
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
# Cost of fit to data, for given delta and theta (Eqn 3)
# 
# data = matrix format data, see 'multiTS' definition
# tps = matrix format timepoints, see 'multiTS' definition
# delta = vector of offsets
# lambda = legacy option, essentially a ridge regression penalty
# traj_fun = trajectory model (i.e. phi(t;theta))
# traj_args = trajectory parameters (i.e. theta)
  
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

loss_traj_arg <- function(traj_args,data,tps,lambda,traj_fun,traj_range){
# Cost of fit to data for given theta, including estimatation of delta (Eqn 3)
# 
# traj_args = trajectory parameters (i.e. theta)
# data = matrix format data, see 'multiTS' definition
# tps = matrix format timepoints, see 'multiTS' definition
# lambda = legacy option, essentially a ridge regression penalty
# traj_fun = trajectory model (i.e. phi(t;theta))
# traj_range = delta range, i.e. interval of allowed deltas
  
  delta <- optimal_offset(data,tps,lambda,traj_fun,traj_args,traj_range)
  
  result <- loss(data,tps,delta,lambda,traj_fun,traj_args)
      
  return(result)
  
}


# Small selection of trajectory models, other functions can be defined and used
line <- function(t,fun_args){fun_args[1]*t + fun_args[2]}
quadratic <- function(t,fun_args){fun_args[1]*t^2 + fun_args[2]*t + fun_args[3]}
cubic <- function(t,fun_args){fun_args[1]*t^3 + fun_args[2]*t^2 + fun_args[3]*t + fun_args[4]}
exp_drop2 <- function(t,fun_args){fun_args[1] - exp(fun_args[2]*t )}
logistic3 <- function(t,fun_args){fun_args[3]+fun_args[1]/(1 + exp(-fun_args[2]*(t)))}


multiTS.matrix <- function(multiTS){
  # Function to add matrix representation of tps and data to multiTS object.
  # Uses data and tps list data.
  #
  # multiTS = multiTS object

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

  
  return(multiTS)
  
}

temporalClustering <- function(data,tps,K=2,lambda=0,phi=exp_drop2,phi_init=c(28,0.0004),max_iter=30,
                               traj_range=c(-50000,50000),fix_first=T,phi_range=c(0,0.005)){
# Main analysis function, performs disease progression score (K=1) or Temporal Clustering (K>1)
#
# Inputs :
#
# data = matrix representation of data (columns = timepoints, rows = individuals)
# tps = matrix representation of timepoints (columns = timepoints, rows = individuals)
# K = number of clusters to return (default: 2). If K=1 disease progression score is applied.
# lambda = legacy option, ridge regression penalty on offsets (default: 0)
# phi = trajectory model, pre-implemented (line/quadratic/cubic/exp_drop2/logistic3), 
#       or user defined (default: exp_drop2)
# phi_init = initial values for phi parameters (theta), used in estimation of parameters
#            (default: c(28,0.0004))
# max_iter = maximum number of iterations, method stops if converged, 
#            or after this many iterations, whichever comes first (default: 30)
# traj_range = range of allowed offsets (delta), used in offset estimation (default: c(-50000,50000))
# fix_first = If TRUE (default) forces all clusters to have the same baseline parameter (theta_1)
# phi_range = if fix_first = TRUE, this range determines the allowed values of theta_1
#
# Output:
#
# result - list containing:
# result[[1]] <- cluster assignment function (c(i))
# result[[2]] <- estimated trajectory parameters (theta)
# result[[3]] <- estimated offsets (delta_i)
# result[[4]] <- model fit outputs
      
  if (K == 1){
    
    k1_fit <- optim(phi_init,loss_traj_arg,data=data,tps=tps,lambda=lambda,traj_fun=phi,traj_range=traj_range)
    
    result <- vector("list",4)
  
    result[[1]] <- rep(1,dim(data)[1])
  
    result[[2]] <- list(k1_fit$par)
  
    result[[4]] <- list(k1_fit)

    result[[3]] <- optimal_offset(data=data,tps=tps,lambda=lambda,traj_fun=phi,traj_arg=k1_fit$par,traj_range=traj_range)
    
  } else {
  
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
  
    result <- vector("list",4)
  
    result[[1]] <- cluster
  
    result[[2]] <- args
  
    result[[4]] <- fits
  
    final_delta <- deltas[[1]]
  
    for (k in 2:K){
    
      final_delta[cluster==k] <- deltas[[k]][cluster==k]

    }
  
    result[[3]] <- final_delta
  
  }
  
  return(result)
  
}

discrimK2 <- function(data,tps,result,traj_range){
# Discrimation score calculation for K = 2
# 
# Inputs :
#
# data = matrix representation of data (columns = timepoints, rows = individuals)
# tps = matrix representation of timepoints (columns = timepoints, rows = individuals)
# result = output from Temporal Clustering
# traj_range = range of allowed offsets (delta), used in offset estimation
#
# Outputs :
#  
# local_loss table
#     Column 1 - local_loss for assigned cluster
#     Column 2 - local_loss for unassigned cluster
#     Column 3 - Discrimination score = absolute difference of previous columns
  
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
      
    if (result[[1]][i] > 0){
      
      local_loss[i,3] <- abs(local_loss[i,result[[1]][i]] - local_loss[i,3-result[[1]][i]])
      
    }
  
  }
  
  return(local_loss)
  
}

discrimScore <- function(data,tps,result,traj_range){
# Discrimation score calculation
# 
# Inputs :
#
# data = matrix representation of data (columns = timepoints, rows = individuals)
# tps = matrix representation of timepoints (columns = timepoints, rows = individuals)
# result = output from Temporal Clustering
# traj_range = range of allowed offsets (delta), used in offset estimation
#
# Outputs :
#  
# Discrimination score (abs difference of local_loss for assigned cluster, and closest other cluster)
  
  num_ts <- dim(tps)[1]
  
  K <- length(unique(result[[1]]))
  
  local_loss <- data.frame(array(dim=c(num_ts,K+1)))
  colnames(local_loss) <- c(paste("k",1:K,sep=""),"discrim_score")
  
  par <- vector("list",K)
  delta <- vector("list",K)
  
  for (k in 1:K){
    
    par[[k]] <- result[[2]][[k]]
    delta[[k]] <- optimal_offset(data,tps,lambda=0,traj_fun=exp_drop2,traj_args=par[[k]],traj_range=traj_range)
    
    for (i in 1:num_ts){
    
      local_loss[i,k] <- localLoss(delta[[k]][i],data[i,],tps[i,],lambda=0,FUN= exp_drop2,fun_args=par[[k]])
      
    }
    
  }
  
  for (i in 1:num_ts){
    
    assigned_loss <- local_loss[i,result[[1]][i]]
    
    unassigned_loss <- local_loss[i,-c(result[[1]][i],K+1)]
    
    local_loss[i,K+1] <- min(unassigned_loss - assigned_loss)
      
  }  
  
  return(local_loss[,K+1])
  
}

multiTS.matrixPlot <- function(sim_matrix){
# LEGACY, REMOVE?

  plot(sim_matrix@mat_tps[1,],sim_matrix@mat_data[1,],type="l",ylim=c(0,30),xlim=c(0,8000))

  for (i in 2:dim(sim_matrix@mat_tps)[1]){
    lines(sim_matrix@mat_tps[i,],sim_matrix@mat_data[i,],col=i)
  }
  
}

simulateExp <- function(tps,num_time_series=421,shift_range=c(-3000,1000),
                        exp_drop_par=array(c(30,0.0005),dim=c(1,2)),num_clusters=2){
# Simulate data from the exponential decline trajectory, main simulation function
#
# Inputs:
#
# tps = list of timepoints for each individual (tau_i)
# num_time_series = number of individuals required (default 421)
# shift_range = range of allowed offsets (delta_i; default c(-3000,1000))
# exp_drop_par = parameters of exponential decline model (default c(30,0.0005))
# num_clusters = Number of clusters to simulate (default: 2)
#
# Output:
#
# sim = simualted data in a multiTS object, including different transformations
  
  
  sim_data <- simulateData(num_time_series = num_time_series,time_of_visits =tps,exp_drop_par = exp_drop_par,num_clusters = num_clusters,shift_range=shift_range,noise = 1.5)
  
  #plot.multiTS(sim_data,yl=c(0,30))
  
  sim <- multiTS.matrix(sim_data)
  
  ind <- which(sim@mat_data<0,arr.ind=T)
  
  sim@mat_data[ind] <- NA
  sim@mat_tps[ind] <- NA
  sim@mat_data_2[ind] <- NA
  
  sim@mat_data <- round(sim@mat_data)
  
  sim@mat_data_2 <- sim@mat_data
  
  ind2 <- which(sim@mat_data>30,arr.ind=T)
  
  sim@mat_data[ind2] <- 30
  
  return(sim)
  
}

simulateData <- function(num_time_series,time_of_visits,exp_drop_par=array(c(30,0.0005),dim=c(1,2)),
                          num_clusters = 2,noise=2,shift_range=c(-3000,1000)){
  # Simulate data from exponential decline model using time-points
  #
  # Inputs:
  #
  # num_time_series = number of individuals required
  # time_of_visits = list of timepoints for each individual (tau_i)
  # exp_drop_par = parameters of exponential decline model (default c(30,0.0005))
  # num_clusters = Number of clusters to simulate (default: 2)  
  # noise = standard deviation of noise (default: 2)
  # shift_range = range of allowed offsets (delta_i; default c(-3000,1000))
  #
  #
  # Outputs :
  #
  # sim_data = multiTS object
  
  new_tps <- time_of_visits
  
  # assign time series to clusters
  
  cluster_id <- mat.or.vec(num_time_series,1)
  
  cluster_id <- sample(num_clusters,num_time_series,replace=TRUE)

  
  # generate data
  #
  
  data <- vector("list", num_time_series)
  
  time_shifts <- numeric(num_time_series)
  
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
