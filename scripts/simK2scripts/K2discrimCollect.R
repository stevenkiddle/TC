setwd(tmp_path)

pos_ctls <- data.frame(iter=NA,real_clusters=NA,est_clusters=NA,discrim_score=NA)

neg_ctls <- data.frame(iter=NA,real_clusters=NA,est_clusters=NA,discrim_score=NA)
  
setwd("theta1")

for (iter in 1:50){
  
  file_name <- paste(sim_names,iter,".RData",sep="")
  
  load(file_name)

  tmp_table <- data.frame(iter=iter,real_clusters=simRes@clusters[,"real_clusters"],est_clusters=simRes@clusters[,"est_clusters"],discrim_score=simRes@loss[,3])
  
  neg_ctls <- rbind(neg_ctls,tmp_table)
  
}

for (sim_num in 2:5){
  
  print(sim_num)
  
  setwd("..")
  
  setwd(paste("theta",sim_num,sep=""))
  
  for (iter in 1:50){
  
    file_name <- paste(sim_names,iter,".RData",sep="")
  
    load(file_name)

    tmp_table <- data.frame(iter=iter + (sim_num-2)*50,real_clusters=simRes@clusters[,"real_clusters"],est_clusters=simRes@clusters[,"est_clusters"],discrim_score=simRes@loss[,3])
  
    pos_ctls <- rbind(pos_ctls,tmp_table)
  
  }
  
}

setwd("..")
setwd("..")
setwd("..")

save(neg_ctls,pos_ctls,file=tmp_file_name)