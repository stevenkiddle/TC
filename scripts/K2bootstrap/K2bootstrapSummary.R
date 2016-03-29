setwd("../../results/K2bootstrap/")

files <- list.files()

num_files_each <- length(files)/2

# AAC

aac_boot <- c()

for (i in 1:num_files_each){

  load(files[i])
  
  for (j in 1:5){
    
    if (aac_cluster_boot[2*j-1,3] > aac_cluster_boot[2*j,3]){
    
      aac_cluster_boot[c((2*j - 1),2*j),2:3] <- aac_cluster_boot[c(2*j,(2*j - 1)),2:3]
      
    }
    
  }
  
  aac_boot <- rbind(aac_boot,aac_cluster_boot)
  
}

print("AAC TC K=2 bootstrap CI for k=1, theta 1")
print(quantile(aac_boot[aac_boot[,1]==1,2],0.025))
print(quantile(aac_boot[aac_boot[,1]==1,2],0.975))

print("AAC TC K=2 bootstrap CI for k=2, theta 1")
print(quantile(aac_boot[aac_boot[,1]==2,2],0.025))
print(quantile(aac_boot[aac_boot[,1]==2,2],0.975))

print("AAC TC K=2 bootstrap CI for k=1, theta 2")
print(quantile(aac_boot[aac_boot[,1]==1,3],0.025))
print(quantile(aac_boot[aac_boot[,1]==1,3],0.975))

print("AAC TC K=2 bootstrap CI for k=2, theta 2")
print(quantile(aac_boot[aac_boot[,1]==2,3],0.025))
print(quantile(aac_boot[aac_boot[,1]==2,3],0.975))


# PAQ

paq_boot <- c()

for (i in 1:num_files_each){

  load(files[i + num_files_each])
  
  for (j in 1:5){
    
    if (paq_cluster_boot[2*j-1,3] > paq_cluster_boot[2*j,3]){
    
      paq_cluster_boot[c((2*j - 1),2*j),2:3] <- paq_cluster_boot[c(2*j,(2*j - 1)),2:3]
      
    }
    
  }
  
  paq_boot <- rbind(paq_boot,paq_cluster_boot)
  
}




print("PAQUID TC K=2 bootstrap CI for k=1, theta 1")
print(quantile(paq_boot[paq_boot[,1]==1,2],0.025))
print(quantile(paq_boot[paq_boot[,1]==1,2],0.975))

print("PAQUID TC K=2 bootstrap CI for k=2, theta 1")
print(quantile(paq_boot[paq_boot[,1]==2,2],0.025))
print(quantile(paq_boot[paq_boot[,1]==2,2],0.975))

print("PAQUID TC K=2 bootstrap CI for k=1, theta 2")
print(quantile(paq_boot[paq_boot[,1]==1,3],0.025))
print(quantile(paq_boot[paq_boot[,1]==1,3],0.975))

print("PAQUID TC K=2 bootstrap CI for k=2, theta 2")
print(quantile(paq_boot[paq_boot[,1]==2,3],0.025))
print(quantile(paq_boot[paq_boot[,1]==2,3],0.975))