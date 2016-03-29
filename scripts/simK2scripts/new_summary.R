setwd("../../results/K2sim/standard")

rb <- c("1_1","1_2","1_3","2_1","2_2","2_3","3_1","3_2","3_3")

thetas <- array(c(29,29,29,29,29,0.0005,0.0003,0.0005,0.0006,0.0008,29,29,29,29,29,0.0005,0.0001,0.0001,0.001,0.001),dim=c(5,4))


for (theta_num in 1:5){

sim_df <- 
data.frame(paq_aac=c(rep("PAQUID",50),rep("AAC",50)),range_bias=NA,label_swap=logical(100),theta1_1=NA,theta1_2=NA,theta2_1=NA,theta2_2=NA)

k <- 1

setwd("paq_sim")

setwd(paste("theta",theta_num,sep=""))

files <- list.files()

for (i in 1:length(files)){
  
  load(files[i])
  

    
    
    
    #if (sum(table(results[[j]]@clusters)[array(c(1,2,1,2),dim=c(2,2))]) < sum(table(results[[j]]@clusters)[array(c(1,2,2,1),dim=c(2,2))])){
    if (dist(rbind(thetas[theta_num,c(2,4)],simRes@params[c(2,4)])) < dist(rbind(thetas[theta_num,c(2,4)],simRes@params[c(4,2)]))) { 
     
      sim_df[k,3] = T
      
      sim_df[k,4:7] <- simRes@params
      
    } else {
    
      sim_df[k,c(6,7,4,5)] <- simRes@params
        
    }
    
    
    k <- k + 1
     
  
  
}


setwd("..")
setwd("..")

setwd("aac_sim")

setwd(paste("theta",theta_num,sep=""))

files <- list.files()

for (i in 1:length(files)){
  
  load(files[i])

    
    
    #if (sum(table(results[[j]]@clusters)[array(c(1,2,1,2),dim=c(2,2))]) < sum(table(results[[j]]@clusters)[array(c(1,2,2,1),dim=c(2,2))])){
    if (dist(rbind(thetas[theta_num,c(2,4)],simRes@params[c(2,4)])) < dist(rbind(thetas[theta_num,c(2,4)],simRes@params[c(4,2)]))) { 
     
      sim_df[k,3] = T
      
      sim_df[k,4:7] <- simRes@params
      
    } else {
    
      sim_df[k,c(6,7,4,5)] <- simRes@params
        
    }
    
    
    k <- k + 1
     
  
  
}


setwd("..")
setwd("..")

file_name <- paste("new_",theta_num,".RData",sep="")

save(sim_df,file=file_name)

#file_name <- paste("k2_",theta_num,"_1_1.pdf",sep="")

#p <- ggplot(sim_df, aes(factor(paq_aac), theta1_1))
#p + geom_boxplot(aes(fill = factor(range_bias))) + xlab("Simulations") + ylab("Theta 1") + scale_fill_discrete(name="Range & Bias",labels=c("Range 1, Bias 1", "Range 1, Bias 2","Range 1, Bias 3","Range 2, Bias 1","Range 2, Bias 2","Range 2, Bias 3","Range 3, Bias 1","Range 3, Bias 2","Range 3, Bias 3"))+ geom_hline(yintercept = 30,color="red",linetype="dashed",size=1)
#ggsave(file_name,width=5,height=2.8)

#file_name <- paste("k2_",theta_num,"_1_2.pdf",sep="")

#p <- ggplot(sim_df, aes(factor(paq_aac), theta1_2))
#p + geom_boxplot(aes(fill = factor(range_bias))) + xlab("Simulations") + ylab("Theta 2") + scale_fill_discrete(name="Range & Bias",labels=c("Range 1, Bias 1", "Range 1, Bias 2","Range 1, Bias 3","Range 2, Bias 1","Range 2, Bias 2","Range 2, Bias 3","Range 3, Bias 1","Range 3, Bias 2","Range 3, Bias 3"))+ geom_hline(yintercept = thetas[theta_num,4],color="red",linetype="dashed",size=1)
#ggsave(file_name,width=5,height=2.8)

#file_name <- paste("k2_",theta_num,"_2_1.pdf",sep="")

#p <- ggplot(sim_df, aes(factor(paq_aac), theta2_1))
#p + geom_boxplot(aes(fill = factor(range_bias))) + xlab("Simulations") + ylab("Theta 1") + scale_fill_discrete(name="Range & Bias",labels=c("Range 1, Bias 1", "Range 1, Bias 2","Range 1, Bias 3","Range 2, Bias 1","Range 2, Bias 2","Range 2, Bias 3","Range 3, Bias 1","Range 3, Bias 2","Range 3, Bias 3"))+ geom_hline(yintercept = 30,color="red",linetype="dashed",size=1)
#ggsave(file_name,width=5,height=2.8)

#file_name <- paste("k2_",theta_num,"_2_2.pdf",sep="")

#p <- ggplot(sim_df, aes(factor(paq_aac), theta2_2))
#p + geom_boxplot(aes(fill = factor(range_bias))) + xlab("Simulations") + ylab("Theta 2") + scale_fill_discrete(name="Range & Bias",labels=c("Range 1, Bias 1", "Range 1, Bias 2","Range 1, Bias 3","Range 2, Bias 1","Range 2, Bias 2","Range 2, Bias 3","Range 3, Bias 1","Range 3, Bias 2","Range 3, Bias 3"))+ geom_hline(yintercept = thetas[theta_num,2],color="red",linetype="dashed",size=1)
#ggsave(file_name,width=5,height=2.8)

}

