# Script to collect and summarise Disease Progression Score 
# simulation performed by TCsimDPS.R on a computing cluster

#install.packages('ggplot2')
library(ggplot2)

setwd("../results/")
setwd("DPSsim")

for (theta_num in 1:5){

# get ready to load PAQUID-style simulations
setwd("paq_sim")

files <- list.files()

# array to store results
# dim 1 = simulation number
# dim 2 = different theta values
# dim 3 = different data transformations, to simulate biases
# dim 4 = different ranges used to estimate offsets/deltas
# dim 5 = theta (theta_1,theta_2)
align_sim <- array(dim=c(length(files),5,3,3,2))

# Load in data one simulation at a time
for (i in 1:length(files)){
  
  load(files[i])
  
  align_sim[i,,,,] <- sim_thetas
  
}

# Data frame to store data in format for plotting
sim_df <- data.frame(paq_aac=c(rep("PAQUID",450),rep("AAC",450)),bias=NA,range=NA,theta1=NA,theta2=NA)

k <- 1

for (iter in 1:50){
  
  for (bias in 1:3){
    
    for (range in 1:3){
      
        # Reformat data into format for plotting
        sim_df[k,2:5] <- c(bias,range,align_sim[iter,theta_num,bias,range,])
        
        k <- k + 1 
      
      }
    
  }
  
}

# as above but for AAC-style simulations
setwd("../aac_sim")

files <- list.files()

align_sim <- array(dim=c(length(files),5,3,3,2))

for (i in 1:length(files)){
  
  load(files[i])
  
  align_sim[i,,,,] <- sim_thetas
  
}



for (iter in 1:50){
  
  for (bias in 1:3){
    
    for (range in 1:3){
      
        sim_df[k,2:5] <- c(bias,range,align_sim[iter,theta_num,bias,range,])
        
        k <- k + 1 
      
      }
    
  }
  
}

sim_df[,"range_bias"] <- paste(sim_df[,"range"],sim_df[,"bias"])

setwd("../")



real_theta2 <- (seq(1,10,2)  * 10^-4)[theta_num] # theta_2 values used in the 5 simulations

# generate plots for Supplementary Table 1
file_name <- paste("align_",theta_num,"_1.pdf",sep="")

p <- ggplot(sim_df, aes(factor(paq_aac), theta1))
p + geom_boxplot(aes(fill = factor(range_bias))) + xlab("Simulations") + ylab("Theta 1") + scale_fill_discrete(name="Range & Bias",labels=c("Range 1, Bias 1", "Range 1, Bias 2","Range 1, Bias 3","Range 2, Bias 1","Range 2, Bias 2","Range 2, Bias 3","Range 3, Bias 1","Range 3, Bias 2","Range 3, Bias 3"))+ geom_hline(yintercept = 29,color="red",linetype="dashed",size=1)
ggsave(file_name,width=5,height=2.8)

# generate plots for Supplementary Table 2
file_name <- paste("align_",theta_num,"_2.pdf",sep="")

p <- ggplot(sim_df, aes(factor(paq_aac), theta2))
p + geom_boxplot(aes(fill = factor(range_bias))) + xlab("Simulations") + ylab("Theta 2") + scale_fill_discrete(name="Range & Bias",labels=c("Range 1, Bias 1", "Range 1, Bias 2","Range 1, Bias 3","Range 2, Bias 1","Range 2, Bias 2","Range 2, Bias 3","Range 3, Bias 1","Range 3, Bias 2","Range 3, Bias 3"))+ geom_hline(yintercept = real_theta2,color="red",linetype="dashed",size=1)
ggsave(file_name,width=5,height=2.8)

}

print("Suppl Figure 1&2 plots created in results/DPSsim")
setwd("..")