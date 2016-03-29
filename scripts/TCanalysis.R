set.seed(1)

# Load scripts
source("TCscripts.R")

# Code to extract and gather data
source("TCload.R")


# Code to generate Table 1 - cohort demographics 
 source("TCdemog.R")


# NOTE: Figure 1 illustrative only and constructed in powerpoint

# Figures 2a - b will be introduced later, when temporal clustering has been performed


setwd("../scripts/")

# Generate summary of DPS simulation and
# Supplementary Figures 1 & 2 (saved into DPSsim folder)
source("simDPSsummary.R")

setwd("../scripts/")

# Generate summary of Temporal Clustering simulation and
# Figure 2 and Supplementary Figures 3 & 4
source("simK2summary.R") # needs scripts and data from cluster to complete


# fit sigmoids and exponential decline models to AAC and PAQUID
#

aac_logistic <- optim(c(30,-0.001,0),loss_traj_arg,data=mmse_aac@mat_data,tps=mmse_aac@mat_tps,lambda=0,traj_fun=logistic3,traj_range=c(-50000,50000))
aac_logistic$value # Loss of logistic AAC model
aac_logistic$par # Parameters of logistic AAC model

aac_exp <- optim(c(28,-0.0004),loss_traj_arg,data=mmse_aac@mat_data,tps=mmse_aac@mat_tps,lambda=0,traj_fun=exp_drop2,traj_range=c(-50000,50000))
aac_exp$value # Loss of exponential decline AAC model
aac_exp$par # Parameters of exponential decline AAC model


paq_logistic <- optim(c(30,-0.001,0),loss_traj_arg,data=mmse_paq@mat_data,tps=mmse_paq@mat_tps,lambda=0,traj_fun=logistic3,traj_range=c(-50000,50000))
paq_logistic$value # Loss of logistic PAQUID model
paq_logistic$par # Parameters of logistic PAQUID model

paq_exp <- optim(c(28,-0.0004),loss_traj_arg,data=mmse_paq@mat_data,tps=mmse_paq@mat_tps,lambda=0,traj_fun=exp_drop2,traj_range=c(-50000,50000))
paq_exp$value # Loss of exponential decline PAQUID model
paq_exp$par # Parameters of exponential decline PAQUID model

save(aac_logistic,aac_exp,paq_logistic,paq_exp,file="align.RData")

#load("align.RData")

# Create AAC disease progression score plots for Figure S5
#

aac_logistic_delta <- optimal_offset(data=mmse_aac@mat_data,tps=mmse_aac@mat_tps,lambda=0,traj_fun=logistic3,traj_arg=aac_logistic$par,traj_range=c(-50000,50000))

pdf("aac_sigmoid.pdf",width=5,height=5)
par(mar=c(4,4,0.5,0.5))
mmse_aac@clusters <- 1:2412
plot.multiTS(mmse_aac,yl=c(0,30),xl=c(0,10000),shift = aac_logistic_delta+20089,xlabel = "Disease progression score in years + 55",ylabel="MMSE")
points(0:10300,logistic3(-20000:-9700,aac_logistic$par))
dev.off()

aac_exp_delta <- optimal_offset(data=mmse_aac@mat_data,tps=mmse_aac@mat_tps,lambda=0,traj_fun=exp_drop2,traj_arg=aac_exp$par,traj_range=c(-50000,50000))

pdf("aac_exp.pdf",width=5,height=5)
par(mar=c(4,4,0.5,0.5))
plot.multiTS(mmse_aac,xl=c(0,8200),yl=c(0,30),shift = aac_exp_delta,xlabel = "Disease progression score in days",ylabel="MMSE")
points(0:8400,exp_drop2(0:8400,aac_exp$par))
dev.off()


# Perform Temporal Clustering on AAC and PAQUID
#

paq_res <- temporalClustering(mmse_paq@mat_data,mmse_paq@mat_tps,K=2,fix_first=F)

print("PAQUID temporal clustering K2: theta estimates")
print(paq_res[[2]])

print("PAQUID temporal clustering K2: cluster sizes")
print(table(paq_res[[1]]))

aac_res <- temporalClustering(mmse_aac@mat_data,mmse_aac@mat_tps,K=2,fix_first=T)

print("AAC temporal clustering K2: theta estimates")
print(aac_res[[2]])

print("AAC temporal clustering K2: cluster sizes")
print(table(aac_res[[1]]))

save(paq_res,aac_res,file="orig_k2.RData")

#load("orig_k2.RData")



# Get the bootstrap intervals (DPS)
setwd("../scripts/DPSbootstrap/")
source("DPSbootstrap.R")





# Get the bootstrap intervals (Temporal Clustering K = 2)
setwd("../K2bootstrap/")
source("K2bootstrapSummary.R")

setwd("..")

# Generate Figure 2
#

pdf("aac.pdf",width=5,height=5)
par(mar=c(4,4,0.5,0.5))
mmse_aac@clusters <- aac_res[[1]]
plot.multiTS(mmse_aac,yl=c(0,30),xlabel = "Years since first visit",ylabel="MMSE",days2years=T)
legend("bottomright",c("AAC k=1","AAC k=2"),lty=1,lwd=5,col=1:2)
dev.off()

pdf("paquid.pdf",width=5,height=5)
par(mar=c(4,4,0.5,0.5))
mmse_paq@clusters <- paq_res[[1]]
plot.multiTS(mmse_paq,yl=c(-4,30),xlabel = "Years since first visit",ylabel="MMSE",days2years=T)
legend("bottomleft",c("PAQUID k=1","PAQUID k=2"),lty=1,lwd=5,col=1:2)
dev.off()

pdf("aac_k2.pdf",width=5,height=5)
par(mar=c(4,4,0.5,1.5))
plot.multiTS(mmse_aac,yl=c(0,30),xl=c(-1000,19200),shift = aac_res[[3]],xlabel = "Disease progression score in years",ylabel="MMSE",recolor = T,days2years=T)
points((-1000:19500)/365.25,exp_drop2(-1000:19500,aac_res[[2]][[1]]),col="black")
points((-1000:4250)/365.25,exp_drop2(-1000:4250,aac_res[[2]][[2]]),col="red")

legend("topright",c("AAC k=1","AAC k=2"),lty=1,lwd=10,col=1:2)

dev.off()

pdf("paquid_k2.pdf",width=5,height=5)
par(mar=c(4,4,0.5,0.5))
plot.multiTS(mmse_paq,yl=c(0,30),xl=c(-1000,13200),shift = paq_res[[3]],xlabel = "Disease progression score in years",ylabel="MMSE",recolor = T,days2years=T)
points((-1000:13000)/365.25,exp_drop2(-1000:13000,paq_res[[2]][[2]]),col="black")
points((-1000:4920)/365.25,exp_drop2(-1000:4920,paq_res[[2]][[1]]),col="red")

legend("topright",c("PAQUID k=1","PAQUID k=2"),lty=1,lwd=10,col=1:2)
dev.off()


# Regression models
#
#

# How baseline MMSE is affected by demographic factors across the cohorts

aac_demog <- data.frame(bl_diag=bl_diags,bl_age=comb_age,gender=comb_gender,bl_mmse=mmse_aac@mat_data[,1],cohort=mmse_aac@cohort)

paq_demog <- data.frame(bl_diag=paq_bl_diag,bl_age=paq_age,gender=paq_gender,bl_mmse=mmse_paq@mat_data[,1],cohort=mmse_paq@cohort)

aac_demog_999 <- aac_demog

aac_demog_999[aac_demog[,2]==999,2] <- NA

#save(aac_demog,aac_demog_999,paq_demog,file="demog_df.RData")

aac_999_df <- data.frame(cbind(aac_res[[1]],aac_demog_999))
names(aac_999_df)[1] <- "k"

aac_999_df[,"k"] <- factor(aac_999_df[,"k"])

aac_999_df[,"gender"] <- factor(aac_999_df[,"gender"])

aac_999_df[,"cohort"] <- factor(aac_999_df[,"cohort"])


paq_df <- data.frame(cbind(paq_res[[1]],paq_demog))
names(paq_df)[1] <- "k"

paq_df[,"k"] <- factor(paq_df[,"k"])

paq_df[,"gender"] <- factor(paq_df[,"gender"],c(1,0))

new <- rbind(paq_demog,aac_demog_999)

new[new[,"cohort"]=="paquid","gender"] <- 2- new[new[,"cohort"]=="paquid","gender"]

new[,"cohort"] <- factor(new[,"cohort"])
new[,"gender"] <- factor(new[,"gender"])

# Linear regression of MMSE against age, gender, cohort. For Table S1
mmse_reg <- glm(bl_mmse ~ bl_age + gender + cohort,data=new)

write.csv(signif(summary(mmse_reg)$coef,2),file="mmse_reg.csv")


# For TV vs LCMM collect relevant data
#

adni_apoe_zero <- mmseADNIts@subj_id[adni_num_apoe == 0]
adni_apoe_one <- mmseADNIts@subj_id[adni_num_apoe == 1]
adni_apoe_two <- mmseADNIts@subj_id[adni_num_apoe == 2]

adni_ind <- which(mmse_aac@cohort == "adni")

zero_ind <- numeric(length(adni_apoe_zero))

one_ind <- numeric(length(adni_apoe_one))

two_ind <- numeric(length(adni_apoe_two))

for (i in 1:length(adni_apoe_zero)){
  
  zero_ind[i] <- adni_ind[which(mmse_aac@subj_id[adni_ind] == adni_apoe_zero[i])]
  
}

for (i in 1:length(adni_apoe_one)){
  
  one_ind[i] <- adni_ind[which(mmse_aac@subj_id[adni_ind] == adni_apoe_one[i])]
  
}

for (i in 1:length(adni_apoe_two)){
  
  two_ind[i] <- adni_ind[which(mmse_aac@subj_id[adni_ind] == adni_apoe_two[i])]
  
}

table(aac_res[[1]][c(zero_ind,one_ind,two_ind)],c(rep(0,length(adni_apoe_zero)),rep(1,length(adni_apoe_one)),rep(2,length(adni_apoe_two))))


all_apoe <- array(dim=c(2412))
all_educ <- array(dim=c(2412))

all_apoe[mmse_aac@cohort == "adni"] <- adni_num_apoe # only put in APOE data for ADNI
all_educ[mmse_aac@cohort == "adni"] <- adni_years_educ

setwd("../data/adni_sep_2015/")

csf <- read.csv("UPENNBIOMK2.csv")

tau <- array(dim=c(mmseADNIts@numTS))

for (i in 1:mmseADNIts@numTS){

  inds <- which(csf[,"RID"] == mmseADNIts@subj_id[i])
  
  if (length(inds)>0){
  
    bl <- which(csf[inds,"VISCODE"] == "bl")
    
    if (length(bl)>0){
     
      tau[i] <- csf[inds[bl[1]],"TTAU"]
       
    }
      
  }
  
}

all_tau <- array(dim=c(2412))

all_tau[mmse_aac@cohort == "adni"] <- tau 

setwd("../../results/")

# APOE logistic regression
aac_demog <- data.frame(bl_diag=bl_diags,lv_diag=lv_diags,apoe=factor(all_apoe),educ=all_educ,tau=all_tau,
                        bl_age=comb_age,gender=comb_gender,bl_mmse=mmse_aac@mat_data[,1],cohort=mmse_aac@cohort)

aac_demog_999 <- aac_demog

aac_demog_999[aac_demog[,"bl_age"]==999,"bl_age"] <- NA

aac_999_df <- data.frame(cbind(aac_res[[1]],aac_demog_999))
names(aac_999_df)[1] <- "k"

aac_999_df[,"k"] <- factor(aac_999_df[,"k"])

aac_999_df[,"gender"] <- factor(aac_999_df[,"gender"])

aac_999_df[,"cohort"] <- factor(aac_999_df[,"cohort"])

# Tables comparing Temporal Clustering to LCMM

# Comparison of temporal clustering to LCMMs
setwd("../scripts/")
source("TCvsLCMM.R")


# Distribution of ADNI and AIBL diagnoses through clusters
adni_aibl <- mmse_aac@cohort=="adni" | mmse_aac@cohort=="aibl"
table(aac_res[[1]][adni_aibl],lv_diags[adni_aibl])

# Check that similiar is true with the filter
aac_discrim <- discrimScore(mmse_aac@mat_data,mmse_aac@mat_tps,aac_res,c(-50000,50000))
discrim10 <- aac_discrim>10

adni_aibl_10 <- (adni_aibl & discrim10)
table(aac_res[[1]][adni_aibl_10],lv_diags[adni_aibl_10])
