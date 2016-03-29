# script to convert multiTS to LCMM format

TCvsLCMM <- data.frame(method=
    rep(c("Temporal Clustering","LCMM study time","LCMM age","normalised LCMM study time","normalised LCMM age"),each=5),
  dataset=rep(c(rep("AAC",4),"PAQUID"),5),variable = rep(c("MMSE first visit","1 APOE e4","2 APOE e4","CSF tau","MMSE first visit"),5),
  N=NA,abs_beta=NA,se=NA,abs_z=NA,pvalue=NA)

install.packages("NormPsy")
library(NormPsy)

aac_lcmm <- data.frame(array(dim=c(length(which(!is.na(mmse_aac@mat_data))),4)))

colnames(aac_lcmm) <- c("ID","MMSE","TIME","AGE")

k <- 1

for (i in 1:mmse_aac@numTS){

  for (j in 1:length(which(!is.na(mmse_aac@mat_data[i,])))){
    
    aac_lcmm[k,] <- c(i,mmse_aac@mat_data[i,j],mmse_aac@mat_tps[i,j],aac_demog_999[i,"bl_age"]*365.25 + mmse_aac@mat_tps[i,j])
    
    k <- k + 1
    
  }
  
}

aac_lcmm[,"norm"] <- normMMSE(aac_lcmm[,2])

aac_lcmm[,"AGE50"] <- aac_lcmm[,"AGE"] - 50*365.25

lcmm_cl <- array(dim=c(2412,4))

mspl_raw_zero <- lcmm(MMSE ~ TIME,random=~ TIME,mixture=~ TIME, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[,1] <- mspl_raw_zero$pprob[,"class"]

mspl_norm_zero <- lcmm(norm ~ TIME,random=~ TIME,mixture=~ TIME, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[,2] <- mspl_norm_zero$pprob[,"class"]



aac_lcmm <- aac_lcmm[complete.cases(aac_lcmm),]

# datnew[,"AGE50"]<-datnew[,"TIME"] # something not needed any more?

comp.case <- unique(aac_lcmm[,1])

mspl_raw_age <- lcmm(MMSE ~ AGE50,random=~ AGE50,mixture=~ AGE50, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[comp.case,3] <- mspl_raw_age$pprob[,"class"]

mspl_norm_age <- lcmm(norm ~ AGE50,random=~ AGE50,mixture=~ AGE50, subject="ID",ng=2,data=aac_lcmm,link = "splines")
lcmm_cl[comp.case,4] <- mspl_norm_age$pprob[,"class"]

colnames(lcmm_cl) <- c("raw_zero","raw_age","norm_zero","norm_age")

# now compare to baseline MMSE and APOE, using betas

aac_999_df <- data.frame(cbind(aac_res[[1]],lcmm_cl,aac_999_df))
names(aac_999_df)[1] <- "k"

for (i in 1:5){aac_999_df[,i] <- factor(aac_999_df[,i])}

aac_999_df[,"gender"] <- factor(aac_999_df[,"gender"])

aac_999_df[,"cohort"] <- factor(aac_999_df[,"cohort"])

TCvsLCMM <- data.frame(method=
    rep(c("Temporal Clustering","LCMM study time","LCMM age","normalised LCMM study time","normalised LCMM age"),each=5),
  dataset=rep(c(rep("AAC",4),"PAQUID"),5),variable = rep(c("MMSE first visit","1 APOE e4","2 APOE e4","CSF tau","MMSE first visit"),5),
  N=NA,beta=NA,se=NA,z=NA,pvalue=NA)

# AAC MMSE at first visit
k_mmse_fit <- summary(glm(k ~ bl_age + bl_mmse + gender + cohort, data = aac_999_df, family = "binomial"))
raw_zero_mmse_fit <- summary(glm(raw_zero ~ bl_age + bl_mmse + gender + cohort, data = aac_999_df, family = "binomial"))
raw_age_mmse_fit <- summary(glm(raw_age ~ bl_age + bl_mmse + gender + cohort, data = aac_999_df, family = "binomial"))
norm_zero_mmse_fit <- summary(glm(norm_zero ~ bl_age + bl_mmse + gender + cohort, data = aac_999_df, family = "binomial"))
norm_age_mmse_fit <- summary(glm(norm_age ~ bl_age + bl_mmse + gender + cohort, data = aac_999_df, family = "binomial"))

TCvsLCMM[1,4:8] <- c(k_mmse_fit$df.null,k_mmse_fit$coef[3,])
TCvsLCMM[6,4:8] <- c(raw_zero_mmse_fit$df.null,raw_zero_mmse_fit$coef[3,])
TCvsLCMM[11,4:8] <- c(raw_age_mmse_fit$df.null,raw_age_mmse_fit$coef[3,])
TCvsLCMM[16,4:8] <- c(norm_zero_mmse_fit$df.null,norm_zero_mmse_fit$coef[3,])
TCvsLCMM[21,4:8] <- c(norm_age_mmse_fit$df.null,norm_age_mmse_fit$coef[3,])

# AAC APOE
k_apoe_fit <- summary(glm(k ~ bl_age + bl_mmse + apoe + educ + gender, data = aac_999_df, family = "binomial"))
raw_zero_apoe_fit <- summary(glm(raw_zero ~ bl_age + bl_mmse + apoe+ educ + gender, data = aac_999_df, family = "binomial"))
raw_age_apoe_fit <- summary(glm(raw_age ~ bl_age + bl_mmse + apoe+ educ + gender, data = aac_999_df, family = "binomial"))
norm_zero_apoe_fit <- summary(glm(norm_zero ~ bl_age + bl_mmse + apoe+ educ + gender, data = aac_999_df, family = "binomial"))
norm_age_apoe_fit <- summary(glm(norm_age ~ bl_age + bl_mmse + apoe + educ + gender, data = aac_999_df, family = "binomial"))

TCvsLCMM[2:3,4:8] <- cbind(k_apoe_fit$df.null,k_apoe_fit$coef[4:5,])
TCvsLCMM[7:8,4:8] <- cbind(raw_zero_apoe_fit$df.null,raw_zero_apoe_fit$coef[4:5,])
TCvsLCMM[12:13,4:8] <- cbind(raw_age_apoe_fit$df.null,raw_age_apoe_fit$coef[4:5,])
TCvsLCMM[17:18,4:8] <- cbind(norm_zero_apoe_fit$df.null,norm_zero_apoe_fit$coef[4:5,])
TCvsLCMM[22:23,4:8] <- cbind(norm_age_apoe_fit$df.null,norm_age_apoe_fit$coef[4:5,])

# AAC tau
k_tau_fit <- summary(glm(k ~ bl_age + bl_mmse + apoe + educ + gender + tau, data = aac_999_df, family = "binomial"))
raw_zero_tau_fit <- summary(glm(raw_zero ~ bl_age + bl_mmse + apoe+ educ + gender + tau, data = aac_999_df, family = "binomial"))
raw_age_tau_fit <- summary(glm(raw_age ~ bl_age + bl_mmse + apoe+ educ + gender + tau, data = aac_999_df, family = "binomial"))
norm_zero_tau_fit <- summary(glm(norm_zero ~ bl_age + bl_mmse + apoe+ educ + gender + tau, data = aac_999_df, family = "binomial"))
norm_age_tau_fit <- summary(glm(norm_age ~ bl_age + bl_mmse + apoe + educ + gender + tau, data = aac_999_df, family = "binomial"))

TCvsLCMM[4,4:8] <- c(k_tau_fit$df.null,k_tau_fit$coef[8,])
TCvsLCMM[9,4:8] <- c(raw_zero_tau_fit$df.null,raw_zero_tau_fit$coef[8,])
TCvsLCMM[14,4:8] <- c(raw_age_tau_fit$df.null,raw_age_tau_fit$coef[8,])
TCvsLCMM[19,4:8] <- c(norm_zero_tau_fit$df.null,norm_zero_tau_fit$coef[8,])
TCvsLCMM[24,4:8] <- c(norm_age_tau_fit$df.null,norm_age_tau_fit$coef[8,])





new_paq <- paquid[!is.na(paquid[,"MMSE"]),]

ind <- which(as.numeric(table(new_paq[,1]))>1)

new_paq <- new_paq[new_paq[,1] %in% ind,]

new_paq[,"norm"] <- normMMSE(new_paq[,"MMSE"])



final_subjects <- unique(new_paq[,1])
first_tp <- numeric(421)

for (i in 1:421){
  
  tmp <- which(new_paq[,1] == final_subjects[i])
  
  first_tp[i] <- min(new_paq[tmp,"age"])
  
  new_paq[tmp,"TIME"] <- (new_paq[tmp,"age"] - first_tp[i]) * 365.25
  
  new_paq[tmp,"AGE65"] <- (new_paq[tmp,"age"] - 65) * 365.25
  
}



paq_cl <- array(dim=c(421,4))

paq_raw_zero <- lcmm(MMSE ~ TIME,random=~ TIME,mixture=~ TIME, subject="ID",ng=2,data=new_paq,link = "splines")
paq_cl[,1] <- paq_raw_zero$pprob[,"class"]

paq_norm_zero <- lcmm(norm ~ TIME,random=~ TIME,mixture=~ TIME, subject="ID",ng=2,data=new_paq,link = "splines")
paq_cl[,2] <- paq_norm_zero$pprob[,"class"]

paq_raw_age <- lcmm(MMSE ~ AGE65,random=~ AGE65,mixture=~ AGE65, subject="ID",ng=2,data=new_paq,link = "splines")
paq_cl[,3] <- paq_raw_age$pprob[,"class"]

paq_norm_age <- lcmm(norm ~ AGE65,random=~ AGE65,mixture=~ AGE65, subject="ID",ng=2,data=new_paq,link = "splines")
paq_cl[,4] <- paq_norm_age$pprob[,"class"]


paq_df <- data.frame(cbind(paq_res[[1]],paq_cl,paq_demog))

names(paq_df)[1:5] <- c("k","raw_zero","raw_age","norm_zero","norm_age")


for (i in 1:5){paq_df[,i] <- factor(paq_df[,i])}

paq_df[,"gender"] <- factor(paq_df[,"gender"],c(1,0))

paq_tc_lcmm <- array(dim=c(4,5))

paq_k_mmse_fit <- summary(glm(k ~ bl_age + bl_mmse + gender, data = paq_df, family = "binomial"))
paq_raw_zero_mmse_fit <- summary(glm(raw_zero ~ bl_age + bl_mmse + gender, data = paq_df, family = "binomial"))
paq_raw_age_mmse_fit <- summary(glm(raw_age ~ bl_age + bl_mmse + gender, data = paq_df, family = "binomial"))
paq_norm_zero_mmse_fit <- summary(glm(norm_zero ~ bl_age + bl_mmse + gender, data = paq_df, family = "binomial"))
paq_norm_age_mmse_fit <- summary(glm(norm_age ~ bl_age + bl_mmse + gender, data = paq_df, family = "binomial"))

TCvsLCMM[5,4:8] <- c(paq_k_mmse_fit$df.null,paq_k_mmse_fit$coef[3,])
TCvsLCMM[10,4:8] <- c(paq_raw_zero_mmse_fit$df.null,paq_raw_zero_mmse_fit$coef[3,])
TCvsLCMM[15,4:8] <- c(paq_raw_age_mmse_fit$df.null,paq_raw_age_mmse_fit$coef[3,])
TCvsLCMM[20,4:8] <- c(paq_norm_zero_mmse_fit$df.null,paq_norm_zero_mmse_fit$coef[3,])
TCvsLCMM[25,4:8] <- c(paq_norm_age_mmse_fit$df.null,paq_norm_age_mmse_fit$coef[3,])

TCvsLCMM[,"N"] <- TCvsLCMM[,"N"] + 1

TCvsLCMM[,5:8] <- signif(TCvsLCMM[,5:8],2)

setwd("../results/")

write.csv(TCvsLCMM[1:10,],file="TCvsLCMM_main.csv",row.names = F)
write.csv(TCvsLCMM[c(5*(1:5)-4,5*(1:5)-3,5*(1:5)-2,5*(1:5)-1,5*(1:5)),],file="TCvsLCMM_supp.csv",row.names = F)





#plotting
#pred <- predictY(lin_raw_age,newdata=datnew,var.time="AGE50",draws=T)
#plot(pred$times[["AGE50"]],pred$pred[,"Ypred_50_class2"],ylim=c(0,100))
#points(pred$times[["TIME"]],pred$pred[,"Ypred_50_class1"],col="red")