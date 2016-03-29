mmseTS <- vector("list",5)

camdTS <- mmse_both_cpath

mmseTS[[1]] <- mmse_aac
mmseTS[[2]] <- mmseADNIts
mmseTS[[3]] <- mmseAIBLts
mmseTS[[4]] <- mmse_both_cpath
mmseTS[[5]] <- mmse_paq

mmse_demog <- data.frame(combined=NA,adni=NA,aibl=NA,camd=NA,paquid=NA)

rownames(mmse_demog) <- "Num ind any data"

mmse_demog["Num ind any data",2:5] <- c(1229,861,5358,500)

mmse_demog["Num ind any data","combined"] <- sum(mmse_demog["Num ind any data",2:4],na.rm=T)

for (i in 1:5){
  
  mmse_demog["Num ind > 1 TP",i] <- mmseTS[[i]]@numTS

  tmp_lengths <- unlist(lapply(mmseTS[[i]]@tps,length))
  
  mmse_demog["Median num TP",i] <- median(tmp_lengths)

  mmse_demog["LQ num TP",i] <- as.numeric(quantile(unlist(lapply(mmseTS[[i]]@tps,length)),0.25))

  mmse_demog["UQ TP",i] <- as.numeric(quantile(unlist(lapply(mmseTS[[i]]@tps,length)),0.75))
  
  last_tp <- c()
  
  for (j in 1:length(tmp_lengths)){
   
    last_tp[j] <- mmseTS[[i]]@tps[[j]][tmp_lengths[j]]
    
  }
  
  mmse_demog["Median years",i] <- signif(median(last_tp) / 365.25,2)

  mmse_demog["LQ years",i] <- signif(as.numeric(quantile(last_tp,0.25))/ 365.25,2)

  mmse_demog["UQ years",i] <- signif(as.numeric(quantile(last_tp,0.75))/ 365.25,2)
  
  
}

setwd("../data")

# ADNI

setwd("adni_sep_2015/")

diag <- read.csv("DXSUM_PDXCONV_ADNIALL.csv",header=T,as.is=T)

adni_bl_diag <- array(dim=c(length(mmseADNIts@subj_id),1))
adni_lv_diag <- array(dim=c(length(mmseADNIts@subj_id),1))
adni_bl_date <- list()

for (i in 1:length(mmseADNIts@subj_id)){
  
  ind <- which(diag[,"RID"] == mmseADNIts@subj_id[i])
  
  ind2 <- which(diag[ind,"VISCODE"]=="bl")
  
  adni_bl_diag[i] <- diag[ind[ind2],"DXCURREN"]
  
  adni_bl_date[[i]] <- as.Date(diag[ind[ind2],"EXAMDATE"])
  
  ind3 <- order(as.Date(diag[ind[!is.na(diag[ind,"DXCURREN"])],"EXAMDATE"]),decreasing = T)[1]
  
  adni_lv_diag[i] <- diag[ind[!is.na(diag[ind,"DXCURREN"])][ind3],"DXCURREN"]
  
}

mmse_demog[c("Num CTL","Num MCI","Num AD"),"adni"] <- table(adni_bl_diag)

demog <- read.csv("PTDEMOG.csv",header=T,as.is=T)

adni_df <- data.frame(age=NA,education=NA,gender=NA)

adni_years_educ <- numeric(length(mmseADNIts@subj_id))

for (i in 1:length(mmseADNIts@subj_id)){
  
  ind <- which(demog[,"RID"] == mmseADNIts@subj_id[i])
  
  dob <- as.Date(paste(demog[ind[1],"PTDOBYY"],demog[ind[1],"PTDOBMM"],1,sep="-"))
  
  adni_years_educ[i] <- demog[ind[1],"PTEDUCAT"]
  
  adni_df[i,] <- c(signif(as.numeric(adni_bl_date[[1]] - dob)/365.25,2),demog[ind[1],"PTEDUCAT"],demog[ind[1],"PTGENDER"])
  
}

mmse_demog[c("Num Male","Num Female"),"adni"] <- table(adni_df[,"gender"])

mmse_demog["Median bl age","adni"] <- median(adni_df[,"age"])

mmse_demog["LQ bl age","adni"] <- as.numeric(quantile(adni_df[,"age"],0.25))

mmse_demog["UQ bl age","adni"] <- as.numeric(quantile(adni_df[,"age"],0.75))

apoe <- read.csv("APOERES.csv",header=T,as.is=T)

adni_num_apoe <- array(dim=c(length(mmseADNIts@subj_id),1))

for (i in 1:length(mmseADNIts@subj_id)){
  
  ind <- which(apoe[,"RID"] == mmseADNIts@subj_id[i])
  
  adni_num_apoe[i] <- (apoe[ind[1],"APGEN1"] == 4) + (apoe[ind[1],"APGEN2"] == 4)
  
}

mmse_demog[c("Num 0 E4","Num 1 E4","Num 2 E4"),"adni"] <- table(adni_num_apoe)

setwd("..")



# AIBL

setwd("aibl_28apr_2015/")

#aibl_bl_date[[i]] <- as.Date(diag[ind[ind2],"EXAMDATE"])

visits <- read.csv("aibl_mmse_28-Apr-2015.csv",header=T,as.is=T)

aibl_bl_date <- list()

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(visits[,"RID"] == mmseAIBLts@subj_id[i])
  
  ind2 <- which(visits[ind,"VISCODE"]=="bl")
  
  aibl_bl_date[[i]] <- as.Date(visits[ind[ind2],"EXAMDATE"],"%m/%d/%Y")
  
}


diag <- read.csv("aibl_pdxconv_28-Apr-2015.csv",header=T,as.is=T)

aibl_bl_diag <- array(dim=c(length(mmseAIBLts@subj_id),1))
aibl_lv_diag <- array(dim=c(length(mmseAIBLts@subj_id),1))

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(diag[,"RID"] == mmseAIBLts@subj_id[i])
  
  ind2 <- which(diag[ind,"VISCODE"]=="bl")
  
  aibl_bl_diag[i] <- diag[ind[ind2],"DXCURREN"]
  
  inds <- which(diag[ind,"DXCURREN"]!= -4)
  
  ind3 <- order(diag[ind[inds],"VISCODE"],decreasing = T)[1]
  
  aibl_lv_diag[i] <- diag[ind[inds[ind3]],"DXCURREN"]
    
}

mmse_demog[c("Num CTL","Num MCI","Num AD"),"aibl"] <- table(aibl_bl_diag)

demog <- read.csv("aibl_ptdemog_28-Apr-2015.csv",header=T,as.is=T)

aibl_df <- data.frame(age=NA,gender=NA)

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(demog[,"RID"] == mmseAIBLts@subj_id[i])
  
  dob <- as.Date(paste(substr(demog[ind[1],"PTDOB"],2,5),"07",1,sep="-"))
  
  aibl_df[i,] <- c(signif(as.numeric(aibl_bl_date[[1]] - dob)/365.25,2),demog[ind[1],"PTGENDER"])
  
}

mmse_demog[c("Num Male","Num Female"),"aibl"] <- table(aibl_df[,"gender"])

mmse_demog["Median bl age","aibl"] <- median(aibl_df[,"age"])

mmse_demog["LQ bl age","aibl"] <- as.numeric(quantile(aibl_df[,"age"],0.25))

mmse_demog["UQ bl age","aibl"] <- as.numeric(quantile(aibl_df[,"age"],0.75))



apoe <- read.csv("aibl_apoeres_28-Apr-2015.csv",header=T,as.is=T)

aibl_num_apoe <- array(dim=c(length(mmseAIBLts@subj_id),1))

for (i in 1:length(mmseAIBLts@subj_id)){
  
  ind <- which(apoe[,"RID"] == mmseAIBLts@subj_id[i])
  
  aibl_num_apoe[i] <- (apoe[ind[1],"APGEN1"] == 4) + (apoe[ind[1],"APGEN2"] == 4)
  
}

mmse_demog[c("Num 0 E4","Num 1 E4","Num 2 E4"),"aibl"] <- table(aibl_num_apoe)

setwd("..")





# CAMD-1013  (Tombaugh and McIntyre (1992) for MMSE cutoffs, except no better than MCI due to recruitment of MCI and AD)

setwd("camd_june_2015/")

camd_bl_mmse <- array(dim=c(length(mmse_both_cpath@subj_id),1))
camd_lv_mmse <- array(dim=c(length(mmse_both_cpath@subj_id),1))

for (i in 1:length(mmse_both_cpath@subj_id)){
  
  camd_bl_mmse[i] <- mmse_both_cpath@data[[i]][1]
  camd_lv_mmse[i] <- mmse_both_cpath@data[[i]][length(mmse_both_cpath@data[[i]])]
  
}

camd_diag <- (camd_bl_mmse<18)+2
camd_lv_diag <- (camd_lv_mmse<18)+2

mmse_demog[c("Num MCI","Num AD"),"camd"] <- table(camd_diag)


demog <- read.csv("dm.csv",header=T,as.is=T)

camd_df <- data.frame(age=NA,gender=NA)

for (i in 1:length(camdTS@subj_id)){
  
  ind <- which(demog[,"USUBJID"] == camdTS@subj_id[i])
  
  camd_df[i,] <- c(demog[ind[1],"AGE"],as.numeric(demog[ind[1],"SEX"]=="F")+1)
  
}

mmse_demog[c("Num Male","Num Female"),"camd"] <- table(camd_df[,"gender"])

mmse_demog["Median bl age","camd"] <- median(camd_df[,"age"])

mmse_demog["LQ bl age","camd"] <- as.numeric(quantile(camd_df[,"age"],0.25))

mmse_demog["UQ bl age","camd"] <- as.numeric(quantile(camd_df[,"age"],0.75))


mmse_demog[c("Num CTL","Num MCI","Num AD"),"combined"] <- c(sum(mmse_demog["Num CTL",2:4],na.rm=T),sum(mmse_demog["Num MCI",2:4],na.rm=T),sum(mmse_demog["Num AD",2:4],na.rm=T))

mmse_demog[c("Num Male","Num Female"),"combined"] <- c(sum(mmse_demog["Num Male",2:4],na.rm=T),sum(mmse_demog["Num Female",2:4],na.rm=T))

comb_age <- c(camd_df[,"age"],adni_df[,"age"],aibl_df[,"age"])

comb_gender <- c(camd_df[,"gender"],adni_df[,"gender"],aibl_df[,"gender"])

mmse_demog["Median bl age","combined"] <- median(comb_age)

mmse_demog["LQ bl age","combined"] <- as.numeric(quantile(comb_age,0.25))

mmse_demog["UQ bl age","combined"] <- as.numeric(quantile(comb_age,0.75))



# PAQUID

mmse_paq <- multiTS.matrix(mmse_paq)

paq_bl_nc <- which(mmse_paq@mat_data[,1]>23)
paq_bl_mci <- which((mmse_paq@mat_data[,1]<24)&(mmse_paq@mat_data[,1]>17))
paq_bl_sci <- which(mmse_paq@mat_data[,1]<18)

paq_bl_diag <- mat.or.vec(421,1)
paq_bl_diag[paq_bl_nc] <- 1
paq_bl_diag[paq_bl_mci] <- 2
paq_bl_diag[paq_bl_sci] <- 3

mmse_demog[c("Num CTL","Num MCI"),"paquid"] <- table(paq_bl_diag)

mmse_demog[c("Num Male","Num Female"),"paquid"] <- table(paq_gender)

mmse_demog["Median bl age","paquid"] <- round(median(paq_age))

mmse_demog["LQ bl age","paquid"] <- round(as.numeric(quantile(paq_age,0.25)))

mmse_demog["UQ bl age","paquid"] <- round(as.numeric(quantile(paq_age,0.75)))

setwd("../..")

setwd("results/")

write.csv(mmse_demog[-1,],file="mmse_demog.csv")

bl_diags <- c(camd_diag,adni_bl_diag,aibl_bl_diag)
lv_diags <- c(camd_lv_diag,adni_lv_diag,aibl_lv_diag)

paquidTS <- mmse_paq

#table(bl_diags,lv_diags)

num_tp <- c()
last_tp <- c()
last_mmse <- c()

for (i in 1:421){

  num_tp[i] <- length(paquidTS@tps[[i]])
  
  last_mmse[i] <- paquidTS@data[[i]][num_tp[i]]
  
  last_tp[i] <- paquidTS@tps[[i]][num_tp[i]]
  
}

paq_lv_nc <- which(last_mmse>23)
paq_lv_mci <- which((last_mmse<24)&(last_mmse>17))
paq_lv_sci <- which(last_mmse<18)

paq_lv_diag <- mat.or.vec(421,1)
paq_lv_diag[paq_lv_nc] <- 1
paq_lv_diag[paq_lv_mci] <- 2
paq_lv_diag[paq_lv_sci] <- 3

#table(paq_bl_diag,paq_lv_diag)

print("Table 1 (demographics) stored in results/")

