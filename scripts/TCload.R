print("starting to load data")
print("starting to load CAMD - this is the longest step, please be patient")

setwd("../data/")

# Load C-PATH CAMD data
#

setwd("camd_june_2015/")


cpath_cog <- read.csv("qs.csv",header=T,as.is=T) # load CAMD questionaire data
tests <- unique(cpath_cog[,"QSCAT"])

mmse <- cpath_cog[cpath_cog[,"QSCAT"]=="MMSE",] # take only MMSE

mmse1013 <- mmse[mmse[,"STUDYID"]==1013,] # extract MMSE data for C-1013

mmse1014 <- mmse[mmse[,"STUDYID"]==1014,] # extract MMSE data for C-1014

# These are the longest running clinical trials in CAMD

# Total scores not given, need to be calculated by summing scores of subtests


# calculate the total MMSE for C-1013 and gather in lists
#

people1013 <- unique(mmse1013[,"USUBJID"])

final_people <- c()

mmse1013_total <- array(dim=c(length(people1013),9))

mmse1013_days <- array(dim=c(length(people1013),9))

for (i in 1:length(people1013)){
 
  ind <- which(mmse1013[,"USUBJID"] == people1013[i])
  
  mmse_individual <- mmse1013[ind,]
  
  days <- sort(unique(mmse_individual[,"QSDY"]))
  mmse1013_days[i,1:length(days)] <- days
  scores <- array(dim=c(8,length(days)))
    
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Time Score","QSDY"])
  scores[1,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Time Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Place Score","QSDY"])
  scores[2,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Place Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Registration Score","QSDY"])
  scores[3,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Registration Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Attention and Calc Score","QSDY"])
  scores[4,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Attention and Calc Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Recall Score","QSDY"])
  scores[5,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Recall Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Naming Score","QSDY"])
  scores[6,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Naming Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Repetition Score","QSDY"])
  scores[7,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Repetition Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Drawing Score","QSDY"])
  scores[8,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Drawing Score","QSSTRESN"][tmp_order]
  
  mmse1013_total[i,1:length(days)] <- apply(scores,2,sum) # total scores for each time point and individual

}

# create lists to store timepoints (tps; tau_i) and MMSE scores (data; x_i(tau_i))
tps <- list()
data <- list()

# index to allow us to drop individuals if have only 1 timepoint 
k <- 1

# Gather tps and data
for (i in 1:length(people1013)){

  tmp <- mmse1013_days[i,][1:length(which(!is.na(mmse1013_days[i,])))]
  
  if (length(tmp) > 1){
   
    tps[[k]] <- tmp - tmp[1]
    data[[k]] <- mmse1013_total[i,][1:length(which(!is.na(mmse1013_days[i,])))]
    
    final_people[k] <- people1013[i]
    
    k <- k + 1
    
  }

}

# create multiTS object for C-1013 MMSE
mmse1013ts <- new("multiTS", numTS = length(tps) , tps = tps , data = data, clusters = rep(3,length(tps)), subj_id = final_people , cohort = "camd-1013")



# calculate the total MMSE for C-1014 and gather in lists
#

people1014 <- unique(mmse1014[,"USUBJID"])

mmse1014_total <- array(dim=c(length(people1014),9))

mmse1014_days <- array(dim=c(length(people1014),9))

final_people <- c()

for (i in 1:length(people1014)){
 
  ind <- which(mmse1014[,"USUBJID"] == people1014[i])
  
  mmse_individual <- mmse1014[ind,]
  
  days <- sort(unique(mmse_individual[,"QSDY"]))
  mmse1014_days[i,1:length(days)] <- days
  scores <- array(dim=c(8,length(days)))
    
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Time Score","QSDY"])
  scores[1,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Time Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Place Score","QSDY"])
  scores[2,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Orientation to Place Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Registration Score","QSDY"])
  scores[3,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Registration Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Attention and Calc Score","QSDY"])
  scores[4,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Attention and Calc Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Recall Score","QSDY"])
  scores[5,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Recall Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Naming Score","QSDY"])
  scores[6,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Naming Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Repetition Score","QSDY"])
  scores[7,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Repetition Score","QSSTRESN"][tmp_order]
  
  tmp_order <- order(mmse_individual[mmse_individual[,"QSTEST"]=="Drawing Score","QSDY"])
  scores[8,]  <-  mmse_individual[mmse_individual[,"QSTEST"]=="Drawing Score","QSSTRESN"][tmp_order]
  
  mmse1014_total[i,1:length(days)] <- apply(scores,2,sum)

}

# create lists to store timepoints (tps; tau_i) and MMSE scores (data; x_i(tau_i))
tps <- list()
data <- list()

# index to allow us to drop individuals if have only 1 timepoint 
k <- 1

# Gather tps and data
for (i in 1:length(people1014)){

  tmp <- mmse1014_days[i,][1:length(which(!is.na(mmse1014_days[i,])))]
  
  if (length(tmp) > 1){
   
    tps[[k]] <- tmp - tmp[1]
    data[[k]] <- mmse1014_total[i,][1:length(which(!is.na(mmse1014_days[i,])))]
    
    final_people[k] <- people1014[i]
    
    k <- k + 1
    
  }

}

# create multiTS object for C-1014 MMSE
mmse1014ts <- new("multiTS", numTS = length(tps) , tps = tps , data = data, clusters = rep(4,length(tps)), subj_id = final_people, cohort = "camd-1014")

setwd("..")

print("loaded CAMD")


# Load and gather ADNI MMSE
#

setwd("adni_sep_2015")

MMSE <- read.table("MMSE.csv",sep=",",header=T)

MMSE <- MMSE[!(as.character(MMSE[,"EXAMDATE"])==""),] # Need a date, so have to exclude scores without a date

unq_rid <- unique(MMSE[,"RID"]) # participant IDs

final_people <- c()

num_ts <- length(unq_rid)

# create lists to store timepoints (time_of_visits; tau_i) and MMSE scores (data; x_i(tau_i))
time_of_visits <- vector("list")
data <- vector("list")

# index to allow us to drop individuals if have only 1 timepoint 
k <- 1

for (i in 1:num_ts){
  
  ind <- which(MMSE[,"RID"]==unq_rid[i])
  
  # index of missing data for this individual, to allow them to be removed later
  nas <- ( MMSE[ind,"MMSCORE"] == -4 | MMSE[ind,"MMSCORE"] == -1 )
  
  # if greater than one timepoint with non-missing data
  if (length(ind[!nas])>1){
    
    ind <- ind[!nas] # remove timepoints with missing data
    
    # dates in original format
    tmp_dates <- as.character(MMSE[ind,"EXAMDATE"])
    
    # store as Date objects
    tmp_f_dates <- as.Date(tmp_dates) 
    
    # time order for timepoints
    tmp_order <- order(tmp_f_dates)
    
    # call first timepoint time zero
    tmp_days <- c(0)
    
    if(length(ind)>1){
      
      # for all subsequent timepoints, store number of days since first visit
      for (j in 2:length(ind)){
        
        tmp_days <- c(tmp_days,as.numeric(tmp_f_dates[tmp_order[j]] - tmp_f_dates[tmp_order[1]]))
        
      }
      
    }
    
    time_of_visits[[k]] <- tmp_days
    
    # order data based on time order of timepoints
    data[[k]] <- MMSE[ind[tmp_order],"MMSCORE"]
    
    final_people[k] <- unq_rid[i]
    
    k <- k + 1
    
  }
  
}

# create multiTS object for ADNI MMSE
mmseADNIts <- new("multiTS", numTS = length(time_of_visits) , tps = time_of_visits , data = data , clusters = rep(1,length(time_of_visits) ), subj_id = final_people, cohort = "adni")

# remove duplicated time point for one individual
mmseADNIts@tps[[785]] <- mmseADNIts@tps[[785]][1:3] 
mmseADNIts@data[[785]] <- mmseADNIts@data[[785]][1:3] 

setwd("..")

print("loaded ADNI")



# Load and gather AIBL MMSE
#

setwd("aibl_28apr_2015")

MMSE <- read.table("aibl_mmse_28-Apr-2015.csv",sep=",",header=T,as.is=T)

unq_rid <- unique(MMSE[,"RID"])

final_people <- c()

num_ts <- length(unq_rid)

time_of_visits <- vector("list")
data <- vector("list")

k <- 1

for (i in 1:num_ts){
  
  ind <- which(MMSE[,"RID"]==unq_rid[i])
  
  # have to remove rows with missing score or dates
  nas <- ( MMSE[ind,"MMSCORE"] == -4 | MMSE[ind,"MMSCORE"] == -1 | MMSE[ind,"EXAMDATE"] == -4 | MMSE[ind,"EXAMDATE"] == -1 )
  
  if (length(ind[!nas])>1){
    
    ind <- ind[!nas]
    
    # dates in original format
    tmp_dates <- as.character(MMSE[ind,"EXAMDATE"])
    
    # dates loaded as Date objects
    tmp_f_dates <- as.Date(tmp_dates,"%m/%d/%Y")
    
    # time ordering of timepoints
    tmp_order <- order(tmp_f_dates)
    
    tmp_days <- c(0)
    
    if(length(ind)>1){
       
      # for all subsequent timepoints, store number of days since first visit
      for (j in 2:length(ind)){
        
        tmp_days <- c(tmp_days,as.numeric(tmp_f_dates[tmp_order[j]] - tmp_f_dates[tmp_order[1]]))
        
      }
      
    }
    
    time_of_visits[[k]] <- tmp_days
    
    # data added in time order
    data[[k]] <- MMSE[ind[tmp_order],"MMSCORE"]
    
    final_people[k] <- unq_rid[i]
    
    k <- k + 1
    
  }
  
}

# create multiTS object for AIBL MMSE
mmseAIBLts <- new("multiTS", numTS = length(time_of_visits) , tps = time_of_visits , data = data  ,clusters = rep(2,length(time_of_visits) ), subj_id = final_people, cohort = "aibl")

# remove duplicated timepoints in two individuals
mmseAIBLts@tps[[30]] <- mmseAIBLts@tps[[30]][1:3] 
mmseAIBLts@data[[30]] <- mmseAIBLts@data[[30]][1:3]

mmseAIBLts@tps[[177]] <- mmseAIBLts@tps[[177]][1:3] 
mmseAIBLts@data[[177]] <- mmseAIBLts@data[[177]][1:3]

print("loaded AIBL")


# Merge MMSE
mmse_both_cpath <- mergeMultiTS(mmse1013ts,mmse1014ts)
mmse_cpath_adni <- mergeMultiTS(mmse_both_cpath,mmseADNIts)
mmse_aac <- mergeMultiTS(mmse_cpath_adni,mmseAIBLts)

print("merged AAC")

# load PAQUID
install.packages("lcmm")
library(lcmm)
data(paquid)

unqs <- unique(paquid[,"ID"])

tps <- list()
data <- list()

final_people <- c()

k <- 1

# store demographic data
paq_age <- c()
paq_gender <- c()

for (i in 1:length(unqs)){

  ind <- which(paquid[,"ID"]==unqs[i])
  
  # remove individuals with only one timepoint
  if (length(ind) > 1){
  
    nas <- is.na(paquid[ind,"MMSE"])
    
    # remove individuals with only one non-missing timepoint
    if (length(which(!nas)) > 1){
    
      paq_age[k] <- paquid[ind[!nas][1],"age"]
      paq_gender[k] <- paquid[ind[!nas][1],"male"]
      
      data[[k]] <- paquid[ind[!nas],"MMSE"]
      tps[[k]] <- round((paquid[ind[!nas],"age"]-paquid[ind[!nas][1],"age"])*365.25)
  
      final_people[k] <- unqs[i]
    
      k <- k + 1
      
    }
    
  }

}

# create multiTS object for PAQUID MMSE
mmse_paq <- new("multiTS", numTS = k-1 , tps = tps , data = data, clusters = rep(1,length(tps)), subj_id = final_people , cohort = "paquid")

print("loaded PAQUID")

# add matrix format data
mmse_aac <- multiTS.matrix(mmse_aac)
mmse_paq <- multiTS.matrix(mmse_paq)

setwd("..")

# save loaded data
save(mmse_aac,mmse_paq,file="mmse.RData")

setwd("../scripts/")