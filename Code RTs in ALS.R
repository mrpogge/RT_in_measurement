# load necessary packages
library(dplyr)
library(mirt)
library(tidyr)
library(statmod)
library(Rcpp)
library(tictoc)
library(parallel)
library(dplyr)
library(purrr)

# Load all custom functions
sourceCpp('EMupd.cpp')
sourceCpp('EM_NRMupd.cpp')
source("functions_RTs_in_ALS.R")
sourceCpp('elo.cpp')

# Preparing data for analysis
for(J in 1:4){ # loop over domains
  
  # load the data from both the old (I=1) and the new (I=2) scoring rule periods
  I=1
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  I=2
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  
  # count the number of persons
  N_total1=length(unique(logs_t1$user_id))
  
  # Divide the sample in four folds
  if(N_total1%%4==0){
    fold1=rep(c(1:4),N_total1/4)
  }
  if(N_total1%%4!=0){
    fold1=c(rep(c(1:4),(N_total1-N_total1%%4)/4),c(1:(N_total1%%4)))
  }
  set.seed(1085018265)
  fold1=sample(fold1,N_total1)
  
  folds1=list(
    unique(logs_t1$user_id)[which(fold1==1)],
    unique(logs_t1$user_id)[which(fold1==2)],
    unique(logs_t1$user_id)[which(fold1==3)],
    unique(logs_t1$user_id)[which(fold1==4)]
  )
  
  # count the number of persons
  N_total2=length(unique(logs_t2$user_id))
  
  # Divide the sample in four folds
  if(N_total2%%4==0){
    fold2=rep(c(1:4),N_total2/4)
  }
  if(N_total2%%4!=0){
    fold2=c(rep(c(1:4),(N_total2-N_total2%%4)/4),c(1:(N_total2%%4)))
  }
  set.seed(1085018215)
  fold2=sample(fold1,N_total2)
  
  folds2=list(
    unique(logs_t2$user_id)[which(fold2==1)],
    unique(logs_t2$user_id)[which(fold2==2)],
    unique(logs_t2$user_id)[which(fold2==3)],
    unique(logs_t2$user_id)[which(fold2==4)]
  )
  
  # select calibration data
  
  Data=lapply(folds1,FUN=function(X){
    data_calib <- logs_t1 %>%
      filter(user_id %in% X) %>% # select the persons within a fold
      mutate(session = paste0(user_id, created)) %>% # create a session variable
      group_by(session) %>% # within each session
      distinct(item_id, .keep_all = TRUE) %>% # keep only the first response to an item 
      filter(n() >= 5 )%>% # remove sessions with less than 5 responses
      ungroup() %>%
      as.data.frame()
    return(data_calib)
  })
  
  # pairs of data folds which will be used for item calibration
  comb=list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  data_calib1=lapply(comb,FUN=function(X){rbind(Data[[X[1]]],Data[[X[2]]])})
  
  # select items with at least 300 responses in each of the data splits
  item_obs1=lapply(data_calib1,FUN=function(X){
    item_counts <- X %>%
      count(item_id)
    return(item_counts$item_id[item_counts$n>=300])
  })
  
  selected_itemIDs1<-Reduce(intersect,item_obs1)
  
  Data=lapply(folds2,FUN=function(X){
    data_calib <- logs_t2 %>%
      filter(user_id %in% X) %>% # select the persons within a fold
      mutate(session = paste0(user_id, created)) %>% # create a session variable
      group_by(session) %>% # within each session
      distinct(item_id, .keep_all = TRUE) %>% # keep only the first response to an item 
      filter(n() >= 5 )%>% # remove sessions with less than 5 responses
      ungroup() %>%
      as.data.frame()
    return(data_calib)
  })
  
  comb=list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  data_calib2=lapply(comb,FUN=function(X){rbind(Data[[X[1]]],Data[[X[2]]])})
  
  item_obs2=lapply(data_calib2,FUN=function(X){
    item_counts <- X %>%
      count(item_id)
    return(item_counts$item_id[item_counts$n>=300])
  })
  
  selected_itemIDs2<-Reduce(intersect,item_obs2)
  
  selected_itemIDs=intersect(selected_itemIDs1,selected_itemIDs2)
  
  # select relevant items, after which remove persons with less than 3 responses
  data_calib1=lapply(data_calib1,FUN=function(X){
    
    X<-X %>%
      filter(item_id %in% selected_itemIDs) %>%
      group_by(session) %>%
      filter(n() >= 3)%>%
      ungroup() %>%
      as.data.frame()
    
    return(X)
  })
  
  data_calib2=lapply(data_calib2,FUN=function(X){
    
    X<-X %>%
      filter(item_id %in% selected_itemIDs) %>%
      group_by(session) %>%
      filter(n() >= 3)%>%
      ungroup() %>%
      as.data.frame()
    
    return(X)
  })
  
  # Prepare the data for the Elo analysis (Steps 2 and 3)
  
  data_elo_folds1<-lapply(folds1,FUN=function(X){
    data_elo <- logs_t1%>%
      filter(user_id %in% X) %>%   # select only persons within a specific data fold
      filter(item_id %in% selected_itemIDs) %>%# select only relevant items
      group_by(user_id) %>%
      filter(n() >= 20) %>% # select persons with at least 20 responses
      ungroup() %>%
      as.data.frame()
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data_elo$item <- item_recode[as.character(data_elo$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data_elo$user_id)), unique(data_elo$user_id))
    data_elo$person <- pers_recode[as.character(data_elo$user_id)]
    
    data_elo<-data_elo %>%
      select(item,person,created,correct_answered,response_in_milliseconds) %>%
      as.data.frame()
    
    return(data_elo)
  })
  data_elo_folds2<-lapply(folds2,FUN=function(X){
    data_elo <- logs_t2 %>%
      filter(user_id %in% X) %>%   # select only persons within a specific data fold
      filter(item_id %in% selected_itemIDs) %>%# select only relevant items
      group_by(user_id) %>%
      filter(n() >= 20) %>% # select persons with at least 20 responses
      ungroup() %>%
      as.data.frame()
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data_elo$item <- item_recode[as.character(data_elo$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data_elo$user_id)), unique(data_elo$user_id))
    data_elo$person <- pers_recode[as.character(data_elo$user_id)]
    
    data_elo<-data_elo %>%
      select(item,person,created,correct_answered,response_in_milliseconds) %>%
      as.data.frame()
    return(data_elo)
  })
  
  # combine the data folds prepared for Elo analysis into 6 pairs of datasets
  
  comb2=lapply(comb,FUN=function(X){setdiff(c(1:4),X)})
  
  data_elo1=lapply(comb2,FUN=function(X){list(data_elo_folds1[[X[1]]],data_elo_folds1[[X[[2]]]])})
  data_elo2=lapply(comb2,FUN=function(X){list(data_elo_folds2[[X[1]]],data_elo_folds2[[X[[2]]]])})
  
  # Descriptives for the data set (number of responses, persons, items)
  initial1=c(nrow(logs_t1),length(unique(logs_t1$user_id)),length(unique(logs_t1$item_id)))
  initial2=c(nrow(logs_t2),length(unique(logs_t2$user_id)),length(unique(logs_t2$item_id)))
  print(initial1)
  print(initial2)
  
  
  save(initial1,initial2,data_calib1,data_calib2,data_elo1,data_elo2,selected_itemIDs,file=paste('clean_upd_dataset',J,'.RData',sep=''))
  
}

# Preparing data for sensitivity analysis:
# selecting the data for the Elo part with a threshold of 5 seconds and 10 seconds

# 1. Threshold for including students in steps 2 and 3 is 10 seconds
for(J in 1:4){ # loop over domains
  
  comb=list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  
  load(file=paste('clean_upd_dataset',J,'.RData',sep=''))
  
  # load the data from both the old (I=1) and the new (I=2) scoring rule periods
  I=1
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  I=2
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  
  # count the number of persons
  N_total1=length(unique(logs_t1$user_id))
  
  # Divide the sample in four folds
  if(N_total1%%4==0){
    fold1=rep(c(1:4),N_total1/4)
  }
  if(N_total1%%4!=0){
    fold1=c(rep(c(1:4),(N_total1-N_total1%%4)/4),c(1:(N_total1%%4)))
  }
  set.seed(1085018265)
  fold1=sample(fold1,N_total1)
  
  folds1=list(
    unique(logs_t1$user_id)[which(fold1==1)],
    unique(logs_t1$user_id)[which(fold1==2)],
    unique(logs_t1$user_id)[which(fold1==3)],
    unique(logs_t1$user_id)[which(fold1==4)]
  )
  
  # count the number of persons
  N_total2=length(unique(logs_t2$user_id))
  
  # Divide the sample in four folds
  if(N_total2%%4==0){
    fold2=rep(c(1:4),N_total2/4)
  }
  if(N_total2%%4!=0){
    fold2=c(rep(c(1:4),(N_total2-N_total2%%4)/4),c(1:(N_total2%%4)))
  }
  set.seed(1085018215)
  fold2=sample(fold1,N_total2)
  
  folds2=list(
    unique(logs_t2$user_id)[which(fold2==1)],
    unique(logs_t2$user_id)[which(fold2==2)],
    unique(logs_t2$user_id)[which(fold2==3)],
    unique(logs_t2$user_id)[which(fold2==4)]
  )
  
  # Prepare the data for the Elo analysis (Steps 2 and 3)
  
  data_elo_folds1<-lapply(folds1,FUN=function(X){
    data_elo <- logs_t1%>%
      filter(user_id %in% X) %>%   # select only persons within a specific data fold
      filter(item_id %in% selected_itemIDs) %>%# select only relevant items
      group_by(user_id) %>%
      filter(n() >= 10) %>% # select persons with at least 10 responses
      ungroup() %>%
      as.data.frame()
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data_elo$item <- item_recode[as.character(data_elo$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data_elo$user_id)), unique(data_elo$user_id))
    data_elo$person <- pers_recode[as.character(data_elo$user_id)]
    
    data_elo<-data_elo %>%
      select(item,person,created,correct_answered,response_in_milliseconds) %>%
      as.data.frame()
    
    return(data_elo)
  })
  data_elo_folds2<-lapply(folds2,FUN=function(X){
    data_elo <- logs_t2 %>%
      filter(user_id %in% X) %>%   # select only persons within a specific data fold
      filter(item_id %in% selected_itemIDs) %>%# select only relevant items
      group_by(user_id) %>%
      filter(n() >= 10) %>% # select persons with at least 10 responses
      ungroup() %>%
      as.data.frame()
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data_elo$item <- item_recode[as.character(data_elo$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data_elo$user_id)), unique(data_elo$user_id))
    data_elo$person <- pers_recode[as.character(data_elo$user_id)]
    
    data_elo<-data_elo %>%
      select(item,person,created,correct_answered,response_in_milliseconds) %>%
      as.data.frame()
    return(data_elo)
  })
  
  # combine the data folds prepared for Elo analysis into 6 pairs of datasets
  
  comb2=lapply(comb,FUN=function(X){setdiff(c(1:4),X)})
  
  data_elo1=lapply(comb2,FUN=function(X){list(data_elo_folds1[[X[1]]],data_elo_folds1[[X[[2]]]])})
  data_elo2=lapply(comb2,FUN=function(X){list(data_elo_folds2[[X[1]]],data_elo_folds2[[X[[2]]]])})
  
  save(data_elo1,data_elo2,file=paste('clean_upd_dataset',J,'threshold10sec.RData',sep=''))
  
}

# 2. Threshold for including students in steps 2 and 3 is 5 seconds
for(J in 1:4){ # loop over domains
  
  comb=list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  
  load(file=paste('/home/user-006/RTs in ALS/clean_upd_dataset',J,'.RData',sep=''))
  
  # load the data from both the old (I=1) and the new (I=2) scoring rule periods
  I=1
  load(file=paste("/home/user-006/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  I=2
  load(file=paste("/home/user-006/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  
  # count the number of persons
  N_total1=length(unique(logs_t1$user_id))
  
  # Divide the sample in four folds
  if(N_total1%%4==0){
    fold1=rep(c(1:4),N_total1/4)
  }
  if(N_total1%%4!=0){
    fold1=c(rep(c(1:4),(N_total1-N_total1%%4)/4),c(1:(N_total1%%4)))
  }
  set.seed(1085018265)
  fold1=sample(fold1,N_total1)
  
  folds1=list(
    unique(logs_t1$user_id)[which(fold1==1)],
    unique(logs_t1$user_id)[which(fold1==2)],
    unique(logs_t1$user_id)[which(fold1==3)],
    unique(logs_t1$user_id)[which(fold1==4)]
  )
  
  # count the number of persons
  N_total2=length(unique(logs_t2$user_id))
  
  # Divide the sample in four folds
  if(N_total2%%4==0){
    fold2=rep(c(1:4),N_total2/4)
  }
  if(N_total2%%4!=0){
    fold2=c(rep(c(1:4),(N_total2-N_total2%%4)/4),c(1:(N_total2%%4)))
  }
  set.seed(1085018215)
  fold2=sample(fold1,N_total2)
  
  folds2=list(
    unique(logs_t2$user_id)[which(fold2==1)],
    unique(logs_t2$user_id)[which(fold2==2)],
    unique(logs_t2$user_id)[which(fold2==3)],
    unique(logs_t2$user_id)[which(fold2==4)]
  )
  
  # Prepare the data for the Elo analysis (Steps 2 and 3)
  
  data_elo_folds1<-lapply(folds1,FUN=function(X){
    data_elo <- logs_t1%>%
      filter(user_id %in% X) %>%   # select only persons within a specific data fold
      filter(item_id %in% selected_itemIDs) %>%# select only relevant items
      group_by(user_id) %>%
      filter(n() >= 5) %>% # select persons with at least 5 responses
      ungroup() %>%
      as.data.frame()
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data_elo$item <- item_recode[as.character(data_elo$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data_elo$user_id)), unique(data_elo$user_id))
    data_elo$person <- pers_recode[as.character(data_elo$user_id)]
    
    data_elo<-data_elo %>%
      select(item,person,created,correct_answered,response_in_milliseconds) %>%
      as.data.frame()
    
    return(data_elo)
  })
  data_elo_folds2<-lapply(folds2,FUN=function(X){
    data_elo <- logs_t2 %>%
      filter(user_id %in% X) %>%   # select only persons within a specific data fold
      filter(item_id %in% selected_itemIDs) %>%# select only relevant items
      group_by(user_id) %>%
      filter(n() >= 5) %>% # select persons with at least 5 responses
      ungroup() %>%
      as.data.frame()
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data_elo$item <- item_recode[as.character(data_elo$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data_elo$user_id)), unique(data_elo$user_id))
    data_elo$person <- pers_recode[as.character(data_elo$user_id)]
    
    data_elo<-data_elo %>%
      select(item,person,created,correct_answered,response_in_milliseconds) %>%
      as.data.frame()
    return(data_elo)
  })
  
  # combine the data folds prepared for Elo analysis into 6 pairs of datasets
  
  comb2=lapply(comb,FUN=function(X){setdiff(c(1:4),X)})
  
  data_elo1=lapply(comb2,FUN=function(X){list(data_elo_folds1[[X[1]]],data_elo_folds1[[X[[2]]]])})
  data_elo2=lapply(comb2,FUN=function(X){list(data_elo_folds2[[X[1]]],data_elo_folds2[[X[[2]]]])})
  
  save(data_elo1,data_elo2,file=paste('clean_upd_dataset',J,'threshold5sec.RData',sep=''))
  
}

# Data summaries:
# compute the proportion of correct response in each dataset
for(J in 1:4){
  I=1
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  I=2
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  
  print(round(mean(logs_t1$correct_answered),3))
  print(round(mean(logs_t2$correct_answered),3))
}


# Main analysis

# Step 1
for(I in 1:4){
  print(I)
  load(file=paste('clean_upd_dataset',I,'.RData',sep=''))
  for(tt in 1:2){
    print(tt)
    data_list=list()
    # prepare the data
    for(J in 1:6){
      if(tt==1){data=data_calib1[[J]]}
      if(tt==2){data=data_calib2[[J]]}
      
      # Prepare the data for analysis
      deadline=20000
      
      data <- data %>%
        group_by(user_id) %>%
        mutate(weight = 1/n_distinct(session)) %>%
        ungroup() 
      
      # add item ID from 1 to the number of items
      item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
      data$item <- item_recode[as.character(data$item_id)]
      
      # add person ID from 1 to the number of persons
      pers_recode <- setNames(1:length(unique(data$session)), unique(data$session))
      data$observation <- pers_recode[as.character(data$session)]
      
      data<-data %>%
        rename(RA = correct_answered) %>%
        rename(RT = response_in_milliseconds) %>%
        mutate(RT = RT/deadline) %>%
        select(RA,RT,item,observation,weight) %>%
        as.data.frame()
      
      data_list[[J]]=data
    }
    # run the models
    tic()
    res0=mclapply(data_list,EMupd,model=0,mc.cores=6)
    toc()
    save(res0,file=paste('calib',tt,'_upd_dataset',I,'.RData',sep=''))
    
    tic()
    res1=mclapply(data_list,EMupd,model=1,mc.cores=6)
    toc()
    save(res0,res1,file=paste('calib',tt,'_upd_dataset',I,'.RData',sep=''))
    
    tic()
    res2=mclapply(data_list,EMupd,model=2,mc.cores=6)
    toc()
    save(res0,res1,res2,file=paste('calib',tt,'_upd_dataset',I,'.RData',sep=''))
    
    res3=list()
    r=rev(seq(0.1,by=0.2,length=10))
    for(i in 1:10){
      tic()
      res3[[i]]=mclapply(data_list,EMupd,model=3,r=r[i],mc.cores=6)
      toc()
      save(res0,res1,res2,res3,file=paste('calib',tt,'_upd_dataset',I,'.RData',sep=''))
    }
    
    
    tic()
    res4=mclapply(data_list,EMupd,model=4,mc.cores=6)
    toc()
    save(res0,res1,res2,res3,res4,file=paste('calib',tt,'_upd_dataset',I,'.RData',sep=''))
    
    res5=list()
    A= seq(0.5,by=-0.5,length=10)
    for(i in 1:10){
      tic()
      res5[[i]]=mclapply(data_list,EMupd,model=5,A=A[i],mc.cores=6)
      toc()
      save(res0,res1,res2,res3,res4,res5,file=paste('calib',tt,'_upd_dataset',I,'.RData',sep=''))
    }
  }
}  

# Steps 2 and 3

for(dataset in 1:4){
  
  load(file=paste('clean_upd_dataset',dataset,'.RData',sep=''))
  
  deadline=20000
  # preparing data of the old scoring rule for analysis
  load(file=paste('calib1_upd_dataset',dataset,'.RData',sep=''))
  
  n=max(data_elo1[[1]][[1]]$item)
  
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo1[[ncal]][[3]]=res
    
    data<-data_elo1[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[1]]<-data
    
    data<-data_elo1[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[2]]<-data
    
  }
  
  # preparing data of the new scoring rule for analysis
  load(file=paste('calib2_upd_dataset',dataset,'.RData',sep=''))
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo2[[ncal]][[3]]=res
    
    data<-data_elo2[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[1]]<-data
    
    data<-data_elo2[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[2]]<-data
    
  }
  
  res_elo1=list()
  for(ncal in 1:6){
    X=data_elo1[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    # Step 1
    k3<-k5<-c(0,0)
    selected_K<-matrix(ncol=7,nrow=2)
    
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      res_step1=list()
      for(i in 1:6){
        if((i!=4)&(i!=6)){
          tic()
          res_step1[[i]]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,interval=c(0.01,1),maximum=FALSE,tol=0.001))
          toc()
        }
        if((i==4)|(i==6)){
          res_step1[[i]]=matrix(ncol=10,nrow=2)
          for(j in 1:10){
            tic()
            res_step1[[i]][,j]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,k=j,interval=c(0.01,1),maximum=FALSE,tol=0.001))
            toc()
          }
        }
      }  
      k3[J]=which.min(res_step1[[4]][2,])
      k5[J]=which.min(res_step1[[6]][2,])
      
      for(i in 1:3){selected_K[J,i]=res_step1[[i]][[1]]}
      selected_K[J,4]=res_step1[[4]][1,k3[J]]
      selected_K[J,5]=res_step1[[5]][[1]]
      selected_K[J,6]=res_step1[[6]][1,k5[J]]
      selected_K[J,7]=res_step1[[4]][1,10]
      
    } 
    # Step 2
    res_step2=matrix(nrow=2,ncol=7)
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        res_step2[J,i]=elo(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      res_step2[J,4]=elo(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      res_step2[J,6]=elo(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      res_step2[J,7]=elo(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
    }
    
    res_elo1[[ncal]]=list(mse=res_step2,selected_K=selected_K,k3=k3,k5=k5)
    
  }  
  
  mse=rowMeans(sapply(res_elo1,FUN=function(X){X$mse}))
  save(mse,res_elo1,file=paste('res1_elo_upd_dataset',dataset,'.RData',sep=''))
  
  res_elo2=list()
  for(ncal in 1:6){
    X=data_elo2[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    # Step 1
    k3<-k5<-c(0,0)
    selected_K<-matrix(ncol=7,nrow=2)
    
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      res_step1=list()
      for(i in 1:6){
        if((i!=4)&(i!=6)){
          tic()
          res_step1[[i]]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,interval=c(0.01,1),maximum=FALSE,tol=0.001))
          toc()
        }
        if((i==4)|(i==6)){
          res_step1[[i]]=matrix(ncol=10,nrow=2)
          for(j in 1:10){
            tic()
            res_step1[[i]][,j]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,k=j,interval=c(0.01,1),maximum=FALSE,tol=0.001))
            toc()
          }
        }
      }  
      k3[J]=which.min(res_step1[[4]][2,])
      k5[J]=which.min(res_step1[[6]][2,])
      
      for(i in 1:3){selected_K[J,i]=res_step1[[i]][[1]]}
      selected_K[J,4]=res_step1[[4]][1,k3[J]]
      selected_K[J,5]=res_step1[[5]][[1]]
      selected_K[J,6]=res_step1[[6]][1,k5[J]]
      selected_K[J,7]=res_step1[[4]][1,10]
      
    } 
    # Step 2
    res_step2=matrix(nrow=2,ncol=7)
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        res_step2[J,i]=elo(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      res_step2[J,4]=elo(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      res_step2[J,6]=elo(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      res_step2[J,7]=elo(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
    }
    
    res_elo2[[ncal]]=list(mse=res_step2,selected_K=selected_K,k3=k3,k5=k5)
    
  }  
  
  mse=rowMeans(sapply(res_elo2,FUN=function(X){X$mse}))
  print(mse)
  save(mse,res_elo2,file=paste('res2_elo_upd_dataset',dataset,'.RData',sep=''))
}

# Repeat step 3 saving the MSE per response to look at the results stratified by the number of responses
for(dataset in 1:4){
  
  load(file=paste('clean_upd_dataset',dataset,'.RData',sep=''))
  load(file=paste('res2_elo_upd_dataset',dataset,'.RData',sep=''))
  load(file=paste('res1_elo_upd_dataset',dataset,'.RData',sep=''))
  
  deadline=20000
  # preparing data of the old scoring rule for analysis
  load(file=paste('calib1_upd_dataset',dataset,'.RData',sep=''))
  
  n=max(data_elo1[[1]][[1]]$item)
  
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo1[[ncal]][[3]]=res
    
    data<-data_elo1[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[1]]<-data
    
    data<-data_elo1[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[2]]<-data
    
  }
  
  # preparing data of the new scoring rule for analysis
  load(file=paste('calib2_upd_dataset',dataset,'.RData',sep=''))
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo2[[ncal]][[3]]=res
    
    data<-data_elo2[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[1]]<-data
    
    data<-data_elo2[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[2]]<-data
    
  }
  
  res1=list()
  for(ncal in 1:6){
    X=data_elo1[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    selected_K=res_elo1[[ncal]]$selected_K
    k3=res_elo1[[ncal]]$k3
    k5=res_elo1[[ncal]]$k5
    
    mse_per_response=list()
    for(J in 1:2){
      
      data=X[[J]]
      data<-data%>%
        group_by(person)%>%
        mutate(n_resp=c(1:n()))%>%
        ungroup()%>%
        mutate(n_resp=1*(n_resp>5)+1*(n_resp>10)+1*(n_resp>15)+1*(n_resp>20))
      
      n_resp=nrow(data)
      mse_per_response[[J]]=matrix(ncol=7,nrow=n_resp)
      
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        mse_per_response[[J]][,i]=elo2(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      mse_per_response[[J]][,4]=elo2(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      mse_per_response[[J]][,6]=elo2(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      mse_per_response[[J]][,7]=elo2(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
    
      mse_per_response[[J]]=apply(mse_per_response[[J]],2,FUN=function(X){tapply(X,data$n_resp,mean)})
      
    }
    res1[[ncal]]=mse_per_response
    
  }  
  
  lapply(res1,FUN=function(X){lapply(X,FUN=function(X){colMeans(X)[1:7]})})
  res2=list()
  for(ncal in 1:6){
    X=data_elo2[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    selected_K=res_elo2[[ncal]]$selected_K
    k3=res_elo2[[ncal]]$k3
    k5=res_elo2[[ncal]]$k5
    
    mse_per_response=list()
    for(J in 1:2){
      
      data=X[[J]]
      data<-data%>%
        group_by(person)%>%
        mutate(n_resp=c(1:n()))%>%
        ungroup()%>%
        mutate(n_resp=1*(n_resp>5)+1*(n_resp>10)+1*(n_resp>15)+1*(n_resp>20))
      
      n_resp=nrow(data)
      mse_per_response[[J]]=matrix(ncol=7,nrow=n_resp)
      
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        mse_per_response[[J]][,i]=elo2(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      mse_per_response[[J]][,4]=elo2(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      mse_per_response[[J]][,6]=elo2(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      mse_per_response[[J]][,7]=elo2(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
      
      mse_per_response[[J]]=apply(mse_per_response[[J]],2,FUN=function(X){tapply(X,data$n_resp,mean)})
      
    }
    res2[[ncal]]=mse_per_response
    
  }  
  
  
  
  save(res1,res2,file=paste0('results_stratified',dataset,'.RData'))
}

res=list()
j=0
for(dataset in 1:4){
j=j+1
load(file=paste0('results_stratified',dataset,'.RData'))
av_res1=matrix(rowMeans(sapply(res1,FUN=function(X){rowMeans(matrix(unlist(X),ncol=2))})),nrow=5)
av_res2=matrix(rowMeans(sapply(res2,FUN=function(X){rowMeans(matrix(unlist(X),ncol=2))})),nrow=5)
av_res1=apply(av_res1,2,FUN=function(X){X-av_res1[,1]})[,-1]
av_res2=apply(av_res2,2,FUN=function(X){X-av_res2[,1]})[,-1]
res[[j]]=av_res1
j=j+1
res[[j]]=av_res2
}

for(i in 1:8){
  res[[i]]=res[[i]][,c(1,2,6,3,4,5)]
}

for(i in 1:8){
res[[i]]=data.frame(c(res[[i]]),rep(c('1-5','6-10','11-15','16-20','21+'),6),rep(c('SRTM','DSRTM','CRTM P=0.1','CRTM','PCM','NRM'),each=5))
}

res_long=NULL
for(i in 1:8){
  res_long=rbind(res_long,res[[i]])
}

res_long<-cbind(res_long,rep(c('+','-','x',':'),each=2*30),rep(rep(c('Old scoring rule','New scoring rule'),each=30),4))

your_data_frame<-res_long

names(your_data_frame)<-c('MSE','# responses','Model','Domain','Scoring rule')
# Ensure levels respect the data's original order
your_data_frame$`# responses` <- factor(your_data_frame$`# responses`, levels = unique(your_data_frame$`# responses`))
your_data_frame$Model        <- factor(your_data_frame$Model, levels = unique(your_data_frame$Model))
your_data_frame$Domain       <- factor(your_data_frame$Domain, levels = unique(your_data_frame$Domain))
your_data_frame$`Scoring rule` <- factor(your_data_frame$`Scoring rule`, levels = unique(your_data_frame$`Scoring rule`))

# Create the plot
ggplot(your_data_frame, aes(x = `# responses`, y = MSE, color = Model, group = Model)) +
  
  # 1. Model results plotted first (in the background)
  geom_line(linewidth = 0.6) + 
  geom_point(size = 1.2) +
  
  # 2. Reference line plotted second (now layers on top of the data)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  
  # 3. Grid setup with independent Y-axes
  facet_grid(`Scoring rule`~Domain, scales = "free_y") +
  
  # 4. Styling and labels
  theme_bw() + 
  scale_color_brewer(palette = "Set1") + 
  labs(
    x = "Number of Responses",
    y = "MSE difference from the Rasch model",
    color = "Model"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.text.y = element_text(angle = 0), # Keeps domain labels horizontal
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Data summaries:
# computing the proportion of removed responses

D1=NULL
D2=NULL

for(J in 1:4){ # loop over domains
  
  # load the data from both the old (I=1) and the new (I=2) scoring rule periods
  I=1
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  I=2
  load(file=paste("~/research-collaboration/data-ro/data-old-new-scoring-rule/update_domain_",J,"_2025-10-06-t",I,".Rdata",sep=''))
  
  # count the number of persons
  N_total1=length(unique(logs_t1$user_id))
  
  # Divide the sample in four folds
  if(N_total1%%4==0){
    fold1=rep(c(1:4),N_total1/4)
  }
  if(N_total1%%4!=0){
    fold1=c(rep(c(1:4),(N_total1-N_total1%%4)/4),c(1:(N_total1%%4)))
  }
  set.seed(1085018265)
  fold1=sample(fold1,N_total1)
  
  folds1=list(
    unique(logs_t1$user_id)[which(fold1==1)],
    unique(logs_t1$user_id)[which(fold1==2)],
    unique(logs_t1$user_id)[which(fold1==3)],
    unique(logs_t1$user_id)[which(fold1==4)]
  )
  
  # count the number of persons
  N_total2=length(unique(logs_t2$user_id))
  
  # Divide the sample in four folds
  if(N_total2%%4==0){
    fold2=rep(c(1:4),N_total2/4)
  }
  if(N_total2%%4!=0){
    fold2=c(rep(c(1:4),(N_total2-N_total2%%4)/4),c(1:(N_total2%%4)))
  }
  set.seed(1085018215)
  fold2=sample(fold1,N_total2)
  
  folds2=list(
    unique(logs_t2$user_id)[which(fold2==1)],
    unique(logs_t2$user_id)[which(fold2==2)],
    unique(logs_t2$user_id)[which(fold2==3)],
    unique(logs_t2$user_id)[which(fold2==4)]
  )
  
  Data=lapply(folds1,FUN=function(X){
    data_calib <- logs_t1 %>%
      filter(user_id %in% X) %>% # select the persons within a fold
      mutate(session = paste0(user_id, created)) %>% # create a session variable
      group_by(session) %>% # within each session
      distinct(item_id, .keep_all = TRUE) %>% # keep only the first response to an item 
      filter(n() >= 5 )%>% # remove sessions with less than 5 responses
      ungroup() %>%
      as.data.frame()
    return(data_calib)
  })
  
  Data0=lapply(folds1,FUN=function(X){
    data_calib <- logs_t1 %>%
      filter(user_id %in% X) %>% # select the persons within a fold
      mutate(session = paste0(user_id, created)) %>% # create a session variable
      as.data.frame()
    return(data_calib)
  })
  
  descr1=c(sum(sapply(Data0,FUN=function(X){length(unique(X$session))})), sum(sapply(Data,FUN=function(X){length(unique(X$session))})),
           sum(sapply(Data0,FUN=function(X){length(unique(X$user_id))})), sum(sapply(Data,FUN=function(X){length(unique(X$user_id))})))
  
  Data=lapply(folds2,FUN=function(X){
    data_calib <- logs_t2 %>%
      filter(user_id %in% X) %>% # select the persons within a fold
      mutate(session = paste0(user_id, created)) %>% # create a session variable
      group_by(session) %>% # within each session
      distinct(item_id, .keep_all = TRUE) %>% # keep only the first response to an item 
      filter(n() >= 5 )%>% # remove sessions with less than 5 responses
      ungroup() %>%
      as.data.frame()
    return(data_calib)
  })
  
  Data0=lapply(folds2,FUN=function(X){
    data_calib <- logs_t2 %>%
      filter(user_id %in% X) %>% # select the persons within a fold
      mutate(session = paste0(user_id, created)) %>% # create a session variable
      as.data.frame()
    return(data_calib)
  })
  
  descr2=c(sum(sapply(Data0,FUN=function(X){length(unique(X$session))})), sum(sapply(Data,FUN=function(X){length(unique(X$session))})),
           sum(sapply(Data0,FUN=function(X){length(unique(X$user_id))})), sum(sapply(Data,FUN=function(X){length(unique(X$user_id))})))
  
  load(file=paste0("clean_upd_dataset",J,".RData"))
  
  descr1=c(descr1,length(unique(data_calib1[[1]]$session))+length(unique(data_calib1[[6]]$session)),
           length(unique(data_calib1[[1]]$user_id))+length(unique(data_calib1[[6]]$user_id)))
  
  descr2=c(descr2,length(unique(data_calib2[[1]]$session))+length(unique(data_calib2[[6]]$session)),
           length(unique(data_calib2[[1]]$user_id))+length(unique(data_calib2[[6]]$user_id)))
  
  D1=cbind(D1,descr1)
  D2=cbind(D2,descr2)
}

# % removed pseudo-respondents, <5 resp
range(c((D1[1,]-D1[2,])/D1[1,],(D2[1,]-D2[2,])/D2[1,]))
# % removed students, <5 resp per session
range(c((D1[3,]-D1[4,])/D1[3,],(D2[3,]-D2[4,])/D2[3,]))

# additional % removed pseudo-respondents, <3 resp
range(c((D1[2,]-D1[5,])/D1[1,],(D2[2,]-D2[5,])/D2[1,]))
# additional % removed students, <3 resp per session
range(c((D1[4,]-D1[6,])/D1[3,],(D2[4,]-D2[6,])/D2[3,]))

# Data summaries:
# making plots for the number of responses per person

dat=NULL
domains=c('+','-','x',':')
for(j in 1:4){
  load(file=names[j])
  for(I in c(1,6)){
    for(J in c(1,2)){
      data1=data.frame(n_resp=(data_elo1[[I]][[J]]%>%group_by(person)%>%summarise(n_resp=n())%>%select(n_resp))$n_resp,
                domain=domains[j],time=1)
      data2=data.frame(n_resp=(data_elo2[[I]][[J]]%>%group_by(person)%>%summarise(n_resp=n())%>%select(n_resp))$n_resp,
                       domain=domains[j],time=2)
      dat=rbind(dat,data1,data2)
    }
  }
}

dat$domain <- factor(
  dat$domain,
  levels = c("+", "-", "x", ":")
)

ggplot(
  dat,
  aes(
    x = n_resp,
    color = domain,
    linetype = factor(time),
    group = interaction(domain, time)
  )
) +
  stat_ecdf(
    geom = "step",
    linewidth = 0.6,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = c(
      "+" = "orange",  # green
      "-" = "darkblue",  # blue
      "x" = "lightblue",  # pink/purple
      ":" = "green"   # orange
    ),
    name = "Domain"
  ) +
  scale_linetype_manual(
    values = c(
      `1` = "solid",
      `2` = "dashed"
    ),
    labels = c(
      `1` = "Old scoring rule",
      `2` = "New scoring rule"
    ),
    name = "Dataset"
  ) +
  scale_x_continuous(
    breaks = seq(0, 500, by = 50)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00)
  ) +
  coord_cartesian(
    xlim = c(0, 500),
    ylim = c(0, 1)
  ) +
  labs(
    x = "Number of responses",
    y = "Cumulative proportion"
  ) +
  theme_classic() +
  theme(
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    # Horizontal grid only
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.5
    ),
    panel.grid.minor.y = element_blank(),
    
    # No vertical grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    legend.position = "right"
  )


# Sensitivity analysis for steps 2 and 3

for(dataset in 1:4){
  
  load(file=paste('clean_upd_dataset',dataset,'threshold10sec.RData',sep=''))
  
  deadline=20000
  # preparing data of the old scoring rule for analysis
  load(file=paste('calib1_upd_dataset',dataset,'.RData',sep=''))
  
  n=max(data_elo1[[1]][[1]]$item)
  
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo1[[ncal]][[3]]=res
    
    data<-data_elo1[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[1]]<-data
    
    data<-data_elo1[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[2]]<-data
    
  }
  
  # preparing data of the new scoring rule for analysis
  load(file=paste('calib2_upd_dataset',dataset,'.RData',sep=''))
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo2[[ncal]][[3]]=res
    
    data<-data_elo2[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[1]]<-data
    
    data<-data_elo2[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[2]]<-data
    
  }
  
  res_elo1=list()
  for(ncal in 1:6){
    X=data_elo1[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    # Step 1
    k3<-k5<-c(0,0)
    selected_K<-matrix(ncol=7,nrow=2)
    
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      res_step1=list()
      for(i in 1:6){
        if((i!=4)&(i!=6)){
          tic()
          res_step1[[i]]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,interval=c(0.01,1),maximum=FALSE,tol=0.001))
          toc()
        }
        if((i==4)|(i==6)){
          res_step1[[i]]=matrix(ncol=10,nrow=2)
          for(j in 1:10){
            tic()
            res_step1[[i]][,j]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,k=j,interval=c(0.01,1),maximum=FALSE,tol=0.001))
            toc()
          }
        }
      }  
      k3[J]=which.min(res_step1[[4]][2,])
      k5[J]=which.min(res_step1[[6]][2,])
      
      for(i in 1:3){selected_K[J,i]=res_step1[[i]][[1]]}
      selected_K[J,4]=res_step1[[4]][1,k3[J]]
      selected_K[J,5]=res_step1[[5]][[1]]
      selected_K[J,6]=res_step1[[6]][1,k5[J]]
      selected_K[J,7]=res_step1[[4]][1,10]
      
    } 
    # Step 2
    res_step2=matrix(nrow=2,ncol=7)
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        res_step2[J,i]=elo(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      res_step2[J,4]=elo(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      res_step2[J,6]=elo(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      res_step2[J,7]=elo(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
    }
    
    res_elo1[[ncal]]=list(mse=res_step2,selected_K=selected_K,k3=k3,k5=k5)
    
  }  
  
  mse=rowMeans(sapply(res_elo1,FUN=function(X){X$mse}))
  save(mse,res_elo1,file=paste('res1_elo_upd_dataset',dataset,'threshold10sec.RData',sep=''))
  
  res_elo2=list()
  for(ncal in 1:6){
    X=data_elo2[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    # Step 1
    k3<-k5<-c(0,0)
    selected_K<-matrix(ncol=7,nrow=2)
    
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      res_step1=list()
      for(i in 1:6){
        if((i!=4)&(i!=6)){
          tic()
          res_step1[[i]]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,interval=c(0.01,1),maximum=FALSE,tol=0.001))
          toc()
        }
        if((i==4)|(i==6)){
          res_step1[[i]]=matrix(ncol=10,nrow=2)
          for(j in 1:10){
            tic()
            res_step1[[i]][,j]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,k=j,interval=c(0.01,1),maximum=FALSE,tol=0.001))
            toc()
          }
        }
      }  
      k3[J]=which.min(res_step1[[4]][2,])
      k5[J]=which.min(res_step1[[6]][2,])
      
      for(i in 1:3){selected_K[J,i]=res_step1[[i]][[1]]}
      selected_K[J,4]=res_step1[[4]][1,k3[J]]
      selected_K[J,5]=res_step1[[5]][[1]]
      selected_K[J,6]=res_step1[[6]][1,k5[J]]
      selected_K[J,7]=res_step1[[4]][1,10]
      
    } 
    # Step 2
    res_step2=matrix(nrow=2,ncol=7)
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        res_step2[J,i]=elo(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      res_step2[J,4]=elo(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      res_step2[J,6]=elo(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      res_step2[J,7]=elo(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
    }
    
    res_elo2[[ncal]]=list(mse=res_step2,selected_K=selected_K,k3=k3,k5=k5)
    
  }  
  
  mse=rowMeans(sapply(res_elo2,FUN=function(X){X$mse}))
  print(mse)
  save(mse,res_elo2,file=paste('res2_elo_upd_dataset',dataset,'threshold10sec.RData',sep=''))
}

for(dataset in 1:4){
  
  load(file=paste('clean_upd_dataset',dataset,'threshold5sec.RData',sep=''))
  
  deadline=20000
  # preparing data of the old scoring rule for analysis
  load(file=paste('calib1_upd_dataset',dataset,'.RData',sep=''))
  
  n=max(data_elo1[[1]][[1]]$item)
  
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo1[[ncal]][[3]]=res
    
    data<-data_elo1[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[1]]<-data
    
    data<-data_elo1[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo1[[ncal]][[2]]<-data
    
  }
  
  # preparing data of the new scoring rule for analysis
  load(file=paste('calib2_upd_dataset',dataset,'.RData',sep=''))
  for(ncal in 1:6){
    
    res=list(res0[[ncal]]$par[1:n],res1[[ncal]]$par[1:n],res2[[ncal]]$par[1:n],res3[[ncal]],
             matrix(res4[[ncal]]$par[1:(3*n)],nrow=3),res5[[ncal]])
    
    res[[4]]=list()
    for(i in 1:10){
      res[[4]][[i]]=res3[[i]][[ncal]]$par[1:n]
    }
    res[[6]]=list()
    for(i in 1:10){
      res[[6]][[i]]=matrix(res5[[i]][[ncal]]$par[1:(2*n)],nrow=2)
    }
    
    data_elo2[[ncal]][[3]]=res
    
    data<-data_elo2[[ncal]][[1]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[1]]<-data
    
    data<-data_elo2[[ncal]][[2]]
    data <- data %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)%>%
      rename(RT=response_in_milliseconds)%>%
      mutate(RT=RT/deadline)
    data_elo2[[ncal]][[2]]<-data
    
  }
  
  res_elo1=list()
  for(ncal in 1:6){
    X=data_elo1[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    # Step 1
    k3<-k5<-c(0,0)
    selected_K<-matrix(ncol=7,nrow=2)
    
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      res_step1=list()
      for(i in 1:6){
        if((i!=4)&(i!=6)){
          tic()
          res_step1[[i]]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,interval=c(0.01,1),maximum=FALSE,tol=0.001))
          toc()
        }
        if((i==4)|(i==6)){
          res_step1[[i]]=matrix(ncol=10,nrow=2)
          for(j in 1:10){
            tic()
            res_step1[[i]][,j]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,k=j,interval=c(0.01,1),maximum=FALSE,tol=0.001))
            toc()
          }
        }
      }  
      k3[J]=which.min(res_step1[[4]][2,])
      k5[J]=which.min(res_step1[[6]][2,])
      
      for(i in 1:3){selected_K[J,i]=res_step1[[i]][[1]]}
      selected_K[J,4]=res_step1[[4]][1,k3[J]]
      selected_K[J,5]=res_step1[[5]][[1]]
      selected_K[J,6]=res_step1[[6]][1,k5[J]]
      selected_K[J,7]=res_step1[[4]][1,10]
      
    } 
    # Step 2
    res_step2=matrix(nrow=2,ncol=7)
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        res_step2[J,i]=elo(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      res_step2[J,4]=elo(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      res_step2[J,6]=elo(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      res_step2[J,7]=elo(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
    }
    
    res_elo1[[ncal]]=list(mse=res_step2,selected_K=selected_K,k3=k3,k5=k5)
    
  }  
  
  mse=rowMeans(sapply(res_elo1,FUN=function(X){X$mse}))
  save(mse,res_elo1,file=paste('res1_elo_upd_dataset',dataset,'threshold5sec.RData',sep=''))
  
  res_elo2=list()
  for(ncal in 1:6){
    X=data_elo2[[ncal]]
    
    # Calibrated item parameters
    est=X[[3]]
    
    # Step 1
    k3<-k5<-c(0,0)
    selected_K<-matrix(ncol=7,nrow=2)
    
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      res_step1=list()
      for(i in 1:6){
        if((i!=4)&(i!=6)){
          tic()
          res_step1[[i]]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,interval=c(0.01,1),maximum=FALSE,tol=0.001))
          toc()
        }
        if((i==4)|(i==6)){
          res_step1[[i]]=matrix(ncol=10,nrow=2)
          for(j in 1:10){
            tic()
            res_step1[[i]][,j]=unlist(optimize(elo,data=data,position=position,res=est,model=i-1,k=j,interval=c(0.01,1),maximum=FALSE,tol=0.001))
            toc()
          }
        }
      }  
      k3[J]=which.min(res_step1[[4]][2,])
      k5[J]=which.min(res_step1[[6]][2,])
      
      for(i in 1:3){selected_K[J,i]=res_step1[[i]][[1]]}
      selected_K[J,4]=res_step1[[4]][1,k3[J]]
      selected_K[J,5]=res_step1[[5]][[1]]
      selected_K[J,6]=res_step1[[6]][1,k5[J]]
      selected_K[J,7]=res_step1[[4]][1,10]
      
    } 
    # Step 2
    res_step2=matrix(nrow=2,ncol=7)
    for(J in 1:2){
      data=X[[J]]
      
      n_resp=nrow(data)
      N=max(data$person)
      n=max(data$item)
      
      position=c(match(c(1:N),data$person)-1,n_resp)
      
      for(i in c(1,2,3,5)){
        res_step2[J,i]=elo(K=selected_K[3-J,i],data=data,position=position,res=est,model=i-1)
      }
      res_step2[J,4]=elo(K=selected_K[3-J,4],data=data,position=position,res=est,model=3,k=k3[3-J])
      res_step2[J,6]=elo(K=selected_K[3-J,6],data=data,position=position,res=est,model=5,k=k5[3-J])
      res_step2[J,7]=elo(K=selected_K[3-J,7],data=data,position=position,res=est,model=3,k=10)
    }
    
    res_elo2[[ncal]]=list(mse=res_step2,selected_K=selected_K,k3=k3,k5=k5)
    
  }  
  
  mse=rowMeans(sapply(res_elo2,FUN=function(X){X$mse}))
  print(mse)
  save(mse,res_elo2,file=paste('res2_elo_upd_dataset',dataset,'threshold5sec.RData',sep=''))
}

# Simulation study

# The design of the study introduces a condition for which the RTs are non-informative about ability.
# Here we generate RA under the RM, and the RTs we generate by sampling a value from the observed distribution of response times conditional on the item id and response accuracy, but doing so we remove any relationship to ability but preserve the distribution of RTs
# In this scenario we do not have true values for the RT-based model, so we cannot use the same design as before where we only do Steps 2 and 3, and skip Step 1 as we use simple the true values for the item parameters

# Creating input for the simulation (item and person parameters)

eap=function(x,i){
  dd=d[i]
  th=seq(-6,6,length=100)
  P=matrix(ncol=length(dd),nrow=length(th))
  for(i in 1:length(dd)){P[,i]=1/(1+exp(-dd[i]-th))}
  P=dnorm(th,0,sd=s)*apply(P,1,FUN=function(X){prod(X^x*(1-X)^(1-x))})
  P=P/sum(P)
  return(sum(P*th))
}


for(dataset in 1:4){
  tic()
  load(file=paste0('clean_upd_dataset',dataset,'.RData'))
  rm(data_elo1)
  rm(data_elo2)
  gc()
  
  for(I in 1:2){
    
    load(file=paste0('res',I,'_elo_upd_dataset',dataset,'.RData'))
    if(I==1){
      data=rbind(data_calib1[[1]],data_calib1[[6]])# [[1]] contains folds 1 & 2, and [[6]] contains 3 & 4, so combining them gives me everything
      #rm(data_calib1)
    }
    if(I==2){
      data=rbind(data_calib2[[1]],data_calib2[[6]])
      rm(data_calib2)
    }
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data$item <- item_recode[as.character(data$item_id)]
    
    tmp <- split(
      data$response_in_milliseconds,
      list(data$correct_answered + 1, data$item)
    )
    
    n <- max(data$item)
    
    RT_distributions <- lapply(1:2, function(i)
      lapply(seq_len(n), function(j)
        tmp[[paste0(i, ".", j)]]
      )
    )
    rm(tmp)
    
    load(file=paste0('calib',I,'_upd_dataset',dataset,'.RData'))
    est=rowMeans(sapply(res0,FUN=function(X){X[[1]]}))
    s=est[n+1]
    d=est[1:n]
    
    if(I==1){K=mean(sapply(res_elo1,FUN=function(X){X$selected_K[,1]}))}
    if(I==2){K=mean(sapply(res_elo2,FUN=function(X){X$selected_K[,1]}))}
    
    # for persons with >=20 observations we will run elo to estimate the learning trajectories
    data_elo<-data %>%
      group_by(user_id) %>%
      filter(n() >= 20) %>% # select persons with at least 20 responses
      ungroup() %>%
      as.data.frame()
    
    # for persons with <20 responses, we will obtain a single EAP estimate using all their responses
    data_eap<-data %>%
      group_by(user_id) %>%
      filter(n() < 20) %>% # select persons with less than 20 response
      ungroup() %>%
      as.data.frame()

    #rm(data)
    
    # First we deal with the persons with >=20 responses
    
    # recode persons from 1 to #persons
    pers_recode <- setNames(1:length(unique(data_elo$user_id)), unique(data_elo$user_id))
    data_elo$person <- pers_recode[as.character(data_elo$user_id)]
    
    # select relevant variables
    data_elo<-data_elo %>%
      select(item,person,created,correct_answered) %>%
      as.data.frame() %>%
      arrange(person, created)%>%
      rename(RA=correct_answered) %>%
      group_by(person) %>%
      group_modify(~ bind_rows(
        .x,
        .x[nrow(.x):1, ],
        .x
      )) %>%
      ungroup()
    
    data_elo<-as.data.frame(data_elo)
    
    # Run elo (three times for each person)
    th=elo_save(K=K,data=data_elo,res=list(d),model=0)
    
    data_elo$th=th
    
    data_elo <- data_elo %>%
      group_by(person) %>%
      slice((2 * n() / 3 + 1):n()) %>% # for each person select only the last third of their data (i.e., where we go through their responses with elo for the last time)
      ungroup()%>%
      group_by(created, person) %>%
      mutate(avg_th = mean(th, na.rm = TRUE))%>% # compute the average theta per session
      ungroup()%>%
      select(item,person,created,avg_th)%>%
      as.data.frame()
    
    # Second we deal with persons with less than 20 responses
    
    
    data_eap<-data_eap %>%
      group_by(user_id) %>%
      mutate(avg_th=eap(correct_answered,item))%>%
      ungroup() %>%
      as.data.frame()
    
    pers_recode <- setNames(1:length(unique(data_eap$user_id)), unique(data_eap$user_id))
    data_eap$person <- pers_recode[as.character(data_eap$user_id)]
    
    data_eap$person=data_eap$person+max(data_elo$person)
    
    data_eap <- data_eap %>%
      arrange(person, created)%>%
      rename(RA=correct_answered)
    
    data_eap<-data_eap%>%
      select(item,person,created,avg_th)
    
    data=rbind(data_elo,data_eap)
    
    rm(data_elo)
    rm(data_eap)
    
    
    avg_th=data$avg_th
    item=data$item
    
    data_structure<-data%>%
      select(person,created,item)
    rm(data)
    
    save(data_structure,file=paste0('simulation_datastructure_dataset',dataset,'_t',I,'.RData'))
    save(avg_th,item,RT_distributions,s,d,file=paste0('simulation_parameters_dataset',dataset,'_t',I,'.RData'))
    
    rm(item)
    rm(avg_th)
    rm(data_structure)
    rm(RT_distributions)
    
    gc()
  }
  toc()
}


# For CRTM P=0.1 we need to get parameters which would match the estimates under the RM

d=seq(-1,4,by=0.01)
t=rnorm(100000,0,1)
t=t-mean(t)
Y=matrix(ncol=length(d),nrow=length(t))
for(i in 1:length(d)){
  Y[,i]=rbinom(length(t),1,1/(1+exp(-t-d[i])))
}
library(mirt)
my_custom_function <- function(par, Theta,ncat=2) {
  prob <- (exp(Theta+par[1])-1)/((exp(Theta+par[1])-1)+(Theta+par[1])*exp(-0.1*(Theta+par[1])))
  return(cbind(1-prob,prob))
}
colnames(Y)<-c(1:length(d))
custom_item <- createItem(name = "CRTM",
                          par = c(d=0.001),
                          est = c(TRUE),
                          P = my_custom_function)
m <- mirt(Y, model=mirt.model('F1=1-501\n COV=F1*F1'), itemtype = rep("CRTM", 501),
          customItems = list(CRTM = custom_item))

d.=unlist(coef(m))[1:501]
s.=coef(m)[[502]][2] # variance of ability in the CRTM (the abilities in the RM need to be rescaled to generate matching CRTM data)

d1=d
d2=d^2
d3=d^3
d4=d^4

C=coef(lm(d.~d1+d2+d3+d4))

save(d.,d,s.,C,file='transform_CRTM.RData')

# Now generate data and run the whole analysis

load(file='transform_CRTM.RData')

# Run 10 replications
for(rep in 1:10){

for(dataset in 1:4){
  for(I in 1:2){
  tic()
  # First, the old scoring rule datasets. The normative RT model is the SRT    
      
      load(file=paste0('simulation_parameters_dataset',dataset,'_t',I,'.RData'))
      load(file=paste0('simulation_datastructure_dataset',dataset,'_t',I,'.RData'))
      #sample(c(1:length(avg_th)),10000)->P
      data=data.frame(avg_th=avg_th,item=item)
      data$d=d[data$item]
      
      
      if(I==1){
        data<-data%>%
          mutate(sim_s1=gen_SRT(avg_th+d),
               sim_x1=as.integer(1*(sim_s1>0)),
               sim_t=1-abs(sim_s1),
               sim_x2=as.integer(1*(gen_SRT(avg_th+d)>0)))
      }
      if(I==2){
        dCRT=C[1]+C[2]*d+C[3]*d^2+C[4]*d^3+C[5]*d^4
        data$dCRT=dCRT[data$item]
        
        crt_data=gen_hSRT(data$avg_th/s*sqrt(s.)+data$dCRT,p=0.1)
        
        data$sim_x1=crt_data[,1]
        data$sim_t=crt_data[,2]
        data$sim_x3=1*(crt_data[,3]<1/(1+exp(-data$d-data$avg_th)))
        data<-data%>%
          mutate(sim_x2=gen_hSRT(avg_th/s*sqrt(s.)+dCRT,0.1)[,1])
      }
      # Create logical vectors that randomly selects 10,30,50% of rows
      # first get 10%
      include_sim_x2_1<- sample(c(TRUE, FALSE),prob=c(0.1,0.9), size = nrow(data), replace = TRUE)
      include_sim_x2_2<-include_sim_x2_1
      # then extra 20%
      include_sim_x2_2[include_sim_x2_1==FALSE]<-sample(c(TRUE,FALSE),prob=c(0.2,0.7),size=sum(include_sim_x2_1==FALSE),replace=TRUE)
      include_sim_x2_3<-include_sim_x2_2
      # then extra 20%
      include_sim_x2_3[include_sim_x2_2==FALSE]<-sample(c(TRUE,FALSE),prob=c(0.2,0.5),size=sum(include_sim_x2_2==FALSE),replace=TRUE)
      
      # combine all information in one
      include_sim_x2=include_sim_x2_1+include_sim_x2_2+include_sim_x2_3
      # Add this vector to the data
      data <- data %>%
        mutate(include_sim_x2=include_sim_x2) 
      
      rm(include_sim_x2_1)
      rm(include_sim_x2_2)
      rm(include_sim_x2_3)
      rm(include_sim_x2)
      
      data<-data%>%
        rename(RA=sim_x1,RT=sim_t,RA2=sim_x2)
      
      if(I==2){
        data<-data%>%
          rename(RA_independent=sim_x3)
      }
      # in addition to the RTs generated under the SRTM, also sample values from the empirical RA- and item-specific distributions for the scenario where RTs are non-informative
      
      if(I==1){
      data <- data %>%
        group_by(RA, item) %>%
        mutate(
          RT_independent = sample(
            RT_distributions[[RA[1] + 1]][[item[1]]],
            n(),
            replace = TRUE
          )
        ) %>%
        ungroup()
      data=as.data.frame(data)    
      data<-data%>%select(RA,RT,RA2,include_sim_x2,RT_independent)
      }
      if(I==2){
        data <- data %>%
          group_by(RA_independent, item) %>%
          mutate(
            RT_independent = sample(
              RT_distributions[[RA_independent[1] + 1]][[item[1]]],
              n(),
              replace = TRUE
            )
          ) %>%
          ungroup()
        data=as.data.frame(data)    
        data<-data%>%select(RA,RT,RA2,include_sim_x2,RT_independent,RA_independent)
      }
      
      # dividing the data into folds
      
      # count the number of persons
      N_total=length(unique(data_structure$person))
      
      # Divide the sample in four folds
      if(N_total%%4==0){
        fold=rep(c(1:4),N_total/4)
      }
      if(N_total%%4!=0){
        fold=c(rep(c(1:4),(N_total-N_total%%4)/4),c(1:(N_total%%4)))
      }
      
      fold=sample(fold,N_total)
      
      persons=unique(data_structure$person)
      folds=list(
        persons[which(fold==1)],
        persons[which(fold==2)],
        persons[which(fold==3)],
        persons[which(fold==4)]
      )
      
      save(data,folds,file=paste0('simulated_data_dataset',dataset,'_t',I,'rep',rep,'.RData'))
      
      toc()
  }    
}

Design=cbind(rep(c(1:4),2),rep(c(1,2),each=4))

for(DD in 1:8){
dataset=Design[DD,1]
I=Design[DD,2]
  
    tic()
  load(file=paste0('simulated_data_dataset',dataset,'_t',I,'rep',rep,'.RData'))
  load(file=paste0('simulation_datastructure_dataset',dataset,'_t',I,'.RData'))
  
  data<-cbind(data,data_structure)

  if(I==1){
  data_calib<-data%>%
    mutate(session = paste0(person, created)) %>% # create a session variable
    group_by(session) %>% # within each session
    #distinct(item_id, .keep_all = TRUE) %>% # keep only the first response to an item. we do not do it here because the actual duplicates have been removed
    filter(n() >= 3 )%>% # remove sessions with less than 3 responses, since I will not be removing items, I can directly put 3 here
    ungroup() %>%
    select(RA,RT,item,person,created,RT_independent)%>%
    as.data.frame() 
  }
  if(I==2){
    data_calib<-data%>%
      mutate(session = paste0(person, created)) %>% # create a session variable
      group_by(session) %>% # within each session
      #distinct(item_id, .keep_all = TRUE) %>% # keep only the first response to an item. we do not do it here because the actual duplicates have been removed
      filter(n() >= 3 )%>% # remove sessions with less than 3 responses, since I will not be removing items, I can directly put 3 here
      ungroup() %>%
      select(RA,RT,item,person,created,RT_independent,RA_independent)%>%
      as.data.frame() 
  }
  Data=lapply(folds,FUN=function(X){
    Data <- data_calib %>%
      filter(person %in% X) %>% # select the persons within a fold
    return(Data)
  })
  rm(data_calib)
  
  comb=list(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
  
  data_calib=lapply(comb,FUN=function(X){rbind(Data[[X[1]]],Data[[X[2]]])})
  
  rm(Data)
  for(tt in c(1:2)){
    data_list=list()
    # prepare the data
    for(J in 1:6){
      data=data_calib[[J]]
      if(tt==2){
        data$RT=data$RT_independent/20000
        if(I==2){
          data$RA=data$RA_independent
        }
      }
      
      data <- data %>%
        mutate(session=paste0(person,created))%>%
        group_by(person) %>%
        mutate(weight = 1/n_distinct(session)) %>%
        ungroup() 
      
      # add person ID from 1 to the number of persons
      pers_recode <- setNames(1:length(unique(data$session)), unique(data$session))
      data$observation <- pers_recode[as.character(data$session)]
      
      data<-data %>%
        select(RA,RT,item,observation,weight) %>%
        as.data.frame()
      
      data_list[[J]]=data
    }
    # run the models
    res0=NULL
    n=max(data_list[[1]]$item)
    
    if(I==2){
      res1=mclapply(data_list,EMupd,model=3,r=0.1,x=c(rnorm(n,0,0.1),1),mc.cores=6)
    }
    if(I==1){
      res1=mclapply(data_list,EMupd,model=1,x=c(rnorm(n,0,0.1),1),mc.cores=6)
    }
    print(sapply(res1,FUN=function(X){X$iter}))
    res1=lapply(res1,FUN=function(X){X[1]})
    save(res0,res1,file=paste('calib',tt,'_simulated_dataset',dataset,'t',I,'rep',rep,'.RData',sep=''))
    
    if(tt==1){
      res0=mclapply(data_list,EMupd,model=0,x=c(rnorm(n,0,0.1),1),mc.cores=6)
      print(sapply(res0,FUN=function(X){X$iter}))
      res0=lapply(res0,FUN=function(X){X[1]})
      save(res0,res1,file=paste('calib',tt,'_simulated_dataset',dataset,'t',I,'rep',rep,'.RData',sep=''))
    }
    rm(data)
  }
  rm(data_calib)
  
  gc()
  
  # now we transform the data into a long format (RA and RA2 observations) and save those that are needed 
  load(file=paste0('simulated_data_dataset',dataset,'_t',I,'rep',rep,'.RData'))
  load(file=paste0('simulation_datastructure_dataset',dataset,'_t',I,'.RData'))
  
  data<-cbind(data,data_structure)
  
  data_long <- data %>%
    pivot_longer(cols = c(RA, RA2),
                 names_to = "source",
                 values_to = "RA") %>%
    filter(source == "RA" | (source == "RA2" & (include_sim_x2>0))) %>%
    arrange(person, created)
  
  # combine the information in source and include_sim_x2. Those with source=RA get values of at least 4
  data_long<-data_long%>%mutate(include_sim_x2=include_sim_x2+4*(source=='RA'))%>%
  mutate(session = paste0(person, created))%>% # create a session variable
  select(-source)
    
  rm(data)
  
  data_calib_long<-list()
  for(i in 1:3){
    data<-data_long
    if(i>1){
      data<-data%>%
      filter(include_sim_x2>(i-1))
    }
    data<-data%>%
      group_by(session) %>% # within each session
      filter(n() >= 3 )%>% # remove sessions with less than 3 responses, since I will not be removing items, I can directly put 3 here
      ungroup() %>%
      select(RA,RT,item,person,created,include_sim_x2)
    
    data_calib_long[[i]]<-lapply(folds,FUN=function(X){
      data_calib <- data %>%
        filter(person %in% X) %>% # select the persons within a fold
        as.data.frame()
      return(data_calib)
    })
    
    rm(data)
    gc()
  }
   
  for(tt in 3:5){
    
    data_list=lapply(comb,FUN=function(X){
      data=rbind(data_calib_long[[tt-2]][[X[1]]],data_calib_long[[tt-2]][[X[2]]])
      data <- data %>%
        mutate(session=paste0(person,created))%>%
        group_by(person) %>%
        mutate(weight = 1/n_distinct(session)) %>%
        ungroup() 
      
      # add person ID from 1 to the number of persons
      pers_recode <- setNames(1:length(unique(data$session)), unique(data$session))
      data$observation <- pers_recode[as.character(data$session)]
      
      data<-data %>%
        select(RA,RT,item,observation,weight) %>%
        as.data.frame()
    })
    
    res1=NULL
    res0=mclapply(data_list,EMupd,model=0,x=c(rnorm(n,0,0.1),1),mc.cores=6)
    print(sapply(res0,FUN=function(X){X$iter}))
    res0=lapply(res0,FUN=function(X){X[1]})
    save(res0,res1,file=paste('calib',tt,'_simulated_dataset',dataset,'t',I,'rep',rep,'.RData',sep=''))
  }
  rm(data_calib_long)
  gc()
  
  # Prepare the data for the Elo analysis (Steps 2 and 3)
  data_elo <- lapply(folds, FUN = function(X){
      
      data_fold <- data_long %>% 
        filter(person %in% X)
      
      eligible_persons <- data_fold %>%
        filter(include_sim_x2 > 3) %>%
        count(person) %>%
        filter(n >= 20) %>%
        pull(person)
      
      data_elo <- data_fold %>%
        filter(person %in% eligible_persons) %>%
        as.data.frame()
      
      pers_recode <- setNames(
        1:length(unique(data_elo$person)),
        unique(data_elo$person)
      )
      
      data_elo$person <- pers_recode[as.character(data_elo$person)]
      if(I==1){
        data_elo <- data_elo %>%
          select(item, person, created, RA, RT,include_sim_x2,RT_independent) %>%
          as.data.frame()
      }
      if(I==2){
        data_elo <- data_elo %>%
          select(item, person, created, RA, RT,include_sim_x2,RT_independent,RA_independent) %>%
          as.data.frame()
      }
      
      data_elo
    })
    
  rm(data_long)
  
  comb2=lapply(comb,FUN=function(X){setdiff(c(1:4),X)})
  data_elo=lapply(comb2,FUN=function(X){list(data_elo[[X[1]]],data_elo[[X[[2]]]])})
  
  res_elo=list()
  mse=list()
  
  n=max(data_elo[[1]][[1]]$item)
  
  tt=c(3,4,5,1,2)
  for(j in 1:5){ 
    for(ncal in 1:6){
      for(J in 1:2){
        if((j<5)&(j>1)){ # first take +50% data, then only +30% data, then only +10% data, and then only original responses
          data_elo[[ncal]][[J]]<-data_elo[[ncal]][[J]]%>%
            filter(include_sim_x2>(j-1))%>%
            as.data.frame()
        }
        if(j==5){
          data_elo[[ncal]][[J]]<-data_elo[[ncal]][[J]]%>%
            mutate(RT=RT_independent/20000)
          if(I==2){
            data_elo[[ncal]][[J]]$RA<-data_elo[[ncal]][[J]]$RA_independent
          }
        }  
      }
    }
    
    load(file=paste('calib',tt[j],'_simulated_dataset',dataset,'t',I,'rep',rep,'.RData',sep=''))
    
    res_elo[[j]]=list()
    for(ncal in 1:6){
      # get the estimates of the item parameters from the right data split
      if(j==4){
        est=list(res0[[ncal]]$par[1:n])
        est[[2]]=res1[[ncal]]$par[1:n]
        est[[4]]=list()
        est[[4]][[10]]=res1[[ncal]]$par[1:n]
      }
      if(j<4){
        est=list(res0[[ncal]]$par[1:n])
      }
      if(j==5){
        est=list(0,res1[[ncal]]$par[1:n])
        est[[4]]=list()
        est[[4]][[10]]=res1[[ncal]]$par[1:n]
      }
      
      # Step 2
      selected_K<-matrix(ncol=2,nrow=2)
      
      for(J in 1:2){
        data=data_elo[[ncal]][[J]]
        
        n_resp=nrow(data)
        N=max(data$person)
        position=c(match(c(1:N),data$person)-1,n_resp)
        
        if(j>3){
          if(I==1){
            selected_K[J,2]=unlist(optimize(elo,data=data,position=position,res=est,model=1,interval=c(0.01,1),maximum=FALSE,tol=0.001))[[1]]
          }
          if(I==2){
            selected_K[J,2]=unlist(optimize(elo,data=data,position=position,res=est,model=3,k=10,interval=c(0.01,1),maximum=FALSE,tol=0.001))[[1]]
          }
        }
        if(j<5){
          selected_K[J,1]=unlist(optimize(elo,data=data,position=position,res=est,model=0,interval=c(0.01,1),maximum=FALSE,tol=0.001))[[1]]
        }
      }  
      # Step 3
      res_step3=matrix(nrow=2,ncol=2)
      for(J in 1:2){
        data=data_elo[[ncal]][[J]]
        
        n_resp=nrow(data)
        N=max(data$person)
        
        position=c(match(c(1:N),data$person)-1,n_resp)
        
        if(j>3){
          if(I==1){
            res_step3[J,2]=elo(K=selected_K[3-J,2],data=data,position=position,res=est,model=1)
          }
          if(I==2){
            res_step3[J,2]=elo(K=selected_K[3-J,2],data=data,position=position,res=est,model=3,k=10)
          }
        }
        if(j<5){
          res_step3[J,1]=elo(K=selected_K[3-J,1],data=data,position=position,res=est,model=0)
        }
      }
      
      res_elo[[j]][[ncal]]=list(mse=res_step3,selected_K=selected_K)
    }  
    
    mse[[j]]=rowMeans(sapply(res_elo[[j]],FUN=function(X){colMeans(X$mse)}))
    
  }  
  
  
  save(res_elo,mse,file=paste0('res_simulated_dataset',dataset,'t',I,"rep",rep,'.RData'))
  
  rm(data_elo)
  gc()
  
  toc()
  }

}


# Collecting the results of the simulation study
Results=array(dim=c(5,8,10))
for(rep in 1:10){
for(DD in 1:8){
  dataset=Design[DD,1]
  I=Design[DD,2]
  load(file=paste0('res_simulated_dataset',dataset,'t',I,"rep",rep,'.RData'))
  Results[,DD,rep]=c(mse[[4]][1],c(mse[[4]][2],mse[[5]][2],mse[[3]][1],mse[[2]][1],mse[[1]][1])-mse[[4]][1])
}
}

# Naive baselines in the simulation

Naive=array(dim=c(4,2,10))

for(rep in 1:10){

for(dataset in 1:4){
  
  for(I in 1:2){
    Nb=NULL
    load(file=paste0("simulated_data_dataset",dataset,'_t',I,'rep',rep,".RData"))
    
    load(file=paste0('simulation_datastructure_dataset',dataset,'_t',I,'.RData'))
    
    data<-cbind(data,data_structure)
    
    data_elo <- lapply(folds, FUN = function(X){
      
      data_elo <- data %>% 
        filter(person %in% X)
      
      
      pers_recode <- setNames(
        1:length(unique(data_elo$person)),
        unique(data_elo$person)
      )
      
      data_elo$person <- pers_recode[as.character(data_elo$person)]
      if(I==1){
        data_elo <- data_elo %>%
          select(item, person, created, RA) %>%
          as.data.frame()
      }
      if(I==2){
        data_elo <- data_elo %>%
          select(item, person, created, RA) %>%
          as.data.frame()
      }
      
      data_elo
    })
    rm(data)
    
    comb2=lapply(comb,FUN=function(X){setdiff(c(1:4),X)})
    data_elo=lapply(comb2,FUN=function(X){list(data_elo[[X[1]]],data_elo[[X[[2]]]])})
    
    for(i in 1:6){
      for(j in 1:2){
        data=data_elo[[i]][[j]]
        P=mean(data$RA)
        nb=mean(tapply(data$RA,data$person,FUN=function(X){mean((X-P)^2)}))
        Nb=c(Nb,nb)
      }
    }
    
    Naive[dataset,I,rep]=mean(Nb)
    
  }
 
}

}  

# Collecting main results in terms of differences from the Benchmark
RES=NULL
meanRES=NULL
meanRES2=NULL
for(I in 1:4){
load(file=paste0("~/RTs in ALS/res1_elo_upd_dataset",I,".RData"))

cbind(sapply(res_elo1,FUN=function(X){X$mse})[1:7,],sapply(res_elo1,FUN=function(X){X$mse})[8:14,])->res
res[2:7,]=t(apply(res[2:7,],1,FUN=function(X){X-res[1,]}))
meanRES2=cbind(meanRES2,colMeans(matrix(mse,nrow=2)))
RES=cbind(RES,
c((round(apply(res,1,range),5))))
meanRES=cbind(meanRES,
              c((round(apply(res,1,mean),5))))
}

for(I in 1:4){
  load(file=paste0("~/RTs in ALS/res2_elo_upd_dataset",I,".RData"))
  
  cbind(sapply(res_elo2,FUN=function(X){X$mse})[1:7,],sapply(res_elo2,FUN=function(X){X$mse})[8:14,])->res
  res[2:7,]=t(apply(res[2:7,],1,FUN=function(X){X-res[1,]}))
  
  RES=cbind(RES,
            c((round(apply(res,1,range),5))))
  
  meanRES=cbind(meanRES,
                c((round(apply(res,1,mean),5))))
  meanRES2=cbind(meanRES2,colMeans(matrix(mse,nrow=2)))
}

rbind(
c(RES[1:2,])[1:8],
c(RES[3:4,])[1:8],
c(RES[5:6,])[1:8],
c(RES[7:8,])[1:8],
c(RES[9:10,])[1:8],
c(RES[11:12,])[1:8],
c(RES[13:14,])[1:8])

# Collecting the results of the sensitivity analysis

meanRESth5=NULL
meanRESth10=NULL
for(I in 1:4){
load(file=paste0("~/RTs in ALS/res1_elo_upd_dataset",I,"threshold5sec.RData"))

meanRESth5=cbind(meanRESth5,colMeans(matrix(mse,nrow=2)))

load(file=paste0("~/RTs in ALS/res1_elo_upd_dataset",I,"threshold10sec.RData"))

meanRESth10=cbind(meanRESth10,colMeans(matrix(mse,nrow=2)))
}

for(I in 1:4){
  load(file=paste0("~/RTs in ALS/res2_elo_upd_dataset",I,"threshold5sec.RData"))
  
  meanRESth5=cbind(meanRESth5,colMeans(matrix(mse,nrow=2)))
  
  load(file=paste0("~/RTs in ALS/res2_elo_upd_dataset",I,"threshold10sec.RData"))
  
  meanRESth10=cbind(meanRESth10,colMeans(matrix(mse,nrow=2)))
}

meanRESth5=apply(meanRESth5,2,FUN=function(X){c(X[1],X[-1]-X[1])})
meanRESth10=apply(meanRESth10,2,FUN=function(X){c(X[1],X[-1]-X[1])})



for(I in 1:4){
  load(file=paste0("~/RTs in ALS/res2_elo_upd_dataset",I,".RData"))
  
  cbind(sapply(res_elo2,FUN=function(X){X$mse})[1:7,],sapply(res_elo2,FUN=function(X){X$mse})[8:14,])->res
  res[2:7,]=t(apply(res[2:7,],1,FUN=function(X){X-res[1,]}))
  
  RES=cbind(RES,
            c((round(apply(res,1,range),5))))
  
  meanRES=cbind(meanRES,
                c((round(apply(res,1,mean),5))))
  meanRES2=cbind(meanRES2,colMeans(matrix(mse,nrow=2)))
}





### Looking at the correlations between the RTs of correct responses and ability
### Correlations between the RTs of incorrect responses and ability

diag_SRT=list()
for(I in 1:4){
  load(file=paste0("~/clean_upd_dataset",I,".RData"))
  
  diag_SRT[[I]]=list()
  
  for(G in 1:2){
    
    Cors_obs1=0
    Cors_obs0=0
    Cors_sim1=0
    Cors_sim0=0
    
    J=1
    
    
    if(G==1){data=data_calib1[[J]]
    load(file=paste0("~/calib1_upd_dataset",I,".RData"))
    }
    if(G==2){data=data_calib2[[J]]
    load(file=paste0("~/calib2_upd_dataset",I,".RData"))
    }
    # Prepare the data for analysis
    
    deadline=20000
    
    data <- data %>%
      group_by(user_id) %>%
      mutate(weight = 1/n_distinct(session)) %>%
      ungroup() 
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data$item <- item_recode[as.character(data$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data$session)), unique(data$session))
    data$observation <- pers_recode[as.character(data$session)]
    
    data<-data %>%
      rename(RA = correct_answered) %>%
      rename(RT = response_in_milliseconds) %>%
      mutate(RT = RT/deadline) %>%
      select(RA,RT,item,observation,weight) %>%
      as.data.frame()
    
    # Plausible values for ability under the SRT
    theta_pv=rnorm(length(res1[[1]]$EAP),res1[[1]]$EAP,sqrt(res1[[1]]$Post_var))
    
    data$pv=theta_pv[data$observation]
    # Estimates of the item parameters under the SRT
    data$d=res1[[1]]$par[data$item]
    
    # Generate data under the SRT model with the PVs of ability and difficulty estimates as true values
    data$simS=genSRT(data$pv+data$d)
    data$simX=1*(data$simS>0)
    data$simT=1-abs(data$simS)
    
    data_sim=data.frame(RA=data$simX,RT=data$simT,observation=data$observation,item=data$item,weight=data$weight)
    
    # Estimate the Rasch model for the simulated data and save the EAPs of ability
    res0sim=EMupd(data_sim,model=0)
    data_sim$eap=res0sim$EAP[data_sim$observation]
    n_resp=tapply(data,data$observation,nrow)
    
    #data_sim_select=data_sim[n_resp[data$observation]==10,]
    Cors_sim1=tapply(data_sim,data_sim$item,FUN=function(X){cor(X$eap[X$RA==1],X$RT[X$RA==1])})
    Cors_sim0=tapply(data_sim,data_sim$item,FUN=function(X){cor(X$eap[X$RA==0],X$RT[X$RA==0])})
    
    # get EAP estimates of ability under the Rasch model for the observed data
    data$eap=res0[[J]]$EAP[data$observation]
    #data_select=data[n_resp[data$observation]==10,]
    Cors_obs1=tapply(data,data$item,FUN=function(X){cor(X$eap[X$RA==1],X$RT[X$RA==1])})
    Cors_obs0=tapply(data,data$item,FUN=function(X){cor(X$eap[X$RA==0],X$RT[X$RA==0])})
    
    diag_SRT[[I]][[G]]=cbind(Cors_sim1,Cors_obs1,Cors_sim0,Cors_obs0)
    
  }
}

save(diag_SRT,file='res_cors.RData')

diag=list()
for(I in 1:4){
  load(file=paste0("~/clean_upd_dataset",I,".RData"))
  
  diag[[I]]=list()
  
  for(G in 1:2){
    Cors_obs1=0
    Cors_sim1=0
    
    J=1
    
    if(G==1){data=data_calib1[[J]]
    load(file=paste0("~/calib1_upd_dataset",I,".RData"))
    }
    if(G==2){data=data_calib2[[J]]
    load(file=paste0("~/calib2_upd_dataset",I,".RData"))
    }
    # Prepare the data for analysis
    
    deadline=20000
    
    data <- data %>%
      group_by(user_id) %>%
      mutate(weight = 1/n_distinct(session)) %>%
      ungroup() 
    
    # add item ID from 1 to the number of items
    item_recode <- setNames(1:length(selected_itemIDs), selected_itemIDs)
    data$item <- item_recode[as.character(data$item_id)]
    
    # add person ID from 1 to the number of persons
    pers_recode <- setNames(1:length(unique(data$session)), unique(data$session))
    data$observation <- pers_recode[as.character(data$session)]
    
    data<-data %>%
      rename(RA = correct_answered) %>%
      rename(RT = response_in_milliseconds) %>%
      mutate(RT = RT/deadline) %>%
      select(RA,RT,item,observation,weight) %>%
      as.data.frame()
    
    # Plausible values for ability under the SRT
    theta_pv=rnorm(length(res3[[10]][[1]]$EAP),res3[[10]][[1]]$EAP,sqrt(res3[[10]][[1]]$Post_var))
    
    data$pv=theta_pv[data$observation]
    # Estimates of the item parameters under the SRT
    data$d=res3[[10]][[1]]$par[data$item]
    
    # Generate data under the SRT model with the PVs of ability and difficulty estimates as true values
    sim=gen_hSRT(data$pv+data$d,p=0.1)
    data$simX=sim[,1]
    data$simT=sim[,2]
    
    data_sim=data.frame(RA=data$simX,RT=data$simT,observation=data$observation,item=data$item,weight=data$weight)
    
    # Estimate the Rasch model for the simulated data and save the EAPs of ability
    res0sim=EMupd(data_sim,model=0)
    data_sim$eap=res0sim$EAP[data_sim$observation]
    
    n_resp=tapply(data,data$observation,nrow)
    
    #data_sim_select=data_sim[n_resp[data$observation]==10,]
    Cors_sim1=tapply(data_sim,data_sim$item,FUN=function(X){cor(X$eap[X$RA==1],X$RT[X$RA==1])})
    
    
    diag[[I]][[G]]=Cors_sim1
    
  }
}

diag_CRT0.1=diag

save(diag_CRT0.1,diag_SRT,file='observed_expected_cors.RData')

# Descriptive statistics for the datasets


names=c("clean_upd_dataset1.RData","clean_upd_dataset2.RData",
        "clean_upd_dataset3.RData","clean_upd_dataset4.RData")

I1=NULL
I2=NULL
F1=NULL
F2=NULL

for(name in names){
  
  load(file=name)
  
  n_resp_total=mean(sapply(X=data_elo1,FUN=function(X){sapply(X,nrow)}))*2+mean((sapply(data_calib1,nrow)))
  n_resp_calib=mean((sapply(data_calib1,nrow)))
  n_resp_elo=mean(sapply(X=data_elo1,FUN=function(X){sapply(X,nrow)}))
  n_student_calib=mean(sapply(X=data_calib1,FUN=function(X){length(unique(X$user_id))}))
  median_n_resp_per_item=mean(sapply(X=data_calib1,FUN=function(X){median(table(X$item_id))}))
  #av_n_resp_per_item=mean(sapply(X=data_calib1,FUN=function(X){mean(table(X$item_id))}))
  #sd_n_resp_per_item=mean(sapply(X=data_calib1,FUN=function(X){sd(table(X$item_id))}))
  #min_n_resp_per_item=mean(sapply(X=data_calib1,FUN=function(X){min(table(X$item_id))}))
  n_items=length(selected_itemIDs)
  n_student_elo=mean(sapply(data_elo1,FUN=function(X){sapply(X,FUN=function(X){length(unique(X$person))})}))
  median_n_resp_per_pers=mean(sapply(data_elo1,FUN=function(X){sapply(X,FUN=function(X){median(table(X$person))})}))
  n_student_total=n_student_elo*2+n_student_calib
  
  final1=c(n_resp_total=n_resp_total,
           prop_resp_total=n_resp_total/initial1[1],
           n_resp_calib=n_resp_calib,
           n_resp_elo=n_resp_elo,
           n_student_total=n_student_total,
           prop_student_total=n_student_total/initial1[2],
           n_student_calib=n_student_calib,
           n_student_elo=n_student_elo,
           per_item=median_n_resp_per_item,
           per_pers=median_n_resp_per_pers,
           n_items=n_items)
  
  n_resp_total=mean(sapply(X=data_elo2,FUN=function(X){sapply(X,nrow)}))*2+mean((sapply(data_calib2,nrow)))
  n_resp_calib=mean((sapply(data_calib2,nrow)))
  n_resp_elo=mean(sapply(X=data_elo2,FUN=function(X){sapply(X,nrow)}))
  n_student_calib=mean(sapply(X=data_calib2,FUN=function(X){length(unique(X$user_id))}))
  median_n_resp_per_item=mean(sapply(X=data_calib2,FUN=function(X){median(table(X$item_id))}))
  #av_n_resp_per_item=mean(sapply(X=data_calib1,FUN=function(X){mean(table(X$item_id))}))
  #sd_n_resp_per_item=mean(sapply(X=data_calib1,FUN=function(X){sd(table(X$item_id))}))
  #min_n_resp_per_item=mean(sapply(X=data_calib1,FUN=function(X){min(table(X$item_id))}))
  n_items=length(selected_itemIDs)
  n_student_elo=mean(sapply(data_elo2,FUN=function(X){sapply(X,FUN=function(X){length(unique(X$person))})}))
  median_n_resp_per_pers=mean(sapply(data_elo2,FUN=function(X){sapply(X,FUN=function(X){median(table(X$person))})}))
  n_student_total=n_student_elo*2+n_student_calib
  
  final2=c(n_resp_total=n_resp_total,
           prop_resp_total=n_resp_total/initial2[1],
           n_resp_calib=n_resp_calib,
           n_resp_elo=n_resp_elo,
           n_student_total=n_student_total,
           prop_student_total=n_student_total/initial2[2],
           n_student_calib=n_student_calib,
           n_student_elo=n_student_elo,
           per_item=median_n_resp_per_item,
           per_pers=median_n_resp_per_pers,
           n_items=n_items)
  
  F1=rbind(F1,final1)
  F2=rbind(F2,final2)
  
  I1=rbind(I1,c(initial1))
  I2=rbind(I2,c(initial2))
  
}

# Naive baselines

names=c("clean_upd_dataset1.RData","clean_upd_dataset2.RData",
        "clean_upd_dataset3.RData","clean_upd_dataset4.RData")

Nb1<-Nb2<-NULL

for(name in names){
  load(file=name)
  
  for(i in 1:6){
    for(j in 1:2){
      data=data_elo1[[i]][[j]]
      P=mean(data$correct_answered)
      nb1=mean(tapply(data$correct_answered,data$person,FUN=function(X){mean((X-P)^2)}))
      Nb1=c(Nb1,nb1)
      data=data_elo2[[i]][[j]]
      P=mean(data$correct_answered)
      nb2=mean(tapply(data$correct_answered,data$person,FUN=function(X){mean((X-P)^2)}))
      Nb2=c(Nb2,nb2)
    }
  }
  
}

Nb1=matrix(Nb1,nrow=12)
Nb2=matrix(Nb2,nrow=12)

c(colMeans(Nb1),colMeans(Nb2))

# Plots with correlations

srt1=NULL
srt0=NULL
crt0.1=NULL
for(i in 1:4){
  for(j in 1:2){srt1=c(srt1,diag_SRT[[i]][[j]][,1])
  srt0=c(srt0,diag_SRT[[i]][[j]][,3])
  crt0.1=c(crt0.1,diag_CRT0.1[[i]][[j]])
  }
}
operations <- c("+ data", "- data", "× data", "÷ data")
scoring_rules <- c("Old rule data", "New rule data",'Simulated data')

ecdf_data <- purrr::map2_df(1:4, diag_SRT, function(i, op_list) {
  purrr::map2_df(1:2, op_list, function(j, mat) {
    tibble(
      value = mat[, 2],  # or [,4] if you prefer
      Operation = operations[i],
      ScoringRule = scoring_rules[j],
      Source = "Observed"
    )
  })
})
ref_data <- bind_rows(
  tibble(value = srt1, Model = "SRTM simulation"),
  tibble(value = crt0.1, Model = "CRTM simulation")
) %>%
  mutate(Source = "Model-based")

ggplot() +
  # Observed ECDFs
  stat_ecdf(data = ecdf_data, 
            aes(x = value, color = Operation, linetype = ScoringRule), 
            size = 1) +
  
  # Model-based ECDFs
  stat_ecdf(data = ref_data, 
            aes(x = value, color = Model), 
            size = 1.2, linetype = "solid") +
  
  # Axis and labels
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  labs(title = "Ability-RT correlations: Correct responses",
       x = "Value", y = "ECDF") +
  
  # Legends
  scale_color_manual(
    values = c("red", "blue", "darkgreen", "purple", "black", "gray70"),
    name = " "
  ) +
  scale_linetype_manual(
    values = c("dotted", "dashed",'solid'),
    name = " "
  ) +
  
  theme_minimal(base_size = 14)


operations <- c("+ data", "- data", "× data", "÷ data")
scoring_rules <- c("Old rule data", "New rule data",'Simulated data')

ecdf_data <- purrr::map2_df(1:4, diag_SRT, function(i, op_list) {
  purrr::map2_df(1:2, op_list, function(j, mat) {
    tibble(
      value = mat[, 4],  # or [,4] if you prefer
      Operation = operations[i],
      ScoringRule = scoring_rules[j],
      Source = "Observed"
    )
  })
})
ref_data <- bind_rows(
  tibble(value = srt0, Model = "SRTM simulation"),
) %>%
  mutate(Source = "Model-based")

ggplot() +
  # Observed ECDFs
  stat_ecdf(data = ecdf_data, 
            aes(x = value, color = Operation, linetype = ScoringRule), 
            size = 1) +
  
  # Model-based ECDFs
  stat_ecdf(data = ref_data, 
            aes(x = value, color = Model), 
            size = 1.2, linetype = "solid") +
  
  # Axis and labels
  coord_cartesian(xlim = c(-0.4, 0.4)) +
  labs(title = "Ability-RT correlations: Incorrect responses",
       x = "Value", y = "ECDF") +
  
  # Legends
  scale_color_manual(
    values = c("red", "blue", "darkgreen", "purple", "black", "gray70"),
    name = " "
  ) +
  scale_linetype_manual(
    values = c("dotted", "dashed",'solid'),
    name = " "
  ) +
  
  theme_minimal(base_size = 14)

operations <- c("+", "-", "\\times", "\\div")
scoring_rules <- c("Old", "New")

