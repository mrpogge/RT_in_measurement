EMupd=function(data,model,K=60,A=0,r=0,x=NULL,iterEM=100){
  data=data[order(data$observation),]
  
  n_resp=nrow(data)
  N=max(data$observation)
  n=max(data$item)
  
  position=c(match(c(1:N),data$observation)-1,n_resp)
  
  if(model==4){z=c(0,0.5,1,1.5)}
  if(model==5){z=c(A,1,1.5);z=z-z[1]}
  
  if(model==0){
    data$score=data$RA
  }
  if(model==1){
    data$score<-(2*data$RA-1)*(1-data$RT)
  }
  if(model==2){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==3){
    data$score<-(1-data$RA)*(-r)+data$RA*(1-data$RT)
  }
  if(model==4){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==5){
    data$score<-z[2]*data$RA*(data$RT>0.5)+z[3]*data$RA*(data$RT<=0.5)
  }
  
  weight<-(data%>%
             group_by(observation)%>%
             distinct(weight)%>%
             as.data.frame())[,2]
  
  th=gauss.quad(K,'hermite')
  w=log(th[[2]]/sqrt(pi))
  th=th[[1]]*sqrt(2)
  
  if(model>3){
    
    n_par=length(z)-1  
    
    scores_order <- z
    items_order <- c(1:n)  # Optional: define custom item order here
    
    RR=(data %>%
          group_by(item, score) %>%
          summarise(weighted_count = sum(weight), .groups = "drop") %>%
          mutate(score = factor(score, levels = scores_order),
                 item = factor(item, levels = items_order)) %>%
          arrange(item, score) %>%
          pivot_wider(names_from = score, values_from = weighted_count, values_fill = 0) %>%
          arrange(item)%>%
          as.matrix())[,-c(1:2)]
    
    if(is.null(x)){
      x=c(rep(0,n*n_par),1)
      
      res=EMnrm_npar_cpp(
        Position = as.integer(position),
        I = as.integer(data$item - 1),
        S = as.double(data$score),
        W = as.double(weight),
        th = as.double(th),
        x = as.double(x),
        w = as.double(w),
        N=as.integer(N),
        n=as.integer(n),
        K=as.integer(K),    
        max_iterEM=as.integer(1),
        max_iterNR=as.integer(1),
        scores=as.double(z[-1]),
        RR=as.double(c(RR)),
        tol=as.double(10^(-5)),
        n_par=as.integer(n_par))
      
      m=mean(res$EAP)  
      s=sqrt(var(res$EAP)+mean(res$Post_var))  
      
      x=c(rep(m*z[-1],n),s)
    }
    res=EMnrm_npar_cpp(
      Position = as.integer(position),
      I = as.integer(data$item - 1),
      S = as.double(data$score),
      W = as.double(weight),
      th = as.double(th),
      x = as.double(x),
      w = as.double(w),
      N=as.integer(N),
      n=as.integer(n),
      K=as.integer(K),    
      max_iterEM=as.integer(iterEM),
      max_iterNR=as.integer(100),
      scores=as.double(z[-1]),
      RR=as.double(c(RR)),
      tol=as.double(10^(-5)),
      n_par=as.integer(n_par))
  }
  
  if(model<4){
    if(is.null(x)){
      x=c(rep(0,n),1)
      res=EMupd_cpp(
        Position = as.integer(position),
        I = as.integer(data$item - 1),
        S = as.double(data$score),
        Weight = as.double(weight),
        th = as.double(th),
        x = as.double(x),
        w = as.double(w),
        N=as.integer(N),
        n=as.integer(n),
        K=as.integer(K),    
        model=as.integer(model),
        r=as.double(r),
        max_iterEM=as.integer(1),
        max_iterNR=as.integer(1),
        tol=as.double(10^(-5))
      )
      
      m=mean(res$EAP)
      s=sqrt(var(res$EAP)+mean(res$Post_var))
      
      x=c(rep(m,n),s)
    }
    res=EMupd_cpp(
      Position = as.integer(position),
      I = as.integer(data$item - 1),
      S = as.double(data$score),
      Weight  = as.double(weight),
      th = as.double(th),
      x = as.double(x),
      w = as.double(w),
      N=as.integer(N),
      n=as.integer(n),
      K=as.integer(K),    
      model=as.integer(model),
      r=as.double(r),
      max_iterEM=as.integer(iterEM),
      max_iterNR=as.integer(100),
      tol=as.double(10^(-5))
    )
    x=res$par
  }
  return(res)
}
elo=function(K,data,position,res,model,k=1,A=NULL){ 
  n_resp=nrow(data)
  N=max(data$person)
  n=max(data$item)
  
  AA=seq(0.5,by=-0.5,length=10)
  rr=rev(seq(0.1,by=0.2,length=10))
  
  if((model<3)|(model==4)){
    est=res[[model+1]]
  }
  if((model==3)|(model==5)){
    est=res[[model+1]][[k]]
  }
  A=AA[k]
  r=rr[k]
  
  z=c(0,0,0,0)
  if(model==4){z=c(0,0.5,1,1.5)}
  if(model==5){z=c(A,1,1.5);z=z-z[1]}
  if(model==0){
    data$score=data$RA
  }
  if(model==1){
    data$score<-(2*data$RA-1)*(1-data$RT)
  }
  if(model==2){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==3){
    data$score<-(1-data$RA)*(-r)+data$RA*(1-data$RT)
  }
  if(model==4){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==5){
    data$score<-z[2]*data$RA*(data$RT>0.5)+z[3]*data$RA*(data$RT<=0.5)
  }
  
  M=elo_cpp(
    position=as.integer(position),
    item_ids=as.integer(data$item-1),
    X=as.integer(data$RA),
    S=as.double(data$score),
    delta=as.double(c(est)),
    N=as.integer(N),
    K=as.double(K),
    z=as.double(z[-1]),
    r=as.double(r),
    model=as.integer(model)
  )
  return(M[[1]])
}  

elo2=function(K,data,position,res,model,k=1,A=NULL){ 
  n_resp=nrow(data)
  N=max(data$person)
  n=max(data$item)
  
  AA=seq(0.5,by=-0.5,length=10)
  rr=rev(seq(0.1,by=0.2,length=10))
  
  if((model<3)|(model==4)){
    est=res[[model+1]]
  }
  if((model==3)|(model==5)){
    est=res[[model+1]][[k]]
  }
  A=AA[k]
  r=rr[k]
  
  z=c(0,0,0,0)
  if(model==4){z=c(0,0.5,1,1.5)}
  if(model==5){z=c(A,1,1.5);z=z-z[1]}
  if(model==0){
    data$score=data$RA
  }
  if(model==1){
    data$score<-(2*data$RA-1)*(1-data$RT)
  }
  if(model==2){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==3){
    data$score<-(1-data$RA)*(-r)+data$RA*(1-data$RT)
  }
  if(model==4){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==5){
    data$score<-z[2]*data$RA*(data$RT>0.5)+z[3]*data$RA*(data$RT<=0.5)
  }
  
  M=elo_cpp(
    position=as.integer(position),
    item_ids=as.integer(data$item-1),
    X=as.integer(data$RA),
    S=as.double(data$score),
    delta=as.double(c(est)),
    N=as.integer(N),
    K=as.double(K),
    z=as.double(z[-1]),
    r=as.double(r),
    model=as.integer(model)
  )
  return(M[[4]])
}  


elo_save=function(K,data,res,model,k=1){ 
  
  n_resp=nrow(data)
  N=max(data$person)
  n=max(data$item)
  position=c(match(c(1:N),data$person)-1,n_resp)
  
  AA=seq(-1,0.8,by=0.2)
  rr=rev(seq(0.1,by=0.2,length=10))
  
  if((model<3)|(model==4)){
    est=res[[model+1]]
  }
  if((model==3)|(model==5)){
    est=res[[model+1]][[k]]
  }
  A=AA[k]
  r=rr[k]
  
  z=c(0,0,0,0)
  if(model==4){z=c(0,0.5,1,1.5)}
  if(model==5){z=c(A,1,1.5);z=z-z[1]}
  if(model==0){
    data$score=data$RA
  }
  if(model==1){
    data$score<-(2*data$RA-1)*(1-data$RT)
  }
  if(model==2){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==3){
    data$score<-(1-data$RA)*(-r)+data$RA*(1-data$RT)
  }
  if(model==4){
    data$score<-0.5*(data$RA==0)*(data$RT>0.5)+1*(data$RA==1)*(data$RT>0.5)+1.5*(data$RA==1)*(data$RT<=0.5)
  }
  if(model==5){
    data$score<-z[2]*data$RA*(data$RT>0.5)+z[3]*data$RA*(data$RT<=0.5)
  }
  
  M=elo_cpp(
    position=as.integer(position),
    item_ids=as.integer(data$item-1),
    X=as.integer(data$RA),
    S=as.double(data$score),
    delta=as.double(c(est)),
    N=as.integer(N),
    K=as.double(K),
    z=as.double(z[-1]),
    r=as.double(r),
    model=as.integer(model)
  )
  return(M[[3]])
}


# Generate SRT-model data using inverse-probability sampling
gen_SRT=function(a){
  # a is a vector of differences between ability and difficulty
  u=runif(length(a),0,1)
  s=log(u*exp(a)-u*exp(-a)+exp(-a))/a
  return(s)
}

# Generate DSRT-model data with penalty p. First sample accuracy, and then RT from its conditional using inverse probability sampling
gen_hSRT=function(a,p){
  # a is a vector of differences between ability and difficulty
  # p is the penalty for incorrect responses
  u_x=runif(length(a));
  x=1*(u_x<((exp(a)-1)/(a*exp(-p*a)+exp(a)-1)))
  t=x
  t[x==0]=runif(sum(x==0),0,1)
  u=runif(sum(x==1),0,1)
  t[x==1]=1-log(exp(a[x==1])-u*(exp(a[x==1])-1))/a[x==1]
  return(cbind(x,t,u_x))
}

# This function only samples accuracy under the DSRT model
gen_hSRTx=function(a,p){
  rbinom(length(a),1,(exp(a)-1)/(a*exp(-p*a)+exp(a)-1))
}

