
h_t=function(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,raw=F,inf_dat=NULL) {
  
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  surv_dat$lp=as.numeric(as.matrix(surv_dat[,Z_names])%*%beta)
  A=surv_dat$starts[surv_dat$status==1]
  Z_beta=as.numeric(as.matrix(surv_dat[surv_dat$status==1,Z_names])%*%beta)
  X=surv_dat$eventtime[surv_dat$status==1]
  LT=surv_dat$lt[surv_dat$status==1]
  A=A[order(X)]
  Z_beta=Z_beta[order(X)]
  LT=LT[order(X)]
  X=X[order(X)]
  surv_dat=surv_dat[order(surv_dat$eventtime),]
  X_indices=match(unique(surv_dat$eventtime[surv_dat$status==1]),surv_dat$eventtime)
  n_times=sapply(split(surv_dat$eventtime[surv_dat$status==1],surv_dat$eventtime[surv_dat$status==1]),length)
  X_indices=rep(X_indices,n_times)
  n=dim(surv_dat)[1]
  s_i=A+X-LT
  if(raw) {
    numerator=exp(alpha*inf_dat$Rhat[findInterval(s_i,inf_dat$s)]+Z_beta) 
  } else {
    numerator=exp(alpha*exp(theta0+theta1*s_i+theta2*(s_i^2)+theta3*(s_i^3))+Z_beta)
  }
  denominator=vector("numeric")
  for(i in 1:length(X)) {
    A_risk=surv_dat$starts[X_indices[i]:n]
    Z_beta_risk=surv_dat$lp[X_indices[i]:n]
    A_risk=A_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    Z_beta_risk=Z_beta_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    LT_risk=surv_dat$lt[X_indices[i]:n]
    LT_risk=LT_risk[LT_risk < X[i]]
    s_risk=A_risk+X[i]-LT_risk
    s_i=surv_dat$starts+X[i]-surv_dat$lt
    if(raw) {
      denominator[i]=sum(exp(alpha*inf_dat$Rhat[findInterval(s_risk,inf_dat$s)]+Z_beta_risk))
    } else {
      denominator[i]=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2)+theta3*(s_risk^3))+Z_beta_risk))
    }
  } 
  h=numerator/denominator
  X_unique=unique(X)
  
  return(list(X=X_unique,h=h))
}

h_t_Naive=function(beta,alpha,surv_dat,period=1) {
  
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  surv_dat$lp=as.numeric(as.matrix(surv_dat[,Z_names])%*%beta)
  A=surv_dat$starts[surv_dat$status==1]
  Z_beta=as.numeric(as.matrix(surv_dat[surv_dat$status==1,Z_names])%*%beta)
  X=surv_dat$eventtime[surv_dat$status==1]
  LT=surv_dat$lt[surv_dat$status==1]
  A=A[order(X)]
  Z_beta=Z_beta[order(X)]
  LT=LT[order(X)]
  X=X[order(X)]
  surv_dat=surv_dat[order(surv_dat$eventtime),]
  X_indices=match(unique(surv_dat$eventtime[surv_dat$status==1]),surv_dat$eventtime)
  n_times=sapply(split(surv_dat$eventtime[surv_dat$status==1],surv_dat$eventtime[surv_dat$status==1]),length)
  X_indices=rep(X_indices,n_times)
  n=dim(surv_dat)[1]
  s_i=A+X-LT
  if(period==1) {
    numerator=exp(alpha*as.numeric(s_i >= 1 & s_i <= 5 | s_i >= 17 & s_i <= 20)+Z_beta)
  } else {
    numerator=exp(alpha*as.numeric(s_i >= 1 & s_i <= 13)+Z_beta)
  }
  denominator=vector("numeric")
  for(i in 1:length(X)) {
    A_risk=surv_dat$starts[X_indices[i]:n]
    Z_beta_risk=surv_dat$lp[X_indices[i]:n]
    A_risk=A_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    Z_beta_risk=Z_beta_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    LT_risk=surv_dat$lt[X_indices[i]:n]
    LT_risk=LT_risk[LT_risk < X[i]]
    s_risk=A_risk+X[i]-LT_risk
    if(period==1) {
      denominator[i]=sum(exp(alpha*as.numeric(s_risk >= 1 & s_risk <= 5 | s_risk >= 17 & s_risk <= 20)+Z_beta_risk))
    } else {
      denominator[i]=sum(exp(alpha*as.numeric(s_risk >= 1 & s_risk <= 13)+Z_beta_risk))
    }
  } 
  h=numerator/denominator
  X_unique=unique(X)
  
  return(list(X=X_unique,h=h))
}

D_Helper=function(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,f=NULL,db=F,Z_index=1) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  surv_dat$lp=as.numeric(as.matrix(surv_dat[,Z_names])%*%beta)
  A=surv_dat$starts[surv_dat$status==1]
  Z_beta=as.numeric(as.matrix(surv_dat[surv_dat$status==1,Z_names])%*%beta)
  X=surv_dat$eventtime[surv_dat$status==1]
  LT=surv_dat$lt[surv_dat$status==1]
  A=A[order(X)]
  Z_beta=Z_beta[order(X)]
  LT=LT[order(X)]
  X=X[order(X)]
  surv_dat=surv_dat[order(surv_dat$eventtime),]
  X_indices=match(unique(surv_dat$eventtime[surv_dat$status==1]),surv_dat$eventtime)
  n_times=sapply(split(surv_dat$eventtime[surv_dat$status==1],surv_dat$eventtime[surv_dat$status==1]),length)
  X_indices=rep(X_indices,n_times)
  n=dim(surv_dat)[1]
  out=vector("numeric")
  for(i in 1:length(X)) {
    A_risk=surv_dat$starts[X_indices[i]:n]
    Z_beta_risk=surv_dat$lp[X_indices[i]:n]
    A_risk=A_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    Z_beta_risk=Z_beta_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    LT_risk=surv_dat$lt[X_indices[i]:n]
    LT_risk=LT_risk[LT_risk < X[i]]
    s_risk=A_risk+X[i]-LT_risk
    if(db) {
      numerator=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2)+theta3*(s_risk^3))+Z_beta_risk)*(as.numeric(surv_dat[,Z_names[Z_index]][X_indices[i]:n])[surv_dat$lt[X_indices[i]:n] < X[i]]))
    } else{
      numerator=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2)+theta3*(s_risk^3))+Z_beta_risk)*f(s_risk))
    }
    denominator=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2)+theta3*(s_risk^3))+Z_beta_risk))
    out[i]=numerator/denominator
  }
  return(out)
}

Log_L=function(beta,alpha,theta0,theta1,theta2,theta3,phi,surv_dat,inf_dat) {
  
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2+theta3*inf_dat$s^3)
  l_I=sum(log(dnbinom(inf_dat$I,mu=R*inf_dat$I_prev,size=phi)))
  
  haz_result=h_t(beta,alpha,theta0,theta1,theta2,theta3,surv_dat)
  h_vals=haz_result$h
  
  return(l_I+sum(log(h_vals)))
}

Neg_Log_L=function(par,surv_dat,inf_dat) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  result=-Log_L(par[1:length(Z_names)],par[length(Z_names)+1],par[length(Z_names)+2],par[length(Z_names)+3],par[length(Z_names)+4],par[length(Z_names)+5],par[length(Z_names)+6],surv_dat,inf_dat)
  show(c(par,result))
  return(result)
}

D2_Log_L=function(beta,alpha,theta0,theta1,theta2,theta3,phi,surv_dat,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  
  R=exp(theta0+theta1*s+theta2*s^2+theta3*s^3)
  n=dim(surv_dat)[1]
  dtheta0=c(rep(0,n),I-(I+phi)*R*I_prev/(R*I_prev+phi))
  dtheta1=c(rep(0,n),s*(I-(I+phi)*R*I_prev/(R*I_prev+phi)))
  dtheta2=c(rep(0,n),(s^2)*(I-(I+phi)*R*I_prev/(R*I_prev+phi)))
  dtheta3=c(rep(0,n),(s^3)*(I-(I+phi)*R*I_prev/(R*I_prev+phi)))
  dbeta=matrix(0,nrow=length(Z_names),ncol=(n+length(s)))
  dalpha=rep(0,n+length(s))
  
  fill=which(surv_dat$status==1)
  
  s_i=surv_dat$starts[surv_dat$status==1]+surv_dat$eventtime[surv_dat$status==1]-surv_dat$lt[surv_dat$status==1]
  exp_i=exp(theta0+theta1*s_i+theta2*s_i^2+theta3*s_i^3)
  
  for(i in 1:length(Z_names)) {
    dbeta[i,fill]=dbeta[i,fill]+as.numeric(surv_dat[surv_dat$status==1,Z_names[i]])-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,db=T,Z_index=i)
  }
  dtheta0[fill]=dtheta0[fill]+alpha*exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,f=function(s) {alpha*exp(theta0+theta1*s+theta2*(s^2)+theta3*(s^3))})
  dtheta1[fill]=dtheta1[fill]+alpha*s_i*exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,f=function(s) {alpha*s*exp(theta0+theta1*s+theta2*(s^2)+theta3*(s^3))})
  dtheta2[fill]=dtheta2[fill]+alpha*(s_i^2)*exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,f=function(s) {alpha*(s^2)*exp(theta0+theta1*s+theta2*(s^2)+theta3*(s^3))})
  dtheta3[fill]=dtheta3[fill]+alpha*(s_i^3)*exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,f=function(s) {alpha*(s^3)*exp(theta0+theta1*s+theta2*(s^2)+theta3*(s^3))})
  dalpha[fill]=dalpha[fill]+exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,f=function(s) {exp(theta0+theta1*s+theta2*(s^2)+theta3*(s^3))})
  
  dbeta_beta=matrix(0,nrow=length(Z_names),ncol=length(Z_names))
  dbeta_alpha=vector("numeric")
  dbeta_theta0=vector("numeric")
  dbeta_theta1=vector("numeric")
  dbeta_theta2=vector("numeric")
  dbeta_theta3=vector("numeric")
  
  for(i in 1:length(Z_names)) {
    dbeta_alpha[i]=sum(dbeta[i,]*dalpha,na.rm=T)
    dbeta_theta0[i]=sum(dbeta[i,]*dtheta0,na.rm=T)
    dbeta_theta1[i]=sum(dbeta[i,]*dtheta1,na.rm=T)
    dbeta_theta2[i]=sum(dbeta[i,]*dtheta2,na.rm=T)
    dbeta_theta3[i]=sum(dbeta[i,]*dtheta3,na.rm=T)
    for(j in 1:length(Z_names)) {
      dbeta_beta[i,j]=sum(dbeta[i,]*dbeta[j,],na.rm=T)
    }
  }
  
  dalpha_alpha=sum(dalpha^2,na.rm=T)
  dalpha_theta0=sum(dalpha*dtheta0,na.rm=T)
  dalpha_theta1=sum(dalpha*dtheta1,na.rm=T)
  dalpha_theta2=sum(dalpha*dtheta2,na.rm=T)
  dalpha_theta3=sum(dalpha*dtheta3,na.rm=T)
  dtheta0_theta0=sum(dtheta0^2,na.rm=T)
  dtheta0_theta1=sum(dtheta0*dtheta1,na.rm=T)
  dtheta0_theta2=sum(dtheta0*dtheta2,na.rm=T)
  dtheta0_theta3=sum(dtheta0*dtheta3,na.rm=T)
  dtheta1_theta1=sum(dtheta1^2,na.rm=T)
  dtheta1_theta2=sum(dtheta1*dtheta2,na.rm=T)
  dtheta1_theta3=sum(dtheta1*dtheta3,na.rm=T)
  dtheta2_theta2=sum(dtheta2^2,na.rm=T)
  dtheta2_theta3=sum(dtheta2*dtheta3,na.rm=T)
  dtheta3_theta3=sum(dtheta3^2,na.rm=T)
  
  hess=matrix(NA,nrow=5+length(Z_names),ncol=5+length(Z_names))
  
  hess[1:length(Z_names),1:length(Z_names)]=dbeta_beta
  
  hess[1:length(Z_names),(length(Z_names)+1):(length(Z_names)+5)]=matrix(c(dbeta_alpha,dbeta_theta0,dbeta_theta1,dbeta_theta2,dbeta_theta3),nrow=length(Z_names),ncol=5)
  
  hess[(length(Z_names)+1):(length(Z_names)+5),1:length(Z_names)]=matrix(c(dbeta_alpha,dbeta_theta0,dbeta_theta1,dbeta_theta2,dbeta_theta3),nrow=5,ncol=length(Z_names),byrow=T)
  
  hess[(length(Z_names)+1):(length(Z_names)+5),(length(Z_names)+1):(length(Z_names)+5)]=matrix(c(dalpha_alpha,dalpha_theta0,dalpha_theta1,dalpha_theta2,dalpha_theta3,
                                                                                                 dalpha_theta0,dtheta0_theta0,dtheta0_theta1,dtheta0_theta2,dtheta0_theta3,
                                                                                                 dalpha_theta1,dtheta0_theta1,dtheta1_theta1,dtheta1_theta2,dtheta1_theta3,
                                                                                                 dalpha_theta2,dtheta0_theta2,dtheta1_theta2,dtheta2_theta2,dtheta2_theta3,
                                                                                                 dalpha_theta3,dtheta0_theta3,dtheta1_theta3,dtheta2_theta3,dtheta3_theta3),nrow=5,ncol=5)
  
  return(hess)
}

D2_Inf_Log_L=function(theta0,theta1,theta2,theta3,phi,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  
  R=exp(theta0+theta1*s+theta2*s^2+theta3*s^3)
  n=dim(surv_dat)[1]
  dtheta0=c(rep(0,n),I-(I+phi)*R*I_prev/(R*I_prev+phi))
  dtheta1=c(rep(0,n),s*(I-(I+phi)*R*I_prev/(R*I_prev+phi)))
  dtheta2=c(rep(0,n),(s^2)*(I-(I+phi)*R*I_prev/(R*I_prev+phi)))
  dtheta3=c(rep(0,n),(s^3)*(I-(I+phi)*R*I_prev/(R*I_prev+phi)))
  
  dtheta0_theta0=sum(dtheta0^2,na.rm=T)
  dtheta0_theta1=sum(dtheta0*dtheta1,na.rm=T)
  dtheta0_theta2=sum(dtheta0*dtheta2,na.rm=T)
  dtheta0_theta3=sum(dtheta0*dtheta3,na.rm=T)
  dtheta1_theta1=sum(dtheta1^2,na.rm=T)
  dtheta1_theta2=sum(dtheta1*dtheta2,na.rm=T)
  dtheta1_theta3=sum(dtheta1*dtheta3,na.rm=T)
  dtheta2_theta2=sum(dtheta2^2,na.rm=T)
  dtheta2_theta3=sum(dtheta2*dtheta3,na.rm=T)
  dtheta3_theta3=sum(dtheta3^2,na.rm=T)
  hess_theta=matrix(c(dtheta0_theta0,dtheta0_theta1,dtheta0_theta2,dtheta0_theta3,
                      dtheta0_theta1,dtheta1_theta1,dtheta1_theta2,dtheta1_theta3,
                      dtheta0_theta2,dtheta1_theta2,dtheta2_theta2,dtheta2_theta3,
                      dtheta0_theta3,dtheta1_theta3,dtheta2_theta3,dtheta3_theta3),nrow=4)
  
  return(hess_theta)
}

Surv_Log_L=function(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,raw=F,inf_dat=NULL) {
  haz_result=h_t(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,raw,inf_dat)
  h_vals=haz_result$h
  return(sum(log(h_vals)))
}

Surv_Log_L_Naive=function(beta,alpha,surv_dat,period=1) {
  haz_result=h_t_Naive(beta,alpha,surv_dat,period=period)
  h_vals=haz_result$h
  return(sum(log(h_vals)))
}

Inf_Log_L=function(theta0,theta1,theta2,theta3,phi,inf_dat,raw=F) {
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2+theta3*inf_dat$s^3)
  l_I=sum(log(dnbinom(inf_dat$I,mu=R*inf_dat$I_prev,size=phi)))
  return(l_I)
}

Test_Log_L=function(beta,alpha,theta0,theta1,theta2,theta3,surv_dat) {
  
  haz_result=h_t(beta,alpha,theta0,theta1,theta2,theta3,surv_dat)
  h_vals=haz_result$h
  
  return(sum(log(h_vals)))
}

Naive_Neg_Log_L=function(par,surv_dat,period=1,naive=F) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  if(!naive) {
    result=-Surv_Log_L_Naive(par[1:length(Z_names)],par[length(Z_names)+1],surv_dat,period)
  } else {
    result=-Surv_Log_L_Naive(par[1:length(Z_names)],0,surv_dat,period)
    
  }
  return(result)
}

Test_Neg_Log_L=function(par,theta0,theta1,theta2,theta3,surv_dat,raw=F,inf_dat=NULL) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  result=-Surv_Log_L(par[1:length(Z_names)],par[length(Z_names)+1],theta0,theta1,theta2,theta3,surv_dat,raw=raw,inf_dat=inf_dat)
  show(c(par,result))
  return(result)
}

Test_D2_Log_L=function(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  
  R=exp(theta0+theta1*s+theta2*s^2+theta3*s^3)
  n=dim(surv_dat)[1]
  dbeta=matrix(0,nrow=length(Z_names),ncol=n)
  dalpha=rep(0,n)
  
  surv_dat=surv_dat[order(surv_dat$eventtime),]
  fill=which(surv_dat$status==1)
  
  s_i=surv_dat$starts[surv_dat$status==1]+surv_dat$eventtime[surv_dat$status==1]-surv_dat$lt[surv_dat$status==1]
  exp_i=exp(theta0+theta1*s_i+theta2*s_i^2+theta3*s_i^3)
  
  for(i in 1:length(Z_names)) {
    dbeta[i,fill]=dbeta[i,fill]+as.numeric(surv_dat[surv_dat$status==1,Z_names[i]])-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,db=T,Z_index=i)
  }
  dalpha[fill]=dalpha[fill]+exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,f=function(s) {exp(theta0+theta1*s+theta2*(s^2)+theta3*(s^3))})
  
  dbeta_beta=matrix(0,nrow=length(Z_names),ncol=length(Z_names))
  dbeta_alpha=vector("numeric")
  
  for(i in 1:length(Z_names)) {
    dbeta_alpha[i]=sum(dbeta[i,]*dalpha,na.rm=T)
    for(j in 1:length(Z_names)) {
      dbeta_beta[i,j]=sum(dbeta[i,]*dbeta[j,],na.rm=T)
    }
  }
  
  dalpha_alpha=sum(dalpha^2,na.rm=T)
  
  hess=matrix(NA,nrow=3,ncol=3)
  
  hess[1:length(Z_names),1:length(Z_names)]=dbeta_beta
  
  hess[1:length(Z_names),length(Z_names)+1]=dbeta_alpha
  
  hess[(length(Z_names)+1),1:length(Z_names)]=dbeta_alpha
  
  hess[(length(Z_names)+1),(length(Z_names)+1)]=dalpha_alpha
  
  return(hess)
}











