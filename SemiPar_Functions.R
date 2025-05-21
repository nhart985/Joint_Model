
h_t=function(beta,alpha,theta0,theta1,theta2,surv_dat) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  surv_dat$lp=as.numeric(as.matrix(surv_dat[,Z_names])%*%beta)
  A=surv_dat$starts[surv_dat$status==1]
  Z_beta=as.numeric(as.matrix(surv_dat[surv_dat$status==1,Z_names])%*%beta)
  X=surv_dat$eventtime[surv_dat$status==1]
  A=A[order(X)]
  Z_beta=Z_beta[order(X)]
  X=X[order(X)]
  surv_dat=surv_dat[order(surv_dat$eventtime),]
  X_indices=which(surv_dat$status==1)
  n=dim(surv_dat)[1]
  s_i=A+X
  numerator=exp(alpha*exp(theta0+theta1*s_i+theta2*(s_i^2))+Z_beta)
  denominator=vector("numeric")
  for(i in 1:length(X)) {
    A_risk=surv_dat$starts[X_indices[i]:n]
    Z_beta_risk=surv_dat$lp[X_indices[i]:n]
    A_risk=A_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    Z_beta_risk=Z_beta_risk[surv_dat$lt[X_indices[i]:n] < X[i]]
    s_risk=A_risk+X[i]
    denominator[i]=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk))
  }
  h=numerator/denominator
  return(list(X=X,h=h,H=cumsum(h*c(X[1],diff(X)))))
}

D_Helper=function(beta,alpha,theta0,theta1,theta2,surv_dat,f=NULL,db=F,Z_index=1) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  surv_dat$lp=as.numeric(as.matrix(surv_dat[,Z_names])%*%beta)
  A=surv_dat$starts[surv_dat$status==1]
  Z_beta=as.numeric(as.matrix(surv_dat[surv_dat$status==1,Z_names])%*%beta)
  X=surv_dat$eventtime[surv_dat$status==1]
  A=A[order(X)]
  Z_beta=Z_beta[order(X)]
  X=X[order(X)]
  surv_dat=surv_dat[order(surv_dat$eventtime),]
  X_indices=which(surv_dat$status==1)
  n=dim(surv_dat)[1]
  out=vector("numeric")
  for(i in 1:length(X)) {
    A_risk=surv_dat$starts[X_indices[i]:n]
    Z_beta_risk=surv_dat$lp[X_indices[i]:n]
    A_risk=A_risk[surv_dat$lt[X_indices[i]:n] <= X[i]]
    Z_beta_risk=Z_beta_risk[surv_dat$lt[X_indices[i]:n] <= X[i]]
    s_risk=A_risk+X[i]
    if(db) {
      numerator=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk)*(as.numeric(surv_dat[,Z_names[Z_index]][X_indices[i]:n])[surv_dat$lt[X_indices[i]:n] <= X[i]]))
    } else{
      numerator=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk)*f(s_risk))
    }
    denominator=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk))
    out[i]=numerator/denominator
  }
  return(out)
}

D_Helper2=function(beta,alpha,theta0,theta1,theta2,surv_dat,f=NULL,db=F,Z_index=1) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  surv_dat$lp=as.numeric(as.matrix(surv_dat[,Z_names])%*%beta)
  A=surv_dat$starts[surv_dat$status==1]
  Z_beta=as.numeric(as.matrix(surv_dat[surv_dat$status==1,Z_names])%*%beta)
  X=surv_dat$eventtime[surv_dat$status==1]
  A=A[order(X)]
  Z_beta=Z_beta[order(X)]
  X=X[order(X)]
  surv_dat=surv_dat[order(surv_dat$eventtime),]
  X_indices=which(surv_dat$status==1)
  n=dim(surv_dat)[1]
  out=vector("numeric")
  for(i in 1:length(X)) {
    s_i=A[i]+X[i]
    A_risk=surv_dat$starts[X_indices[i]:n]
    Z_beta_risk=surv_dat$lp[X_indices[i]:n]
    A_risk=A_risk[surv_dat$lt[X_indices[i]:n] <= X[i]]
    Z_beta_risk=Z_beta_risk[surv_dat$lt[X_indices[i]:n] <= X[i]]
    s_risk=A_risk+X[i]
    if(db) {
      numerator1=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk))*exp(alpha*exp(theta0+theta1*s_i+theta2*s_i^2)+Z_beta[i])*as.numeric(surv_dat[i,Z_names[Z_index]])
      numerator2=exp(alpha*exp(theta0+theta1*s_i+theta2*s_i^2)+Z_beta[i])*sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk)*(as.numeric(surv_dat[X_indices[i]:n,Z_names[Z_index]]))[surv_dat$lt[X_indices[i]:n] <= X[i]])
    } else{
      numerator1=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk))*exp(alpha*exp(theta0+theta1*s_i+theta2*s_i^2)+Z_beta[i])*f(s_i)
      numerator2=exp(alpha*exp(theta0+theta1*s_i+theta2*s_i^2)+Z_beta[i])*sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk)*f(s_risk))
    }
    denominator=sum(exp(alpha*exp(theta0+theta1*s_risk+theta2*(s_risk^2))+Z_beta_risk))^2
    out[i]=(numerator1-numerator2)/denominator
  }
  out_step=stepfun(X,c(0,cumsum(out*c(X[1],diff(X)))))
  return(out_step(surv_dat$eventtime))
}

Log_L_Inf=function(theta0,theta1,theta2,inf_dat) {
  
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2)
  l_I=sum(log(dpois(inf_dat$I,R*inf_dat$I_prev)))
  
  return(l_I)
}

Neg_Log_L_Inf=function(par,N,inf_dat) {
  return(-Log_L_Inf(par[1],par[2],par[3],inf_dat))
}

Log_L=function(beta,alpha,theta0,theta1,theta2,surv_dat,inf_dat) {
  
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2)
  l_I=sum(log(dpois(inf_dat$I,R*inf_dat$I_prev)))
  
  haz_result=h_t(beta,alpha,theta0,theta1,theta2,surv_dat)
  h_vals=haz_result$h
  H_fun=stepfun(haz_result$X,c(0,haz_result$H))
  S_vals=exp(-H_fun(surv_dat$eventtime))
  
  show(c(beta,alpha,theta0,theta1,theta2,l_I+sum(log(h_vals))+sum(log(S_vals))))
  
  return(l_I+sum(log(h_vals))+sum(log(S_vals)))
}

Neg_Log_L=function(par,surv_dat,inf_dat) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  return(-Log_L(par[1:length(Z_names)],par[length(Z_names)+1],par[length(Z_names)+2],par[length(Z_names)+3],par[length(Z_names)+4],surv_dat,inf_dat))
}

D_Log_L=function(beta,alpha,theta0,theta1,theta2,surv_dat,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  
  dtheta0=sum(I-I_prev*exp(theta0+theta1*s+theta2*s^2))
  dtheta1=sum(s*(I-I_prev*exp(theta0+theta1*s+theta2*s^2)))
  dtheta2=sum((s^2)*(I-I_prev*exp(theta0+theta1*s+theta2*s^2)))
  dbeta=rep(0,length(Z_names))
  dalpha=0
  
  s_i=surv_dat$starts[surv_dat$status==1]+surv_dat$eventtime[surv_dat$status==1]
  exp_i=exp(theta0+theta1*s_i+theta2*s_i^2)
  
  for(i in 1:length(Z_names)) {
    dbeta[i]=dbeta[i]+sum(as.numeric(surv_dat[surv_dat$status==1,Z_names[i]]))-sum(D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,db=T,Z_index=i))
  }
  dtheta0=dtheta0+alpha*sum(exp_i)-sum(D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*exp(theta0+theta1*s+theta2*(s^2))}))
  dtheta1=dtheta1+alpha*sum(s_i*exp_i)-sum(D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*s*exp(theta0+theta1*s+theta2*(s^2))}))
  dtheta2=dtheta2+alpha*sum((s_i^2)*exp_i)-sum(D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*(s^2)*exp(theta0+theta1*s+theta2*(s^2))}))
  dalpha=dalpha+sum(exp_i)-sum(D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {exp(theta0+theta1*s+theta2*(s^2))}))
  
  for(i in 1:length(Z_names)) {
    dbeta[i]=dbeta[i]+sum(D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,db=T,Z_index=i))
  }
  dtheta0=dtheta0+sum(D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*exp(theta0+theta1*s+theta2*(s^2))}))
  dtheta1=dtheta1+sum(D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*s*exp(theta0+theta1*s+theta2*(s^2))}))
  dtheta2=dtheta2+sum(D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*(s^2)*exp(theta0+theta1*s+theta2*(s^2))}))
  dalpha=dalpha+sum(D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {exp(theta0+theta1*s+theta2*(s^2))}))
  
  return(c(dbeta,dalpha,dtheta0,dtheta1,dtheta2))
}

Neg_D_Log_L=function(par,surv_dat,inf_dat) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  return(-D_Log_L(par[1:length(Z_names)],par[length(Z_names)+1],par[length(Z_names)+2],par[length(Z_names)+3],par[length(Z_names)+4],surv_dat,inf_dat))
}

D2_Log_L_Inf=function(theta0,theta1,theta2,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  
  dtheta0_theta0=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2))
  dtheta0_theta1=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s)
  dtheta0_theta2=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^2)
  dtheta1_theta1=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^2)
  dtheta1_theta2=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^3)
  dtheta2_theta2=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^4)
  
  hess=matrix(c(dtheta0_theta0,dtheta0_theta1,dtheta0_theta2,
                dtheta0_theta1,dtheta1_theta1,dtheta1_theta2,
                dtheta0_theta2,dtheta1_theta2,dtheta2_theta2),nrow=3,ncol=3)
  
  return(hess)
}

D2_Log_L=function(beta,alpha,theta0,theta1,theta2,surv_dat,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  
  n=dim(surv_dat)[1]
  dtheta0=c(rep(0,n),I-I_prev*exp(theta0+theta1*s+theta2*s^2))
  dtheta1=c(rep(0,n),s*(I-I_prev*exp(theta0+theta1*s+theta2*s^2)))
  dtheta2=c(rep(0,n),(s^2)*(I-I_prev*exp(theta0+theta1*s+theta2*s^2)))
  dbeta=matrix(0,nrow=length(Z_names),ncol=(n+length(s)))
  dalpha=rep(0,n+length(s))
  
  fill=which(surv_dat$status==1)
  
  s_i=surv_dat$starts[surv_dat$status==1]+surv_dat$eventtime[surv_dat$status==1]
  exp_i=exp(theta0+theta1*s_i+theta2*s_i^2)
  
  for(i in 1:length(Z_names)) {
    dbeta[i,fill]=dbeta[i,fill]+as.numeric(surv_dat[surv_dat$status==1,Z_names[i]])-D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,db=T,Z_index=i)
  }
  dtheta0[fill]=dtheta0[fill]+alpha*exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*exp(theta0+theta1*s+theta2*(s^2))})
  dtheta1[fill]=dtheta1[fill]+alpha*s_i*exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*s*exp(theta0+theta1*s+theta2*(s^2))})
  dtheta2[fill]=dtheta2[fill]+alpha*(s_i^2)*exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*(s^2)*exp(theta0+theta1*s+theta2*(s^2))})
  dalpha[fill]=dalpha[fill]+exp_i-D_Helper(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {exp(theta0+theta1*s+theta2*(s^2))})
  
  for(i in 1:length(Z_names)) {
    dbeta[i,1:n]=dbeta[i,1:n]+D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,db=T,Z_index=i)
  }
  dtheta0[1:n]=dtheta0[1:n]+D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*exp(theta0+theta1*s+theta2*(s^2))})
  dtheta1[1:n]=dtheta1[1:n]+D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*s*exp(theta0+theta1*s+theta2*(s^2))})
  dtheta2[1:n]=dtheta2[1:n]+D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {alpha*(s^2)*exp(theta0+theta1*s+theta2*(s^2))})
  dalpha[1:n]=dalpha[1:n]+D_Helper2(beta,alpha,theta0,theta1,theta2,surv_dat,f=function(s) {exp(theta0+theta1*s+theta2*(s^2))})
  
  dbeta_beta=matrix(0,nrow=length(Z_names),ncol=length(Z_names))
  dbeta_alpha=vector("numeric")
  dbeta_theta0=vector("numeric")
  dbeta_theta1=vector("numeric")
  dbeta_theta2=vector("numeric")
  
  for(i in 1:length(Z_names)) {
    dbeta_alpha[i]=sum(dbeta[i,]*dalpha)
    dbeta_theta0[i]=sum(dbeta[i,]*dtheta0)
    dbeta_theta1[i]=sum(dbeta[i,]*dtheta1)
    dbeta_theta2[i]=sum(dbeta[i,]*dtheta2)
    for(j in 1:length(Z_names)) {
      dbeta_beta[i,j]=sum(dbeta[i,]*dbeta[j,])
    }
  }
  
  dalpha_alpha=sum(dalpha^2)
  dalpha_theta0=sum(dalpha*dtheta0)
  dalpha_theta1=sum(dalpha*dtheta1)
  dalpha_theta2=sum(dalpha*dtheta2)
  dtheta0_theta0=sum(dtheta0^2)
  dtheta0_theta1=sum(dtheta0*dtheta1)
  dtheta0_theta2=sum(dtheta0*dtheta2)
  dtheta1_theta1=sum(dtheta1^2)
  dtheta1_theta2=sum(dtheta1*dtheta2)
  dtheta2_theta2=sum(dtheta2^2)
  
  hess=matrix(NA,nrow=4+length(Z_names),ncol=4+length(Z_names))
  
  hess[1:length(Z_names),1:length(Z_names)]=dbeta_beta
  
  hess[1:length(Z_names),(length(Z_names)+1):(length(Z_names)+4)]=matrix(c(dbeta_alpha,dbeta_theta0,dbeta_theta1,dbeta_theta2),nrow=length(Z_names),ncol=4)
  
  hess[(length(Z_names)+1):(length(Z_names)+4),1:length(Z_names)]=matrix(c(dbeta_alpha,dbeta_theta0,dbeta_theta1,dbeta_theta2),nrow=4,ncol=length(Z_names),byrow=T)
  
  hess[(length(Z_names)+1):(length(Z_names)+4),(length(Z_names)+1):(length(Z_names)+4)]=matrix(c(dalpha_alpha,dalpha_theta0,dalpha_theta1,dalpha_theta2,
                                                                         dalpha_theta0,dtheta0_theta0,dtheta0_theta1,dtheta0_theta2,
                                                                         dalpha_theta1,dtheta0_theta1,dtheta1_theta1,dtheta1_theta2,
                                                                         dalpha_theta2,dtheta0_theta2,dtheta1_theta2,dtheta2_theta2))
  
  return(hess)
}
