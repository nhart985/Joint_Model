h_t=function(q,start,Z,beta,lambda,alpha,theta0,theta1,theta2) {
  s=start+q
  return(lambda*exp(as.numeric(beta%*%Z)+alpha*exp(theta0+theta1*s+theta2*s^2)))
}

S_t=function(t,start,Z,beta,lambda,alpha,theta0,theta1,theta2) {
  H_t=integrate(f=h_t,lower=0,upper=t,start,Z,beta,lambda,alpha,theta0,theta1,theta2)$value
  return(exp(-H_t))
}

D_Helper=function(t,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun) {
  fun=function(q,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun) {
    return(s_fun(start+q)*h_t(q,start,Z,beta,lambda,alpha,theta0,theta1,theta2))
  }
  return(-integrate(f=fun,lower=0,upper=t,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun)$value)
}

Log_L_Inf=function(theta0,theta1,theta2,inf_dat) {
  
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2)
  l_I=sum(log(dpois(inf_dat$I,R*inf_dat$I_prev)))
  
  return(l_I)
}

Neg_Log_L_Inf=function(par,N,inf_dat) {
  return(-Log_L_Inf(par[1],par[2],par[3],inf_dat))
}

Log_L=function(beta,lambda,alpha,theta0,theta1,theta2,surv_dat,inf_dat) {
  
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2)
  l_I=sum(log(dpois(inf_dat$I,R*inf_dat$I_prev)))
  
  result=vector("numeric")
  for(i in 1:length(surv_dat$id)) {
    X=surv_dat$eventtime[surv_dat$id==i]
    Delta=surv_dat$status[surv_dat$id==i]
    start=surv_dat$starts[surv_dat$id==i]
    Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
    Z=as.numeric(surv_dat[,Z_names][surv_dat$id==i])
    
    result[i]=log(h_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2)^Delta)+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2))
  }
  return(l_I+sum(result))
}

Neg_Log_L=function(par,N,surv_dat,inf_dat) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  return(-Log_L(par[1:length(Z_names)],par[length(Z_names)+1],par[length(Z_names)+2],par[length(Z_names)+3],par[length(Z_names)+4],par[length(Z_names)+5],surv_dat,inf_dat))
}

D_Log_L=function(beta,lambda,alpha,theta0,theta1,theta2,surv_dat,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  
  dtheta0=sum(I-I_prev*exp(theta0+theta1*s+theta2*s^2))
  dtheta1=sum(s*(I-I_prev*exp(theta0+theta1*s+theta2*s^2)))
  dtheta2=sum((s^2)*(I-I_prev*exp(theta0+theta1*s+theta2*s^2)))
  dlambda=0
  dbeta=rep(0,length(Z_names))
  dalpha=0
  
  for(i in 1:length(surv_dat$id)) {
    X=surv_dat$eventtime[surv_dat$id==i]
    Delta=surv_dat$status[surv_dat$id==i]
    start=surv_dat$starts[surv_dat$id==i]
    Z=as.numeric(surv_dat[,Z_names][surv_dat$id==i])
    
    dbeta=dbeta+Delta*Z+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2))*Z
    dtheta0=dtheta0+Delta*alpha*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)})
    dtheta1=dtheta1+Delta*alpha*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)*(start+X)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*s})
    dtheta2=dtheta2+Delta*alpha*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)*(start+X)^2+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*s^2})
    dlambda=dlambda+Delta/lambda+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2))/lambda
    dalpha=dalpha+Delta*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2)})
  }
  return(c(dbeta,dlambda,dalpha,dtheta0,dtheta1,dtheta2))
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

D2_Log_L=function(beta,lambda,alpha,theta0,theta1,theta2,surv_dat,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]

  dbeta_dbeta=matrix(0,length(Z_names),length(Z_names))
  dbeta_lambda=rep(0,length(Z_names))
  dbeta_alpha=rep(0,length(Z_names))
  dbeta_theta0=rep(0,length(Z_names))
  dbeta_theta1=rep(0,length(Z_names))
  dbeta_theta2=rep(0,length(Z_names))
  dlambda_lambda=0
  dlambda_alpha=0
  dlambda_theta0=0
  dlambda_theta1=0
  dlambda_theta2=0
  dalpha_alpha=0
  dalpha_theta0=0
  dalpha_theta1=0
  dalpha_theta2=0
  dtheta0_theta0=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2))
  dtheta0_theta1=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s)
  dtheta0_theta2=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^2)
  dtheta1_theta1=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^2)
  dtheta1_theta2=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^3)
  dtheta2_theta2=sum(-I_prev*exp(theta0+theta1*s+theta2*s^2)*s^4) 
  
  for(i in 1:length(surv_dat$id)) {
    X=surv_dat$eventtime[surv_dat$id==i]
    Delta=surv_dat$status[surv_dat$id==i]
    start=surv_dat$starts[surv_dat$id==i]
    Z=as.numeric(surv_dat[,Z_names][surv_dat$id==i])
    
    dbeta_dbeta=dbeta_dbeta+(Z%o%Z)*log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2))
    dbeta_lambda=dbeta_lambda+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2))*Z/lambda
    dbeta_alpha=dbeta_alpha+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2)})*Z
    dbeta_theta0=dbeta_theta0+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)})*Z
    dbeta_theta1=dbeta_theta1+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*s*exp(theta0+theta1*s+theta2*s^2)})*Z
    dbeta_theta2=dbeta_theta2+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*(s^2)*exp(theta0+theta1*s+theta2*s^2)})*Z
    dlambda_lambda=dlambda_lambda-Delta/(lambda^2)
    dlambda_alpha=dlambda_alpha+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2)})/lambda
    dlambda_theta0=dlambda_theta0+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)})/lambda
    dlambda_theta1=dlambda_theta1+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*s*exp(theta0+theta1*s+theta2*s^2)})/lambda
    dlambda_theta2=dlambda_theta2+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*(s^2)*exp(theta0+theta1*s+theta2*s^2)})/lambda
    dalpha_alpha=dalpha_alpha+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2)^2})
    dalpha_theta0=dalpha_theta0+Delta*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2)*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dalpha_theta1=dalpha_theta1+Delta*(start+X)*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2)*s*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dalpha_theta2=dalpha_theta2+Delta*((start+X)^2)*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2)*(s^2)*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dtheta0_theta0=dtheta0_theta0+alpha*Delta*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dtheta0_theta1=dtheta0_theta1+alpha*Delta*(start+X)*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*s*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dtheta0_theta2=dtheta0_theta2+alpha*Delta*((start+X)^2)*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*(s^2)*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dtheta1_theta1=dtheta1_theta1+alpha*Delta*((start+X)^2)*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*(s^2)*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dtheta1_theta2=dtheta1_theta2+alpha*Delta*((start+X)^3)*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*(s^3)*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
    dtheta2_theta2=dtheta2_theta2+alpha*Delta*((start+X)^4)*exp(theta0+theta1*(start+X)+theta2*(start+X)^2)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2)*(s^4)*(1+alpha*exp(theta0+theta1*s+theta2*s^2))})
  }
  
  hess=matrix(NA,nrow=5+length(Z),ncol=5+length(Z))
  
  hess[1:length(Z),1:length(Z)]=dbeta_dbeta
  
  hess[1:length(Z),(length(Z)+1):(length(Z)+5)]=matrix(c(dbeta_lambda,dbeta_alpha,dbeta_theta0,dbeta_theta1,dbeta_theta2),nrow=length(Z_names),ncol=5)
  
  hess[(length(Z)+1):(length(Z)+5),1:length(Z)]=matrix(c(dbeta_lambda,dbeta_alpha,dbeta_theta0,dbeta_theta1,dbeta_theta2),nrow=5,ncol=length(Z_names),byrow=T)
  
  hess[(length(Z)+1):(length(Z)+5),(length(Z)+1):(length(Z)+5)]=matrix(c(dlambda_lambda,dlambda_alpha,dlambda_theta0,dlambda_theta1,dlambda_theta2,
                                                                         dlambda_alpha,dalpha_alpha,dalpha_theta0,dalpha_theta1,dalpha_theta2,
                                                                         dlambda_theta0,dalpha_theta0,dtheta0_theta0,dtheta0_theta1,dtheta0_theta2,
                                                                         dlambda_theta1,dalpha_theta1,dtheta0_theta1,dtheta1_theta1,dtheta1_theta2,
                                                                         dlambda_theta2,dalpha_theta2,dtheta0_theta2,dtheta1_theta2,dtheta2_theta2))
  return(hess)
}

Boot=function(dat,theta,sigma_theta,nboot=100) {
  boots=matrix(0,nboot,3)
  theta_boot=mvrnorm(nboot,theta,sigma_theta)
  for(i in 1:nboot) {
    dat$Rhat_stage1=exp(theta_boot[i,1]+theta_boot[i,2]*dat$s+theta_boot[i,3]*dat$s^2)
    initial=phreg(Surv(tstart,tstop,tdelta)~Z+Rhat_stage1,data=dat,dist="weibull",shape=1)
    beta_stage1=initial$coefficients["Z"]
    alpha_stage1=initial$coefficients["Rhat_stage1"]
    lambda_stage1=1/exp(initial$coefficients["log(scale)"])
    boots[i,]=c(beta_stage1,lambda_stage1,alpha_stage1)
  }
  return(apply(boots,2,sd))
}
