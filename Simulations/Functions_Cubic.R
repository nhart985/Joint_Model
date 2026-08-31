
h_t=function(q,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3) {
  s=start+q
  return(lambda*exp(as.numeric(beta%*%Z)+alpha*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)))
}

h_t_simple=function(q,start,Z,beta,lambda,alpha,method,inf_dat) {
  s=start+q
  if(method=="naive") {
    out=rep(lambda*exp(as.numeric(beta%*%Z)),length(s))
  } else if(method=="period 1") {
    out=lambda*exp(as.numeric(beta%*%Z)+alpha*as.numeric(s >= 15 & s <= 25))
  } else if (method=="period 2") {
    out=lambda*exp(as.numeric(beta%*%Z)+alpha*as.numeric(s >= 10 & s <= 30))
  } else if(method=="raw") {
    out=lambda*exp(as.numeric(beta%*%Z)+alpha*inf_dat$Rhat[findInterval(s,inf_dat$s)])
  }
  return(out)
}

S_t=function(t,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3) {
  H_t=integrate(f=h_t,lower=0,upper=t,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3)$value
  return(exp(-H_t))
}

S_t_simple=function(t,start,Z,beta,lambda,alpha,method,inf_dat) {
  H_t=integrate(f=h_t_simple,lower=0,upper=t,start,Z,beta,lambda,alpha,method,inf_dat)$value
  return(exp(-H_t))
}

D_Helper=function(t,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun) {
  fun=function(q,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun) {
    return(s_fun(start+q)*h_t(q,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3))
  }
  return(-integrate(f=fun,lower=0,upper=t,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun)$value)
}

Log_L_Inf=function(theta0,theta1,theta2,theta3,inf_dat) {
  
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2+theta3*inf_dat$s^3)
  l_I=sum(log(dpois(inf_dat$I,R*inf_dat$I_prev)))
  
  return(l_I)
}

Neg_Log_L_Inf=function(par,inf_dat) {
  return(-Log_L_Inf(par[1],par[2],par[3],par[4],inf_dat))
}

Log_L=function(beta,lambda,alpha,theta0,theta1,theta2,theta3,surv_dat,inf_dat) {
  
  R=exp(theta0+theta1*inf_dat$s+theta2*inf_dat$s^2+theta3*inf_dat$s^3)
  l_I=sum(log(dpois(inf_dat$I,R*inf_dat$I_prev)))
  
  result=vector("numeric")
  for(i in 1:length(surv_dat$id)) {
    X=surv_dat$eventtime[surv_dat$id==i]
    Delta=surv_dat$status[surv_dat$id==i]
    start=surv_dat$starts[surv_dat$id==i]
    Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
    Z=as.numeric(as.matrix(surv_dat[,Z_names])[surv_dat$id==i,])
    
    result[i]=log(h_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3)^Delta)+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3))
  }
  return(l_I+sum(result))
}

Neg_Log_L=function(par,N,surv_dat,inf_dat) {
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  return(-Log_L(par[1:length(Z_names)],par[length(Z_names)+1],par[length(Z_names)+2],par[length(Z_names)+3],par[length(Z_names)+4],par[length(Z_names)+5],par[length(Z_names)+6],surv_dat,inf_dat))
}

D_Log_L=function(beta,lambda,alpha,theta0,theta1,theta2,theta3,surv_dat,inf_dat) {
  
  s=inf_dat$s
  I=inf_dat$I
  I_prev=inf_dat$I_prev
  Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
  
  dtheta0=sum(I-I_prev*exp(theta0+theta1*s+theta2*s^2+theta3*s^3))
  dtheta1=sum(s*(I-I_prev*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)))
  dtheta2=sum((s^2)*(I-I_prev*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)))
  dtheta3=sum((s^3)*(I-I_prev*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)))
  dlambda=0
  dbeta=rep(0,length(Z_names))
  dalpha=0
  
  for(i in 1:length(surv_dat$id)) {
    X=surv_dat$eventtime[surv_dat$id==i]
    Delta=surv_dat$status[surv_dat$id==i]
    start=surv_dat$starts[surv_dat$id==i]
    Z=as.numeric(surv_dat[,Z_names][surv_dat$id==i,])
    
    dbeta=dbeta+Delta*Z+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3))*Z
    dtheta0=dtheta0+Delta*alpha*exp(theta0+theta1*(start+X)+theta2*(start+X)^2+theta3*(start+X)^3)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)})
    dtheta1=dtheta1+Delta*alpha*exp(theta0+theta1*(start+X)+theta2*(start+X)^2+theta3*(start+X)^3)*(start+X)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)*s})
    dtheta2=dtheta2+Delta*alpha*exp(theta0+theta1*(start+X)+theta2*(start+X)^2+theta3*(start+X)^3)*(start+X)^2+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)*s^2})
    dtheta3=dtheta2+Delta*alpha*exp(theta0+theta1*(start+X)+theta2*(start+X)^2+theta3*(start+X)^3)*(start+X)^3+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun=function(s){alpha*exp(theta0+theta1*s+theta2*s^2+theta3*s^3)*s^3})
    dlambda=dlambda+Delta/lambda+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3))/lambda
    dalpha=dalpha+Delta*exp(theta0+theta1*(start+X)+theta2*(start+X)^2+theta3*(start+X)^3)+D_Helper(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3,s_fun=function(s){exp(theta0+theta1*s+theta2*s^2+theta3*s^3)})
  }
  return(c(dbeta,dlambda,dalpha,dtheta0,dtheta1,dtheta2,dtheta3))
}

Log_L_Surv=function(beta,lambda,alpha,theta0,theta1,theta2,theta3,surv_dat) {
  
  result=vector("numeric")
  for(i in 1:length(surv_dat$id)) {
    X=surv_dat$eventtime[surv_dat$id==i]
    Delta=surv_dat$status[surv_dat$id==i]
    start=surv_dat$starts[surv_dat$id==i]
    Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
    Z=as.numeric(as.matrix(surv_dat[,Z_names])[surv_dat$id==i,])
    
    result[i]=log(h_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3)^Delta)+log(S_t(X,start,Z,beta,lambda,alpha,theta0,theta1,theta2,theta3))
  }
  return(sum(result))
}

Log_L_Surv_Simple=function(beta,lambda,alpha,surv_dat,inf_dat,method="naive") {
  
  result=vector("numeric")
  for(i in 1:length(surv_dat$id)) {
    X=surv_dat$eventtime[surv_dat$id==i]
    Delta=surv_dat$status[surv_dat$id==i]
    start=surv_dat$starts[surv_dat$id==i]
    Z_names=names(surv_dat)[substr(names(surv_dat),1,1)=="Z"]
    Z=as.numeric(as.matrix(surv_dat[,Z_names])[surv_dat$id==i,])
    
    result[i]=log(h_t_simple(X,start,Z,beta,lambda,alpha,method,inf_dat)^Delta)+log(S_t_simple(X,start,Z,beta,lambda,alpha,method,inf_dat))
  }
  return(sum(result))
}


