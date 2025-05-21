library(survival)
library(ggplot2)
library(ggpubr)
library(eha)
library(simsurv)
library(rootSolve)

n=1000
I0=10
inf_maxt=52
surv_maxt=12
theta0=-0.06
theta1=0.04
theta2=-0.001
alpha=-1
lambda=0.1
beta=0.05

Gen_Pois=function(inf_maxt,I0,theta0,theta1,theta2) {
  s=1:inf_maxt
  R=exp(theta0+theta1*s+theta2*(s^2))
  I=vector("numeric")
  I[1]=rpois(1,R[1]*I0)
  for(i in 2:inf_maxt) {
    I[i]=rpois(1,R[i]*I[i-1])
  }
  result=data.frame(s,I,R)
  return(result)
}

inf_dat=Gen_Pois(inf_maxt,I0,theta0,theta1,theta2)
inf_dat$I_prev=c(I0,inf_dat$I[1:(inf_maxt-1)])
inf_dat$Rhat=inf_dat$I/inf_dat$I_prev

sum(inf_dat$I==0)

starts=sample(2:(inf_maxt-surv_maxt),n,replace=T)
R=exp(theta0+theta1*inf_dat$s+theta2*(inf_dat$s^2))
Z=rnorm(n,R[starts],0.25)
Xmat=data.frame(id=1:n,starts=starts,Z=Z)
hazard=function(t,x,betas,...) {
  s=x$starts+t
  return(lambda*exp(betas[1]*x$Z+betas[2]*exp(theta0+theta1*s+theta2*s^2)))
}
betas=c(beta,alpha)
names(betas)="alpha"
surv_dat=simsurv(hazard=hazard,x=Xmat,betas=betas,lambda=lambda,inf_dat=inf_dat,maxt=surv_maxt)
surv_dat$starts=starts
surv_dat$Z=Z

inf_dat$lRhat=log(inf_dat$Rhat)
fit=lm(lRhat~s+I(s^2),data=inf_dat)
theta0_smooth=coef(fit)[1]
theta1_smooth=coef(fit)[2]
theta2_smooth=coef(fit)[3]
inf_dat$Rhat_smooth=exp(theta0_smooth+theta1_smooth*inf_dat$s+theta2_smooth*inf_dat$s^2)
smooth_inf_ses=summary(fit)$coefficients[,2]

opt_inf=optim(par=c(theta0_smooth,theta1_smooth,theta2_smooth),fn=Neg_Log_L_Inf,inf_dat=inf_dat,control=list(maxit=1000))
theta0_stage1=opt_inf$par[1]
theta1_stage1=opt_inf$par[2]
theta2_stage1=opt_inf$par[3]
inf_dat$Rhat_stage1=exp(theta0_stage1+theta1_stage1*inf_dat$s+theta2_stage1*inf_dat$s^2)
stage1_inf_ses=sqrt(diag(solve(-D2_Log_L_Inf(theta0_stage1,theta1_stage1,theta2_stage1,inf_dat=inf_dat))))
stage1_inf_ses

indices=unlist(lapply(1:n,function(i) {return(starts[i]:(starts[i]+surv_maxt-1))}))
cov_dat=data.frame(id=rep(1:n,each=surv_maxt),
                   t=rep(0:(surv_maxt-1),n),
                   I=inf_dat$I[indices],
                   Rhat=inf_dat$Rhat[indices],
                   Rhat_smooth=inf_dat$Rhat_smooth[indices],
                   Rhat_stage1=inf_dat$Rhat_stage1[indices])
temp=tmerge(surv_dat,surv_dat,id=id,tdelta=event(eventtime,status))
dat=tmerge(temp,cov_dat,id=id,I=tdc(t,I),
           Rhat=tdc(t,Rhat),
           Rhat_smooth=tdc(t,Rhat_smooth),
           Rhat_stage1=tdc(t,Rhat_stage1))

initial=phreg(Surv(tstart,tstop,tdelta)~Z+Rhat,data=dat,dist="weibull",shape=1)
beta_raw=initial$coefficients["Z"]
alpha_raw=initial$coefficients["Rhat"]
lambda_raw=1/exp(initial$coefficients["log(scale)"])
raw_ses=c(rev(sqrt(diag(initial$var)))[c("Z","log(scale)","Rhat")],NA,NA,NA)
raw_ses["log(scale)"]=lambda*raw_ses["log(scale)"]

initial=phreg(Surv(tstart,tstop,tdelta)~Z+Rhat_smooth,data=dat,dist="weibull",shape=1)
beta_smooth=initial$coefficients["Z"]
alpha_smooth=initial$coefficients["Rhat_smooth"]
lambda_smooth=1/exp(initial$coefficients["log(scale)"])
smooth_ses=c(rev(sqrt(diag(initial$var)))[c("Z","log(scale)","Rhat_smooth")],smooth_inf_ses)
smooth_ses["log(scale)"]=lambda*smooth_ses["log(scale)"]

initial=phreg(Surv(tstart,tstop,tdelta)~Z+Rhat_stage1,data=dat,dist="weibull",shape=1)
beta_stage1=initial$coefficients["Z"]
alpha_stage1=initial$coefficients["Rhat_stage1"]
lambda_stage1=1/exp(initial$coefficients["log(scale)"])
stage1_ses=c(rev(sqrt(diag(initial$var)))[c("Z","log(scale)","Rhat_stage1")],stage1_inf_ses)
stage1_ses["log(scale)"]=lambda*stage1_ses["log(scale)"]

opt=optim(par=c(beta_stage1,lambda_stage1,alpha_stage1,theta0_stage1,theta1_stage1,theta2_stage1),fn=Neg_Log_L,surv_dat=surv_dat,inf_dat=inf_dat,control=list(maxit=2000))
opt$par
opt$convergence
joint_ses=sqrt(diag(solve(-D2_Log_L(opt$par[1:length(beta)],opt$par[length(beta)+1],opt$par[length(beta)+2],opt$par[length(beta)+3],opt$par[length(beta)+4],opt$par[length(beta)+5],surv_dat=surv_dat,inf_dat=inf_dat))))

alpha_result_table=matrix(c(beta_raw,lambda_raw,alpha_raw,NA,NA,NA,NA,
                           beta_smooth,lambda_smooth,alpha_smooth,theta0_smooth,theta1_smooth,theta2_smooth,NA,
                           beta_stage1,lambda_stage1,alpha_stage1,theta0_stage1,theta1_stage1,theta2_stage1,NA,
                           opt$par,opt$convergence),nrow=4,ncol=7,byrow=T)
alpha_result_table

alpha_ses_table=matrix(c(raw_ses,smooth_ses,stage1_ses,joint_ses),nrow=4,ncol=6,byrow=T)
alpha_ses_table

