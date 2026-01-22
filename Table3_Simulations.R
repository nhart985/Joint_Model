library(survival)
library(eha)
library(simsurv)
library(MASS)
source("Functions.R")

k=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

n=10000
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
  Rbias=exp(theta0+theta1*s+1.5*theta2*(s^2)+0.00001*s^3)
  I=vector("numeric")
  I[1]=rpois(1,Rbias[1]*I0)
  for(i in 2:inf_maxt) {
    I[i]=rpois(1,Rbias[i]*I[i-1])
  }
  result=data.frame(s,I,R)
  return(result)
}

inf_dat=Gen_Pois(inf_maxt,I0,theta0,theta1,theta2)
while(sum(inf_dat$I==0)>0) {
  inf_dat=Gen_Pois(inf_maxt,I0,theta0,theta1,theta2)
}
inf_dat$I_prev=c(I0,inf_dat$I[1:(inf_maxt-1)])
inf_dat$Rhat=inf_dat$I/inf_dat$I_prev

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

n_test=500
starts_test=sample(2:(inf_maxt-surv_maxt),n_test,replace=T)
R_test=exp(theta0+theta1*inf_dat$s+theta2*(inf_dat$s^2))
Z_test=rnorm(n,R[starts],0.25)
Xmat_test=data.frame(id=1:n_test,starts=starts_test,Z=Z_test)
test=simsurv(hazard=hazard,x=Xmat_test,betas=betas,lambda=lambda,inf_dat=inf_dat,maxt=surv_maxt)
test$starts=starts_test
test$Z=Z_test

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
sigma=solve(-D2_Log_L_Inf(theta0_stage1,theta1_stage1,theta2_stage1,inf_dat=inf_dat))
stage1_inf_ses=sqrt(diag(sigma))

indices=unlist(lapply(1:n,function(i) {return(starts[i]:(starts[i]+surv_maxt-1))}))
cov_dat=data.frame(id=rep(1:n,each=surv_maxt),
                   t=rep(0:(surv_maxt-1),n),
                   I=inf_dat$I[indices],
                   s=inf_dat$s[indices],
                   Rhat=inf_dat$Rhat[indices],
                   Rhat_smooth=inf_dat$Rhat_smooth[indices],
                   Rhat_stage1=inf_dat$Rhat_stage1[indices])
temp=tmerge(surv_dat,surv_dat,id=id,tdelta=event(eventtime,status))
dat=tmerge(temp,cov_dat,id=id,I=tdc(t,I),s=tdc(t,s),
           Rhat=tdc(t,Rhat),
           Rhat_smooth=tdc(t,Rhat_smooth),
           Rhat_stage1=tdc(t,Rhat_stage1))

indices=unlist(lapply(1:n_test,function(i) {return(starts_test[i]:(starts_test[i]+surv_maxt-1))}))
cov_dat=data.frame(id=rep(1:n_test,each=surv_maxt),
                   t=rep(0:(surv_maxt-1),n_test),
                   I=inf_dat$I[indices],
                   s=inf_dat$s[indices],
                   Rhat=inf_dat$Rhat[indices],
                   Rhat_smooth=inf_dat$Rhat_smooth[indices],
                   Rhat_stage1=inf_dat$Rhat_stage1[indices])
temp=tmerge(test,test,id=id,tdelta=event(eventtime,status))
test_dat=tmerge(temp,cov_dat,id=id,I=tdc(t,I),s=tdc(t,s),
                Rhat=tdc(t,Rhat),
                Rhat_smooth=tdc(t,Rhat_smooth),
                Rhat_stage1=tdc(t,Rhat_stage1))

initial=phreg(Surv(tstart,tstop,tdelta)~Z,data=dat,dist="weibull",shape=1)
beta_naive=initial$coefficients["Z"]
lambda_naive=1/exp(initial$coefficients["log(scale)"])
naive_ses=c(rev(sqrt(diag(initial$var)))[c("Z","log(scale)")],NA)
naive_ses["log(scale)"]=lambda*naive_ses["log(scale)"]

dat$period=as.numeric(dat$s >= 15 & dat$s <= 25)
initial=phreg(Surv(tstart,tstop,tdelta)~Z+period,data=dat,dist="weibull",shape=1)
beta_period1=initial$coefficients["Z"]
alpha_period1=initial$coefficients["period"]
lambda_period1=1/exp(initial$coefficients["log(scale)"])
period1_ses=c(rev(sqrt(diag(initial$var)))[c("Z","log(scale)","period")])
period1_ses["log(scale)"]=lambda*period1_ses["log(scale)"]

dat$period=as.numeric(dat$s >= 10 & dat$s <= 30)
initial=phreg(Surv(tstart,tstop,tdelta)~Z+period,data=dat,dist="weibull",shape=1)
beta_period2=initial$coefficients["Z"]
alpha_period2=initial$coefficients["period"]
lambda_period2=1/exp(initial$coefficients["log(scale)"])
period2_ses=c(rev(sqrt(diag(initial$var)))[c("Z","log(scale)","period")])
period2_ses["log(scale)"]=lambda*period2_ses["log(scale)"]

alpha_result_naive=matrix(c(beta_naive,lambda_naive,NA,
                            beta_period1,lambda_period1,alpha_period1,
                            beta_period2,lambda_period2,alpha_period2),
                          nrow=3,ncol=3,byrow=T)
save(alpha_result_naive,file=paste0("alpha_result_naive_",k,".Rdata"))

alpha_ses_naive=matrix(c(naive_ses,period1_ses,period2_ses),nrow=3,ncol=3,byrow=T)
save(alpha_ses_naive,file=paste0("alpha_ses_naive_",k,".Rdata"))




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
robust_stage1_ses=sqrt(stage1_ses^2+c(Boot(dat,opt_inf$par,sigma)^2,0,0,0))

opt=optim(par=c(beta_stage1,lambda_stage1,alpha_stage1,theta0_stage1,theta1_stage1,theta2_stage1),fn=Neg_Log_L,surv_dat=surv_dat,inf_dat=inf_dat,control=list(maxit=1000))
joint_ses=sqrt(diag(solve(-D2_Log_L(opt$par[1:length(beta)],opt$par[length(beta)+1],opt$par[length(beta)+2],opt$par[length(beta)+3],opt$par[length(beta)+4],opt$par[length(beta)+5],surv_dat=surv_dat,inf_dat=inf_dat))))
beta_joint=opt$par[1]
alpha_joint=opt$par[3]
theta0_joint=opt$par[4]
theta1_joint=opt$par[5]
theta2_joint=opt$par[6]

alpha_result=matrix(c(beta_raw,lambda_raw,alpha_raw,NA,NA,NA,NA,
                      beta_smooth,lambda_smooth,alpha_smooth,theta0_smooth,theta1_smooth,theta2_smooth,NA,
                      beta_stage1,lambda_stage1,alpha_stage1,theta0_stage1,theta1_stage1,theta2_stage1,NA,
                      opt$par,opt$convergence),nrow=4,ncol=7,byrow=T)
save(alpha_result,file=paste0("alpha_result_",k,".Rdata"))

alpha_ses=matrix(c(raw_ses,smooth_ses,stage1_ses,robust_stage1_ses,joint_ses),nrow=5,ncol=6,byrow=T)
save(alpha_ses,file=paste0("alpha_ses_",k,".Rdata"))