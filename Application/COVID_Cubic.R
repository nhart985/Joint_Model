library(survival)
library(haven)
library(ggplot2)
library(ggpubr)
library(eha)
library(risksetROC)
library(xtable)
library(data.table)
library(risksetROC)
library(dplyr)

source("//Documents//Joint Model//Real_Functions_Cubic.R")

###########
#COVID Data
###########
covid=read.csv("//Downloads//time_series_covid19_confirmed_US.csv")
dates=names(covid)[substr(names(covid),1,1)=="X"]
dates=substr(dates,start=2,stop=9)
dates=as.Date(dates,format="%m.%d.%y")
names(covid)[substr(names(covid),1,1)=="X"]=as.character(dates)
dates=dates[seq(1,length(dates),by=7)]
dates=dates[dates >= as.Date("2020-03-01") & dates <= as.Date("2020-08-01")]
us_counts=diff(apply(covid[,as.character(dates)],2,sum,na.rm=T))
I0=us_counts[1]
Rhat=us_counts[2:length(us_counts)]/us_counts[1:(length(us_counts)-1)]
dates=dates[3:length(dates)]
us_counts=us_counts[2:length(us_counts)]
s=1:length(Rhat)
logR=log(Rhat)
inf_dat=data.frame(Rhat=Rhat,s=s,I=us_counts)
inf_dat$I_prev=c(I0,inf_dat$I[1:(length(inf_dat$I)-1)])

fig_dat=data.frame(lRhat=logR,s=s,dates=format(dates,"%m-%d-%Y"))
g=ggplot(fig_dat)+geom_point(aes(x=s,y=lRhat),size=1.75)
g=g+stat_smooth(aes(x=s,y=lRhat))
g=g+theme_classic()+xlab("Calendar Time")+ylab(bquote(log(hat(R)[s])))
g=g+theme(text=element_text(size=24),axis.text.x=element_text(angle=90),axis.title.x=element_text(vjust=0.75))
g=g+scale_x_continuous(breaks=seq(1,length(fig_dat$dates),by=2),labels=fig_dat$dates[seq(1,length(fig_dat$dates),by=2)])

########
#Tx Data
########
df=as.data.frame(read_sas("//Downloads//SAS Dataset 202406//Liver//liver_data.sas7bdat"))
dat=df[!is.na(df$INIT_DATE) & !is.na(df$END_DATE),]
dat=dat[dat$END_DATE > as.Date(dates[1]) & dat$INIT_DATE < as.Date(dates[length(dates)]),]
dat$INIT_DATE=as.Date(dat$INIT_DATE)
dat$TX_DATE=as.Date(dat$TX_DATE)
dat$END_DATE=as.Date(dat$END_DATE)

##################
#Outcome Variables
##################
dat$time=as.numeric(dat$END_DATE-dat$INIT_DATE)/7
dat$tx_time=as.numeric(dat$TX_DATE-dat$INIT_DATE)/7
dat$statusDD=as.numeric(!is.na(dat$TX_DATE) & dat$tx_time==dat$time & dat$REM_CD %in% c(2,4,18))

#####
#Plot
#####
tx=vector("numeric")
pweeks=vector("numeric")
for(i in 1:length(dates)) {
  tx[i]=sum(dat$END_DATE > as.Date(dates[i]) & dat$END_DATE < (as.Date(dates[i])+7) & dat$statusDD==1)
  origin=pmax(dat$INIT_DATE,as.Date(dates[i]))
  pweeks[i]=sum(as.numeric(dat$END_DATE[dat$INIT_DATE < (as.Date(dates[i])+7) & dat$END_DATE > as.Date(dates[i])]-origin[dat$INIT_DATE < (as.Date(dates[i])+7) & dat$END_DATE > as.Date(dates[i])])/7)
}

tr=1000*tx/pweeks
fig_dat=data.frame(tr=tr,s=s,dates=format(dates,"%m-%d-%Y"))
g2=ggplot(fig_dat)+geom_point(aes(x=s,y=tr),size=1.75)
g2=g2+stat_smooth(aes(x=s,y=tr))
g2=g2+theme_classic()+xlab("Calendar Time")+ylab("Transplant Rate")
g2=g2+theme(text=element_text(size=24),axis.text.x=element_text(angle=90),axis.title.x=element_text(vjust=0.75))
g2=g2+scale_x_continuous(breaks=seq(1,length(fig_dat$dates),by=2),labels=fig_dat$dates[seq(1,length(fig_dat$dates),by=2)])

ggarrange(g,g2,nrow=1)
ggsave(file="//Documents//Joint Model//SideBySide.pdf",width=18,height=8)

#########################
#Censoring and Truncation
#########################
dat$statusDD[dat$END_DATE > as.Date(dates[length(dates)])]=0
dat$time[dat$END_DATE > as.Date(dates[length(dates)])]=as.numeric(as.Date(dates[length(dates)])-dat$INIT_DATE[dat$END_DATE > as.Date(dates[length(dates)])])/7
dat$lt=pmax(as.numeric(as.Date(dates[1])-dat$INIT_DATE)/7,0)
dat=dat[dat$time > 0,]

####################
#Predictor Variables
####################
dat$DIABETES=as.numeric(dat$DIAB!=1)
dat$DIABETES[dat$DIAB==998 & is.na(dat$DIAB)]=NA

dat$BMI=as.numeric(cut(as.numeric(dat$INIT_BMI_CALC),c(0,18.5,25,30,Inf),labels=1:4))

dat$PREV=as.numeric(dat$NUM_PREV_TX > 0)
dat$PREV[is.na(dat$NUM_PREV_TX)]=NA

dat$MELD=pmin(pmax(dat$INIT_MELD_PELD_LAB_SCORE,6),40)

dat$id=1:dim(dat)[1]

surv_dat=data.frame(eventtime=dat$time,status=dat$statusDD,lt=dat$lt,
                    Z1=dat$MELD,Z2=dat$GENDER=="F")
dat=dat[complete.cases(surv_dat),]
surv_dat=surv_dat[complete.cases(surv_dat),]
surv_dat$id=1:dim(surv_dat)[1]

f=function(x) {
  out=1
  if(x >= dates[1]) {
    out=max((1:length(dates))[x >= dates])
  } 
  return(out)
}
surv_dat$starts=sapply(split(dat$INIT_DATE,dat$id),f)

#Models
inf_dat$lRhat=log(inf_dat$Rhat)
fit=lm(lRhat~s+I(s^2)+I(s^3),data=inf_dat)
theta0_smooth=coef(fit)[1]
theta1_smooth=coef(fit)[2]
theta2_smooth=coef(fit)[3]
theta3_smooth=coef(fit)[4]
inf_dat$Rhat_smooth=exp(theta0_smooth+theta1_smooth*inf_dat$s+theta2_smooth*inf_dat$s^2+theta3_smooth*inf_dat$s^3)
phi_smooth=summary(glm.nb(I~offset(log(Rhat_smooth)+log(I_prev)),data=inf_dat))$theta
sigma=solve(D2_Inf_Log_L(theta0_smooth,theta1_smooth,theta2_smooth,theta3_smooth,phi_smooth,inf_dat))
smooth_inf_ses=sqrt(diag(sigma))
Vcov_smooth=sigma

inf_dat$s2=inf_dat$s^2
inf_dat$s3=inf_dat$s^3
fit=glm.nb(I~s+s2+s3+offset(log(I_prev)),data=inf_dat)
theta0_stage1=coef(fit)[1]
theta1_stage1=coef(fit)[2]
theta2_stage1=coef(fit)[3]
theta3_stage1=coef(fit)[4]
inf_dat$Rhat_stage1=exp(theta0_stage1+theta1_stage1*inf_dat$s+theta2_stage1*inf_dat$s^2+theta3_stage1*inf_dat$s^3)
phi_stage1=summary(fit)$theta
sigma=solve(D2_Inf_Log_L(theta0_stage1,theta1_stage1,theta2_stage1,theta3_stage1,phi_stage1,inf_dat))
stage1_inf_ses=sqrt(diag(sigma))
Vcov_stage1=sigma

n=dim(surv_dat)[1]
indices=unlist(lapply(1:n,function(i) {return((surv_dat$starts[i]):(surv_dat$starts[i]+floor(surv_dat$eventtime[i])))}))
cov_dat=data.frame(id=rep(1:n,floor(surv_dat$eventtime)+1),
                   t=indices-rep(surv_dat$starts,floor(surv_dat$eventtime)+1),
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
dat$tstart=dat$tstart+dat$lt
dat$tstop=dat$tstop+dat$lt

surv_dat$study_time=surv_dat$eventtime-surv_dat$lt+0.5
n=dim(surv_dat)[1]
indices=unlist(lapply(1:n,function(i) {return((surv_dat$starts[i]):(surv_dat$starts[i]+floor(surv_dat$study_time[i])))}))
cov_dat=data.frame(id=rep(1:n,floor(surv_dat$study_time)+1),
                   t=indices-rep(surv_dat$starts,floor(surv_dat$study_time)+1),
                   I=inf_dat$I[indices],
                   s=inf_dat$s[indices],
                   Rhat=inf_dat$Rhat[indices],
                   Rhat_smooth=inf_dat$Rhat_smooth[indices],
                   Rhat_stage1=inf_dat$Rhat_stage1[indices])
temp=tmerge(surv_dat,surv_dat,id=id,tdelta=event(study_time,status))
dat=tmerge(temp,cov_dat,id=id,I=tdc(t,I),s=tdc(t,s),
           Rhat=tdc(t,Rhat),
           Rhat_smooth=tdc(t,Rhat_smooth),
           Rhat_stage1=tdc(t,Rhat_stage1))
dat$tstart=dat$tstart+dat$lt
dat$tstop=dat$tstop+dat$lt

fit=coxph(Surv(tstart,tstop,tdelta)~Z1+Z2,data=dat,method="breslow")
beta_naive=coef(fit)[1:2]
naive_ses=c(summary(fit)$coefficients[,"se(coef)"],NA)
beta_naive_initial=coef(fit)[1:2]
alpha_naive_initial=coef(fit)[3]
opt=optim(par=beta_naive_initial,fn=Naive_Neg_Log_L,surv_dat=surv_dat,period=1,naive=T,hessian=T)
beta_naive=opt$par[1:2]
naive_ses=c(sqrt(diag(solve(opt$hessian))),NA)

dates[s >= 1 & s <= 5  | s >= 17 & s <= 20]
dat$covid1=as.numeric(dat$s >= 1 & dat$s <= 5 | dat$s >= 17 & dat$s <= 20)
fit=coxph(Surv(tstart,tstop,tdelta)~Z1+Z2+covid1,data=dat,method="breslow")
beta_period1_initial=coef(fit)[1:2]
alpha_period1_initial=coef(fit)[3]
opt=optim(par=c(beta_period1_initial,alpha_period1_initial),fn=Naive_Neg_Log_L,surv_dat=surv_dat,period=1,hessian=T)
beta_period1=opt$par[1:2]
alpha_period1=opt$par[3]
period1_ses=sqrt(diag(solve(opt$hessian)))

dates[s >= 1 & s <= 13]
dat$covid2=as.numeric(dat$s >= 1 & dat$s <= 13)
fit=coxph(Surv(tstart,tstop,tdelta)~Z1+Z2+covid2,data=dat,method="breslow")
beta_period2_initial=coef(fit)[1:2]
alpha_period2_initial=coef(fit)[3]
opt=optim(par=c(beta_period2_initial,alpha_period2_initial),fn=Naive_Neg_Log_L,surv_dat=surv_dat,period=2,hessian=T)
beta_period2=opt$par[1:2]
alpha_period2=opt$par[3]
period2_ses=sqrt(diag(solve(opt$hessian)))

fit=coxph(Surv(tstart,tstop,tdelta)~Z1+Z2+Rhat,data=dat,method="breslow")
beta_raw_initial=coef(fit)[1:2]
alpha_raw_initial=coef(fit)[3]
raw_ses=c(summary(fit)$coefficients[,"se(coef)"],NA,NA,NA)
opt=optim(par=c(beta_raw_initial,alpha_raw_initial),fn=Test_Neg_Log_L,theta0=NA,theta1=NA,theta2=NA,theta3=NA,surv_dat=surv_dat,raw=T,inf_dat=inf_dat,hessian=T)
beta_raw=opt$par[1:2]
alpha_raw=opt$par[3]
raw_ses=c(sqrt(diag(solve(opt$hessian))),NA,NA,NA)

fit=coxph(Surv(tstart,tstop,tdelta)~Z1+Z2+Rhat_smooth,data=dat,method="breslow")
beta_smooth_initial=coef(fit)[1:2]
alpha_smooth_initial=coef(fit)[3]
opt=optim(par=c(beta_smooth_initial,alpha_smooth_initial),fn=Test_Neg_Log_L,theta0=theta0_smooth,theta1=theta1_smooth,theta2=theta2_smooth,theta3=theta3_smooth,surv_dat=surv_dat)
beta_smooth=opt$par[1:2]
alpha_smooth=opt$par[3]
temp_var=diag(solve(Test_D2_Log_L(beta_smooth,alpha_smooth,theta0_smooth,theta1_smooth,theta2_smooth,theta3_smooth,surv_dat=surv_dat,inf_dat=inf_dat)))
smooth_ses=c(sqrt(temp_var),smooth_inf_ses)

fit=coxph(Surv(tstart,tstop,tdelta)~Z1+Z2+Rhat_stage1,data=dat,method="breslow")
beta_stage1_initial=coef(fit)[1:2]
alpha_stage1_initial=coef(fit)[3]
opt=optim(par=c(beta_stage1_initial,alpha_stage1_initial),fn=Test_Neg_Log_L,theta0=theta0_stage1,theta1=theta1_stage1,theta2=theta2_stage1,theta3=theta3_stage1,surv_dat=surv_dat)
beta_stage1=opt$par[1:2]
alpha_stage1=opt$par[3]
temp_var=diag(solve(Test_D2_Log_L(beta_stage1,alpha_stage1,theta0_stage1,theta1_stage1,theta2_stage1,theta3_stage1,surv_dat=surv_dat,inf_dat=inf_dat)))
stage1_ses=c(sqrt(temp_var),stage1_inf_ses)

opt=optim(par=c(beta_stage1_initial,alpha_stage1_initial,theta0_stage1,theta1_stage1,theta2_stage1,theta3_stage1,phi_stage1),fn=Neg_Log_L,surv_dat=surv_dat,inf_dat=inf_dat)
save(opt,file="//Documents//Joint Model//Cubic//joint_opt.Rdata")
result=opt$par[1:7]
phi_joint=opt$par[8]
Vcov_joint=solve(D2_Log_L(opt$par[1:length(beta_stage1)],opt$par[length(beta_stage1)+1],opt$par[length(beta_stage1)+2],opt$par[length(beta_stage1)+3],opt$par[length(beta_stage1)+4],opt$par[length(beta_stage1)+5],opt$par[length(beta_stage1)+6],surv_dat=surv_dat,inf_dat=inf_dat))
joint_ses=sqrt(diag(Vcov_joint))
Vcov_joint=Vcov_joint[4:7,4:7]

Vcov_list=list(Vcov_smooth,Vcov_stage1,Vcov_joint)

alpha_result_table_naive=matrix(c(beta_naive,NA,
                                  beta_period1,alpha_period1,
                                  beta_period2,alpha_period2),
                                nrow=3,ncol=3,byrow=T)
alpha_result_table_naive
write.csv(alpha_result_table_naive,"//Documents//Joint Model//Cubic//alpha_result_table_naive.csv",row.names=F)

alpha_ses_table_naive=matrix(c(naive_ses,period1_ses,period2_ses),nrow=3,ncol=3,byrow=T)
write.csv(alpha_ses_table_naive,"//Documents//Joint Model//Cubic//alpha_ses_table_naive.csv",row.names=F)

alpha_result_table=matrix(c(beta_raw,alpha_raw,NA,NA,NA,NA,
                            beta_smooth,alpha_smooth,theta0_smooth,theta1_smooth,theta2_smooth,theta3_smooth,
                            beta_stage1,alpha_stage1,theta0_stage1,theta1_stage1,theta2_stage1,theta3_stage1,
                            result),nrow=4,ncol=7,byrow=T)
write.csv(alpha_result_table,"//Documents//Joint Model//Cubic//alpha_result_table.csv",row.names=F)

alpha_ses_table=matrix(c(raw_ses,smooth_ses,stage1_ses,joint_ses),nrow=4,ncol=7,byrow=T)
write.csv(alpha_ses_table,"//Documents//Joint Model//Cubic//alpha_ses_table.csv",row.names=F)

Vcov_list=list(Vcov_smooth,Vcov_stage1,Vcov_joint)
save(Vcov_list,file="//Documents//Joint Model//Cubic//Vcov_list.Rdata")

alpha_result_table_naive=read.csv("//Documents//Joint Model//Cubic//alpha_result_table_naive.csv")
alpha_ses_table_naive=read.csv("//Documents//Joint Model//Cubic//alpha_ses_table_naive.csv")
alpha_result_table=read.csv("//Documents//Joint Model//Cubic//alpha_result_table.csv")
alpha_ses_table=read.csv("//Documents//Joint Model//Cubic//alpha_ses_table.csv")

load("//Documents//Joint Model//Cubic//Vcov_list.Rdata")

get_tab=function(est,ses) {
  hr=format(round(exp(est),2),nsmall=2)[1:3]
  lower=format(round(exp(est-1.96*ses),2),nsmall=2)[1:3]
  upper=format(round(exp(est+1.96*ses),2),nsmall=2)[1:3]
  return(paste0(hr," (",lower,", ",upper,")"))
}

surv_sub_naive=as.data.frame(t(as.matrix(data.frame(naive=get_tab(alpha_result_table_naive[1,],alpha_ses_table_naive[1,]),
                                                    period1=get_tab(alpha_result_table_naive[2,],alpha_ses_table_naive[2,]),
                                                    period2=get_tab(alpha_result_table_naive[3,],alpha_ses_table_naive[3,])))))
names(surv_sub_naive)=c("MELD","Sex","COVID")
row.names(surv_sub_naive)=c("Hazard-Only","Period1","Period2")
print(xtable(surv_sub_naive))

surv_sub=as.data.frame(t(as.matrix(data.frame(raw=get_tab(alpha_result_table[1,],alpha_ses_table[1,]),
                                              smooth=get_tab(alpha_result_table[2,],alpha_ses_table[2,]),
                                              stage1=get_tab(alpha_result_table[3,],alpha_ses_table[3,]),
                                              joint=get_tab(alpha_result_table[4,],alpha_ses_table[4,])))))
names(surv_sub)=c("MELD","Sex","R_s")
row.names(surv_sub)=c("Raw","Smooth","Two-Stage","Joint")
print(xtable(surv_sub))

#Goodness of Fit
run_fit=function(i,raw=F) {
  beta=as.numeric(alpha_result_table[i,1:2])
  alpha=as.numeric(alpha_result_table[i,3])
  theta0=as.numeric(alpha_result_table[i,4])
  theta1=as.numeric(alpha_result_table[i,5])
  theta2=as.numeric(alpha_result_table[i,6])
  theta3=as.numeric(alpha_result_table[i,7])
  prof=Surv_Log_L(beta,alpha,theta0,theta1,theta2,theta3,surv_dat,raw=raw,inf_dat=inf_dat)
  prof_AIC=2*sum(!is.na(alpha_result_table[i,]))-2*prof
  if(i==3 | i==4) {
    phi=ifelse(i==3,phi_stage1,phi_joint)
    inf=Inf_Log_L(theta0,theta1,theta2,theta3,phi,inf_dat,raw=raw)
    inf_AIC=2*5-2*inf
  } else {
    inf=NA
    inf_AIC=NA
  }
  
  return(list(prof=prof,prof_AIC=prof_AIC,
              inf=inf,inf_AIC=inf_AIC))
}
fit1=run_fit(1,raw=T)
fit2=run_fit(2)
fit3=run_fit(3)
fit4=run_fit(4)
fit_stats=data.table(inf=c(fit1$inf,fit2$inf,fit3$inf,fit4$inf),
                     inf_AIC=c(fit1$inf_AIC,fit2$inf_AIC,fit3$inf_AIC,fit4$inf_AIC),
                     prof=c(fit1$prof,fit2$prof,fit3$prof,fit4$prof),
                     prof_AIC=c(fit1$prof_AIC,fit2$prof_AIC,fit3$prof_AIC,fit4$prof_AIC))

run_fit_Naive=function(i) {
  beta=as.numeric(alpha_result_table_naive[i,1:2])
  alpha=as.numeric(alpha_result_table_naive[i,3])
  if(i==1) {
    alpha=0 
  }
  prof=Surv_Log_L_Naive(beta,alpha,surv_dat,period=i-1)
  prof_AIC=2*sum(!is.na(alpha_result_table_naive[i,]))-2*prof
  return(list(prof=prof,prof_AIC=prof_AIC))
}
fit1=run_fit_Naive(1)
fit2=run_fit_Naive(2)
fit3=run_fit_Naive(3)
fit_stats_naive=data.table(inf=NA,
                           inf_AIC=NA,
                           prof=c(fit1$prof,fit2$prof,fit3$prof),
                           prof_AIC=c(fit1$prof_AIC,fit2$prof_AIC,fit3$prof_AIC))

fit_stats=rbind(fit_stats_naive,fit_stats)
row.names(fit_stats)=c("Hazard","Period 1","Period 2","Raw","Smooth","Stage 1","Joint")
print(xtable(fit_stats))

