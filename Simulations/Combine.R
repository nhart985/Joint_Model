
load("alpha_result_1.Rdata")
raw=alpha_result[1,]
smooth=alpha_result[2,]
stage1=alpha_result[3,]
joint=alpha_result[4,]
for(i in 2:2500) {
    if(file.exists(paste0("alpha_result_",i,".Rdata"))) {
	load(paste0("alpha_result_",i,".Rdata"))
    raw=rbind(raw,alpha_result[1,])
    smooth=rbind(smooth,alpha_result[2,])
    stage1=rbind(stage1,alpha_result[3,])
    joint=rbind(joint,alpha_result[4,])
  }
}

alpha_result_table_raw=as.data.frame(raw)
names(alpha_result_table_raw)=c("beta","lambda","alpha","theta0","theta1","theta2","convergence")
write.csv(alpha_result_table_raw,"alpha_result_table_raw.csv",row.names=F)

alpha_result_table_smooth=as.data.frame(smooth)
names(alpha_result_table_smooth)=c("beta","lambda","alpha","theta0","theta1","theta2","convergence")
write.csv(alpha_result_table_smooth,"alpha_result_table_smooth.csv",row.names=F)

alpha_result_table_stage1=as.data.frame(stage1)
names(alpha_result_table_stage1)=c("beta","lambda","alpha","theta0","theta1","theta2","convergence")
write.csv(alpha_result_table_stage1,"alpha_result_table_stage1.csv",row.names=F)

alpha_result_table_joint=as.data.frame(joint)
names(alpha_result_table_joint)=c("beta","lambda","alpha","theta0","theta1","theta2","convergence")
write.csv(alpha_result_table_joint,"alpha_result_table_joint.csv",row.names=F)

load("alpha_ses_1.Rdata")
raw=alpha_ses[1,]
smooth=alpha_ses[2,]
stage1=alpha_ses[3,]
robust_stage1=alpha_ses[4,]
joint=alpha_ses[5,]
for(i in 2:2500) {
    if(file.exists(paste0("alpha_ses_",i,".Rdata"))) {
	load(paste0("alpha_ses_",i,".Rdata"))
    raw=rbind(raw,alpha_ses[1,])
    smooth=rbind(smooth,alpha_ses[2,])
    stage1=rbind(stage1,alpha_ses[3,])
    robust_stage1=rbind(robust_stage1,alpha_ses[4,])
    joint=rbind(joint,alpha_ses[5,])
  }
}

alpha_ses_table_raw=as.data.frame(raw)
names(alpha_ses_table_raw)=c("beta","lambda","alpha","theta0","theta1","theta2")
write.csv(alpha_ses_table_raw,"alpha_ses_table_raw.csv",row.names=F)

alpha_ses_table_smooth=as.data.frame(smooth)
names(alpha_ses_table_smooth)=c("beta","lambda","alpha","theta0","theta1","theta2")
write.csv(alpha_ses_table_smooth,"alpha_ses_table_smooth.csv",row.names=F)

alpha_ses_table_stage1=as.data.frame(stage1)
names(alpha_ses_table_stage1)=c("beta","lambda","alpha","theta0","theta1","theta2")
write.csv(alpha_ses_table_stage1,"alpha_ses_table_stage1.csv",row.names=F)

alpha_ses_table_robust_stage1=as.data.frame(robust_stage1)
names(alpha_ses_table_robust_stage1)=c("beta","lambda","alpha","theta0","theta1","theta2")
write.csv(alpha_ses_table_robust_stage1,"alpha_ses_table_robust_stage1.csv",row.names=F)

alpha_ses_table_joint=as.data.frame(joint)
names(alpha_ses_table_joint)=c("beta","lambda","alpha","theta0","theta1","theta2")
write.csv(alpha_ses_table_joint,"alpha_ses_table_joint.csv",row.names=F)





load("alpha_result_naive_1.Rdata")
naive=alpha_result[1,]
period1=alpha_result[2,]
period2=alpha_result[3,]
for(i in 2:2500) {
  if(file.exists(paste0("alpha_result_naive_",i,".Rdata"))) {
    load(paste0("alpha_result_",i,".Rdata"))
    naive=rbind(naive,alpha_result[1,])
    period1=rbind(period1,alpha_result[2,])
    period2=rbind(period2,alpha_result[3,])
  }
}

alpha_result_table_naive=as.data.frame(naive)
names(alpha_result_table_naive)=c("beta","lambda","alpha")
write.csv(alpha_result_table_naive,"alpha_result_table_naive.csv",row.names=F)

alpha_result_table_period1=as.data.frame(period1)
names(alpha_result_table_period1)=c("beta","lambda","alpha")
write.csv(alpha_result_table_period1,"alpha_result_table_period1.csv",row.names=F)

alpha_result_table_period2=as.data.frame(period2)
names(alpha_result_table_period2)=c("beta","lambda","alpha")
write.csv(alpha_result_table_period2,"alpha_result_table_period2.csv",row.names=F)

load("alpha_ses_naive_1.Rdata")
naive=alpha_ses[1,]
period1=alpha_ses[2,]
period2=alpha_ses[3,]
for(i in 2:2500) {
  if(file.exists(paste0("alpha_ses_naive_",i,".Rdata"))) {
    load(paste0("alpha_ses_",i,".Rdata"))
    naive=rbind(naive,alpha_ses[1,])
    period1=rbind(period1,alpha_ses[2,])
    period2=rbind(period2,alpha_ses[3,])
  }
}

alpha_ses_table_naive=as.data.frame(naive)
names(alpha_ses_table_naive)=c("beta","lambda","alpha")
write.csv(alpha_ses_table_naive,"alpha_ses_table_naive.csv",row.names=F)

alpha_ses_table_period1=as.data.frame(period1)
names(alpha_ses_table_period1)=c("beta","lambda","alpha")
write.csv(alpha_ses_table_period1,"alpha_ses_table_period1.csv",row.names=F)

alpha_ses_table_period2=as.data.frame(period2)
names(alpha_ses_table_period2)=c("beta","lambda","alpha")
write.csv(alpha_ses_table_period2,"alpha_ses_table_period2.csv",row.names=F)
