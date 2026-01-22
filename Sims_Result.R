library(ggplot2)
library(ggpubr)
library(xtable)

mse=function(x) {return(1000*mean((x-mean(x))^2))}
mape=function(x,v) {return(100*mean(abs((x-v)/v)))}

get_results=function(spath,fun=mape) {
  path1=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_raw.csv")
  path2=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_smooth.csv")
  path3=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_stage1.csv")
  path4=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_joint.csv")
  alpha_result_table_raw=read.csv(path1)
  alpha_result_table_smooth=read.csv(path2)
  alpha_result_table_stage1=read.csv(path3)
  alpha_result_table_joint=read.csv(path4)
  alpha_result_table_joint=alpha_result_table_joint[alpha_result_table_joint$convergence==0,]
  true=c(0.05,0.1,-1,-0.06,0.04,-0.001)
  out=data.frame(method=c("raw","smooth","stage1","joint"))
  for(i in 1:length(true)) {
    out=cbind(out,c(fun(alpha_result_table_raw[,i],v=true[i]),
                    fun(alpha_result_table_smooth[,i],v=true[i]),
                    fun(alpha_result_table_stage1[,i],v=true[i]),
                    fun(alpha_result_table_joint[,i],v=true[i])))
  }
  names(out)=c("method","beta","lambda","alpha","theta0","theta1","theta2")
  return(out)
}

get_ses=function(spath) {
  path1=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_raw.csv")
  path2=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_smooth.csv")
  path3=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_stage1.csv")
  path4=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_joint.csv")
  alpha_ses_table_raw=read.csv(path1)
  alpha_ses_table_smooth=read.csv(path2)
  alpha_ses_table_stage1=read.csv(path3)
  alpha_ses_table_joint=read.csv(path4)
  f=function(x) {
    return(mean(x))
  }
  out=as.data.frame(rbind(apply(alpha_ses_table_raw,2,f),
                          rbind(apply(alpha_ses_table_smooth,2,f),
                                rbind(apply(alpha_ses_table_stage1,2,f),
                                      apply(alpha_ses_table_joint,2,f)))))
  out=out[,!names(out)%in% c("convergence")]
  out$method=c("raw","smooth","stage1","joint")
  return(out)
}

get_coverage=function(spath) {
  path1=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_raw.csv")
  path2=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_smooth.csv")
  path3=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_stage1.csv")
  path4=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_joint.csv")
  alpha_result_table_raw=read.csv(path1)
  alpha_result_table_smooth=read.csv(path2)
  alpha_result_table_stage1=read.csv(path3)
  alpha_result_table_joint=read.csv(path4)
  path1=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_raw.csv")
  path2=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_smooth.csv")
  path3=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_stage1.csv")
  path4=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_robust_stage1.csv")
  path5=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_joint.csv")
  alpha_ses_table_raw=read.csv(path1)
  alpha_ses_table_smooth=read.csv(path2)
  alpha_ses_table_stage1=read.csv(path3)
  alpha_ses_table_robust_stage1=read.csv(path4)
  alpha_ses_table_joint=read.csv(path5)
  cp=matrix(NA,nrow=5,ncol=6)
  true=c(0.05,0.1,-1,-0.06,0.04,-0.001)
  cp_fun=function(result,se,i) {
    lower=result[,i]-1.96*se[,i]
    upper=result[,i]+1.96*se[,i]
    return(mean(lower < true[i] & upper > true[i],na.rm=T))
  }
  for(i in 1:6) {
    cp[1,i]=cp_fun(alpha_result_table_raw,alpha_ses_table_raw,i)
    cp[2,i]=cp_fun(alpha_result_table_smooth,alpha_ses_table_smooth,i)
    cp[3,i]=cp_fun(alpha_result_table_stage1,alpha_ses_table_stage1,i)
    cp[4,i]=cp_fun(alpha_result_table_stage1,alpha_ses_table_robust_stage1,i)
    cp[5,i]=cp_fun(alpha_result_table_joint,alpha_ses_table_joint,i)
  }
  cp=as.data.frame(cp)
  names(cp)=c("beta","lambda","alpha","theta0","theta1","theta2")
  cp$method=c("raw","smooth","stage1","robust stage1","joint")
  return(cp)
}

get_results_naive=function(spath,fun=mape) {
  path1=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_naive.csv")
  path2=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_period1.csv")
  path3=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_period2.csv")
  alpha_result_table_naive=read.csv(path1)
  alpha_result_table_period1=read.csv(path2)
  alpha_result_table_period2=read.csv(path3)
  true=c(0.05,0.1)
  out=data.frame(method=c("naive","period1","period2"))
  for(i in 1:length(true)) {
    out=cbind(out,c(fun(alpha_result_table_naive[,i],v=true[i]),
                    fun(alpha_result_table_period1[,i],v=true[i]),
                    fun(alpha_result_table_period2[,i],v=true[i])))
  }
  names(out)=c("method","beta","lambda")
  return(out)
}

get_coverage_naive=function(spath) {
  path1=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_naive.csv")
  path2=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_period1.csv")
  path3=paste0("//Joint Model//Simulation Results//",spath,"alpha_result_table_period2.csv")
  alpha_result_table_naive=read.csv(path1)
  alpha_result_table_period1=read.csv(path2)
  alpha_result_table_period2=read.csv(path3)
  path1=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_naive.csv")
  path2=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_period1.csv")
  path3=paste0("//Joint Model//Simulation Results//",spath,"alpha_ses_table_period2.csv")
  alpha_ses_table_naive=read.csv(path1)
  alpha_ses_table_period1=read.csv(path2)
  alpha_ses_table_period2=read.csv(path3)
  cp=matrix(NA,nrow=3,ncol=2)
  true=c(0.05,0.1)
  cp_fun=function(result,se,i) {
    lower=result[,i]-1.96*se[,i]
    upper=result[,i]+1.96*se[,i]
    return(mean(lower < true[i] & upper > true[i],na.rm=T))
  }
  for(i in 1:2) {
    cp[1,i]=cp_fun(alpha_result_table_naive,alpha_ses_table_naive,i)
    cp[2,i]=cp_fun(alpha_result_table_period1,alpha_ses_table_period1,i)
    cp[3,i]=cp_fun(alpha_result_table_period2,alpha_ses_table_period2,i)
  }
  cp=as.data.frame(cp)
  names(cp)=c("beta","lambda")
  cp$method=c("naive","period1","period2")
  return(cp)
}

######
#Table
######
get_results("I0=10//Sample Size//n=1000//")
get_results("I0=10//Sample Size//n=5000//")
get_results("I0=10//Sample Size//n=10000//")

get_results_naive("I0=10//Sample Size//n=1000//")
get_results_naive("I0=10//Sample Size//n=5000//")
get_results_naive("I0=10//Sample Size//n=10000//")

get_coverage("I0=10//Sample Size//n=1000//")
get_coverage("I0=10//Sample Size//n=5000//")
get_coverage("I0=10//Sample Size//n=10000//")

get_coverage_naive("I0=10//Sample Size//n=1000//")
get_coverage_naive("I0=10//Sample Size//n=5000//")
get_coverage_naive("I0=10//Sample Size//n=10000//")

print(xtable(get_results("I0=10//Sample Size//n=1000//"),digits=1),include.rownames=F)
print(xtable(get_results("I0=10//Sample Size//n=5000//"),digits=1),include.rownames=F)
print(xtable(get_results("I0=10//Sample Size//n=10000//"),digits=1),include.rownames=F)

print(xtable(get_results_naive("I0=10//Sample Size//n=1000//"),digits=1),include.rownames=F)
print(xtable(get_results_naive("I0=10//Sample Size//n=5000//"),digits=1),include.rownames=F)
print(xtable(get_results_naive("I0=10//Sample Size//n=10000//"),digits=1),include.rownames=F)

print(xtable(get_coverage("I0=10//Sample Size//n=1000//"),digits=3),include.rownames=F)
print(xtable(get_coverage("I0=10//Sample Size//n=5000//"),digits=3),include.rownames=F)
print(xtable(get_coverage("I0=10//Sample Size//n=10000//"),digits=3),include.rownames=F)

print(xtable(get_coverage_naive("I0=10//Sample Size//n=1000//"),digits=3),include.rownames=F)
print(xtable(get_coverage_naive("I0=10//Sample Size//n=5000//"),digits=3),include.rownames=F)
print(xtable(get_coverage_naive("I0=10//Sample Size//n=10000//"),digits=3),include.rownames=F)

bias=function(result,vals=c(0.05,0.1,-1,-0.06,0.04,-0.001)) {
  result$beta=100*abs(result$beta-vals[1])/abs(vals[1])
  result$lambda=100*abs(result$lambda-vals[2])/abs(vals[2])
  result$alpha=100*abs(result$alpha-vals[3])/abs(vals[3])
  result$theta0=100*abs(result$theta0-vals[4])/abs(vals[4])
  result$theta1=100*abs(result$theta1-vals[5])/abs(vals[5])
  result$theta2=100*abs(result$theta2-vals[6])/abs(vals[6])
  return(result)
}
result=bias(get_results("Bias//",fun=mean))
print(xtable(result,digits=1),include.rownames=F)

bias_naive=function(result,vals=c(0.05,0.1)) {
  result$beta=100*abs(result$beta-vals[1])/abs(vals[1])
  result$lambda=100*abs(result$lambda-vals[2])/abs(vals[2])
  return(result)
}
result=bias_naive(get_results_naive("Bias//",fun=mean))
print(xtable(result,digits=1),include.rownames=F)







