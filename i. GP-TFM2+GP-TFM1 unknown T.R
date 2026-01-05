# Simulating the data
set.seed(199999)
#install.packages("mvtnorm")
library(mvtnorm)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(LaplacesDemon)
library(cowplot)
library(magick)
library(grid)
library(rmutil)
library(coda)
rgenpois_custom <- function(n, lambda, psi) {
  y <- integer(n)
  for (i in 1:n) {
    repeat {
      U <- runif(1)
      X <- 0
      p0 <- exp(-lambda)
      F <- p0
      while (U > F) {
        X <- X + 1
        p0 <- (lambda + psi * X) * p0 / X
        if (p0 <= 0) break
        F <- F + p0
      }
      if (p0 > 0) {
        y[i] <- X
        break
      }
    }
  }
  return(y)
}



constraint_con_delta <- function(star){
  return(star[1]>0 & star[2]>0 & (star[3]>0 & star[3]<1))
}

constraint_time_star <- function(time_star){
  return((time_star[1] > 5 & time_star[1] < 50) & # 45 
           ((time_star[2] > (time_star[1] + 50)) & time_star[2] < 135))
}

nosim  = 1
result = NULL
DIC=matrix(0,nosim,1)
S1	 = matrix(0,nosim,5) # omega
S2	 = matrix(0,nosim,5) # alpha
S3	 = matrix(0,nosim,5) # beta
S4	 = matrix(0,nosim,5) # psi
S5	 = matrix(0,nosim,5) # con 1
S6   = matrix(0,nosim,5) # con 2 
S7   = matrix(0,nosim,5) # delta 1
S8	 = matrix(0,nosim,5) # T1 
S9	 = matrix(0,nosim,5) # T2

r1	= NULL
r2	= NULL
r3  = NULL
r4  = NULL
r5  = NULL
r6	= NULL
r7  = NULL
r8  = NULL
r9	= NULL

# Likelihood for 2B 2A
liko<-function(nob,om,alpha,beta,yt,psi,con,delta,time){
  liko=0
  Pt[,1]=c(rep(0,round(time[1])-1),1,rep(1,nob-round(time[1]))) # 2B
  Pt[,2]=c(rep(0,round(time[2])-1),1,rep(0,nob-round(time[2]))) # 2A
  for (t in 2:nob)
  { 
    
    if(t<time[1]){
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
    }else if(t>=time[1] & t<time[2]){
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con[1]*((delta[1])^(t-time[1]))*Pt[t,1] # 2B
    }else{
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con[2]*Pt[t,2] # 2A
    }
    liko= liko+log((lam_o[t]))+(yt[t]-1)*log((lam_o[t])+psi*yt[t])-(lam_o[t])-psi*yt[t]
  }
  return(liko)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

df <- read.csv("C:/Users/ACER/Desktop/Dipolog crimes.csv", stringsAsFactors = FALSE)
df <- df[,c(1,20)]
par(mfrow=c(1,1))
plot(df$COMPREHENSIVE.DANGEROUS.DRUGS.ACT.OF.2002....RA.9165, type = "l")
drug_data <- df$COMPREHENSIVE.DANGEROUS.DRUGS.ACT.OF.2002....RA.9165

for (isi in 1:nosim){
  
  yt <- drug_data
  nob=length(drug_data)
  
  #****************************************************************New Part
  # *******************Random walk Metropolis
  # Setting up starting values
  om	=	0.2
  alpha	=	0.1
  beta	=	0.1
  con=c(0.2, 0.2)
  delta= 0.1
  time=c(1, 55)
  psi=0.1
  a=3
  b=1
  lam_o	=	NULL
  lam_n	=	NULL
  star=NULL
  star1=NULL
  yt_pred=NULL
  m=NULL
  v=NULL
  res=NULL
  Pt=matrix(0, nrow=length(yt), ncol = 2)
  Ptt=matrix(0, nrow=length(yt), ncol = 2)
  time_star=c()
  lam_o[1]	=	2
  lam_n[1] =	2
  count	=	0
  count1=0
  count2=0
  count3=0
  
  #Hyperparameters
  a_1=2
  a_2=2
  b_1=10 # 6
  b_2=1
  c_1=15 # 15
  c_2=15
  
  ### set stepsize 
  step_om	=	0.2
  step_al	= 	0.1 #0.01
  step_beta  =	0.1
  step_c    =1.3
  step_delta=0.5
  step_psi= 0.05 
  step_time=c(5, 5)
  
  
  M 	= 	20000
  ind	=   	8000
  draws = matrix(0,M, length(paste0("S",1:9)))
  d=matrix(0,M,1)
  resid=matrix(0,M,nob)
  resid_square=matrix(0,M,nob)
  yt_predict=matrix(0,M,nob)
  lower_bound = matrix(0,M,nob)
  upper_bound = matrix(0,M,nob)
  #result  = matrix(0,M-ind,3)
  
  for (i in 1:M){
    
    lik7=0
    lik8=0
    
    lognor1=0
    lognor2=0
    lognor3=0
    lognor4=0
    lognor5=0
    lognor6=0
    
    old=c(om,alpha,beta)
    old1=c(con,delta)
    old2=psi
    
    # 
    if(i  <= ind)
    {
      repeat 
      {
        star[1]	= om +  rnorm(1,0,1)*step_om
        star[2]	=  alpha +  rnorm(1,0,1)*step_al
        star[3]		=  beta  +   rnorm(1,0,1)*step_beta
        if(star[1]>0& star[2]>=0& star[3]>=0&((star[2]+star[3]) <1) ) {break}
      }
    }else
      repeat
      {
        star=rmvnorm(1,IK_mean,IK_cov)
        if(star[1]>0& star[2]>=0& star[3]>=0&((star[2]+star[3]) <1) ) {break}
        
        lognor2=dmvnorm(star,IK_mean,IK_cov,log=T)
        lognor1=dmvnorm(old,IK_mean,IK_cov,log=T)
      }
    
    
    lik1=liko(nob,om,alpha,beta,yt,psi,con,delta,time)
    lik2=liko(nob,star[1],star[2],star[3],yt,psi,con,delta,time)
    
    lik = lognor1-lognor2+lik2- lik1
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      om      	=  star[1]
      alpha 	=  star[2]
      beta    	=  star[3]
      count	= count+1
    }
    
    
    # if(i  <= ind)
    # {
    repeat 
    {
      star2 = psi +rnorm(1,0,1)*step_psi
      if((star2>0 & star2 < 1)){break}
    }
    # }
    # else
    #   repeat
    #   {
    #     star2=rnorm(1,IK_mean2,IK_cov2)
    #     if((star2>1)){break}
    #     
    #     lognor6=dnorm(star2,IK_mean2,IK_cov2,log=T)
    #     lognor5=dnorm(old2,IK_mean2,IK_cov2,log=T)
    #   }
    
    lik3=liko(nob,om,alpha,beta,yt,psi,con,delta,time)
    lik4=liko(nob,om,alpha,beta,yt,star2,con,delta,time)
    
    lik = lognor5-lognor6+ lik4-lik3
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      #con     =  con_star
      psi 	=  star2
      count1	= count1+1
    }
    
    if(i  <= ind)
    {
      repeat 
      {
        star1[1]	= con[1] + rnorm(1,0,1)*step_c
        star1[2]	= con[2] + rnorm(1,0,1)*step_c
        star1[3]	= delta[1] +rnorm(1,0,1)*step_delta
        if(constraint_con_delta(star1)){break}
      }
    }else
      repeat
      {
        star1=rmvnorm(1,IK_mean1,IK_cov1)
        if(constraint_con_delta(star1)){break}
        lognor4=dmvnorm(star1,IK_mean1,IK_cov1,log=T)
        lognor3=dmvnorm(old1,IK_mean1,IK_cov1,log=T)
      }
    
    
    lik5=liko(nob,om,alpha,beta,yt,psi,con,delta,time) + 
      sum(dgamma(con, shape = b_1, rate = b_2, log = TRUE)) + 
      sum(dbeta(delta, shape1 = a_1, shape2 = a_2, log = TRUE))
    lik6=liko(nob,om,alpha,beta,yt,psi,star1[1:2],star1[3],time) + 
      sum(dgamma(star1[1:2], shape = b_1, rate = b_2, log = TRUE)) + 
      sum(dbeta(star1[3], shape1 = a_1, shape2 = a_2, log = TRUE))
    
    
    lik = lognor3-lognor4 +lik6-lik5
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      con     =  star1[1:2]
      delta 	=  star1[3]
      count2	= count2+1
    }
    
    repeat
    { 
      
      time_star[1]= round(time[1]+rnorm(1,0,1)*step_time[1])
      time_star[2]= round((time_star[1]+50)+rnorm(1,0,1)*step_time[2])
      if(constraint_time_star(time_star)){break} 
    }
    
    lik9=liko(nob,om,alpha,beta,yt,psi,con,delta,time)
    lik10=liko(nob,om,alpha,beta,yt,psi,con,delta,time_star)
    
    lik = lik10-lik9
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      time[1]     =  time_star[1]
      time[2]     =  time_star[2]
      count3	= count3+1
    }
    
    draws[i,] = c(om ,alpha,beta,psi,con,delta, time)
    
    
    # if(isi == 1){
    #### In-sample 
    if(i  > ind){
      for (t in 2:nob)
      {
        Pt[,1]=c(rep(0,round(time[1])-1),1,rep(1,nob-round(time[1])))
        Pt[,2]=c(rep(0,round(time[2])-1),1,rep(0,nob-round(time[2])))
        
        if(t<time[1]){
          lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
        }else if(t>=time[1] & t<time[2]){
          lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con[1]*((delta[1])^(t-time[1]))*Pt[t,1] # 2B
        }else{
          lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con[2]*Pt[t,2] # 2A
        }
        
        m[t]=lam_o[t]/(1-psi)
        v[t]=m[t]/(1-psi)^2
        res[t]=(yt[t]-m[t])/sqrt(v[t])
        yt_pred[t] = rgenpois_custom(1, lambda = lam_o[t], psi = psi)
      }
      
      resid[i,]=res
      resid_square[i,]=res^2
      yt_predict[i,]=yt_pred
      
      lower_bound[i, ] <- mean(quantile(yt_pred[-1], probs = 0.025))
      upper_bound[i, ] <- mean(quantile(yt_pred[-1], probs = 0.975))
    }
    # }
    
    
    
    
    # for DIC
    lik11 <- liko(nob,om,alpha,beta,yt,psi,con,delta,time)
    
    d[i,]=c(lik11)
    
    
    
    # print(i)
    if(i%%2000 == 0)
    { cat("************  SIMU AT: ",isi,"\n")
      cat(i,"\n")
      cat("-------------------------------------------------------------\n")
      cat("om, alpha,beta",om, alpha,beta,"\n")         
      cat("accept. rate",100*count /i,"\n")
      cat("-------------------------------------------------------------\n")
      cat("psi",psi, "\n")         
      cat("accept. rate",100*count1 /i,"\n")
      cat("-------------------------------------------------------------\n")
      cat("c ",con, "\n")   
      cat("delta",delta, "\n") 
      cat("accept. rate",100*count2 /i,"\n")
      cat("-------------------------------------------------------------\n")
      cat("T1-T2",time, "\n")         
      cat("accept. rate",100*count3 /i,"\n")
      cat("-------------------------------------------------------------\n")
    }
    if(i==ind)
    {
      IK_mean=c(mean(draws[1001:ind,1]),mean(draws[1001:ind,2]),mean(draws[1001:ind,3]))
      IK_cov=cov(draws[1001:ind,1:3])
      IK_mean2=mean(draws[1001:ind,4]) # r
      IK_cov2=var(draws[1001:ind,4]) # r
      IK_mean1=colMeans(draws[1001:ind,5:7]) # con & delta
      IK_cov1=cov(draws[1001:ind,5:7]) # con & delta
      
      sum1=count/i # ind
      sum2=log(sum1/(1-sum1))-log(0.4/0.6)
      step_om=exp(log(step_om)+sum2)
      step_al=exp(log(step_al)+sum2)
      step_beta=exp(log(step_beta)+sum2)
      
      sum3=count2/i
      sum4=log(sum3/(1-sum3))-log(0.4/0.6)
      step_c=exp(log(step_c)+sum4)
      
      
      sum5=count1/i
      sum6=log(sum5/(1-sum5))-log(0.4/0.6)
      step_psi=exp(log(step_psi)+sum6)
      
      # sum7=count3/i
      # sum8=log(sum7/(1-sum7))-log(0.4/0.6)
      # step_time=exp(log(step_time)+sum8)
    }
    
  }
  
  
  ############################# DIC
  
  MCMC=(ind+1):M
  om_mean=mean(draws[MCMC,1])
  alpha_mean=mean(draws[MCMC,2])
  beta_mean=mean(draws[MCMC,3])
  psi_mean=mean(draws[MCMC,4])
  con_mean=colMeans(draws[MCMC,5:6])
  delta_mean=mean(draws[MCMC,7])
  time_mode=c(getmode(draws[MCMC,8]), getmode(draws[MCMC,9]))
  
  liko12 <- liko(nob,om_mean,alpha_mean,beta_mean,yt,psi_mean,con_mean,delta_mean,time_mode)
  
  DIC[isi,1]=2*(mean(d[MCMC,1]))-liko12
  
  ############################# PLOT
  
  names		 = c(list(bquote(alpha[0])),list(bquote(alpha[1])),list(expression(beta)), 
              list(expression(psi)),list(bquote(omega[0][2])),list(bquote(omega[0][1])), 
              list(bquote(delta[1])))
  if (isi==1)
  {
    par(mai=c(0.6,0.4,0.5,0.4),mfrow=c(4,4))
    for (i in 1:length(names)){
      ts.plot(draws[8001:M,i],xlab="iterations",ylab="",main=names[i])
      # abline(v=0,col=4)
      #abline(h=true[i],col=2,lwd=2)
      acf(draws[8001:M,i],main=names[i])
      # hist(draws[8001:M,i],prob=T,main="",xlab="")
      #abline(v=true[i],col=2,lwd=2)
    }
  }
  
  ######################################################################
  #############################
  MCMC=(ind+1):M
  S1[isi,1]=mean(draws[MCMC,1])
  S1[isi,2]=median(draws[MCMC,1])
  S1[isi,3]=sd(draws[MCMC,1])
  S1[isi,4]=quantile(draws[MCMC,1],0.025)
  S1[isi,5]=quantile(draws[MCMC,1],0.975)
  S2[isi,1]=mean(draws[MCMC,2])
  S2[isi,2]=median(draws[MCMC,2])
  S2[isi,3]=sd(draws[MCMC,2])
  S2[isi,4]=quantile(draws[MCMC,2],0.025)
  S2[isi,5]=quantile(draws[MCMC,2],0.975)
  S3[isi,1]=mean(draws[MCMC,3])
  S3[isi,2]=median(draws[MCMC,3])
  S3[isi,3]=sd(draws[MCMC,3])
  S3[isi,4]=quantile(draws[MCMC,3],0.025)
  S3[isi,5]=quantile(draws[MCMC,3],0.975)
  S4[isi,1]=mean(draws[MCMC,4])
  S4[isi,2]=median(draws[MCMC,4])
  S4[isi,3]=sd(draws[MCMC,4])
  S4[isi,4]=quantile(draws[MCMC,4],0.025)
  S4[isi,5]=quantile(draws[MCMC,4],0.975)
  S5[isi,1]=mean(draws[MCMC,5])
  S5[isi,2]=median(draws[MCMC,5])
  S5[isi,3]=sd(draws[MCMC,5])
  S5[isi,4]=quantile(draws[MCMC,5],0.025)
  S5[isi,5]=quantile(draws[MCMC,5],0.975)
  S6[isi,1]=mean(draws[MCMC,6])
  S6[isi,2]=median(draws[MCMC,6])
  S6[isi,3]=sd(draws[MCMC,6])
  S6[isi,4]=quantile(draws[MCMC,6],0.025)
  S6[isi,5]=quantile(draws[MCMC,6],0.975)
  S7[isi,1]=mean(draws[MCMC,7])
  S7[isi,2]=median(draws[MCMC,7])
  S7[isi,3]=sd(draws[MCMC,7])
  S7[isi,4]=quantile(draws[MCMC,7],0.025)
  S7[isi,5]=quantile(draws[MCMC,7],0.975)
  S8[isi,1]=getmode(draws[MCMC,8])
  S9[isi,1]=getmode(draws[MCMC,9])
  Mean_resid=apply(resid[MCMC,-1:-4],2,mean)
  Mean_resid_square=apply(resid_square[MCMC,-1:-4],2,mean)
  Mean_yt_predict=apply(yt_predict[MCMC,-1:-4],2,mean)
  mean_lower <- apply(yt_predict[MCMC,-1:-4], 2, quantile, probs = 0.025)
  mean_upper <- apply(yt_predict[MCMC,-1:-4], 2, quantile, probs = 0.975)
  # MCMCfile<-cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9)
  # colnames(MCMCfile)<-c("mean_omega","median_omega","std_omega","P025_omega","P975_omega",
  #                       "mean_alpha","median_alpha","std_alpha","P025_alpha","P975_alpha",
  #                       "mean_beta","median_beta","std_beta","P025_beta","P975_beta",
  #                       "mean_r","median_r","std_r","P025_r","P975_r",
  #                       "mean_con1","median_con1","std_con1","P025_con1","P975_con1",
  #                       "mean_con2","median_con2","std_con2","P025_con2","P975_con2",
  #                       "mean_delta","median_delta","std_delta","P025_delta","P975_delta",
  #                       "mode_T1","NAN","NAN","NAN","NAN","mode_T2","NAN","NAN","NAN","NAN")
  # write.csv(MCMCfile,"D:/NB2A2B_1.csv")
  
  
}
{
  r1=apply(S1,2,mean)
  r2=apply(S2,2,mean)
  r3=apply(S3,2,mean)
  r4=apply(S4,2,mean)
  r5=apply(S5,2,mean)
  r6=apply(S6,2,mean)
  r7=apply(S7,2,mean)
  r8=apply(S8,2,getmode)
  r9=apply(S9,2,getmode)
  
  result = round(rbind(r1, r2, r3, r4, r5, r6, r7), 4)
  colnames(result) <- c("mean", "median", "std", "P025", "P975")
  rownames(result) <- c("omega", "alpha", "beta", "psi", "c1", "c2", "delta")
  result
  result1 = round(rbind(r8, r9), 4)
  colnames(result1) <- c("mean", "median", "std", "P025", "P975")
  rownames(result1) <- c("T1", "T2")
  print(result)
  print(result1)
  
  D=apply(DIC,2,mean)
  print(D)
  c(acceptrate=count/M,acceptrate1=count1/M,acceptrate2=count2/M,acceptrate3=count3/M)
  print(c(acceptrate = count/M, acceptrate1 = count1/M, acceptrate2 = count2/M, acceptrate3 = count3/M))
  mssr <- sum(resid_square[MCMC,-1:-4])/((nob-1-4)*nrow(resid_square[MCMC,-1:-4])) 
  print(mssr)
  
  Geweke.Diagnostic(draws[MCMC,1:length(names)]) %>%  pnorm(., 0, 1) %>% round(4)
  sapply(1:length(names), function(j) {
    ess <- effectiveSize(draws[MCMC, j])
    ineff <- length(draws[MCMC, j]) / ess
    return(ineff)
  }) %>% round(4)
  print(round(pnorm(Geweke.Diagnostic(draws[MCMC,1:7])), 4))
  
  ineff_factors <- sapply(1:7, function(j) {
    ess <- effectiveSize(draws[MCMC, j])
    ineff <- length(draws[MCMC, j]) / ess
    return(ineff)
  })
  print(round(ineff_factors, 4))
} 
  
  
  
  # Traceplot and ACF plot
  png("traceplot_acf_GP_2B_2A.png", width = 18, height = 10, units = "in", res = 300)
  par(mai=c(0.3,0.3,0.5,0.3),mfrow=c(4,4))
  for (i in 1:length(names)){
    par(cex.main = 2, cex.axis = 1.5)
    ts.plot(draws[8001:M,i],xlab="iterations",ylab="",main=names[i])
    acf(draws[8001:M,i],main=names[i])
  }
  dev.off()
  
  
  newdata <- data.frame(yt[-1:-4], Mean_yt_predict, mean_lower, mean_upper)
  newdata$year=seq(from=1,to=nrow(newdata))
  names(newdata)<- c("Observed","Predicted", "lower", "upper", "year")
  saveRDS(newdata, "in_sample_GP_2B_2A.rds")
  intervention <- data.frame(xint = c(r8[1]-4, r9[1]-4),  # 39 86
                             yint = c("June 2015", "February 2020"),
                             type = "Intervention")
  
  
  
  plot= newdata %>% pivot_longer(cols = c("Observed","Predicted"),names_to="var",values_to="val") %>%
    ggplot(aes(x=year,y=val, color=var,group=var,linetype=var))+
    # geom_line(size=0.75, lineend = "round")+
    geom_line(size=1, lineend = "round")+
    # geom_vline(data = intervention, aes(xintercept = xint, color = type, linetype = type)) +
    # annotate("text", x = intervention$xint, y = c(114, 105, 55, 55), 
    #          label = intervention$yint, hjust = 1.1, size = 5) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.05, colour = NA) +
    geom_vline(data = intervention, aes(xintercept = xint, color = type, linetype = type), size = 1.3, alpha = 1) +
    annotate("text", x = intervention$xint, y = c(9, 13),
             label = intervention$yint, hjust = 1.1, size = 5) +
    scale_color_manual(
      name   = " ",
      values = c("Observed" = "black",
                 "Predicted" = "blue",
                 "Intervention" = "red"),
      breaks = c("Observed", "Predicted", "Intervention"),
      labels = c("Observed", "Predicted", "Intervention line")
    ) +
    scale_linetype_manual(
      name   = " ",
      values = c("Observed" = 2,
                 "Predicted" = 1,
                 "Intervention" = 5),
      breaks = c("Observed", "Predicted", "Intervention"),
      labels = c("Observed", "Predicted", "Intervention line")
    ) +
    theme_bw()+
    theme(panel.grid =element_blank())+
    xlab("Year")+
    ylab("Number of cases")+
    ggtitle("Predicted Drug Cases in Dipolog City")+
    # theme(panel.grid    = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold",margin=margin(0,0,30,0)))+
    theme(axis.text.x  = element_text(size=19,colour = "black"), 
          axis.title.x = element_text(size=19))+
    theme(axis.text.y  = element_text(size=19 ,colour="black"),
          axis.title.y = element_text(size=19))+
    theme(axis.line = element_line(size=0.5, colour = "black"))+
    theme(legend.title =  element_blank())+
    theme(legend.text = element_text(size=19))+
    theme(legend.position = c(0.15,0.85))+
    guides(color=guide_legend(override.aes=list(fill=NA))) + 
    scale_x_continuous(breaks = c(1,  24-4,  48-4, 72-4,  96-4,  120-4),
                       labels = c("2013","2015","2017",  "2019", "2021",  "2023"))
  
  plot

# png("in-sample_GP_2B_2A.png", width = 300, height = 200, unit="mm", res = 300) #350
plot
# dev.off()


# Standardized Residual Plot
# png("residual_plot_GP_2B_2A.png", width = 300, height = 100, units = "mm", res = 300)
par(mai=c(0.3,0.4,0.6,0.3),mfrow=c(1,3))
# layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), widths = c(2, 2))
par(cex.main = 1, cex.main = 1.6, cex.lab = 1.6, cex.axis = 1.6)
plot(ts(Mean_resid),main="Standardized Residual Plot", lwd = 3,  xaxt='n', ylab = "")
axis(1, at = c(1,  24-4,  48-4, 72-4,  96-4,  120-4), labels = c("2013","2015","2017",  "2019", "2021",  "2023"))
abline(h=0,col="red", lwd = 3)
acf(Mean_resid,main="ACF of Standardized Residual", lwd = 3)
acf(Mean_resid_square,main="ACF of Standardized Squared Residual", lwd = 3)
# dev.off()


