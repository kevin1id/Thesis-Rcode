# Simulating the data
set.seed(199999)
#install.packages("mvtnorm")
library(mvtnorm)
nosim  = 1
result = NULL
S1	 = matrix(0,nosim,5)
S2	 = matrix(0,nosim,5)
S3	 = matrix(0,nosim,5)
S4	 = matrix(0,nosim,5)
S5	 = matrix(0,nosim,5)
S6	 = matrix(0,nosim,5)
S7	 = matrix(0,nosim,5)
Mean_resid=matrix(0,1,140)
Mean_resid_square=matrix(0,1,140)
Mean_yt_predict=matrix(0,1,140)
Median_yt_predict=matrix(0,1,140)
Q_yt_predict=matrix(0,1,140)
L_yt_predict=matrix(0,1,140)
r1	 = NULL
r2	 = NULL
r3     	= NULL
r4     	= NULL
r5     	= NULL
r6     	= NULL
r7     	= NULL

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

liko<-function(nob,yt,om,alpha,beta,r,con,delta,time){
  liko=0
  time=round(time)
  Pt=c(rep(0,time-1),1,rep(1,nob-time))
  for (t in 2:nob)
  { 
    if(t<time){
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
      liko= liko+lgamma(yt[t]+r)-lgamma(r) -lgamma(yt[t]+1)- r*log(1+lam_o[t]) + yt[t]*log(lam_o[t])-yt[t]*log(1+lam_o[t])  
    }else{
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt[t]
      liko= liko+lgamma(yt[t]+r)-lgamma(r) -lgamma(yt[t]+1)- r*log(1+lam_o[t]) + yt[t]*log(lam_o[t])-yt[t]*log(1+lam_o[t])   
    }
  }
  return(liko)
}



data <- read.csv("C:/Users/ACER/Desktop/Dipolog crimes.csv", header=TRUE)

for (isi in 1:nosim)
{
  yt<-data$COMPREHENSIVE.DANGEROUS.DRUGS.ACT.OF.2002....RA.9165
  
  nob=140
  

  
  #****************************************************************New Part
  # *******************Random walk Metropolis
  # Setting up starting values
  om	=	0.2
  alpha	=	0.1
  beta	=	0.1
  con=0.2
  delta=0.1
  time=39
  r=1.1
  a=9.5
  b=1
  lam_o	=	NULL
  lam_n	=	NULL
  star=NULL
  star1=NULL
  yt_pred=NULL
  m=NULL
  v=NULL
  res=NULL
  res_square=NULL
  lam_o[1]	=	2
  lam_n[1] =	2
  count	=	0
  count1=0
  count2=0
  count3=0
  
  #Hyperparameters
  a_1=6
  a_2=10
  b_1=7
  b_2=1
  c_1=7
  c_2=1
  
  ### set stepsize 
  step_om	=	0.001
  step_al	= 	0.01
  step_beta  =	0.01
  step_c    =1.0
  step_delta=0.1
  step_r=0.025
  #step_time=0.9
  
  
  
  
  M 	= 	20000
  ind	=   	8000
  draws = matrix(0,M,6)
  d=matrix(0,M,1)
  resid=matrix(0,M,140)
  resid_square=matrix(0,M,140)
  yt_predict=matrix(0,M,140)
  #result  = matrix(0,M-ind,3)
  
  for (i in 1:M){
    lik1=0
    lik2=0
    lik3=0
    lik4=0
    lik5=0
    lik6=0
    lik7=0
    lik8=0
    lik9=0
    lik10=0
    lognor1=0
    lognor2=0
    lognor3=0
    lognor4=0
    lognor5=0
    lognor6=0
    
    old=c(om,alpha,beta)
    old1=c(con,delta)
    old2=r
    
 
    
  
    if(i  <= ind)
    {
      repeat 
      {
        star[1]	= om +  rnorm(1,0,1)*step_om
        star[2]	=  alpha +  rnorm(1,0,1)*step_al
        star[3]		=  beta  +   rnorm(1,0,1)*step_beta
        if(star[1]>0& star[2]>0& star[3]>0&( r*(star[2])^2+( r*star[2]+star[3])^2 <1) ) {break}
      }
    }else
      repeat
      {
        star=rmvnorm(1,IK_mean,IK_cov)
        if(star[1]>0& star[2]>0& star[3]>0&( r*(star[2])^2+( r*star[2]+star[3])^2 <1) ) {break}
        
        lognor2=dmvnorm(star,IK_mean,IK_cov,log=T)
        lognor1=dmvnorm(old,IK_mean,IK_cov,log=T)
      }
    
   lik1=liko(nob,yt,om,alpha,beta,r,con,delta,time)
   lik2=liko(nob,yt,star[1],star[2],star[3],r,con,delta,time)
    
    lik = lognor1-lognor2+lik2- lik1
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      om      	=  star[1]
      alpha 	=  star[2]
      beta    	=  star[3]
      count	= count+1
    }
    
   
    lik3=lik3
    
    if(i  <= ind)
    {
      repeat 
      {
        star2	= r +rnorm(1,0,1)*step_r
        if((star2>1)){break}
      }
    }else
      repeat
      {
        star2=rnorm(1,IK_mean2,IK_cov2)
        if((star2>1)){break}
        
        lognor6=dnorm(star2,IK_mean2,IK_cov2,log=T)
        lognor5=dnorm(old2,IK_mean2,IK_cov2,log=T)
      }
    
    
    lik3=liko(nob,yt,om,alpha,beta,r,con,delta,time)+(c_1-1)*log(r)-c_2*r
    lik4=liko(nob,yt,om,alpha,beta,star2,con,delta,time)+ (c_1-1)*log(star2)-c_2*star2
    
    lik = lognor5-lognor6+lik4-lik3
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      r 	=  star2
      count1	= count1+1
    }
    

     
    if(i  <= ind)
    {
      repeat 
      {
        star1[1]	= con + rnorm(1,0,1.0)*step_c
        star1[2]	=  delta +rnorm(1,0,1)*step_delta
        if((star1[1]>0)&(star1[2]> 0 & star1[2]<1)){break} 
      }
    }else
      repeat
      {
        star1=rmvnorm(1,IK_mean1,IK_cov1)
        if((star1[1]>0)&(star1[2]> 0 & star1[2]<1)){break} 
        
        lognor4=dmvnorm(star1,IK_mean1,IK_cov1,log=T)
        lognor3=dmvnorm(old1,IK_mean1,IK_cov1,log=T)
      }
  
    
    
    lik5=liko(nob,yt,om,alpha,beta,r,con,delta,time)+(b_1-1)*log(con)-b_2*con+ (a_1-1)*log(delta)+(a_2-1)*log(1-delta)
    lik6=liko(nob,yt,om,alpha,beta,r,star1[1],star1[2],time)+(b_1-1)*log(star1[1])-b_2*star1[1]+ (a_1-1)*log(star1[2])+(a_2-1)*log(1-star1[2])
    
    
    lik = lognor3-lognor4+ lik6-lik5
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      con     =  star1[1]
      delta 	=  star1[2]
      count2	= count2+1
    }
    
    
    
    draws[i,] = c(om ,alpha,beta,r,con,delta)
    
    Pt_t=c(rep(0,time-1),1,rep(1,nob-time))
    
    for (t in 2:nob)
    { 
      
      if(t<time){
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
        m[t]=r*lam_o[t]
        v[t]=r*(lam_o[t]+1)*(lam_o[t])
      }else{
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt_t[t]
        m[t]=r*lam_o[t]
        v[t]=r*(lam_o[t]+1)*(lam_o[t])
      }
      res[t]=(yt[t]-m[t])/sqrt(v[t])
    }
    
    resid[i,]=res
    resid_square[i,]=res^2
    
    for (t in 2:nob)
    { 
      
      if(t<time){
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
        yt_pred[t]	  = rnbinom(1, mu = ceiling(r)*(lam_o[t]), size = ceiling(r))
      }else{
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt_t[t]
        yt_pred[t]	  = rnbinom(1, mu = ceiling(r)*(lam_o[t]), size = ceiling(r))
      }
    }
    
    yt_predict[i,]=yt_pred
    
    for (t in 2:nob)
    { 
  
      if(t<time){
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
        lik7= lik7-2*lgamma(yt[t]+r)+2*lgamma(r) +2*lgamma(yt[t]+1)+2*r*log(1+lam_o[t]) -2*yt[t]*log(lam_o[t])+2*yt[t]*log(1+lam_o[t])   
      }else{
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt_t[t]
        lik7= lik7-2*lgamma(yt[t]+r)+2*lgamma(r) +2*lgamma(yt[t]+1)+2*r*log(1+lam_o[t]) -2*yt[t]*log(lam_o[t])+2*yt[t]*log(1+lam_o[t])  
      }
    }
    
    d[i,]=c(lik7)
    
    if(i%%2000 == 0)
    {cat("************  SIMU AT: ",isi,"\n")
      cat(i,"\n")
      cat("om, alpha,beta",om, alpha,beta,"\n")         
      cat("accept. rate",100*count /i,"\n")
      cat("r",r, "\n")         
      cat("accept. rate",100*count1 /i,"\n")
      cat("c","delta",con,delta, "\n")         
      cat("accept. rate",100*count2 /i,"\n")
  #    cat("T",time, "\n")         
   #   cat("accept. rate",100*count3 /i,"\n")
    }
    if(i==ind)
    {
      IK_mean=c(mean(draws[1001:ind,1]),mean(draws[1001:ind,2]),mean(draws[1001:ind,3]))
      IK_cov=cov(draws[1001:ind,1:3])
      IK_mean1=c(mean(draws[1001:ind,5]),mean(draws[1001:ind,6]))
      IK_cov1=cov(draws[1001:ind,5:6])
      IK_mean2=mean(draws[1001:ind,4])
      IK_cov2=var(draws[1001:ind,4])
      
      sum2=count2/i
      sum1=log(sum2/(1-sum2))-log(0.4/0.6)
      
      step_c=exp(log(step_c)+sum1)
      
    }
  }
  
  MCMC=(ind+1):M
  om_mean=mean(draws[MCMC,1])
  alpha_mean=mean(draws[MCMC,2])
  beta_mean=mean(draws[MCMC,3])
  r_mean=mean(draws[MCMC,4])
  con_mean=mean(draws[MCMC,5])
  delta_mean=mean(draws[MCMC,6])
  #time_mode=getmode(draws[MCMC,7])
  
  for (t in 2:nob)
  { 
    Pttt=c(rep(0,time-1),1,rep(1,nob-time))
    if(t<time){
      lam_o[t] = om_mean+ alpha_mean*yt[t-1] + beta_mean*lam_o[t-1]
      lik8= lik8-2*lgamma(yt[t]+r_mean)+2*lgamma(r_mean) +2*lgamma(yt[t]+1)+2*r_mean*log(1+lam_o[t]) -2*yt[t]*log(lam_o[t])+2*yt[t]*log(1+lam_o[t])   
    }else{
      lam_o[t] = om_mean+ alpha_mean*yt[t-1] + beta_mean*lam_o[t-1]+con_mean*((delta_mean)^(t-time))*Pttt[t]
      lik8= lik8-2*lgamma(yt[t]+r_mean)+2*lgamma(r_mean) +2*lgamma(yt[t]+1)+2*r_mean*log(1+lam_o[t]) -2*yt[t]*log(lam_o[t])+2*yt[t]*log(1+lam_o[t])   
    }
  }
  DIC=2*(mean(d[MCMC,1]))-lik8
  
  ############################# PLOT
  names		 = c(list(expression(omega)),list(expression(alpha)),list(expression(beta)), "r","c" ,list(expression(delta)))
  if (isi==1) 
  {
    par(mai=c(0.6,0.4,0.5,0.4),mfrow=c(3,3))
    for (i in 1:6){
      ts.plot(draws[8001:M,i],xlab="iterations",ylab="",main=names[i])
     # abline(v=0,col=4)
      #abline(h=true[i],col=2,lwd=2)
      acf(draws[8001:M,i],main=names[i])
      hist(draws[8001:M,i],prob=T,main="",xlab="")
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
 #S7[isi,1]=getmode(draws[MCMC,7])
 Mean_resid=apply(resid[MCMC,-1],2,mean)
 Mean_resid_square=apply(resid_square[MCMC,-1],2,mean)
 Mean_yt_predict=apply(yt_predict[MCMC,-1],2,mean)
 Median_yt_predict=apply(yt_predict[MCMC,-1],2,median)
 Q_yt_predict=apply(yt_predict[MCMC,-1],2,quantile,probs=c(0.975))
 L_yt_predict=apply(yt_predict[MCMC,-1],2,quantile,probs=c(0.025))
}
#r0=true
r1=apply(S1,2,mean)
r2=apply(S2,2,mean)
r3=apply(S3,2,mean)
r4=apply(S4,2,mean)
r5=apply(S5,2,mean)
r6=apply(S6,2,mean)
#r7=apply(S7,2,mean)
result=round(rbind(r1,r2,r3,r4,r5,r6),4)
colnames(result) <- c("mean","median","std","P025","P975")
rownames(result) <- names
result
acceptrate=count/M
acceptrate1=count1/M
acceptrate2=count2/M
acceptrate3=count3/M

acceptrate
acceptrate1
acceptrate2
acceptrate3
DIC

MSE=sum(Mean_resid_square)/139
MSE

#jpeg("C:/Users/user/Documents/Dissertation/transfer models final-final/Related_files/Data Examples/Data1/Unknown/new/DS1_resid_unknown.jpeg", height=100,width = 250, units = 'mm', res = 300)
par(mai=c(0.8,0.4,0.6,0.4),mfrow=c(1,2))
plot(ts(Mean_resid),main="Standardized Residual Plot")
abline(h=0,col="red")
#hist(Mean_resid,main="Histogram of Standardized Residual")
acf(Mean_resid,main="ACF of Standardized Residual")
#acf(Mean_resid_square,main="ACF of Standardized Squared Residual")

#dev.off()



library(tidyr)
library(dplyr)
library(ggplot2)
newdata=data.frame(yt[-1],Mean_yt_predict,L_yt_predict,Q_yt_predict)
names(newdata)<- c("ture","pred","lower","upper")
newdata$year=seq(from=1,to=nrow(newdata))
ds1_ture<- data.frame(newdata$ture,newdata$lower,newdata$upper,newdata$year )        

ds1_pred<- data.frame(newdata$pred,newdata$lower,newdata$upper,newdata$year )            



DS1pred_y_res<- rbind(as.matrix(ds1_ture),as.matrix(ds1_pred))
DS1pred_y_res<-as.data.frame(DS1pred_y_res)
DS1pred_y_res$type=as.factor(rep(c(1,2), each=139))
DS1pred_y_res$value=as.factor(rep(c("True","Pred"), each=139))
names(DS1pred_y_res) <- c("data","lower bound","upper bound","year","type","value")


newdata<-data.frame(yt[-1],Mean_yt_predict)
newdata$year=seq(from=1,to=nrow(newdata))
names(newdata)<- c("Observed","Predicted","year")

plot= newdata %>% pivot_longer(cols = 'Observed':'Predicted',names_to="var",values_to="val") %>%
  ggplot(aes(x=year,y=val, color=var,group=var,linetype=var))+
  geom_line(size=0.75)+
  scale_color_manual(name="",values=c("blue","red"))+
  scale_linetype_manual(name="",values = c(2,1))+
  theme_bw()+
  theme(panel.grid =element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Illegal Drugs")+
  theme(plot.title = element_text(hjust = 0.5,size = 15, face="bold",margin=margin(0,0,30,0)))+
  theme(axis.text.x  = element_text(size=10,colour = "black"))+
  theme(axis.text.y  = element_text(size=10 ,colour="black"))+
  theme(axis.line = element_line(size=0.5, colour = "black"))+
  theme(legend.title =  element_blank())+
  theme(legend.text = element_text(size=16))+
  theme(legend.position = c(0.9,0.9))+
  scale_x_continuous(breaks = c(1, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144),
                     labels = c("2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","2023","2024","2025"))


###########################################HINGARCH-NB

#jpeg("C:/Users/user/Documents/Dissertation/transfer models final-final/Related_files/Data Examples/Data1/Unknown/new/DS1_pred_unknown.jpeg", height=200,width = 350, units = 'mm', res = 300)


plot+geom_vline(xintercept = time, linetype="dashed", color = "black",size=1)+
  annotate("text", x = time, y = max(DS1pred_y_res$data, na.rm = TRUE), 
           label = paste("Time =", time), vjust = -0.5, size=5, color="black")

#dev.off()