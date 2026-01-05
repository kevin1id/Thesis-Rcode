# Simulating the data
set.seed(199999)
#install.packages("mvtnorm")
library(mvtnorm)
library(HMMpa)
nosim  = 1
result = NULL
S1	 = matrix(0,nosim,5)
S2	 = matrix(0,nosim,5)
S3	 = matrix(0,nosim,5)
S4	 = matrix(0,nosim,5)
S5	 = matrix(0,nosim,5)
S6   = matrix(0,nosim,5)
S7   = matrix(0,nosim,5)
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
r6  = NULL
r7  = NULL

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

liko<-function(nob,yt,om,alpha,beta,psi,con,delta,time){
  liko=0
  time=round(time)
  Pt=c(rep(0,time-1),1,rep(1,nob-time))
  for (t in 2:nob)
  { 
    
    if(t<time){
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
     
    }else{
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt[t]
     
    }
    liko= liko+log((lam_o[t]))+(yt[t]-1)*log((lam_o[t])+psi*yt[t])-(lam_o[t])-psi*yt[t]
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
  psi=0.1
  time=39
  lam_o	=	NULL
  lam_n	=	NULL
  star=NULL
  star1=NULL
  yt_pred=NULL
  m=NULL
  v=NULL
  res=NULL
  lam_o[1]	=	2
  lam_n[1] =	2
  count	=	0
  count1=0
  count2=0
  count3=0
  
  ### set stepsize 
  step_om	=	0.15
  step_al	= 	0.15
  step_beta  =	0.00001
  step_c    =10.0
  step_delta =0.1
  step_psi=0.05
 # step_time=1.0
  
  ###Hyperparameters####
  
  a_1=6
  a_2=10
  b_1=7
  b_2=1
  
  
  M 	= 	20000
  ind	=   	8000
  draws = matrix(0,M,6)
  d=matrix(0,M,1)
  resid=matrix(0,M,140)
  resid_square=matrix(0,M,140)
  yt_predict=matrix(0,M,140)
  #result  = matrix(0,M-ind,3)
  
  for (i in 1:M){
   
    lik9=0
    lik10=0
    lognor1=0
    lognor2=0
    lognor3=0
    lognor4=0
    
    old=c(om,alpha,beta)
    old1=c(con,delta)
    
    
    
  
    if(i  <= ind)
    {
      repeat 
      {
        star[1]	= om +  rnorm(1,0,1)*step_om
        star[2]	=  alpha +  rnorm(1,0,1)*step_al
        star[3]		=  beta  +   rnorm(1,0,1)*step_beta
        if ( (star[1]>0)&(star[2] >0 )&( star[3] >0) & (star[2]+ star[3]<1) ) {break}
      }
    }else
      repeat
      {
        star=rmvnorm(1,IK_mean,IK_cov)
        if ( (star[1]>0)&(star[2] >0 )&( star[3] >0) & (star[2]+ star[3]<1) ) {break}
        
        lognor2=dmvnorm(star,IK_mean,IK_cov,log=T)
        lognor1=dmvnorm(old,IK_mean,IK_cov,log=T)
      }
    
    lik1<-liko(nob,yt,om,alpha,beta,psi,con,delta,time)
    lik2<-liko(nob,yt,star[1],star[2],star[3],psi,con,delta,time)
    
    lik = lognor1-lognor2+lik2- lik1
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      om      	=  star[1]
      alpha 	=  star[2]
      beta    	=  star[3]
      count	= count+1
    }
    
  
    
    repeat
    {
      #con_star = con + rnorm(1,0,0.2)*step_c
      psi_star = psi +rnorm(1,0,1)*step_psi
      if((psi_star>0& psi_star<1)){break} 
    }
    
    lik3<-liko(nob,yt,om,alpha,beta,psi,con,delta,time)
    lik4<-liko(nob,yt,om,alpha,beta,psi_star,con,delta,time)
  
    lik = lik4-lik3
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      #con     =  con_star
      psi 	=  psi_star
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
    
    lik5<-liko(nob,yt,om,alpha,beta,psi,con,delta,time)+ (a_1-1)*log(delta)+(a_2-1)*log(1-delta)+(b_1-1)*log(con)-b_2*con
    lik6<-liko(nob,yt,om,alpha,beta,psi,star1[1],star1[2],time)+ (a_1-1)*log(star1[2])+(a_2-1)*log(1-star1[2])+(b_1-1)*log(star1[1])-b_2*star1[1]
    
    lik = lognor3-lognor4+ lik6-lik5
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      con     =  star1[1]
      delta 	=  star1[2]
      count2	= count2+1
    }
    
    
   
    
    
    draws[i,] = c(om ,alpha,beta,psi,con,delta)
    
    Pt_t=c(rep(0,time-1),1,rep(1,nob-time))
    
    for (t in 2:nob)
    { 
      
      if(t<time){
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
        m[t]=lam_o[t]/(1-psi)
        v[t]=lam_o[t]/(1-psi)^3
      }else{
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt_t[t]
        m[t]=lam_o[t]/(1-psi)
        v[t]=lam_o[t]/(1-psi)^3
      }
      res[t]=(yt[t]-m[t])/sqrt(v[t])
    }
    
    resid[i,]=res
    resid_square[i,]=res^2
    
    
    for (t in 2:nob)
    { 
      
      if(t<time){
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
        yt_pred[t]	  = rgenpois(1,lam_o[t],psi)
      }else{
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt_t[t]
        yt_pred[t]	  = rgenpois(1,lam_o[t],psi)
      }
    }
    
    yt_predict[i,]=yt_pred
    
    for (t in 2:nob)
    { 
     
      if(t<time){
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]
        lik9= lik9-2*log((lam_o[t]))-2*(yt[t]-1)*log((lam_o[t])+psi*yt[t])+2*(lam_o[t])+2*psi*yt[t]+2*log(factorial(yt[t]))
      }else{
        lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*((delta)^(t-time))*Pt_t[t]
        lik9= lik9-2*log((lam_o[t]))-2*(yt[t]-1)*log((lam_o[t])+psi*yt[t])+2*(lam_o[t])+2*psi*yt[t]+2*log(factorial(yt[t]))
      }
    }
    
    d[i,]=c(lik9)
    
    
    if(i%%2000 == 0)
    {cat("************  SIMU AT: ",isi,"\n")
      cat(i,"\n")
      cat("om, alpha,beta",om, alpha,beta,"\n")         
      cat("accept. rate",100*count /i,"\n")
      cat("psi",psi, "\n")         
      cat("accept. rate",100*count1 /i,"\n")
      cat("c","delta",con,delta, "\n")         
      cat("accept. rate",100*count2 /i,"\n")
      cat("Time",time, "\n")         
      cat("accept. rate",100*count3 /i,"\n")
     
    }
    if(i==ind)
    {
      IK_mean=c(mean(draws[1001:ind,1]),mean(draws[1001:ind,2]),mean(draws[1001:ind,3]))
      IK_cov=cov(draws[1001:ind,1:3])
      IK_mean1=c(mean(draws[1001:ind,5]),mean(draws[1001:ind,6]))
      IK_cov1=cov(draws[1001:ind,5:6])
    }
  }
  
  MCMC=(ind+1):M
  om_mean=mean(draws[MCMC,1])
  alpha_mean=mean(draws[MCMC,2])
  beta_mean=mean(draws[MCMC,3])
  psi_mean=mean(draws[MCMC,4])
  con_mean=mean(draws[MCMC,5])
  delta_mean=mean(draws[MCMC,6])
#  time_mode=getmode(draws[MCMC,7])
  
  for (t in 2:nob)
  { 
    Pttt=c(rep(0,time-1),1,rep(1,nob-time))
    if(t<time){
      lam_o[t] = om_mean+ alpha_mean*yt[t-1] + beta_mean*lam_o[t-1]
      lik10= lik10-2*log((lam_o[t]))-2*(yt[t]-1)*log((lam_o[t])+psi_mean*yt[t])+2*(lam_o[t])+2*psi_mean*yt[t]+2*log(factorial(yt[t]))
    }else{
      lam_o[t] = om_mean+ alpha_mean*yt[t-1] + beta_mean*lam_o[t-1]+con_mean*((delta_mean)^(t-time))*Pttt[t]
      lik10= lik10-2*log((lam_o[t]))-2*(yt[t]-1)*log((lam_o[t])+psi_mean*yt[t])+2*(lam_o[t])+2*psi_mean*yt[t]+2*log(factorial(yt[t]))
    }
  }
  
  DIC=2*(mean(d[MCMC,1]))-lik10
  
  
  ############################# PLOT
  names		 = c(list(expression(omega)),list(expression(alpha)),list(expression(beta)), list(expression(psi)),"c" ,list(expression(delta)))
  if (isi==1) 
  {
    par(mai=c(0.6,0.4,0.5,0.4),mfrow=c(3,3))
    for (i in 1:6){
      ts.plot(draws[8001:M,i],xlab="iterations",ylab="",main=names[i])
      #abline(v=0,col=4)
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
# S7[isi,1]=getmode(draws[MCMC,7])
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
#acceptrate3=count3/M
acceptrate
acceptrate1
acceptrate2
#acceptrate3
DIC

MSE=sum(Mean_resid_square)/139
MSE

newdata=data.frame(yt[-1],Mean_yt_predict,L_yt_predict,Q_yt_predict)
names(newdata)<- c("ture","pred","lower","upper")


#####################################
#############################?B?z????
newdata$year=seq(from=1,to=nrow(newdata))
ds1_ture<- data.frame(newdata$ture,newdata$lower,newdata$upper,newdata$year )        

ds1_pred<- data.frame(newdata$pred,newdata$lower,newdata$upper,newdata$year )            



DS1pred_y_res<- rbind(as.matrix(ds1_ture),as.matrix(ds1_pred))
DS1pred_y_res<-as.data.frame(DS1pred_y_res)
DS1pred_y_res$type=as.factor(rep(c(1,2), each=139))
DS1pred_y_res$value=as.factor(rep(c("True","Pred"), each=139))
names(DS1pred_y_res) <- c("data","lower bound","upper bound","year","type","value")


###########################################HINGARCH-NB
library(ggplot2)

plot=ggplot(DS1pred_y_res, aes(x = year, y = data)) + 
  #geom_ribbon(aes(ymin = DS1pred_y_res$`lower bound`, ymax = DS1pred_y_res$`upper bound`), fill = "grey87")+
  geom_line(aes(color = value, linetype = value),size=0.8) + 
  scale_color_manual(values = c("red", "blue")) +       #yt  
  scale_x_continuous(breaks = c(1,12,24,36,48,60,72,84,96,108,120,132,140),
                     labels = c("2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","2023","2024","2025"))+
  theme_bw()+
  theme(panel.grid =element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Illegal Drugs")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x  = element_text(size=15,colour = "black"))+
  theme(axis.text.y  = element_text(size=15 ,colour="black"))+
  theme(axis.line = element_line(size=0.5, colour = "black"))+
  theme(legend.title =  element_blank()) +
geom_vline(xintercept = time, linetype="dashed", color = "black", size=1) +
  annotate("text", x = time, y = max(DS1pred_y_res$data, na.rm = TRUE), 
           label = paste("Time =", time), vjust = -0.5, size=5, color="black")





#jpeg("data30.jpeg", width = 35, height = 18.5, units = 'cm', res = 300)

plot
