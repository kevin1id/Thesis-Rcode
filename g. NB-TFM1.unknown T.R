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
S6   = matrix(0,nosim,5)
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

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

liko<-function(nob,yt,om,alpha,beta,r,con,time){
  liko=0
  Pt=c(rep(0,round(time)-1),1,rep(0,nob-round(time)))
  for (t in 2:nob)
  { 
    lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*Pt[t]
    
    liko= liko+lgamma(yt[t]+r)-lgamma(r) -lgamma(yt[t]+1)- r*log(1+lam_o[t]) + yt[t]*log(lam_o[t])-yt[t]*log(1+lam_o[t])
    
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
  time=70
  r=1.1
  lam_o	=	NULL
  lam_n	=	NULL
  star=NULL
  star1=NULL
  yt_pred=NULL
  m=NULL
  v=NULL
  res=NULL
  res_square=NULL
  lam=NULL
  lam[1]=2
  lam_o[1]	=	2
  lam_n[1] =	2
  count	=	0
  count1=0
  count2=0
  count3=0
  
  #Hyperparameters
  b_1=7
  b_2=1
  c_1=7
  c_2=1
  
  ### set stepsize 
  step_om	=	0.001
  step_al	= 	0.01
  step_beta  =	0.05
  step_c    =1.0
  step_r=0.05
  step_time=0.5
  
  
  
  
  M 	= 	20000
  ind	=   	8000
  draws = matrix(0,M,6)
  d=matrix(0,M,1)
  resid=matrix(0,M,140)
  resid_square=matrix(0,M,140)
  yt_predict=matrix(0,M,140)
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
    old1=r
    old2=con
    
 
    
  
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
    
  lik1=liko(nob,yt,om,alpha,beta,r,con,time)
  lik2=liko(nob,yt,star[1],star[2],star[3],r,con,time)
    
    lik = lognor1-lognor2+lik2- lik1
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      om      	=  star[1]
      alpha 	=  star[2]
      beta    	=  star[3]
      count	= count+1
    }
    
   
    
    if(i  <= ind)
    {
      repeat 
      {
        star1	= r +rnorm(1,0,1)*step_r
        if((star1>1)){break}
      }
    }else
      repeat
      {
        star1=rnorm(1,IK_mean1,IK_cov1)
        if((star1>1)){break}
        
        lognor4=dnorm(star1,IK_mean1,IK_cov1,log=T)
        lognor3=dnorm(old1,IK_mean1,IK_cov1,log=T)
      }
    
   
    
    lik3=liko(nob,yt,om,alpha,beta,r,con,time)+(c_1-1)*log(r)-c_2*r
    lik4=liko(nob,yt,om,alpha,beta,star1,con,time)+(c_1-1)*log(star1)-c_2*star1
    
    lik = lognor3-lognor4+lik4-lik3
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      #con     =  con_star
      r 	=  star1
      count1	= count1+1
    }
    
   
  
    
    if(i  <= ind)
    {
      repeat 
      {
        star2	= con +rnorm(1,0,1)*step_c
        if((star2>0)){break}
      }
    }else
      repeat
      {
        star2=rnorm(1,IK_mean2,IK_cov2)
        if((star2>0)){break}
        
        lognor6=dnorm(star2,IK_mean2,IK_cov2,log=T)
        lognor5=dnorm(old2,IK_mean2,IK_cov2,log=T)
      }
    

    
    
    lik5=liko(nob,yt,om,alpha,beta,r,con,time)+(b_1-1)*log(con)-b_2*con
    lik6=liko(nob,yt,om,alpha,beta,r,star2,time)+(b_1-1)*log(star2)-b_2*star2
    
    
    lik = lognor5-lognor6+lik6-lik5
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      con     =  star2
     
      count2	= count2+1
    }
    
  
    
 
    
    
    Pa=round(0.2*nob)
    Pb=round(0.8*nob)
    
    repeat
    {
      time_star= time+rnorm(1,0,3)*step_time
      if(time_star>Pa &time_star<Pb){break} 
    }
    
   
    
    lik9=liko(nob,yt,om,alpha,beta,r,con,time)
    lik10=liko(nob,yt,om,alpha,beta,r,con,time_star)
    
  
    lik = lik10-lik9
    u 	=	 runif(1)
    if (log(u)<lik) 
    {
      time    =  round(time_star)
    
      count3	= count3+1
    }
    
    
    draws[i,] = c(om ,alpha,beta,r,con,time)
    
    Pt_t=c(rep(0,time-1),1,rep(0,nob-time))
    
    for (t in 2:nob)
    { 
      
      lam[t] = om+alpha*yt[t-1]+beta*lam[t-1]+con*Pt_t[t]
      m[t]=r*lam[t]
      v[t]=r*(lam[t]+1)*(lam[t])
      res[t]=(yt[t]-m[t])/sqrt(v[t])
    }
    
    resid[i,]=res
    resid_square[i,]=res^2
    
    
    for (t in 2:nob)
    { 
    
      lam[t] = om+alpha*yt[t-1]+beta*lam[t-1]+con*Pt_t[t]
      yt_pred[t]	  = rnbinom(1, mu = ceiling(r)*(lam[t]), size = ceiling(r))
    }
    
    yt_predict[i,]=yt_pred
    
    
    for (t in 2:nob)
    { 
      lam_o[t] = om+ alpha*yt[t-1] + beta*lam_o[t-1]+con*Pt_t[t]
    
        lik7= lik7-2*lgamma(yt[t]+r)+2*lgamma(r) +2*lgamma(yt[t]+1)+2*r*log(1+lam_o[t]) -2*yt[t]*log(lam_o[t])+2*yt[t]*log(1+lam_o[t])  
      
    }
    
    d[i,]=c(lik7)
    
    if(i%%2000 == 0)
    {cat("************  SIMU AT: ",isi,"\n")
      cat(i,"\n")
      cat("om, alpha,beta",om, alpha,beta,"\n")         
      cat("accept. rate",100*count /i,"\n")
      cat("r",r, "\n")         
      cat("accept. rate",100*count1 /i,"\n")
      cat("c",con, "\n")         
      cat("accept. rate",100*count2 /i,"\n")
      
    }
    if(i==ind)
    {
      IK_mean=c(mean(draws[1001:ind,1]),mean(draws[1001:ind,2]),mean(draws[1001:ind,3]))
      IK_cov=cov(draws[1001:ind,1:3])
      IK_mean1=mean(draws[1001:ind,4])
      IK_cov1=var(draws[1001:ind,4])
      IK_mean2=mean(draws[1001:ind,5])
      IK_cov2=var(draws[1001:ind,5])
      
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
 time_mode=getmode(draws[MCMC,6])
  
  for (t in 2:nob)
  { 
    Pttt=c(rep(0,time_mode-1),1,rep(0,nob-time_mode))
    lam_o[t] = om_mean+ alpha_mean*yt[t-1] + beta_mean*lam_o[t-1]+con_mean*Pttt[t]
  
      lik8= lik8-2*lgamma(yt[t]+r_mean)+2*lgamma(r_mean) +2*lgamma(yt[t]+1)+2*r_mean*log(1+lam_o[t]) -2*yt[t]*log(lam_o[t])+2*yt[t]*log(1+lam_o[t])  
    
  }
  
  DIC=2*(mean(d[MCMC,1]))-lik8
  
  ############################# PLOT
#  jpeg("C:/Users/user/Documents/Dissertation/transfer models final/NB/Data-Examples/Data1-NB/DS1_acf_unknown.jpeg", height=350,width = 500, units = 'mm', res = 300)
  
  names		 = c(list(bquote(alpha[0])),list(bquote(alpha[1])),list(bquote(beta[1])), "r",list(bquote(omega[0])), "T")
  if (isi==1) 
  {
    par(mai=c(0.6,0.4,0.5,0.4),mfrow=c(3,4))
    for (i in 1:5){
      ts.plot(draws[8001:M,i],xlab="iterations",ylab="",main=names[i])
      #abline(v=0,col=4)
      #abline(h=true[i],col=2,lwd=2)
      acf(draws[8001:M,i],main=names[i])
      #hist(draws[5001:M,i],prob=T,main="",xlab="")
      #abline(v=true[i],col=2,lwd=2)
    }
  }
 # dev.off()
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
 S6[isi,1]=getmode(draws[MCMC,6])
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

#Mean_resid
#jpeg("C:\\Users\\user\\Documents\\latex\\DS1_resid.jpeg",height=200,width = 350,units='mm',res=300)
#jpeg("C:/Users/user/Documents/Dissertation/transfer models final/NB/Data-Examples/Data1-NB/DS1_resid_unknown.jpeg", height=100,width = 250, units = 'mm', res = 300)
par(mai=c(0.8,0.4,0.6,0.4),mfrow=c(1,2))
plot(ts(Mean_resid),main="Standardized Residual Plot")
abline(h=0,col="red")
#hist(Mean_resid,main="Histogram of Standardized Residual")
acf(Mean_resid,main="ACF of Standardized Residual")
#acf(Mean_resid_square,main="ACF of Standardized Squared Residual")

#dev.off()

par(mai=c(0.9,0.4,0.5,0.4),mfrow=c(1,1))
x = ts(yt, frequency = 12,start = c(2013,1))
#y = ts(yt2, frequency = 12,start = c(2013,2))
z= ts(Mean_yt_predict, frequency = 12,start = c(2013,2))
ts.plot(x, z, gpars = list(col = c("blue", "red")))
#ts.plot(x, y, gpars = list(col = c("black", "red")))
#plot(ts(yt))
#plot(ts(Mean_yt_predict))
#plot(ts(Median_yt_predict))
#plot(ts(Q_yt_predict))
#Mean_yt_predict

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
# Extract the dynamically computed "Time" value from the result matrix
Time_value <- as.numeric(result["T", "mean"])  # Ensure it's numeric

# Ensure Time_value is within the year sequence range
if (Time_value < min(newdata$year) | Time_value > max(newdata$year)) {
  warning("Time_value is outside the range of the x-axis!")
}
library(ggplot2)
# Update the ggplot to include a vertical line and label dynamically
plot <- ggplot(DS1pred_y_res, aes(x = year, y = data)) + 
  geom_line(aes(color = value, linetype = value), linewidth = 0.8) + 
  scale_color_manual(values = c("red", "blue")) +      
  scale_x_continuous(breaks = c(1,12,24,36,48,60,72,84,96,108,120,132,140),
                     labels = c("2013","2014","2015","2016","2017","2018","2019",
                                "2020","2021","2022","2023","2024","2025")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggtitle("Illegal Drugs") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face="bold", margin=margin(0,0,30,0))) +
  theme(axis.text.x  = element_text(size=15, colour = "black")) +
  theme(axis.text.y  = element_text(size=15, colour="black")) +
  theme(axis.line = element_line(linewidth=0.5, colour="black")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=12)) +
  
  # Add a vertical line at the dynamically computed Time_value
  geom_vline(xintercept = time, linetype="dashed", color="black", linewidth=1) +
  
  # Add a label for the computed Time value at the top of the plot
  annotate("text", x = Time_value, y = max(DS1pred_y_res$data, na.rm=TRUE), 
           label = paste("Time =", round(Time_value)), vjust = -0.5, hjust = 1, 
           size = 6, color = "red", fontface = "bold")
plot



# Save the plot to a file
jpeg("C:/Users/ACER/Desktop/Data New/NB_RW_IK_model.jpeg", height = 200, width = 350, units = 'mm', res = 300)
print(plot)  # Ensure the plot is rendered inside jpeg()
dev.off()  # Close the device to save the image

