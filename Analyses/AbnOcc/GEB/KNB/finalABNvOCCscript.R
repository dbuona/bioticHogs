rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
graphics.off()

#load libraries
library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(tidybayes)
  setwd("~/Documents/git/bioticHogs/Analyses/AbnOcc/GEB/KNB/")# setworkinf directory

df<-read.csv("ABNvOCC_data.csv") ## read in environmentally matched data

####calcutate how many points are in each graphical quadrant###
df$quad<-NA
df$quad[which(df$H.A>0 &df$H.O>0)]<-1
df$quad[which(df$H.A<=0&df$H.O>0)]<-4  
df$quad[which(df$H.A<0&df$H.O<0)]<-3 
df$quad[which(df$H.A>0&df$H.O<=0)]<-2
df<-filter(df,!is.na(quad))
############################################################

###calculate descriptive statistics
df$dif<-abs(df$H.O-df$H.A) # absolute value difference between occurence-based homogenization index and abundance-based
  round(mean(df$dif,na.rm=TRUE),3) #mean 
  round(sd(df$dif,na.rm=TRUE),3) #sd

df$number<-ifelse(df$dif>=.5,1,0) # calculate how many rows have > 50% difference between abn and occ based homogenization index
  round(table(df$number)[2]/(table(df$number)[1]+table(df$number)[2]),2) # convert this number to percentage of data


cor.test(df$H.A,df$H.O ,method="pearson") #calculate pearson's correlation coeficient

countit<-df %>%group_by(quad) %>% count() #count how many points are in each graphical quadrant
  countit$perc<-round(countit$n/ sum(countit$n),3) ## convert this number to percentage of data

#### compare to spatially derive matches#########
##same workflow as above, so follow comments in previous section  
df2<-read.csv("ABNvOCC_spatial.csv") # read in sptaitally matched data
  
df2$quad<-NA
df2$quad[which(df2$H.A>0 &df2$H.O>0)]<-1
df2$quad[which(df2$H.A<=0&df2$H.O>0)]<-4  
df2$quad[which(df2$H.A<0&df2$H.O<0)]<-3 
df2$quad[which(df2$H.A>0&df2$H.O<=0)]<-2
df2<-filter(df2,!is.na(quad))

df2$dif<-abs(df2$H.O-df2$H.A)
round(mean(df2$dif,na.rm=TRUE),3)
round(sd(df2$dif,na.rm=TRUE),3)
df2$number<-ifelse(df2$dif>=.5,1,0)
round(table(df2$number)[2]/(table(df2$number)[1]+table(df2$number)[2]),2)


cor.test(df2$H.A,df2$H.O ,method="pearson")

countit2<-df2 %>%group_by(quad) %>% count()    
countit2$perc<-round(countit2$n/ sum(countit2$n),3)

####Make FIGURE 2###### 
heat.time<-df
  heat.time$H.A<-round(heat.time$H.A,1) #round homogenization values
  heat.time$H.O<-round(heat.time$H.O,1) #round homogenization values
  heat.test<-heat.time %>% group_by(H.O,H.A) %>% count() # count how many point are in each gridcell of heatmap
  heat.test$frequency<-heat.test$n/sum(heat.test$n) #convert to %
  heat.test$frequency<-round(heat.test$frequency,2) # round to useful number
  
  #Plot FIGURE 2
  ggplot()+
    geom_point(data=df,aes(x=H.O,y=H.A),size=0.01)+
    geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=frequency),alpha=.8)+
    ggthemes::theme_few(base_size = 14)+
    scale_fill_distiller(palette = "Blues",direction=1 )+
    geom_hline(yintercept=0,color="darkgray")+geom_vline(xintercept=0,color="darkgray")+
    ylab("abundance")+xlab("occurence")+
    #geom_text(data=heat.test,aes(x=H.O,y=H.A,label=frequency3),size=2,color="navy")+
    #geom_smooth(method="lm",color="yellow")
    geom_abline(slope = 1,intercept = 0,color="black",linetype="dotdash",size=1)

  
  
  
  
  
  
###################################################################################
 ###Part III: Bromus tectorum case study
bf<-read.csv("BROMUS_data.csv") ## read in bromus data
  bfsmall<-sample_frac(bf,.1,replace = FALSE) # sample 10% of data for computational tractability
  
###Figure 3a ######
ggplot()+
    geom_point(data=bfsmall,aes(x=H.O,y=H.A,color=domers),size=.1)+
    #geom_tile(data=heat.test,aes(x=H.O,y=H.A,fill=frequency),alpha=.8)+
    ggthemes::theme_few()+scale_color_manual(values=(c("lightblue1","navy")))+
    #scale_fill_distiller(palette = "PuRd",direction=1 )+
    geom_hline(yintercept=0,color="darkgray")+geom_vline(xintercept=0,color="darkgray")+
    ylab("abundance")+xlab("occurence")+
    #geom_text(data=heat.test,aes(x=H.O,y=H.A,label=frequency3),size=2,color="navy")+
    #geom_smooth(method="lm",color="yellow")
    geom_abline(slope = 1,intercept = 0,color="black",linetype="dotdash",size=1)   
#####################################################################################

mod1<-brms::brm(quad~domers*cover_diff+(1|mm(site1,site2)),
                       data=bfsmall,
                       family=categorical(link="logit")) ### model for effect of invasive cover difference in figure 3b
  
mod2<-brms::brm(quad~domers*logdist+(1|mm(site1,site2)),
                       data=bfsmall,
                       family=categorical(link="logit"))### model for effect of distance between invaded plots  in figure 3c
  
###############################################################################################################################
##posterior predictors for plot 3b
d3b<-data.frame(domers=rep(c("mono-specific","hetero-specific"),3),cover_diff=rep(unique(envgoo$cover_diff),each=2))

  p3b<-epred_draws(mo1,newdata = d3b,ndraws=1000,re_formula = ~NULL)

  ggplot(p3b,aes(reorder(cover_diff,-.epred),.epred))+stat_pointinterval(aes(color=domers),.width=c(.5,.95))+
  facet_wrap(~factor(.category, levels=c(2,1,3,4)),ncol=2)+
  ggthemes::theme_few()+scale_color_manual(values=(c("lightblue1","navy")))+
  xlab("invader abundance")+ylab("predicted likelihood")
#############################################################################################      
##posterior predictors for plot 3c  
d3c<-data.frame(domers=rep(c("mono-specific","hetero-specific"),2),logdist=rep(c(0,6),each=2))
  
  p3c-epred_draws(mod2,newdata = d3c,ndraws=1000,re_formula = ~NULL)
  p3c$grouper<-paste( p3c$.draw, p3c$domers)
  
ggplot(p3c,aes(logdist,.epred))+geom_line(aes(color=domers,group=grouper),size=0.004)+
    geom_smooth(aes(color=domers),method="glm")+
    facet_wrap(~factor(.category, levels=c(2,1,3,4)))+
    ggthemes::theme_few()+scale_color_manual(values=(c("lightblue1","navy")))+
  xlab("distance log(km)")+ylab("predicted likelihood")+theme(legend.position = "top")
######################################################################################################  
