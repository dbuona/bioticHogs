

invaded<-filter(use.list,status==1)
uninvaded<-filter(use.list,status!=1)

invaded<-dplyr::select(invaded,Plot,pair)
uninvaded<-dplyr::select(uninvaded,Plot,pair)

uninvaded1<-uninvaded
uninvaded2<-uninvaded  
colnames(uninvaded1)<-c("uninvaded.plot.1","site1")
colnames(uninvaded2)<-c("uninvaded.plot.2","site2") 

invaded1<-invaded
invaded2<-invaded  
colnames(invaded1)<-c("invaded.plot.1","site1")
colnames(invaded2)<-c("invaded.plot.2","site2") 

env<-left_join(env,uninvaded1)
env<-left_join(env,uninvaded2)

env<-left_join(env,invaded1)
env<-left_join(env,invaded2)


##now get dominate invader
d<-d1 %>%
  group_by(Plot) %>%
  filter(Year==max(Year)) %>%
  ungroup() 


d.inv<-filter(d,NativeStatus=="I")

d.inv<-dplyr::select(d.inv,Plot,AcceptedTaxonName,PctCov_100)  

dom<-d.inv %>% group_by(Plot) %>% slice(which.max(PctCov_100))

dom1<-dom
dom2<-dom

colnames(dom1)<-c("invaded.plot.1","dom.invader.1","cov.dom.1")
colnames(dom2)<-c("invaded.plot.2","dom.invader.2","cov.dom.2")

env<-left_join(env,dom1)
env<-left_join(env,dom2)
