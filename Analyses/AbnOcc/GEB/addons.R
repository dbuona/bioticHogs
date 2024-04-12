use.list<-read.csv("Analyses/AbnOcc/GEB/matchedplots.csv")

invaded<-filter(use.list,status==1)
uninvaded<-filter(use.list,status!=1)

invaded<-dplyr::select(invaded,Plot,subclass)
uninvaded<-dplyr::select(uninvaded,Plot,subclass)

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
###add dataset

dset<-dplyr::select(d,Plot,Dataset)
dset<-distinct(dset)

dsetu1<-dset
dsetu2<-dset

dseti1<-dset
dseti2<-dset

colnames(dsetu1)<-c("uninvaded.plot.1","dataset.u.plot.1")
colnames(dsetu2)<-c("uninvaded.plot.2","dataset.u.plot.2")

colnames(dseti1)<-c("invaded.plot.1","dataset.i.plot.1")
colnames(dseti2)<-c("invaded.plot.2","dataset.i.plot.2")

env<-left_join(env,dsetu1)
env<-left_join(env,dsetu2)

env<-left_join(env,dseti1)
env<-left_join(env,dseti2)


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