#Richness, PD, MPD

#should change to referencing dat6 (it is the same data I just added the MPD data to the dat4 dataframe)

options(contrasts=c("contr.helmert","contr.poly"))
options(contrasts=c("contr.treatment","contr.poly"))
options("contrasts")



##### Abundance of isolates per plant species per site #####

databun<-dat6
temp<-dat6%>%
  dplyr::select(OTU0:OTU9)%>%
  mutate(abun=rowSums(.))%>%
  dplyr::select(abun)
databun<-cbind(databun,temp)

#Per host plant per site, abundance and richness
m1<-databun%>%
  group_by(Site,HostPlant)%>%
  summarise(mean=mean(abun),se=std.error(abun),n=n())
m1

m1<-databun%>%
  group_by(Site,HostPlant)%>%
  summarise(mean=mean(SR),se=std.error(SR),n=n())
m1

#Per site abundance
m1<-databun%>%
  group_by(Site)%>%
  summarise(sum=sum(abun),n=n())
m1

#Per species abundance
m1<-databun%>%
  group_by(HostPlant)%>%
  summarise(sum=sum(abun),n=n())
m1

##### Calculating richness by site #####
dat6turtlecove<-dat6%>%
  filter(Site=="Turtle Cove")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6turtlecove)>0)) # 20 taxa in TC

dat6cerf<-dat6%>%
  filter(Site=="CERF")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6cerf)>0)) # 33 taxa in cerf

dat6lumcon<-dat6%>%
  filter(Site=="LUMCON")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6lumcon)>0)) # 35 taxa in lumcon




##### Calculating richness by plant species #####
dat6junroe<-dat6%>%
  filter(HostPlant=="J. roemerianus")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6junroe)>0)) # 16 taxa in junroe

dat6phraus<-dat6%>%
  filter(HostPlant=="P. australis")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6phraus)>0)) # 30 taxa in phrag

dat6saglan<-dat6%>%
  filter(HostPlant=="S. lancifolia")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6saglan)>0)) # 8 taxa in sag

dat6spapat<-dat6%>%
  filter(HostPlant=="S. patens")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6spapat)>0)) # 20 taxa in patens

dat6spaalt<-dat6%>%
  filter(HostPlant=="S. alterniflora")%>%
  dplyr::select(OTU0:OTU9)
length(which(colSums(dat6spaalt)>0)) # 27 taxa in alterniflora


#abundance
rowSums(dat6[,11:66])



##### Richness #####
richnessmean<-dat6%>%
  group_by(HostPlant,Site)%>%
  summarise(mean=mean(SR),se=std.error(SR),count=n())

ggplot(richnessmean,aes(x=Site,y=mean,color=Site,group=HostPlant))+
  labs(x = "",y="Richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) +
  #scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(HostPlant),strip.position = "bottom")

#This is one I used for ms, but the final fig is below in the same pdf as MPD
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/richness2.pdf",width=6.6,height=3.5)
ggplot(richnessmean,aes(x=HostPlant,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="Richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",strip.placement = "outside",panel.spacing=unit(0,"cm"),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.background = element_rect(fill = NA, color = "gray50"))+
  #scale_color_viridis(discrete = T,option="rocket")+
  #scale_color_manual(values = c("#46F884FF", "#3E9BFEFF","#E1DD37FF","#30123BFF","#F05B12FF"))+
  #scale_color_manual(values = c("#EEBAB4", "#F05039","#3D65A5","#12285c","#7ca1cc"))+
  scale_color_manual(values = scales::viridis_pal(option = "viridis")(5))+
  geom_point(size=2.75)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.35,size=1) +
  facet_wrap(vars(Site),strip.position = "bottom")
#dev.off()

#I based the color palette on https://www.datylon.com/blog/data-visualization-for-colorblind-readers
#tubro 6 is pretty good
#scales::viridis_pal(option = "viridis")(6)

#Averaging within each Site
richnessmean<-dat4%>%
  # ungroup()%>%
  group_by(Site)%>%
  summarise(mean=mean(SR),se=std.error(SR))

ggplot(richnessmean,aes(x=Site,y=mean,color=Site,group=Site))+
  labs(x = "",y="Richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) 

#Averaging within host plant
richnessmean<-dat4%>%
  group_by(HostPlant)%>%
  summarise(mean=mean(SR),se=std.error(SR))

ggplot(richnessmean,aes(x=HostPlant,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="Richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) 


#Stats

# m1<-gls(SR~HostPlant+Site,data=dat6,na.action=na.omit)
# anova(m1,type="marginal")
# 
# m2<-gls(SR~HostPlant*Site,data=phragspartina,na.action=na.omit)
# anova(m2,type="marginal")

m1<-lme(SR~HostPlant+Site,random=~1|Year,data=dat6,na.action=na.omit)
anova(m1,type="marginal")
hist(resid(m1))
boxplot(resid(m1)~dat6$HostPlantSite)

#with heterogeneous variances, not significant
#m1a<-lme(SR~HostPlant+Site,random=~1|Year,weights=varIdent(form=~1|Site*HostPlant),data=dat6,na.action=na.omit)
#anova(m1,m1a)
#anova(m1,type="marginal")


m2<-lme(SR~HostPlant*Site,random=~1|Year,data=phragspartina,na.action=na.omit)
anova(m2,type="marginal")

#with heterogeneous variances, not significant
#m2a<-lme(SR~HostPlant*Site,random=~1|Year,weights=varIdent(form=~1|Site*HostPlant),data=phragspartina,na.action=na.omit)
#anova(m2,m2a)
#anova(m2,type="marginal")






##### MPD #####

## an MPD that is > 0 is overdispersion, and MPD < 0 is phylogenetic clustering
## mpd.obs.z is what we want it compares the observed to the randomized data, if > 0 overdispersion, if < 0 clustered

mpdmean<-dat6%>%
  group_by(Site,HostPlant)%>%
  filter(is.na(mpd.obs.z.weighted)==F)%>%
  summarise(mean=mean(mpd.obs.z.weighted,na.rm=T),se=std.error(mpd.obs.z.weighted),n=n())%>%
  ungroup()%>%
  mutate(HostPlant=recode(HostPlant,"Phragmites australis"="PhrAus","Sagittaria lancifolia"="SagLan","Spartina alterniflora"="SpaAlt","Spartina patens"="SpaPat","Juncus roemerianus"="JunRoe"))
mpdmean

#This is one I used for ms
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/mpd2.pdf",width=6.6,height=3.5)
ggplot(mpdmean,aes(x=HostPlant,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="MPD") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",strip.placement = "outside",panel.spacing=unit(0,"cm"),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.background = element_rect(fill = NA, color = "gray50"))+
  scale_color_manual(values = c("#EEBAB4", "#F05039","#3D65A5","#12285c","#7ca1cc"))+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_point(size=2.75)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.35,size=1) +
  facet_wrap(vars(Site),strip.position = "bottom")
dev.off()

ggplot(mpdmean,aes(x=Site,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="MPD") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) +
  facet_wrap(~HostPlant)#,strip.position = "bottom")

ggplot(mpdmean,aes(x=HostPlant,y=mean,color=HostPlant,group=Site))+
  labs(x = "",y="MPD") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(.5,"cm"),strip.placement = "outside")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) +
  facet_wrap(~Site)#,strip.position = "bottom")

#with year
mpdmean<-dat6%>%
  group_by(Site,HostPlant,Year)%>%
  filter(is.na(mpd.obs.z.weighted)==F)%>%
  summarise(mean=mean(mpd.obs.z.weighted,na.rm=T),se=std.error(mpd.obs.z.weighted),n=n())%>%
  ungroup()%>%
  mutate(HostPlant=recode(HostPlant,"Phragmites australis"="PhrAus","Sagittaria lancifolia"="SagLan","Spartina alterniflora"="SpaAlt","Spartina patens"="SpaPat","Juncus roemerianus"="JunRoe"))
mpdmean

ggplot(mpdmean,aes(x=Site,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="MPD") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
#  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) +
  facet_wrap(~HostPlant)#,strip.position = "bottom")

options(contrasts=c("contr.helmert","contr.poly"))
options(contrasts=c("contr.treatment","contr.poly"))

# m1<-gls(mpd.obs.z.weighted~HostPlant+Site,data=dat6,na.action=na.omit)
# anova(m1,type="marginal")
 
# m2<-gls(mpd.obs.z.weighted~HostPlant*Site,data=phragspartina,na.action=na.omit)
# anova(m2,type="marginal")

m1<-lme(mpd.obs.z.weighted~HostPlant+Site,random=~1|Year,data=dat6,na.action=na.omit)
#m2<-lme(mpd.obs.z.weighted~HostPlant+Site,random=~1|Year,weights=varIdent(form=~1|Site),data=dat6,na.action=na.omit)
#m3<-lme(mpd.obs.z.weighted~HostPlant+Site,random=~1|Year,weights=varIdent(form=~1|HostPlant),data=dat6,na.action=na.omit)
#m4<-lme(mpd.obs.z.weighted~HostPlant+Site,random=~1|Year,weights=varIdent(form=~1|HostPlantSite),data=dat6,na.action=na.omit)
#anova(m1,m2,m3,m4)
anova(m1,type="marginal")
hist(resid(m1))
plot(fitted(m1),resid(m1))


#remove the space in HostPlantSite and make it a factor
#dat6$HostPlantSite<-gsub("\\s+", "", dat6$HostPlantSite)
#dat6$HostPlantSite<-factor(dat6$HostPlantSite)

#For t-tests
#m2<-lme(mpd.obs.z.weighted~HostPlantSite,random=~1|Year,data=dat6,na.action=na.omit)#,weights=varIdent(form=~1|HostPlantSite)
m2<-gls(mpd.obs.z.weighted~HostPlantSite,data=dat6,na.action=na.omit)#,weights=varIdent(form=~1|HostPlantSite)

#this does not work, it gives the wrong means
#summary(glht(m2, linfct = mcp(HostPlantSite=c("Phragmitesaustralis_TurtleCove=0","Sagittarialancifolia_TurtleCove=0","Spartinapatens_TurtleCove=0","Phragmitesaustralis_CERF=0","Spartinaalterniflora_CERF=0","Spartinapatens_CERF=0","Juncusroemerianus_LUMCON=0","Phragmitesaustralis_LUMCON=0","Spartinaalterniflora_LUMCON=0","Spartinapatens_LUMCON=0"))))

m2em<-as.data.frame(summary(emmeans(m2,~HostPlantSite)))
m2.r<-ref_grid(m2)
m2.s<-emmeans(m2.r,"HostPlantSite")
test(m2.s,adjust="fdr")#dunnett (what glht uses) or fdr

m2<-lme(mpd.obs.z.weighted~HostPlant*Site,random=~1|Year,data=phragspartina,na.action=na.omit)
anova(m2,type="marginal")

#tyring heterogeneou variances, not significant
#m2a<-lme(mpd.obs.z.weighted~HostPlant*Site,random=~1|Year,weights=varIdent(form=~1|HostPlant*Site),data=phragspartina,na.action=na.omit)
#anova(m2,m2a)

# m2<-lme(mpd.obs.z.weighted~Site,random=~1|HostPlant,data=dat6,na.action=na.omit)
# anova(m2)


##### Plot for mdp and richness together #####

richfig<-ggplot(richnessmean,aes(x=HostPlant,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="Richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",strip.placement = "outside",panel.spacing=unit(0,"cm"),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.background = element_rect(fill = NA, color = "gray50"))+
  scale_color_manual(values = c("#3D65A5", "#F05039","#EEBAB4","#7ca1cc","#12285c"))+#c("#EEBAB4", "#F05039","#3D65A5","#12285c","#7ca1cc")
  geom_point(size=2.75)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.35,size=1) +
  facet_wrap(vars(Site),strip.position = "bottom")

mpdfig<-ggplot(mpdmean,aes(x=HostPlant,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="MPD") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",strip.placement = "outside",panel.spacing=unit(0,"cm"),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),panel.background = element_rect(fill = NA, color = "gray50"))+
  scale_color_manual(values = c("#3D65A5", "#F05039","#EEBAB4","#7ca1cc","#12285c"))+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_point(size=2.75)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.35,size=1) +
  facet_wrap(vars(Site),strip.position = "bottom")


pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/richnessmpd2.pdf")
plot_grid(richfig, mpdfig, nrow = 2)
dev.off()







#####MPD trials #####
#### below are the trials that I did for individual (abundance weighted vs not and grouped by species x site, site, and species)

#MPD by plant individual 
#notes on results of unweighted vs weighted. the overall pattern is very similar, but the weighted results tend to have slightly smaller error bars so that some (junroe) might be significantly positive. However, it probably makes the most sense to use unweighted because (based on the labeling file in google drive) it sounds like they only isolated one "morphotype" per plate (of 5 root pieces) therefore it is unlikely to get an abundance >2 (there were two plates of 5 roots each for each sample).

phydist <- cophenetic(tree)
ses.mpd.result.ind <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=FALSE, runs=999) #takes 5 min with 999
ses.mpd.result.ind$PlantIndividualYear<-rownames(ses.mpd.result.ind)
ses.mpd.result.ind

mpdmean<-ses.mpd.result.ind%>%
  full_join(dat4)%>%
  group_by(Site,HostPlant)%>%
  filter(is.na(mpd.obs.z)==F)%>%
  summarise(mean=mean(mpd.obs.z,na.rm=T),se=std.error(mpd.obs.z),n=n())%>%
  ungroup()%>%
  mutate(HostPlant=recode(HostPlant,"Phragmites australis"="PhrAus","Sagittaria lancifolia"="SagLan","Spartina alterniflora"="SpaAlt","Spartina patens"="SpaPat","Juncus roemerianus"="JunRoe"))
mpdmean

ses.mpd.result.ind.abun <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=T, runs=999) #takes 5 min with 999
ses.mpd.result.ind.abun$PlantIndividualYear<-rownames(ses.mpd.result.ind.abun)
ses.mpd.result.ind.abun

mpdmean<-ses.mpd.result.ind.abun%>%
  full_join(dat4)%>%
  group_by(Site,HostPlant)%>%
  filter(is.na(mpd.obs.z)==F)%>%
  summarise(mean=mean(mpd.obs.z,na.rm=T),se=std.error(mpd.obs.z),n=n())%>%
  ungroup()%>%
  mutate(HostPlant=recode(HostPlant,"Phragmites australis"="PhrAus","Sagittaria lancifolia"="SagLan","Spartina alterniflora"="SpaAlt","Spartina patens"="SpaPat","Juncus roemerianus"="JunRoe"))
mpdmean

ggplot(mpdmean,aes(x=Site,y=mean,color=HostPlant,group=HostPlant))+
  labs(x = "",y="MPD") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) +
  facet_wrap(~HostPlant)#,strip.position = "bottom")

ggplot(mpdmean,aes(x=HostPlant,y=mean,color=HostPlant,group=Site))+
  labs(x = "",y="MPD") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(.5,"cm"),strip.placement = "outside")+
  #geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) +
  facet_wrap(~Site)#,strip.position = "bottom")

dm1<-gls(mpd.obs.z~HostPlant+Site,data=dat6,na.action=na.omit)
m1<-gls(mpd.obs.z~Site,data=dat6,na.action=na.omit)
anova(m1)

m2<-lme(mpd.obs.z~Site,random=~1|HostPlant,data=dat6,na.action=na.omit)
anova(m2)


#Summing by species and site
# Phrag TC is sig, Phrag LUMCON nearly sig; alterniflora CERF nearly sig, alternaflora LUMCON sig; patens CERF nearly sig
phydist <- cophenetic(tree)

dat.comm.sum<-dat6%>%
  group_by(HostPlant,Site)%>%
  summarise(across(OTU0:OTU9,sum))%>%
  unite(HostPlantSite,HostPlant,Site)
dat.comm.sum$HostPlantSite<-apply(as.matrix(dat.comm.sum$HostPlantSite),1,function(x){gsub(" ", "", x)})
dat.comm.sum<-data.frame(dat.comm.sum)
row.names(dat.comm.sum)<-dat.comm.sum$HostPlantSite
dat.comm.sum$HostPlantSite<-NULL

#ses.mntd.result <- ses.mntd(dat.comm.sum, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=999)
#ses.mntd.result

ses.mpd.result.sum <- ses.mpd(dat.comm.sum, phydist, null.model="taxa.labels", abundance.weighted=F, runs=999)
ses.mpd.result.sum


#Summing by species
#spartina alterniflora nearly sig
dat.comm.sum.sp<-dat6%>%
  group_by(HostPlant)%>%
  summarise(across(OTU0:OTU9,sum))
dat.comm.sum.sp$HostPlant<-apply(as.matrix(dat.comm.sum.sp$HostPlant),1,function(x){gsub(" ", "", x)})
dat.comm.sum.sp<-data.frame(dat.comm.sum.sp)
row.names(dat.comm.sum.sp)<-dat.comm.sum.sp$HostPlant
dat.comm.sum.sp$HostPlant<-NULL

ses.mpd.result.sum.sp <- ses.mpd(dat.comm.sum.sp, phydist, null.model="taxa.labels",abundance.weighted=FALSE, runs=999)
ses.mpd.result.sum.sp



#Summing by site
#no significance
dat.comm.sum.site<-dat6%>%
  group_by(Site)%>%
  summarise(across(OTU0:OTU9,sum))
dat.comm.sum.site$Site<-apply(as.matrix(dat.comm.sum.site$Site),1,function(x){gsub(" ", "", x)})
dat.comm.sum.site<-data.frame(dat.comm.sum.site)
row.names(dat.comm.sum.site)<-dat.comm.sum.site$Site
dat.comm.sum.site$Site<-NULL

ses.mpd.result.sum.site <- ses.mpd(dat.comm.sum.site, phydist, null.model="taxa.labels",abundance.weighted=F, runs=999)
ses.mpd.result.sum.site





##### Tree plotting #####
# par(mfrow=c(2,2))
# for (i in names(traits)) {
#   plot(phy, show.tip.label=FALSE, main=i)
#   tiplabels(pch=22, col=traits[,i]+1, bg=traits[,i]+1, cex=1.5)
# }

dat.comm.sum.spt<-as.data.frame(t(dat.comm.sum.sp))
dat.comm.sumt<-as.data.frame(t(dat.comm.sum))

plot(tree, show.tip.label=FALSE)
colort<-ifelse(dat.comm.sumt$Spartinaalterniflora_LUMCON >0,2,1)
tiplabels(pch=22, col=colort, bg=colort, cex=1.5)

plot(tree, show.tip.label=T)
