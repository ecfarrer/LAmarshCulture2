#Ordination 


# Because some species were only collected in some years, I can't put year in TC models (phrag only from 2018, patens only 2017), CF (patens only in 2017), LU (patens only 2017)


###### Year - does year matter ######
#upshot - yes year sometimes matters, esp for phrag, not as much for spartina alterniflora

sample_data(datp)$Year<-factor(sample_data(datp)$Year)

#have multiple years: phragCF, phragLU, SaCF, SaLU
datpyear<-datp%>%
  subset_samples(Site!="Turtle Cove")%>%
  subset_samples(HostPlant!="Juncus roemerianus")%>%
  subset_samples(HostPlant!="Spartina patens")%>%
  subset_samples(HostPlant!="Sagittaria lancifolia")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracyear<-subset_dist(datpyear, unifracp)
sample_data(datpyear)$Year<-factor(sample_data(datpyear)$Year)
sample_data(datpyear)$HostPlantYear<-paste(sample_data(datpyear)$HostPlant,sample_data(datpyear)$Year)
sample_data(datpyear)$HostPlantSiteYear<-paste(sample_data(datpyear)$HostPlant,sample_data(datpyear)$Site,sample_data(datpyear)$Year)

mynmdsyearj <- ordinate(datpyear, "CAP",distance(datpyear, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Site+Year))#
anova(mynmdsyearj,by="margin",permutations = how(nperm=9999))#blocks=sample_data(datpTC)$HostPlant,

mynmdsyearu <- ordinate(datpyear, "CAP",distance=unifracyear,formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsyearu,by="margin",permutations = how(nperm=99999))

plot_ordination(datpyear, mynmdsyearj, type="samples", color="HostPlantSiteYear",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlantSiteYear),level=.95)


#subsetting
datpyear<-datp%>%
  subset_samples(Site=="CERF"&HostPlant=="Spartina alterniflora")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracyear<-subset_dist(datpyear, unifracp)

mynmdsyearj <- ordinate(datpyear, "CAP",distance(datpyear, method = "jaccard", binary = TRUE),formula=as.formula(~Year))#
anova(mynmdsyearj,by="margin",permutations = how(nperm=9999))

mynmdsyearu <- ordinate(datpyear, "CAP",distance=unifracyear,formula=as.formula(~Year))
anova(mynmdsyearu,by="margin",permutations = how(nperm=9999))

plot_ordination(datpyear, mynmdsyearj, type="samples", color="Year",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Year),level=.95)





#calculating MDP (mean pairwise distance) matrix for use in ordination

MPDdist<-comdist(data.frame(t(otu_table(datp))), cophenetic(phy_tree(datp)), abundance.weighted=F)
MPDdistweighted<-comdist(data.frame(t(otu_table(datp))), cophenetic(phy_tree(datp)), abundance.weighted=T)

#unifracp<-unifrac(otus,tree) #doesnt work with subsetting b/c sample names are not the same
unifracp<-unifrac(t(otu_table(datp)),tree)


##### all #####
datp
MPDdist #(unweighted)
MPDdistweighted

#take out singletons. 
datp2<-datp%>%
  filter_taxa(function(x) sum(x>0) >1, prune=T)
datp3<-prune_samples(sample_sums(datp2)>=1, datp2)


#bray curtis
mynmdsALLb <- ordinate(datp3, "CAP",distance(datp3, method = "bray", binary = F),formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsALLb,by="margin",permutations = how(nperm=9999))

#Getting variance explained
m1<-ordinate(datp3, "CAP",distance(datp3, method = "bray", binary = F),formula=as.formula(~HostPlant+Condition(Year+Site)))
summary(m1)
m1<-ordinate(datp3, "CAP",distance(datp3, method = "bray", binary = F),formula=as.formula(~Site+Condition(Year+HostPlant)))
summary(m1)

###plots
plot_ordination(datp3, mynmdsALLb, type="samples", color="Site",shape="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  scale_shape_manual(values=c(16,3,15,17,5))+ 
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(x=CAP1,y=CAP2,fill=Site),level=.95,inherit.aes = F)

plot1<-plot_ordination(datp3, mynmdsALLb, type="samples", color="HostPlant",shape="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  scale_shape_manual(values=c(5,3,16,15,17))+ 
  geom_point(size = 2)+
  #stat_ellipse(geom = "polygon", type="t", alpha=0, aes(x=CAP1,y=CAP2,color=Site),level=.95,inherit.aes = F)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(x=CAP1,y=CAP2,fill=HostPlant),level=.95,inherit.aes = F)

#Plotting convex hulls rather than ellipses, USE THIS FOR MANUSCRIPT
site_scores1 <- data.frame(cbind(sample_data(datp3),vegan::scores(mynmdsALLb)$sites,labels=rownames(vegan::scores(mynmdsALLb)$sites)))

hull1 <- site_scores1 %>%
  group_by(HostPlant) %>%
  slice(chull(CAP1,CAP2))

plot1<-
  ggplot(site_scores1)+
  theme_classic()+#  theme(legend.position = "none")
  xlab("CAP1 [5.2%]") +  # 
  ylab("CAP2 [3.6%]") +  # 
  #scale_color_manual(values = scales::viridis_pal(option = "turbo")(6))+
  #scale_fill_manual(values = scales::viridis_pal(option = "turbo")(6))+
  #scale_color_manual(values = scales::viridis_pal(option = "viridis")(5))+
  #scale_fill_manual(values = scales::viridis_pal(option = "viridis")(5))+
  scale_color_manual(values = c("#3D65A5", "#F05039","#EEBAB4","#7ca1cc","#12285c"))+
  scale_fill_manual(values = c("#3D65A5", "#F05039","#EEBAB4","#7ca1cc","#12285c"))+
  scale_shape_manual(values=c(5,3,16))+ 
  geom_point(aes(x=CAP1, y=CAP2,color=HostPlant,shape=Site),size = 2)+
  geom_polygon(data=hull1,aes(x=CAP1,y=CAP2, fill=HostPlant,colour = HostPlant),alpha=.2)

# p1=plot_ordination(datp3, mynmdsALLb, type="samples",shape="HostPlant",color="HostPlant",axes=c(1,2))
# 
# ggplot(p1$data, aes(x=CAP1,y=CAP2,color="HostPlant",shape="HostPlant")) +
#   stat_ellipse(data=test,geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)+
#   geom_point(data=test)
# 
# test<-plot_ordination(datp3, mynmdsALLb, type="samples",color=NULL,shape=NULL,axes=c(1,2),justDF = T)
  

#MPD
mynmdsALLm <- ordinate(datp, "CAP",distance=MPDdistweighted,formula=as.formula(~HostPlant+Site+Condition(Year)))
#mynmdsALLm <- ordinate(datp, "CAP",distance=MPDdistweighted,formula=as.formula(~HostPlant+Site))
anova(mynmdsALLm,by="margin",permutations = how(nperm=9999))
#anova(mynmdsALLm,by="margin",permutations = how(blocks=sample_data(datp)$Year,nperm=9999)) #adding permutations within year does not change the anova results

#percent explained
m1 <- ordinate(datp, "CAP",distance=MPDdistweighted,formula=as.formula(~HostPlant+Condition(Year+Site)))
summary(m1)
m1 <- ordinate(datp, "CAP",distance=MPDdistweighted,formula=as.formula(~Site+Condition(Year+HostPlant)))
summary(m1)

plot2<-plot_ordination(datp, mynmdsALLm, type="samples", color="HostPlant",shape="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  scale_shape_manual(values=c(5,3,16,15,17))+ 
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(x=CAP1,y=CAP2,fill=HostPlant),level=.95,inherit.aes = F)

#with edges on the ellipse
plot_ordination(datp, mynmdsALLm, type="samples", color="HostPlant",shape="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  scale_shape_manual(values=c(5,3,16,15,17))+ 
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(x=CAP1,y=CAP2,fill=HostPlant, color=HostPlant),level=.95,inherit.aes = F)

plot_ordination(datp, mynmdsALLm, type="samples", color="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlant),level=.95)

#Using convex hulls USE THIS FOR MS
site_scores2 <- data.frame(cbind(sample_data(datp),vegan::scores(mynmdsALLm)$sites,labels=rownames(vegan::scores(mynmdsALLm)$sites)))

hull2 <- site_scores2 %>%
  group_by(HostPlant) %>%
  slice(chull(CAP1,CAP2))

plot2<-ggplot(site_scores2)+
  theme_classic()+#  theme(legend.position = "none")
  xlab("CAP1 [7.6%]") +  # 
  ylab("CAP2 [1.7%]") +  # 
  scale_color_manual(values = c("#3D65A5", "#F05039","#EEBAB4","#7ca1cc","#12285c"))+
  scale_fill_manual(values = c("#3D65A5", "#F05039","#EEBAB4","#7ca1cc","#12285c"))+
  scale_shape_manual(values=c(5,3,16))+ 
  geom_point(aes(x=CAP1, y=CAP2,color=HostPlant,shape=Site),size = 2)+
  geom_polygon(data=hull2,aes(x=CAP1,y=CAP2, fill=HostPlant,colour = HostPlant),alpha=.2)



#jaccard
mynmdsALLj <- ordinate(datp, "CAP",distance(datp, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsALLj,by="margin",permutations = how(nperm=9999))

plot_ordination(datp, mynmdsALLj, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)

# figure for comp and mdp
#I could try flipping over the y axis
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/ordinationcompmpdhulls2.pdf",height = 3.5,width = 10.5) #was ,height = 3.5,width = 11
plot_grid(plot1, plot2, nrow = 1, labels=c("A. Composition","B. MPD"),label_fontface = "plain",hjust=c(-.5,-1.),vjust=1)
dev.off()




##### Phrag and patens #####

datpphragspar<-datp%>%
  subset_samples(HostPlant%in%c("P. australis","S. patens"))%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
colSums(otu_table(datpphragspar))
#don't need to takeout any samples, they all have >0 species in them
unifracphragspar<-subset_dist(datpphragspar, unifracp)
MPDdistphragspar<-subset_dist(datpphragspar, MPDdistweighted)

#Bray
mynmdsphragsparb <- ordinate(datpphragspar, "CAP",distance(datpphragspar, method = "bray", binary = F),formula=as.formula(~HostPlant+Site+HostPlant*Site+Condition(Year)))
anova(mynmdsphragsparb,by="margin",permutations = how(nperm=9999))
mynmdsphragsparb <- ordinate(datpphragspar, "CAP",distance(datpphragspar, method = "bray", binary = F),formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsphragsparb,by="margin",permutations = how(nperm=9999))
#variance explained
m1 <- ordinate(datpphragspar, "CAP",distance(datpphragspar, method = "bray", binary = F),formula=as.formula(~HostPlant*Site+Condition(HostPlant+Site+Year)))
summary(m1)
m1 <- ordinate(datpphragspar, "CAP",distance(datpphragspar, method = "bray", binary = F),formula=as.formula(~HostPlant+Condition(Site+Year)))
summary(m1)
m1 <- ordinate(datpphragspar, "CAP",distance(datpphragspar, method = "bray", binary = F),formula=as.formula(~Site+Condition(HostPlant+Year)))
summary(m1)

mynmdsphragsparm <- ordinate(datpphragspar, "CAP",distance=MPDdistphragspar,formula=as.formula(~HostPlant+Site+HostPlant*Site+Condition(Year)))
anova(mynmdsphragsparm,by="margin",permutations = how(nperm=9999))
mynmdsphragsparm <- ordinate(datpphragspar, "CAP",distance=MPDdistphragspar,formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsphragsparm,by="margin",permutations = how(nperm=9999))
# variance explained
m1 <- ordinate(datpphragspar, "CAP",distance=MPDdistphragspar,formula=as.formula(~HostPlant*Site+Condition(Year+Site+HostPlant)))
summary(m1)
m1 <- ordinate(datpphragspar, "CAP",distance=MPDdistphragspar,formula=as.formula(~HostPlant+Condition(Year+Site)))
summary(m1)
m1 <- ordinate(datpphragspar, "CAP",distance=MPDdistphragspar,formula=as.formula(~Site+Condition(Year+HostPlant)))
summary(m1)


mynmdsphragsparj <- ordinate(datpphragspar, "CAP",distance(datpphragspar, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Site+HostPlant*Site+Condition(Year)))
anova(mynmdsphragsparj,by="margin",permutations = how(nperm=9999))
mynmdsphragsparj <- ordinate(datpphragspar, "CAP",distance(datpphragspar, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsphragsparj,by="margin",permutations = how(nperm=9999))

mynmdsphragmsparu <- ordinate(datpphragspar, "CAP",distance=unifracphragspar,formula=as.formula(~HostPlant+Site+HostPlant*Site+Condition(Year)))
anova(mynmdsphragmsparu,by="margin",permutations = how(nperm=9999))
mynmdsphragmsparu <- ordinate(datpphragspar, "CAP",distance=unifracphragspar,formula=as.formula(~HostPlant+Site+Condition(Year)))
anova(mynmdsphragmsparu,by="margin",permutations = how(nperm=9999))



#for plotting
mynmdsphragsparb <- ordinate(datpphragspar, "CAP",distance(datpphragmitesspartina, method = "bray", binary = F),formula=as.formula(~HostPlantSite+Condition(Year)))
mynmdsphragsparm <- ordinate(datpphragspar, "CAP",distance=MPDdistphragspar,formula=as.formula(~HostPlantSite+Condition(Year)))

plot_ordination(datpphragspar, mynmdsphragsparb, type="samples", color="HostPlantSite",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlantSite),level=.95)
plot_ordination(datpphragspar, mynmdsphragsparm, type="samples", color="HostPlantSite",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlantSite),level=.95)





##### By Site #####
#Does host plant affect composition within each site?
#blocks=sample_data(datpTC)$HostPlant

#TURTLE COVE

datpTC<-datp%>%
  subset_samples(Site=="Turtle Cove")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracTC<-subset_dist(datpTC, unifracp)
MDPdistTC<-subset_dist(datpTC, MDPdist)

mynmdsTCj <- ordinate(datpTC, "CAP",distance(datpTC, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsTCj,by="margin",permutations = how(nperm=9999))

mynmdsTCu <- ordinate(datpTC, "CAP",distance=unifracTC,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsTCu,by="terms",permutations = how(nperm=9999))

mynmdsTCm <- ordinate(datpTC, "CAP",distance=MDPdistTC,formula=as.formula(~HostPlant))
anova(mynmdsTCm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpTC, mynmdsTCu, type="samples", color="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  #scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  #scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlant),level=.95)


#CERF

datpCF<-datp%>%
  subset_samples(Site=="CERF")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracCF<-subset_dist(datpCF, unifracp)
MDPdistTC<-subset_dist(datpCF, MDPdist)

mynmdsCFj <- ordinate(datpCF, "CAP",distance(datpCF, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsCFj,by="terms",permutations = how(nperm=9999))

mynmdsCFu <- ordinate(datpCF, "CAP",distance=unifracCF,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsCFu,by="terms",permutations = how(nperm=9999))#

mynmdsCFm <- ordinate(datpCF, "CAP",distance=MDPdistTC,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsCFm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpCF, mynmdsCFm, type="samples", color="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlant),level=.95)


#LUMCON

datpLU<-datp%>%
  subset_samples(Site=="LUMCON")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracLU<-subset_dist(datpLU, unifracp)
MDPdistLU<-subset_dist(datpLU, MDPdist)

mynmdsLUj <- ordinate(datpLU, "CAP",distance(datpLU, method = "jaccard", binary = TRUE),formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsLUj,by="terms",permutations = how(nperm=9999))

mynmdsLUu <- ordinate(datpLU, "CAP",distance=unifracLU,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsLUu,by="terms",permutations = how(nperm=9999))

mynmdsLUm <- ordinate(datpLU, "CAP",distance=MDPdistLU,formula=as.formula(~HostPlant+Condition(Year)))
anova(mynmdsLUm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpLU, mynmdsLUm, type="samples", color="HostPlant",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=HostPlant),level=.95)





##### By Host #####
#Does site affect composition within each site?


#Phragmites 

datpPa<-datp%>%
  subset_samples(HostPlant=="Phragmites australis")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracPa<-subset_dist(datpPa, unifracp)
MDPdistPa<-subset_dist(datpPa, MDPdist)

mynmdsPaj <- ordinate(datpPa, "CAP",distance(datpPa, method = "jaccard", binary = TRUE),formula=as.formula(~Site+Condition(Year)))
anova(mynmdsPaj,by="terms",permutations = how(nperm=9999))

mynmdsPau <- ordinate(datpPa, "CAP",distance=unifracPa,formula=as.formula(~Site+Condition(Year)))#
anova(mynmdsPau,by="terms",permutations = how(nperm=9999))

mynmdsPam <- ordinate(datpPa, "CAP",distance=MDPdistPa,formula=as.formula(~Site+Condition(Year)))#+Condition(Year)
anova(mynmdsPam,by="terms",permutations = how(nperm=9999))

plot_ordination(datpPa, mynmdsPaj, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  #scale_color_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  #scale_fill_manual(values = c("#0047b3", "#99c2ff","#2d862d","#79d279","#b30000","#ff8080"),labels = c("Fresh Native", "Fresh Phragmites","Brackish Native","Brackish Phragmites","Saline Native","Saline Phragmites"),name = "Marsh class/Invasion")+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)


#Spartina alterniflora 

datpSa<-datp%>%
  subset_samples(HostPlant=="Spartina alterniflora")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracSa<-subset_dist(datpSa, unifracp)
MDPdistSa<-subset_dist(datpSa, MDPdist)

mynmdsSaj <- ordinate(datpSa, "CAP",distance(datpSa, method = "jaccard", binary = TRUE),formula=as.formula(~Site+Condition(Year)))
anova(mynmdsSaj,by="margin",permutations = how(nperm=9999))

mynmdsSau <- ordinate(datpSa, "CAP",distance=unifracSa,formula=as.formula(~Site+Condition(Year)))
anova(mynmdsSau,by="margin",permutations = how(nperm=9999))

mynmdsSam <- ordinate(datpSa, "CAP",distance=MDPdistSa,formula=as.formula(~Site+Condition(Year)))
anova(mynmdsSam,by="margin",permutations = how(nperm=9999))

plot_ordination(datpSa, mynmdsSam, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)


#Spartina patens 
#only collected in 2017

datpSp<-datp%>%
  subset_samples(HostPlant=="Spartina patens")%>%
  filter_taxa(function(x) sum(x>0) >0, prune=T)
unifracSp<-subset_dist(datpSp, unifracp)
MDPdistSp<-subset_dist(datpSp, MDPdist)

mynmdsSpj <- ordinate(datpSp, "CAP",distance(datpSp, method = "jaccard", binary = TRUE),formula=as.formula(~Site))
anova(mynmdsSpj,by="terms",permutations = how(nperm=9999))

mynmdsSpu <- ordinate(datpSp, "CAP",distance=unifracSp,formula=as.formula(~Site))
anova(mynmdsSpu,by="terms",permutations = how(nperm=9999))

mynmdsSpm <- ordinate(datpSp, "CAP",distance=MDPdistSp,formula=as.formula(~Site))
anova(mynmdsSpm,by="terms",permutations = how(nperm=9999))

plot_ordination(datpSp, mynmdsSpm, type="samples", color="Site",axes=c(1,2))+
  theme_classic()+#  theme(legend.position = "none")
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon", type="t", alpha=0.2, aes(fill=Site),level=.95)




