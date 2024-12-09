##### Loading packages #####

library(tidyverse)
#library(picante)
library(plotrix)
library(phyloseq)
library(cowplot)
library(nlme)
library(multcomp)
#library(emmeans)
#library(FUNGuildR)
#devtools::install_github("brendanf/FUNGuildR")

#library(FUNGuildR)

#install.packages("remotes")
#remotes::install_github("jfq3/QsRutils") #distance matrix subsetting function

#library(QsRutils)
#library(remotes)
#library(ape)
#library(phytools)


save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture2/workspace.Rdata")  # 

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture2/workspace.Rdata")  # 


##### Reading in data #####

##### T-BAS data, 97% similarity ####

#OTUold is the otu name that was input for the isolate name in T-BAS, every isolate has a unique OTUold. The numbers of OTUold are from Nelle's original T-BAS run and then we added A, B C, etc to them to make them unique.
#Query.sequence is also unique to each isolate. it is the OTUold plus the species name that Nelle's original T-BAS run came up with
#OTU is the new analysis OTU - what I want
#I think I want "Genusspecies" as for the taxonomy. "taxon assignment" is similar to genusspecies but sometimes has multiple taxa listed
#anything with CERVNAZT is from the input data (so from Nelle's data), not what I want

OTUreport<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/FarrerTBAScleaned5/tbas21_archiveWG5WGW7J_/assignments_report_addvoucherWG5WGW7J.csv",stringsAsFactors = T)
cbind(OTUreport$otu,OTUreport$Query.sequence)
OTUreport2<-OTUreport%>%
  dplyr::select(Query.sequence,Phylum,Taxon.assignment,Genusspecies,otu,OTU_CERVNAZT:Juncus_roemerianus_CERVNAZT,Trophic.Mode,Guild)%>%
  rename(OTUold=OTU_CERVNAZT,HostPlant=HostPlant_CERVNAZT,Site=Site_CERVNAZT,TurtleCove=TurtleCove_CERVNAZT,LUMCON=LUMCON_CERVNAZT,CERF=CERF_CERVNAZT,Spartina_patens=Spartina_patens_CERVNAZT,Spartina_alterniflora=Spartina_alterniflora_CERVNAZT,Phragmites_australis=Phragmites_australis_CERVNAZT,Sagittaria_lancifolia=Sagittaria_lancifolia_CERVNAZT,Juncus_roemerianus=Juncus_roemerianus_CERVNAZT,OTU=otu)%>%  mutate(Site = factor(Site, levels = c("Turtle Cove", "CERF", "LUMCON")))

head(OTUreport2)


##### Culture data ####

#Culture IDs are unique to the sequence
#Otus 18 and 46 did not have an ITS region, remove
culturefile<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Copy of Mastersheet - July 29, 2020, 11_55 AM.csv",stringsAsFactors = T)
culturefile2<-culturefile%>%
  filter(is.na(CultureID)==F,OTU!="OTU18",OTU!="OTU46")%>%
  dplyr::select(-Sequence)%>%
  unite(Query.sequence,OTU,TBAS_Species,remove=F)%>%
  dplyr::select(Query.sequence,Year,CultureID,PlantIndividual,TBAS_Species)%>%
  separate(PlantIndividual,c("Plant","Individual"),sep=" ")%>%
  unite(PlantIndividual,Plant,Individual,remove=T)
head(culturefile2)

length(culturefile2$CultureID)

dat<-full_join(OTUreport2,culturefile2,by="Query.sequence")
head(dat)
dim(dat)

#Collapse identical OTUs into one row. Looking at abundances of OTUs in the whole dataset. there are many many 1's (singletons). also this is needed to create the community dataset below
dat2<-dat%>%
  group_by(HostPlant,Site,Year,PlantIndividual,OTU,Genusspecies)%>%
  summarise(abundance=n())
as.data.frame(dat2)
data.frame(dat2$HostPlant,dat2$PlantIndividual,dat2$Year,dat2$OTU,dat2$abundance)

#make wide
dat3<-dat2%>%
  ungroup()%>%
  unite(PlantIndividualYear,PlantIndividual,Year,remove=F)%>%
  unite(HostPlantSite,HostPlant,Site,remove=F)%>%
  dplyr::select(-Genusspecies)%>%
  spread(OTU,abundance,fill=0)
dat.comm<-data.frame(dat3[,7:66])
row.names(dat.comm)<-dat3$PlantIndividualYear


##### Phylogenetic tree #####
#note - the uncleaned tree in this directory contains the genus species names attached to the OTU numbers
tree<-read.newick("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/FarrerTBAScleaned3/Farrertreecleaned.nwk")

plot(tree)





##### Faith's Phylogenetic distance #####
PD<-pd(dat.comm,tree)

dat4<-data.frame(dat3[,1:6],PD,dat3[,7:66])





##### MPD by plantindividual #####
phydist <- cophenetic(tree)
ses.mpd.result.notweighted <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=FALSE, runs=999) #takes 5 min with 999
ses.mpd.result.notweighted
ses.mpd.result.notweighted$PlantIndividualYear<-rownames(ses.mpd.result.notweighted)
ses.mpd.result.notweighted1<-ses.mpd.result.notweighted%>%
  select(PlantIndividualYear,mpd.obs.z)%>%
  rename(mpd.obs.z.notweighted=mpd.obs.z)

ses.mpd.result.weighted <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=TRUE, runs=999) #takes 5 min with 999
ses.mpd.result.weighted
ses.mpd.result.weighted$PlantIndividualYear<-rownames(ses.mpd.result.weighted)
ses.mpd.result.weighted1<-ses.mpd.result.weighted%>%
  select(PlantIndividualYear,mpd.obs.z)%>%
  rename(mpd.obs.z.weighted=mpd.obs.z)

dat5<-dat4%>%
  full_join(ses.mpd.result.notweighted1)%>%
  full_join(ses.mpd.result.weighted1)%>%
  mutate(HostPlant=recode(HostPlant,"Spartina alterniflora"= "S. alterniflora","Spartina patens"="S. patens","Phragmites australis"="P. australis","Juncus roemerianus"="J. roemerianus","Sagittaria lancifolia"="S. lancifolia"))
  
  
dat6<-data.frame(dat5[,1:8],dat5[,69:70],dat5[,9:68])
head(dat6)


dplyr::mutate(Species=dplyr::recode(Species,"Phragmites australis (L)"="Phragmites australis"))





##### phyloseq object ####
otus<-dat6[,11:70]
otus2<-t(otus)
sampleotus<-dat6[,c(1:10)]
taxonomyotus<-as.matrix(data.frame(Kingdom=row.names(otus2),Phylum=row.names(otus2),Class=row.names(otus2),Order=row.names(otus2),Class=row.names(otus2),Family=row.names(otus2),Genus=row.names(otus2),Species=row.names(otus2)))
rownames(taxonomyotus)<-row.names(otus2)

datp <- merge_phyloseq(otu_table(otus2,taxa_are_rows = T), tax_table(taxonomyotus), sample_data(sampleotus),tree)

#calculate unifrac distances
unifracp<-unifrac(otus,tree)



#### Summary of files #####
datp #phyloseq object
head(dat6) #big data frame, wide data format, dat3 plus PD and MPD data
dat7 #dat6 plus funguild pathogens/symbiotrophs
dat2 #long dataformat 




###### testing blasting to unite ######
#unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAMarsh/Survey/Stats/Gradient/QIIME2/sh_general_release_dynamic_s_04.02.2020.fasta" #

unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAMarsh/Culturing/FiguresStats/LAmarshCulture/sh_general_release_dynamic_s_29.11.2022.fasta"

#to look at all or random sequences, not just the 60 otu rep set. can get seuqences in report_placed_qrydata.csv
test<-read.csv("test.csv",header=T)
taxatest<-assignTaxonomy(test, unite.ref, multithread = TRUE, minBoot=50, tryRC = TRUE,outputBootstraps=T) #was minBoot=70
taxatestonly<-taxatest$tax;rownames(taxatestonly)<-1:dim(taxatest$tax)[1]
taxatestboot<-taxatest$boot;rownames(taxatestboot)<-1:dim(taxatest$tax)[1]
taxatestonly
taxatestboot

#I got this data from assignments_report_nodupsCERVNAZD.csv
my60otus<-read.csv("otus60.csv",header=T);rownames(my60otus)<-my60otus$abundance

taxa<-assignTaxonomy(my60otus, unite.ref, multithread = TRUE, minBoot=70, tryRC = TRUE,outputBootstraps=T) #was minBoot=70
taxaonly<-data.frame(taxa$tax);rownames(taxaonly)<-1:dim(taxa$tax)[1]
taxaboot<-data.frame(taxa$boot);rownames(taxaboot)<-1:dim(taxa$tax)[1]
taxaonly
taxaboot
rownames(taxaonly)<-rownames(my60otus)
genusspecies<-data.frame(otu=rownames(my60otus),genusspecies=paste(gsub("^.*?__","",taxaonly[,"Genus"]),gsub("^.*?__","",taxaonly[,"Species"])),genusboot=taxaboot[,"Genus"],speciesboot=taxaboot[,"Species"])

sort(unique(genusspecies$genusspecies))
sort(unique(taxaonly$Phylum))

taxaonly%>%group_by(Phylum)%>%tally()




##### FunGuild #####

#The weird thing about FUNGuildR is that for some taxa, like when I have Fusarium equiseti it matches to Nectriaceae rather than Fusarium in the database. When I used the FUNGUILD in the terminal, it worked better and matched to the lowest taxonomic level. 

#using FUNGuildR
# genussp<-paste(gsub("^.*?__","",taxaonly$Genus),gsub("^.*?__","",taxaonly$Species))
# taxaonly$Taxonomy<-paste(gsub("^.*?__","",taxaonly$Kingdom),gsub("^.*?__","",taxaonly$Phylum),gsub("^.*?__","",taxaonly$Class),gsub("^.*?__","",taxaonly$Order),gsub("^.*?__","",taxaonly$Family),gsub("^.*?__","",taxaonly$Genus),genussp,sep=";")

#Download the up-to-date database
#fung <- get_funguild_db()

#save it to my computer for reproducable data
#saveRDS(fung, "funguild.rds")

#load it back into workspace
#fung <- loadRDS("funguild.rds")

# fung_guilds <- funguild_assign(taxaonly, db = fung)

#temp3 <- funguild_query("Buergenerula*", "taxon", db = fung)


#Doing it in terminal - use this
taxaonly2<-taxaonly
taxaonly2$taxonomy<-paste(taxaonly2$Kingdom,taxaonly2$Phylum,taxaonly2$Class,taxaonly2$Order,taxaonly2$Family,taxaonly2$Genus,taxaonly2$Species,sep=";")

write.csv(taxaonly2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture/FUNGuild/taxaonly.csv")

#open the file, make the first column name OTU ID, save as a .txt file (add the .txt extension). then run below:

python FUNGuild.py taxa -otu taxaonly.txt -format tsv -column taxonomy -classifier unite
python FUNGuild.py guild -taxa taxaonly.taxa.txt         

guilds<-read.delim("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture/FUNGuild/taxaonly.taxa.guilds.txt")
View(guilds)

guilds%>%filter(trophicMode==("Pathotroph"))
guilds%>%filter(trophicMode==("Symbiotroph"))

View(funguild_query("*Saprotroph*", "trophicMode", db = guilds))

#merge with OTU table (of datp) with guilds from above
otutomerge<-data.frame(otu_table(datp))
otutomerge$OTU<-rownames(otutomerge)
head(guilds)
guilds2<-guilds%>%
  full_join(otutomerge)
head(guilds2)
ind<-which(guilds2$guild=="Plant Pathogen")
plantpathogen<-colSums(guilds2[ind,18:89])
plantpathogentaxa<-colSums(guilds2[ind,18:89]>0)

ind<-funguild_query("*Plant Pathogen*", "guild", db = guilds2)$OTU#trophicMode
ind2<-which(guilds2$OTU%in%ind)
plantpathogenbroad<-colSums(guilds2[ind2,18:89])
plantpathogenbroadtaxa<-colSums(guilds2[ind2,18:89]>0)

ind<-funguild_query("*Plant Pathogen*", "guild", db = guilds2)#trophicMode
ind2<-ind[which(ind$confidenceRanking%in%c("Probable","Highly Probable")),"OTU"]
ind3<-which(guilds2$OTU%in%ind2)
plantpathogenbroadprobablehp<-colSums(guilds2[ind3,18:89])
plantpathogenbroadprobablehptaxa<-colSums(guilds2[ind3,18:89]>0)

ind<-funguild_query("*Symbiotroph*", "trophicMode", db = guilds2)$OTU
ind2<-which(guilds2$OTU%in%ind)
symbiotroph<-colSums(guilds2[ind2,18:89])
symbiotrophtaxa<-colSums(guilds2[ind2,18:89]>0)

ind<-funguild_query("*Symbiotroph*", "trophicMode", db = guilds2)
ind2<-ind[which(ind$confidenceRanking%in%c("Probable","Highly Probable")),"OTU"]
ind3<-which(guilds2$OTU%in%ind2)
symbiotrophprobablehp<-colSums(guilds2[ind3,18:89])
symbiotrophprobablehptaxa<-colSums(guilds2[ind3,18:89]>0)

totalabundance<-colSums(guilds2[,18:89])

temp<-data.frame(sample_data(datp))
temp$rownames<-rownames(temp)
temp<-temp%>%
  dplyr::select(rownames,PlantIndividualYear)
temp2<-data.frame(temp,plantpathogen,plantpathogenbroad,plantpathogenbroadprobablehp,symbiotroph,symbiotrophprobablehp,totalabundance,plantpathogentaxa,plantpathogenbroadtaxa,plantpathogenbroadprobablehptaxa,symbiotrophtaxa,symbiotrophprobablehptaxa)

dat7<-dat6%>%
  full_join(temp2)

phragspartina<-dat6%>%
  filter(HostPlant%in%c('Phragmites australis','Spartina patens'))
phragspartina<-dat7%>%
  filter(HostPlant%in%c('P. australis','S. patens'))



##### git hub token stuff #####
install.packages("gitcreds")
library(gitcreds)
gitcreds_set()
 #first when it asks to enter password or token I put my computer password
 #then do gitcreds_set() again and select 2, then paste my token
#Note: use usename (email) and token, when RStudio wants the github password
