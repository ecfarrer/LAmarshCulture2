##### Loading packages #####

library(tidyverse)
library(picante) # for faith's phylogenetic distance
library(plotrix)
library(phyloseq)
library(cowplot)
library(nlme)
library(multcomp)
library(dada2)
library(emmeans) #for t-tests for MPD
#library(FUNGuildR)
#devtools::install_github("brendanf/FUNGuildR")

#library(FUNGuildR)

#install.packages("remotes")
#remotes::install_github("jfq3/QsRutils") #distance matrix subsetting function

library(QsRutils)
#library(remotes)
#library(ape)
library(phytools) #to read the tree in
library(viridis)#for colorblind friendly colors


save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture2/workspace.Rdata")  # 

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/LAmarshCulture2/workspace.Rdata")  # 


##### Reading in data #####

##### T-BAS data, 97% similarity, there are 56 OTUs OTU0-OTU55 ####

#OTUold is the otu name that was input for the isolate name in T-BAS, every isolate has a unique OTUold. The numbers of OTUold are from Nelle's original T-BAS run and then we added A, B C, etc to them to make them unique.
#Query.sequence is also unique to each isolate. it is the OTUold plus the species name that Nelle's original T-BAS run came up with
#OTU is the new analysis OTU - what I want
#Before I considered "GenusSpecies" for the taxonomy [this used to be called "Genusspecies" with no capital S], but now many of the entries are blank. "taxon assignment" is similar to genusspecies but sometimes has multiple taxa listed. I will still include it now so that my dataframes are the same size
#any genus/species info with WG5WGW7J is from the input data (so from Nelle's data), not what I want

OTUreport<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/FarrerTBAScleaned5/tbas21_archiveWG5WGW7J_/assignments_report_addvoucherWG5WGW7J.csv",stringsAsFactors = T)
cbind(OTUreport$otu,OTUreport$Query.sequence)
OTUreport2<-OTUreport%>%
  dplyr::select(Query.sequence,Phylum,Taxon.assignment,GenusSpecies,otu,OTU_WG5WGW7J:Juncus_roemerianus_WG5WGW7J,Trophic.Mode,Guild)%>%
  rename(OTUold=OTU_WG5WGW7J,HostPlant=HostPlant_WG5WGW7J,Site=Site_WG5WGW7J,TurtleCove=TurtleCove_WG5WGW7J,LUMCON=LUMCON_WG5WGW7J,CERF=CERF_WG5WGW7J,Spartina_patens=Spartina_patens_WG5WGW7J,Spartina_alterniflora=Spartina_alterniflora_WG5WGW7J,Phragmites_australis=Phragmites_australis_WG5WGW7J,Sagittaria_lancifolia=Sagittaria_lancifolia_WG5WGW7J,Juncus_roemerianus=Juncus_roemerianus_WG5WGW7J,OTU=otu)%>%  mutate(Site = factor(Site, levels = c("Turtle Cove", "CERF", "LUMCON")))

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
  group_by(HostPlant,Site,Year,PlantIndividual,OTU,GenusSpecies)%>%
  summarise(abundance=n())
as.data.frame(dat2)
data.frame(dat2$HostPlant,dat2$PlantIndividual,dat2$Year,dat2$OTU,dat2$abundance)

#make wide
dat3<-dat2%>%
  ungroup()%>%
  unite(PlantIndividualYear,PlantIndividual,Year,remove=F)%>%
  unite(HostPlantSite,HostPlant,Site,remove=F)%>%
  dplyr::select(-GenusSpecies)%>%
  spread(OTU,abundance,fill=0)
dat.comm<-data.frame(dat3[,7:62])
row.names(dat.comm)<-dat3$PlantIndividualYear


##### Phylogenetic tree #####
#note - the uncleaned tree in this directory contains the genus species names (query sequence) attached to the OTU numbers
tree<-read.newick("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/FiguresStats/FarrerTBAScleaned5/Farrertreecleaned.nwk")

plot(tree)





##### Faith's Phylogenetic distance #####
PD<-pd(dat.comm,tree)

dat4<-data.frame(dat3[,1:6],PD,dat3[,7:62])





##### MPD by plantindividual #####
phydist <- cophenetic(tree)
ses.mpd.result.notweighted <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=FALSE, runs=999) #takes 5 min with 999
ses.mpd.result.notweighted
ses.mpd.result.notweighted$PlantIndividualYear<-rownames(ses.mpd.result.notweighted)
ses.mpd.result.notweighted1<-ses.mpd.result.notweighted%>%
  dplyr::select(PlantIndividualYear,mpd.obs.z)%>%
  rename(mpd.obs.z.notweighted=mpd.obs.z)

ses.mpd.result.weighted <- ses.mpd(dat.comm, phydist, null.model="taxa.labels",abundance.weighted=TRUE, runs=999) #takes 5 min with 999
ses.mpd.result.weighted
ses.mpd.result.weighted$PlantIndividualYear<-rownames(ses.mpd.result.weighted)
ses.mpd.result.weighted1<-ses.mpd.result.weighted%>%
  dplyr::select(PlantIndividualYear,mpd.obs.z)%>%
  rename(mpd.obs.z.weighted=mpd.obs.z)

dat5<-dat4%>%
  full_join(ses.mpd.result.notweighted1)%>%
  full_join(ses.mpd.result.weighted1)%>%
  mutate(HostPlant=recode(HostPlant,"Spartina alterniflora"= "S. alterniflora","Spartina patens"="S. patens","Phragmites australis"="P. australis","Juncus roemerianus"="J. roemerianus","Sagittaria lancifolia"="S. lancifolia"))
  
  
dat6<-data.frame(dat5[,1:8],dat5[,65:66],dat5[,9:64])
head(dat6)
dim(dat6)




##### phyloseq object ####
otus<-dat6[,11:66]
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

#Exporing dat6 for EDI archiving
#add sample ID so it can link to the MPD distance matrix
dat6withids<-data.frame(sampleID=rownames(sample_data(datp)),dat6)
write.csv(dat6withids,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/Revision2/EDIdata/RichnessMDPAbund.csv",row.names = F)



###### Blasting to unite ######
#unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAMarsh/Survey/Stats/Gradient/QIIME2/sh_general_release_dynamic_s_04.02.2020.fasta" #

unite.ref <- "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAMarsh/Culturing/FiguresStats/LAmarshCulture2/sh_general_release_dynamic_s_04.04.2024.fasta" #sh_general_release_dynamic_s_29.11.2022.fasta

#to look at all or random sequences, not just the 56 otu rep set. can get sequences in report_placed_qrydata.csv

#I got this data from assignments_report_nodupsWG5WGW7J.csv
my56otus<-read.csv("otus56.csv",header=T);rownames(my56otus)<-my56otus$abundance

taxa<-assignTaxonomy(my56otus, unite.ref, multithread = TRUE, minBoot=70, tryRC = TRUE,outputBootstraps=T) #was minBoot=70
taxaonly<-data.frame(taxa$tax);rownames(taxaonly)<-1:dim(taxa$tax)[1]
taxaboot<-data.frame(taxa$boot);rownames(taxaboot)<-1:dim(taxa$tax)[1]
taxaonly
taxaboot
rownames(taxaonly)<-rownames(my56otus)
genusspecies<-data.frame(otu=rownames(my56otus),genusspecies=paste(gsub("^.*?__","",taxaonly[,"Genus"]),gsub("^.*?__","",taxaonly[,"Species"])),genusboot=taxaboot[,"Genus"],speciesboot=taxaboot[,"Species"])

sort(unique(genusspecies$genusspecies))
sort(unique(taxaonly$Phylum))
sort(unique(taxaonly$Order))
sort(unique(taxaonly$Family))
sort(unique(taxaonly$Genus))
sort(unique(taxaonly$Species)) #there are two spartinae species, so this count is 29 so I need to add 1 so 30 species

taxaonly%>%group_by(Phylum)%>%tally()

phragspartina<-dat6%>%
  filter(HostPlant%in%c('P. australis','S. patens'))

#Exporting data for EDI archiving
TaxaAndSeqs<-data.frame(OTU=rownames(taxaonly),taxaonly,sequence=rownames(taxa$tax))
write.csv(TaxaAndSeqs,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Manuscripts/Revision2/EDIdata/TaxaAndSeqs.csv",row.names = F)



##### Reading in isolate data #####

#I just added a "CultureID" column in the csv file by duplicating "Culture ID code". I didn't see any typos. nothing was changed.
isolatefile2017<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Claire_project_fungi.csv",stringsAsFactors = T)
isolatefile2017b<-isolatefile2017%>%
  dplyr::select(PlantIndividual,Species,Site)
isolatefile2017b$Year<-2017
head(isolatefile2017b)
levels(isolatefile2017b$Species)
levels(isolatefile2017b$Site)
levels(isolatefile2017b$PlantIndividual)

#for 2018 I cleaned up the "culture ID code" column in the .csv file as much as I could and called it cultureID. Some of the things I cleaned were Saglan changed to SagLan, JuneRoe changed to JunRoe. there maybe be some cultureIDs that may have been renamed between this voucher file and the master file, for example i'm not sure if we changed the name (added a -1 or -2) when we had the same name repeated in 2017 and 2018. I'm leaving it as is for now. i decided not touse this since it wasn't the smaples that were sent for sequecing, this includes weird pcr fails.
#isolatefile2018<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/Vouchers_2018.csv",stringsAsFactors = T)
#head(isolatefile2018)
#isolatefile2018b<-isolatefile2018%>%
#  dplyr::select(CultureID,Plant,Site)
#isolatefile2018b$Year<-2018
#colnames(isolatefile2018b)[2]<-"Species"
#isolatefiletot<-rbind(isolatefile2017b,isolatefile2018b)

head(culturefile2)

isolatefile2018<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Culturing/isolates2018.csv",stringsAsFactors = T)
head(isolatefile2018)
levels(isolatefile2018$Species)
levels(isolatefile2018$Site)
levels(isolatefile2018$PlantIndividual)
isolatefile2018b<-isolatefile2018%>%
  dplyr::select(PlantIndividual, Species, Site, Year)

isolatefiletot2<-rbind(isolatefile2017b,isolatefile2018b)
dim(isolatefiletot2)
isolatefiletot2$PlantIndividualYear<-paste(isolatefiletot2$PlantIndividual,isolatefiletot2$Year,sep="_")
isolatefiletot2$Site<-factor(isolatefiletot2$Site,levels=c("Turtle cove","CERF","Lumcon"))

#by individual
mi<-isolatefiletot2%>%
  group_by(PlantIndividualYear,Site,Species)%>%
  summarise(isolates=n())
mi
mi2<-mi%>%
  group_by(Site,Species)%>%
  summarise(mean=mean(isolates),se=std.error(isolates),n=n())
mi2

#totals by site
mi<-isolatefiletot2%>%
  group_by(Site)%>%
  summarise(isolates=n())
mi


#totals by species
mi<-isolatefiletot2%>%
  group_by(Species)%>%
  summarise(isolates=n())
mi


##### FunGuild #####

#Note I did not redo this for the 97% clustering since the 99% clustering didn't work

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

phragspartina<-dat7%>%
  filter(HostPlant%in%c('P. australis','S. patens'))



##### git hub token stuff #####
install.packages("gitcreds")
library(gitcreds)
gitcreds_set()
 #first when it asks to enter password or token I put my computer password
 #then do gitcreds_set() again and select 2, then paste my token
#Note: use usename (email) and token, when RStudio wants the github password
