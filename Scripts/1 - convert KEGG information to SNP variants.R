#install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ape", quietly = TRUE))
  install.packages("ape")

BiocManager::install("KEGGgraph")
BiocManager::install("rtracklayer")


library(ape)
library("KEGGgraph")
library("rtracklayer")
#import KGML file
#make sure you are in your project directory
#change file name to your xml file you downloads
TerpenoidKGML <- "data/ec00906.xml"
Terpenoid <- parseKGML(TerpenoidKGML)


#get a list of enzymes 
#identify the number of nodes in the pathway
#nodes can be compounds reactions or enzymes
numberofnodes<-length(Terpenoid@nodes)

#run a for loop to get a list of enzyme names
i<-1
EnzymeList<-list()
for (i in 1:numberofnodes) {
  if(Terpenoid@nodes[[i]]@type %in% "enzyme"){
    EnzymeList<-c(EnzymeList,as.character(Terpenoid@nodes[[i]]@name))
    # EnzymeList[length(EnzymeList)+1,]<-Enzymetobeadded
  }
}


#remove EC: identifies for easier searches
EnzymeList<-sub('...', '', EnzymeList)

#the list of enzymes is ready for use now!

#read in annotation data



SUNFLOWERANNOTATION<-readGFF("data/SUNFLOWER_ANNOTATION_Ha412HOv2.0-20181130.gff3")

#get the number of total annotated genes
SUNFLOWERTOTALGENES<- length(SUNFLOWERANNOTATION@listData$ec_number)
#Get indexes of enzymes that match out ec_numbers


SUNFLOWERENZYMES<-data.frame()

#this loop make take ~2hours  to run its checking every gene
for (i in 1:SUNFLOWERTOTALGENES) {
  #if statement checking that the length of the ec # is not equal to 0 
  if(length(SUNFLOWERANNOTATION@listData$ec_number[i][[1]])!=0){
    #nested if statement checking if the ec_number is in our enzyme list
    if(any(SUNFLOWERANNOTATION@listData$ec_number[i][[1]]%in%EnzymeList)){
      #check if there is more than one enzyme annotations if so add an index for each enzyme
      if(length(SUNFLOWERANNOTATION@listData$ec_number[i][[1]])>1){
        for (j in 1:length(SUNFLOWERANNOTATION@listData$ec_number[i][[1]])) {
          newrow<- c(i,SUNFLOWERANNOTATION@listData$ec_number[i][[1]][j],SUNFLOWERANNOTATION@listData$ID[i][[1]])
          SUNFLOWERENZYMES<-rbind(SUNFLOWERENZYMES,newrow)
        }
      }else{
        newrow<- c(i,SUNFLOWERANNOTATION@listData$ec_number[i][[1]],SUNFLOWERANNOTATION@listData$ID[i][[1]])
        SUNFLOWERENZYMES<-rbind(SUNFLOWERENZYMES,newrow)
      }
      
      print(i)
      print(SUNFLOWERANNOTATION@listData$ec_number[i][[1]])
    }
  }
}




#here we are using a loop to get gene model indexes 



i<-1 
#use for loop to get genemodel indexes
SUNFLOWERGENEMODELS<-data.frame()

for (i in 1:SUNFLOWERTOTALGENES) {
  #if any parent gene matches previously identifed mrna
  if(any(SUNFLOWERANNOTATION@listData$Parent[i][[1]]%in%SUNFLOWERENZYMES[,3])){
    
    #get index, parent gene,IDs
    newrow<- c(i,SUNFLOWERANNOTATION@listData$Parent[i][[1]],SUNFLOWERANNOTATION@listData$ID[i][[1]])
    SUNFLOWERGENEMODELS<-rbind(SUNFLOWERGENEMODELS,newrow)
    print(i)
    print(SUNFLOWERANNOTATION@listData$ID[i][[1]])
  }
}





#combine Genemodels with parent Genes & add EC for each model

SUNFLOWERGENEMODELS_1<-data.frame()

for (i in 1:nrow(SUNFLOWERGENEMODELS)){
  
  newrow<-c(SUNFLOWERGENEMODELS[i,1],SUNFLOWERENZYMES[match(SUNFLOWERGENEMODELS[i,2],SUNFLOWERENZYMES[,3]),2],SUNFLOWERGENEMODELS[i,3])
  SUNFLOWERGENEMODELS_1<-rbind(SUNFLOWERGENEMODELS_1,newrow)
  print(i)
}



#bind enzymes and genemodels

SUNFLOWERENZYMES_GENEMODELS<- rbind(SUNFLOWERENZYMES,setNames(SUNFLOWERGENEMODELS_1,names(SUNFLOWERENZYMES)))



#read in the annotation as a differnt type of data frame that is easier to access
SUNFLOWERANNOTATION<-read.gff("data/SUNFLOWER_ANNOTATION_Ha412HOv2.0-20181130.gff3")



#subset the annotation based of the enzymes & gene models you want
SUNFLOWERANNOTATION_Subset<- SUNFLOWERANNOTATION[as.numeric(SUNFLOWERENZYMES_GENEMODELS[,1]),]

#there should be no gene,rRNA,ncRNA or tRNA 
summary(SUNFLOWERANNOTATION_Subset$type)



#import map file to prune SNPS of interest

#we are making an assumption that the SNPlist and SNPvariants are in the same order
SNPLIST<-read.table("data/SNPLIST_XRQv1_412_261_filtered.map", )
SNPVARIANTS<- read.table("data/SNP_VARIANTS_XRQv1_412_261_filtered.tped")





#create data output files
dir.create("Output")
#create a directory for enzyme specific outputs
dir.create("Output/Enzymes")
#move to Enzymes folder 
setwd("Output/Enzymes")


#create directory 1 for each Enzyme type

EnzymeTypes<-unique(SUNFLOWERENZYMES_GENEMODELS[,2])

#create a directory for genetypes

for (i in 1:length(EnzymeTypes)) {
  dir.create(path=paste(EnzymeTypes[i]))
}

#for loop to subset genes based on type and write csv into file with EC_type as name

#list folders present in the output


Directorylist<-list.dirs(full.names = F,recursive = F)

#the Direstory list is also equal to the Enzyme Types
#loops across Enzyme types and genemodel types 
for (i in 1:length(EnzymeTypes)) {
  print(EnzymeTypes[i])
  
  
  
  
  
  #get the Enzyme
  Enzymetypeofinterest<-EnzymeTypes[i]
  
  #get the indicies of the annotation based on your Enzyme type
  
  SUNFLOWERENZYMES_GENEMODELS_ofinterest<-subset(SUNFLOWERENZYMES_GENEMODELS, SUNFLOWERENZYMES_GENEMODELS[,2] %in% Enzymetypeofinterest)
  
  #subset the annotation file based on the indicies of enzymes
  SUNFLOWERANNOTATION_Enzyme_Subset<- SUNFLOWERANNOTATION[as.numeric(SUNFLOWERENZYMES_GENEMODELS_ofinterest[,1]),]
  
  
  #for loop to cycle across chromosomes @@@@@
  
  j<-1
  Genemodel_SNPLIST<-data.frame()
  Genemodel_VARIANTLIST<-data.frame()
  for (j in 1:length(unique(SUNFLOWERANNOTATION_Enzyme_Subset[,1]))) {
    
    
    #subset by chromosome
    chromosome<-unique(SUNFLOWERANNOTATION_Enzyme_Subset[,1])[j]
    #use the annotation to subset based on the chromosome
    chromosomesetofinterest<-subset(SUNFLOWERANNOTATION_Enzyme_Subset, SUNFLOWERANNOTATION_Enzyme_Subset[,1]%in%chromosome)
    #subset the SNPLIST based on chromsome(SNPLIST is a smaller data file than SNPvariant which is why we are using it)
    chromosomeSNPLISTofinterest<-subset(SNPLIST,SNPLIST[,1]%in%chromosome)
    
    
    #get SNPS that lie with in the range defined by the annotation, on the chromosome of interest
    #create a list of possible positions
    Possiblepositions<-c()
    
    for (L in 1:nrow(chromosomesetofinterest)) {
      newlist<-chromosomesetofinterest[L,4]:chromosomesetofinterest[L,5]
      Possiblepositions<-c(Possiblepositions,newlist)
    }
    #subset the SNP list if the position is %in% one of the possible positions
    subset_Intergenic_SNPLIST<- subset(chromosomeSNPLISTofinterest,chromosomeSNPLISTofinterest[,4]%in% Possiblepositions)
    
    
    #get SNP variants
    #we are using the row indexes from out SNPlist
    
    subset_Intergenic_VARIANTLIST<- SNPVARIANTS[as.numeric(row.names(subset_Intergenic_SNPLIST)),]
    #make sure the names of SNPS are equal before proceeding 
    #this will need to be an if else statement before adding it to the dataframe for writing 
    print(chromosome)
    Genemodel_SNPLIST<-rbind(Genemodel_SNPLIST,subset_Intergenic_SNPLIST)
    Genemodel_VARIANTLIST<-rbind(Genemodel_VARIANTLIST,subset_Intergenic_VARIANTLIST)
  }
  #write .tped & .map file so they can be used in other analyses
  write.table(Genemodel_SNPLIST, paste0(Enzymetypeofinterest,"/",Enzymetypeofinterest,".map"), col.names = F,row.names = F,quote = F, sep = "\t")
  write.table(Genemodel_VARIANTLIST, paste0(Enzymetypeofinterest,"/",Enzymetypeofinterest,".tped"), col.names = F,row.names = F,quote = F)
  
  
  
  
  
  
  #create directory 1 for each genemodel type in each enzyme files 
  Genetypes<-unique(SUNFLOWERANNOTATION_Enzyme_Subset[,3])
  
  #create a directory for genetypes
  
  for (i in 1:length(Genetypes)) {
    dir.create(path=paste0(Enzymetypeofinterest,"/",Genetypes[i]))
  }
  
  #this for loop cycles over gene model types
  i<-1
  for (i in 1:length(Genetypes)) {
    print(Genetypes[i])
    
    #get the genetype
    Genetypeofinterest<-Genetypes[i]
    
    
    
    #subset the annotation based on your gene type
    EnzymeGenemodel_Genesetofinterest<-subset(SUNFLOWERANNOTATION_Enzyme_Subset,SUNFLOWERANNOTATION_Enzyme_Subset[,3]%in%Genetypeofinterest)
    
    
    #for loop to cycle across chromosomes 
    
    j<-1
    Genemodel_SNPLIST<-data.frame()
    Genemodel_VARIANTLIST<-data.frame()
    
    #fix for assessing chromosome level SNPs
    chromosomes<-c(unique(as.character(SNPLIST[,1])), unique(as.character(EnzymeGenemodel_Genesetofinterest[,1])))
    chromosomes_present<-chromosomes[duplicated(chromosomes)]
    
    
    
    
    
    for (j in 1:length(unique(chromosomes_present))) {
      
      
      #subset by chromosome
      chromosome<-chromosomes_present[j]
      #use the annotation to subset based on the chromosome
      chromosomesetofinterest<-subset(EnzymeGenemodel_Genesetofinterest, EnzymeGenemodel_Genesetofinterest[,1]%in%chromosome)
      #subset the SNPLIST based on chromsome(SNPLIST is a smaller data file than SNPvariant which is why we are using it)
      chromosomeSNPLISTofinterest<-subset(SNPLIST,SNPLIST[,1]%in%chromosome)
      
      
      #get SNPS that lie with in the range defined by the annotation, on the chromosome of interest
      #create a list of possible positions
      Possiblepositions<-c()
      
      for (L in 1:nrow(chromosomesetofinterest)) {
        newlist<-chromosomesetofinterest[L,4]:chromosomesetofinterest[L,5]
        Possiblepositions<-c(Possiblepositions,newlist)
      }
      #subset the SNP list if the position is %in% one of the possible positions
      subset_Intergenic_SNPLIST<- subset(chromosomeSNPLISTofinterest,chromosomeSNPLISTofinterest[,4]%in% Possiblepositions)
      
      
      #get SNP variants
      #we are using the row indexes from out SNPlist
      
      subset_Intergenic_VARIANTLIST<- SNPVARIANTS[as.numeric(row.names(subset_Intergenic_SNPLIST)),]
      #make sure the names of SNPS are equal before proceeding 
      #this will need to be an if else statement before adding it to the dataframe for writing 
      print(chromosome)
      Genemodel_SNPLIST<-rbind(Genemodel_SNPLIST,subset_Intergenic_SNPLIST)
      Genemodel_VARIANTLIST<-rbind(Genemodel_VARIANTLIST,subset_Intergenic_VARIANTLIST)
    }
    #write .tped & .map file so they can be used in other analyses
    write.table(Genemodel_SNPLIST, paste0(Enzymetypeofinterest,"/",Genetypeofinterest,"/",Genetypeofinterest,"_",Enzymetypeofinterest,".map"), col.names = F,row.names = F,quote = F, sep = "\t")
    write.table(Genemodel_VARIANTLIST, paste0(Enzymetypeofinterest,"/",Genetypeofinterest,"/",Genetypeofinterest,"_",Enzymetypeofinterest,".tped"), col.names = F,row.names = F,quote = F)
    
  }
  
  
}




#set working directory to output
setwd("..")



#for loop to cycle across chromosomes to create a global list of variants for the pathway

j<-1
Genemodel_SNPLIST<-data.frame()
Genemodel_VARIANTLIST<-data.frame()

#fix for assessing chromosome level SNPs
chromosomes<-c(unique(as.character(SNPLIST[,1])), unique(as.character(SUNFLOWERANNOTATION_Subset[,1])))
chromosomes_present<-chromosomes[duplicated(chromosomes)]


for (j in 1:length(unique(chromosomes_present))) {
  
  
  
  #subset by chromosome
  chromosome<-chromosomes_present[j]
  #use the annotation to subset based on the chromosome
  chromosomesetofinterest<-subset(SUNFLOWERANNOTATION_Subset, SUNFLOWERANNOTATION_Subset[,1]%in%chromosome)
  #subset the SNPLIST based on chromsome(SNPLIST is a smaller data file than SNPvariant which is why we are using it)
  chromosomeSNPLISTofinterest<-subset(SNPLIST,SNPLIST[,1]%in%chromosome)
  
  
  #get SNPS that lie with in the range defined by the annotation, on the chromosome of interest
  #create a list of possible positions
  Possiblepositions<-c()
  
  for (L in 1:nrow(chromosomesetofinterest)) {
    newlist<-chromosomesetofinterest[L,4]:chromosomesetofinterest[L,5]
    Possiblepositions<-c(Possiblepositions,newlist)
  }
  #subset the SNP list if the position is %in% one of the possible positions
  subset_Intergenic_SNPLIST<- subset(chromosomeSNPLISTofinterest,chromosomeSNPLISTofinterest[,4]%in% Possiblepositions)
  
  
  #get SNP variants
  #we are using the row indexes from out SNPlist
  
  subset_Intergenic_VARIANTLIST<- SNPVARIANTS[as.numeric(row.names(subset_Intergenic_SNPLIST)),]
  #make sure the names of SNPS are equal before proceeding 
  #this will need to be an if else statement before adding it to the dataframe for writing 
  print(chromosome)
  Genemodel_SNPLIST<-rbind(Genemodel_SNPLIST,subset_Intergenic_SNPLIST)
  Genemodel_VARIANTLIST<-rbind(Genemodel_VARIANTLIST,subset_Intergenic_VARIANTLIST)
}
#write .tped & .map file so they can be used in other analyses
write.table(Genemodel_SNPLIST, paste0("Global_SNPlist.map"), col.names = F,row.names = F,quote = F, sep = "\t")
write.table(Genemodel_VARIANTLIST, paste0("Global_VARIANTLIST.tped"), col.names = F,row.names = F,quote = F)











#create directory 1 for each genemodel type

Genetypes<-unique(SUNFLOWERANNOTATION_Subset[,3])

#create a directory for genetypes

for (i in 1:length(Genetypes)) {
  dir.create(path=paste(Genetypes[i]))
}

#for loop to subset genes based on type and write csv into file with EC_type as name

#list folders present in the output


Directorylist<-list.dirs(full.names = F,recursive = F)







#set working directory to first of list based on type
#MAKE SURE YOU ARE IN THE OUTPUT DIRECTORY
#setwd("..")
#check the number of genetypes with the summary function and make the for loop only sequence through genetypes with genes
getwd()
#this for loop cycles over gene model types
i<-1
for (i in 1:length(Genetypes)) {
  print(Genetypes[i])
  
  #get the genetype
  Genetypeofinterest<-Genetypes[i]
  
  
  
  #subset the annotation based on your gene type
  Genesetofinterest<-subset(SUNFLOWERANNOTATION_Subset,SUNFLOWERANNOTATION_Subset[,3]%in%Genetypeofinterest)
  
  #for loop to cycle across chromosomes 
  
  j<-1
  Genemodel_SNPLIST<-data.frame()
  Genemodel_VARIANTLIST<-data.frame()
  
  #fix for assessing chromosome level SNPs
  chromosomes<-c(unique(as.character(SNPLIST[,1])), unique(as.character(Genesetofinterest[,1])))
  chromosomes_present<-chromosomes[duplicated(chromosomes)]
  
  for (j in 1:length(unique(chromosomes_present))) {
    
    
    
    #subset by chromosome
    chromosome<-chromosomes_present[j]
    #use the annotation to subset based on the chromosome
    chromosomesetofinterest<-subset(Genesetofinterest, Genesetofinterest[,1]%in%chromosome)
    #subset the SNPLIST based on chromsome(SNPLIST is a smaller data file than SNPvariant which is why we are using it)
    chromosomeSNPLISTofinterest<-subset(SNPLIST,SNPLIST[,1]%in%chromosome)
    
    
    #get SNPS that lie with in the range defined by the annotation, on the chromosome of interest
    #create a list of possible positions
    Possiblepositions<-c()
    
    for (L in 1:nrow(chromosomesetofinterest)) {
      newlist<-chromosomesetofinterest[L,4]:chromosomesetofinterest[L,5]
      Possiblepositions<-c(Possiblepositions,newlist)
    }
    #subset the SNP list if the position is %in% one of the possible positions
    subset_Intergenic_SNPLIST<- subset(chromosomeSNPLISTofinterest,chromosomeSNPLISTofinterest[,4]%in% Possiblepositions)
    
    
    #get SNP variants
    #we are using the row indexes from out SNPlist
    
    subset_Intergenic_VARIANTLIST<- SNPVARIANTS[as.numeric(row.names(subset_Intergenic_SNPLIST)),]
    #make sure the names of SNPS are equal before proceeding 
    #this will need to be an if else statement before adding it to the dataframe for writing 
    print(chromosome)
    Genemodel_SNPLIST<-rbind(Genemodel_SNPLIST,subset_Intergenic_SNPLIST)
    Genemodel_VARIANTLIST<-rbind(Genemodel_VARIANTLIST,subset_Intergenic_VARIANTLIST)
  }
  #write .tped & .map file so they can be used in other analyses
  write.table(Genemodel_SNPLIST, paste0(Genetypeofinterest,"/",Genetypeofinterest,".map"), col.names = F,row.names = F,quote = F, sep = "\t")
  write.table(Genemodel_VARIANTLIST, paste0(Genetypeofinterest,"/",Genetypeofinterest,".tped"), col.names = F,row.names = F,quote = F)
  
}





#ANYTHING PAST HERE IS USELESS
getwd()




