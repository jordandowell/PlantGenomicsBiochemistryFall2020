#install and library required packaged


if(!require("gdsfmt")){
  BiocManager::install("gdsfmt")
  stopifnot(require("gdsfmt"))
}
library(gdsfmt)
# Loading required package: gdsfmt
if(!require("SNPRelate")){
  BiocManager::install("SNPRelate")
  stopifnot(require("SNPRelate"))
}
if(!require("GGally")){
  install.packages("GGally")
  stopifnot(require("GGally"))
}
library("GGally")
library(SNPRelate)
library(ggplot2)
library(SNPRelate)
library(gdsfmt)

# starting off recoding  our tped, & tfam files into bed files via plink


#give your computer permission to run external programs
# this may require you to do some googling on your part if this does not work.
#kewords: command line grant access folder binary
system("chmod -R 755 ../PlantGenomicsBiochemistryFall2020") # grant permissions (to avoid access denied errors when tyring to run the external binary)


#change NEWSNPs to whatever file you are interested 
#this needs to be the precursor characters for the .tped and .tfam files
getwd()
NEWSNPS<-"Global_SNPlist"
#the following is using the plink software to recode the files into more usable forms
system(paste("data/Software/plink --tfile",paste0("data/datasetsScript2/", NEWSNPS),"--covar",paste0("data/datasetsScript2/", NEWSNPS,"_COVARIATES.txt") , "--make-bed", "--allow-no-sex --allow-extra-chr", "--out ",paste0("data/datasetsScript2/", NEWSNPS)))



#Investigate the contents of the PLINK files .fam and .bim files
#Count number of individuals in the study

#read in fam file to get the names of the individuals in the dataset
FAM<-read.table(file=paste0("data/datasetsScript2/",NEWSNPS,".fam"),sep=" ", header=FALSE)

#view top entries
#column 6 referees to a phenotype measurement that will be used in another tutorial
head(FAM)

#1st dimensions is number of individuals
dim(FAM)


#Read through the map information on the number of SNPs

map<- read.table(file=paste0("data/datasetsScript2/",NEWSNPS,".bim"),sep=" ", header=FALSE)
#number of SNPs
dim(map)


#read in file information on breeding groups of samples
BreedingGroupInfo<-read.csv(file="data/datasetsScript2/SAM_Metadata.csv",header=TRUE)

#get information on the number of sample individuals from each breeding group
table(BreedingGroupInfo$GROUP)


#check to make sure the order of individuals is the same as the PLINK file
#following line translated is:
#what is the sum of the breeding group column PPN that is not also in FAM column V1
#if they match the number should be 0 
sum(BreedingGroupInfo$PPN!=FAM$V1)


#Convert bed bim and fam files into a .gds file this helps speed up calculations in large datasets
#name files to be used
bedfile<-paste0("data/datasetsScript2/",NEWSNPS,".bed")
bimfile<-paste0("data/datasetsScript2/",NEWSNPS,".bim")
famfile<-paste0("data/datasetsScript2/",NEWSNPS,".fam")
#place the above files in the function and name the output
snpgdsBED2GDS(bedfile,famfile,bimfile,cvt.chr="char", out.gdsfn = paste0("data/datasetsScript2/",NEWSNPS,".gds"))

#open the .gds file
genofile <- snpgdsOpen(paste0("data/datasetsScript2/",NEWSNPS,".gds"))
#get a visual of the file structure
head(genofile)


#if you need to reopen the gds file for some reason make sure to close the gds file un comment the following line
#closefn.gds(genofile)


#we can see the sample id or snp id
#remove head to see the full list. the head file shows the first 6 entries
head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))


#extract genotype for a specific SNP this can be changed to any SNP.id you would like to see
g <- snpgdsGetGeno(genofile, snp.id="Ha412HOChr01:59890719")
# Genotype matrix: 261 samples X 1 SNPs
# view a histogram of the alleles at that SNP 0 & 2 are homozygous, 1 is missing or heterozygous
hist(g)







#perform a PCA for samples using all pre Selected SNPs. 
#autosome is set to false because we are using character strings for our chromosome names
pca <- snpgdsPCA(genofile,autosome.only = F)



#incorporate the population membership in the PCA plot. This must be categorical in the following example

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
population<-BreedingGroupInfo$BREED
SpecificBreedinggroup<-BreedingGroupInfo$GROUP
Core12<-BreedingGroupInfo$CORE12


#create a data frame to use for plotting 
#each line below is representative of a column

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(population)[match(pca$sample.id, sample.id)],
                  Specificpop = factor(SpecificBreedinggroup)[match(pca$sample.id, sample.id)],
                  Core = factor(Core12,levels = c("NON-CORE","CORE"))[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca2$eigenvect[,3],# the third eigenvector
                  EV4 = pca2$eigenvect[,4],# the fourth eigenvector
                  stringsAsFactors = FALSE)
#view the resulting data frame
head(tab)


#calculate variance explained
#round percent variation explained to 100ths
pc.percent <- round(pca$varprop*100,2)
head((pc.percent))


#make a single PCA plot
# a warning will follow about using a discrete variable for size ignore it! 
#traditionally you want a continuous variable related to size
#if you want to change the color to the more specific population description change color= pop to color = Specificpop
SinglePCA_Plot<-ggplot(tab,
                       aes(x=EV1,y=EV2,color= pop, size=Core))+
  xlab(paste0("PC 1 (",pc.percent[1],"% Variance Explained)"))+
  ylab(paste0("PC 2 (",pc.percent[2],"% Variance Explained)"))+
  geom_point()+theme_linedraw()

#note if you want to view the plot highlight just the plot name and run it the plot can be exported through the GUI


#Now make scatter plots of the top 4 PCs with proportional variance explained included


#pairwise scatter plots of the first 4 components
PairwisePCA_plots<-ggpairs(tab, columns = c("EV1","EV2","EV3","EV4"),
                           aes(color=pop),
                           columnLabels = c(paste0("PC 1 (",pc.percent[1],"% VE)"),
                                            paste0("PC 2 (",pc.percent[2],"% VE)"),
                                            paste0("PC 3 (",pc.percent[3],"% VE)"),
                                            paste0("PC 4 (",pc.percent[4],"% VE)")),
                           axisLabels = 'none',
                           diag = list(continuous=wrap("densityDiag", alpha=0.2)),
                           lower = list(continuous="blank"),
                           upper = list(continuous='points')) +theme_linedraw()
#the following lines add the size changes for the core 12 in the plot
#subset locations to add size for just scatter plots

xpositions<-c(1,1,1,2,2,3)
ypositions<-c(2,3,4,3,4,4)
#for loop to add size 
for (i in 1:length(xpositions)) {
  print(xpositions[i])
  print(ypositions[i])
  
  PairwisePCA_plots4<-PairwisePCA_plots[xpositions[i],ypositions[i]]
  PairwisePCA_plots[xpositions[i],ypositions[i]]<-PairwisePCA_plots4+aes(size=Core)
}

#end of plot


########################Refining the SNP list based on linkage disequilibrium & autocorrelation


#often autocorrelation of SNPs or linkage disequilibrium can alter PCA results
#we can retry our analysis with a reduced set of SNPs based on a threshold for autocorrelation

#identify subsets of SNPs based on LD threshold of 0.2

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)

#the marker dataset may be significantly reduced
snpset.id <- unlist(snpset)

#Now perform a PCA with the subset of SNPs

pca2 <- snpgdsPCA(genofile, snp.id=snpset.id, autosome.only = F)
#create that dataframe for plotting 
tab2 <- data.frame(sample.id = pca2$sample.id,
                   Core = factor(Core12,levels = c("NON-CORE","CORE"))[match(pca2$sample.id, sample.id)],
                   pop = factor(population)[match(pca2$sample.id, sample.id)],
                   Specificpop = factor(SpecificBreedinggroup)[match(pca$sample.id, sample.id)],
                   EV1 = pca2$eigenvect[,1],    # the first eigenvector
                   EV2 = pca2$eigenvect[,2],    # the second eigenvector
                   EV3 = pca2$eigenvect[,3],# the third eigenvector
                     EV4 = pca2$eigenvect[,4],# the fourth eigenvector
                   stringsAsFactors = FALSE)

#round percent variation explained to 100ths
pc.percent <- round(pca2$varprop*100,2)
head(round(pc.percent, 2))

#make a single PCA plot
SinglePCA_Plot<-ggplot(tab2,
       aes(x=EV1,y=EV2,color= pop, size=Core))+
  xlab(paste0("PC 1 (",pc.percent[1],"% Variance Explained)"))+
  ylab(paste0("PC 2 (",pc.percent[2],"% Variance Explained)"))+
  geom_point()+theme_linedraw()



#pairwise scatter plots of the first 4 components
PairwisePCA_plots<-ggpairs(tab2, columns = c("EV1","EV2","EV3","EV4"),
               aes(color=pop),
        columnLabels = c(paste0("PC 1 (",pc.percent[1],"% VE)"),
                         paste0("PC 2 (",pc.percent[2],"% VE)"),
                         paste0("PC 3 (",pc.percent[3],"% VE)"),
                         paste0("PC 4 (",pc.percent[4],"% VE)")),
        axisLabels = 'none',
       diag = list(continuous=wrap("densityDiag", alpha=0.2)),
        lower = list(continuous="blank"),
        upper = list(continuous='points')) +theme_linedraw()
#the following lines add the size changes for the core 12 in the plot
#subset locations to add size for just scatter plots

xpositions<-c(1,1,1,2,2,3)
ypositions<-c(2,3,4,3,4,4)
#for loop to add size 
for (i in 1:length(xpositions)) {
  print(xpositions[i])
  print(ypositions[i])
  
  PairwisePCA_plots4<-PairwisePCA_plots[xpositions[i],ypositions[i]]
  PairwisePCA_plots[xpositions[i],ypositions[i]]<-PairwisePCA_plots4+aes(size=Core)
}
#end of plot


#######now lets compare our results to the population structure of H. annuus. 
#Population structure data taken from Huber et al. 2018 (https://www.nature.com/articles/s41477-018-0329-0)
 

HannuusPopstructure<-read.table(file=paste0("data/datasetsScript2/",NEWSNPS,".cov"),sep=" ", header=TRUE)
#take a look the data set should be very similar to what we have been using so far
View(HannuusPopstructure)



# 
#create that dataframe for plotting adding in the meta data we have been using 
tab3 <- data.frame(sample.id = HannuusPopstructure$FID,
                   Core = factor(Core12,levels = c("NON-CORE","CORE"))[match(HannuusPopstructure$FID, sample.id)],
                   pop = factor(population)[match(HannuusPopstructure$FID, sample.id)],
                   Specificpop = factor(SpecificBreedinggroup)[match(HannuusPopstructure$FID, sample.id)],
                   COV1 =HannuusPopstructure[,3],    # the first eigenvector
                   COV2 = HannuusPopstructure[,4],    # the second eigenvector
                   COV3 = HannuusPopstructure[,5],# the third eigenvector
                   COV4 = HannuusPopstructure[,6],# the fourth eigenvector
                   stringsAsFactors = FALSE)




#make plots.
#make a single PCA plot change color or size based on columns of interest
SinglePCA_Plot<-ggplot(tab3,
                       aes(x=COV1,y=COV2,color= pop, size=Core))+
  xlab(paste0("Population Structure PC1 16.9% VE"))+
  ylab(paste0("Population Structure PC2 14.5 VE"))+
  geom_point()+theme_linedraw()

#the following lines add the size changes for the core 12 in the plot
#pairwise scatter plots of the first 4 components
PairwisePCA_plots<-ggpairs(tab3, columns = c("COV1","COV2","COV3","COV4"),
                           aes(color=pop),
                           columnLabels = c(paste0("Population Structure PC1 16.9% VE"),
                                            paste0("Population Structure PC2 14.5 VE"),
                                            paste0("Population Structure PC3"),
                                            paste0("Population Structure PC4")),
                           axisLabels = 'none',
                           diag = list(continuous=wrap("densityDiag", alpha=0.2)),
                           lower = list(continuous="blank"),
                           upper = list(continuous='points')) +theme_linedraw()
#subset locations to add size for just scatter plots

xpositions<-c(1,1,1,2,2,3)
ypositions<-c(2,3,4,3,4,4)
#for loop to add size 
for (i in 1:length(xpositions)) {
  print(xpositions[i])
  print(ypositions[i])
  
  PairwisePCA_plots4<-PairwisePCA_plots[xpositions[i],ypositions[i]]
  PairwisePCA_plots[xpositions[i],ypositions[i]]<-PairwisePCA_plots4+aes(size=Core)
}

#end of plot


#following line can be used to test for correlations between different columnts and datasets
# use ?cor.test to get more parameters and a better description of the function.
# tab$EV4 translated is dataset "tab" column "EV4"...... "$" is used to access columns within the dataset 

cor.test(tab$EV4,tab3$COV2)

