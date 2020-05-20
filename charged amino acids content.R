# load the libraries required in the script

library(readr)
library(stringr)

# define the functions used in the script

composition<-function(dataset,aa){
  # takes a data frame which include a column 'Sequence' and a list of amino acids 
  # to count the percentage of within the sequence
  # returns a vector with amino acids composition
  sapply(aa,function(a,s){round(str_count(s,a)/str_length(s)*100,2)},dataset$Sequence)}

# import the dataset of all human proteins

human_proteins<-read_csv("human_proteins.csv")

# subset the dataset in two groups, depending on whether the 'Subcellular location' annotation 
# includes the word membrane or not

Membrane<-human_proteins[grep("membrane",human_proteins$`Subcellular location`,ignore.case=TRUE),]
Soluble<-human_proteins[grep("membrane",human_proteins$`Subcellular location`,ignore.case=TRUE,invert=TRUE),]

# list of the charged amino acids

charged<-c("R","H","K","D","E")

# compute vectors with percentage of charged amino acids of all proteins 
# included in each subset

membrane_proteins_charge<-apply(composition(Membrane,charged),1,sum)
soluble_proteins_charge<-apply(composition(Soluble,charged),1,sum)

# plot the distribution of the two vectors and compute their statistics

boxplot(membrane_proteins_charge,soluble_proteins_charge,names=c("Membrane proteins","Soluble proteins"),ylab="charged amino acids (%)",outline=FALSE)
summary(membrane_proteins_charge)
summary(soluble_proteins_charge)

# perform a Mann–Whitney U-test to verify whether the difference of charged amino acids 
# percentage between the two groups is significant

wilcox.test(membrane_proteins_charge,soluble_proteins_charge,conf.int=TRUE)
