#Perry, KI. 2021. Using a functional trait approach to study insect community ecology. Virtual North
#Central Branch Meeting, Entomological Society of America, June 22.

#https://cran.r-project.org/web/packages/FD/FD.pdf

#Set working directory
setwd("~/The Ohio State University/Thesis/Data/R/PNR/Tornado")

#Load datasets, one with species abundances and one with species functional traits
t <-read.csv("Beetle2015_Traits.csv", row.names=1)
a <-read.csv("Beetle2015_Abund.csv", row.names=1)

#Convert abundance data to relative abundance
install.packages("vegan")
library(vegan)
citation("vegan")

#Check total abundance in each sample (i.e. row)
apply(a, 1, sum)

#Convert to relative abundance
ra <- decostand(a, method = "total")

#Check the totals for each row in the new dataset
#All rows should equal 1
apply(ra, 1, sum)

#Check the structure of the dataset
str(t)
str(ra)

#Plot the traits to see if any are redundant
plot(t)
cor(t, method = c("pearson"))

#head width, eye width, and elytra length appear correlated with other traits
#remove redundant traits
t2 <- t[,-2]
t2 <- t2[,-5]
t2 <- t2[,-3]

#Plot the traits again to see if any more need to be removed
plot(t2)
cor(t2)

install.packages("FD")
library(FD)
citation("FD")
#function dbFD

####################################################################################
#community-weighted mean trait values
#works best on continuous traits

#remove dispersal ability (categorical trait) from dataset
t3 <- t2[,-5]

#calculate CWMs for continuous variables
cwm <- dbFD(t3[colnames(a), ], a, w.abun = T, stand.x = T)$CWM
cwm

#Error: Error in dbFD(t3[colnames(a), ], a, w.abun = T, stand.x = T) : 
#Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method
#There may be additional redundancies within the dataset
#sqrt is the default correction method;others include "cailliez" and "lingoes" -- see FD package and Laliberte and Legendre 2010
#stand.x standardizes all traits to mean 0, but traits must be numeric
cwm <- dbFD(t3[colnames(a), ], a, w.abun = T, stand.x = T, corr = c("cailliez"))$CWM
cwm
write.csv(cwm, file = "BeetleCWM.csv")

#CWM values for each site can be used in other analyses to test for treatment effects


#####################################################################################
#Functional diversity indices
#If dataset has a categorical variable, then calculate a gower dissimilarity matrix
#Indices will only calculate value if community has at least 3 species
#Communities that have fewer than 3 species, output will be NA

##Gower dissimilarity matrix
td<-as.matrix(gowdis(t2))

#Functional richness
fric <- dbFD(td[colnames(a), ], a, w.abun = T, stand.x = T)$FRic
fric

#Functional evenness
feve <- dbFD(td[colnames(a), ], a, w.abun = T, stand.x = T)$FEve
feve

#Functional divergence
fdiv <- dbFD(td[colnames(a), ], a, w.abun = T, stand.x = T)$FDiv
fdiv

#Functional dispersion
fdis <- dbFD(td[colnames(a), ], a, w.abun = T, stand.x = T)$FDis
fdis

outputs <- as.data.frame(cbind(fric, feve, fdiv, fdis))
plot(outputs)
write.csv(outputs, file = "FDiversity_Indices.csv")

#Functional diversity indices can be used in other analyses to test for treatment effects

