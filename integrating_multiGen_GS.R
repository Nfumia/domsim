library(AlphaSimR)
library(tidyverse)
library(sommer)
rm(list=ls())
# SET-UP FOUNDER POPULATION
founderHap <- runMacs2(nInd=100,nChr=11,segSites=100)
# New global simulation parameters from founder haplotypes
SP <- SimParam$new(founderHap)
SP$restrSegSites(minQtlPerChr=10,minSnpPerChr=10,overlap=FALSE)
SP$addTraitA(nQtlPerChr=10,mean=0,var=1) # Additive trait genetic architecture
SP$setSexes("no") #all individuals are hermaphrodites
SP$addSnpChip(nSnpPerChr=10) # Observed SNPs per chromosome 
SP$setTrackPed(TRUE) #keeps pedigree information in slot SP@pedigree 
SP$setTrackRec(TRUE) #keeps recomb. records of all individuals in slot of "SP"
# New founder pop
founders <- newPop(founderHap,simParam=SP)
# Initial founder phenotypes
founders <- setPheno(pop=founders,h2=0.5,reps=2)
# loop length / number of cycles of selection
nCycles<-10
# very simple container for each cycles sim output
simOutput<-list(founders)
cycle <- 1

gen.mat <- pullSnpGeno(simOutput[[cycle]]) #pull genotypes from offspringPop object
A <- sommer::A.mat(gen.mat) #create A genomic relationship matrix
Y <- as.data.frame(cbind(simOutput[[cycle]]@id,simOutput[[cycle]]@pheno)) #create dataframe of traits and id for model
colnames(Y) <- c("id","trait") #rename columns something meaningful
Y$trait <- as.numeric(Y$trait) #make response (trait) numeric
#Y <- Y %>%
#    mutate(trait =  replace(trait, sample(row_number(),  
#           size = ceiling(0.3 * n()), replace = FALSE), NA) )

#run univariate mixed model with sommer 
##specifying random covariates of id and relationship
##specifying covariance structure by units
ans_gs <- mmer(trait~1,
               random=~vs(id,Gu=A),
               rcov=~units,
               data=Y) 
#summary(ans_gs)$varcomp #check variance composition
#plot(ans_gs) #check normality with visual
gebv.pb <- as.data.frame(ans_gs$U) #pull genomic estimated breeding values from mixed model
rownames(gebv.pb) <- factor(Y$id,levels = rownames(A)) #place id as rownames of gebv
#how to organize by descending gebv and select top (slice(n:n)) individuals
chosenParents <-tibble(id=Y$id, 
                       gebv=gebv.pb$trait) %>% 
  arrange(desc(gebv.pb$trait)) %>% 
  slice(1:5)
chosenParents <- simOutput[[cycle]][chosenParents$id] #subset genomic selected chosenParents from POP to be used in cross
# make crosses 
offspringPop<-randCross(pop=chosenParents, 
                        nCrosses=5,nProgeny = 10)
# phenotype  new offspring
offspringPop<-setPheno(pop = offspringPop,h2=0.5,reps=2)
# add new offspring to simOutput list
simOutput[[cycle+1]]<-offspringPop

###################################################################
#SECOND CYCLE#
###################################################################

gen.mat.wrk <- pullSnpGeno(simOutput[[cycle+1]]) #pull genotypes from offspringPop object
x <- rbind(gen.mat,gen.mat.wrk)
A <- sommer::A.mat(x) #create A genomic relationship matrix
Y <- as.data.frame(cbind(simOutput[[cycle]]@id,simOutput[[cycle]]@pheno)) #create dataframe of traits and id for model
Y <- rbind(Y,cbind(simOutput[[cycle+1]]@id,simOutput[[cycle+1]]@pheno))

colnames(Y) <- c("id","trait") #rename columns something meaningful
Y$trait <- as.numeric(Y$trait) #make response (trait) numeric
#Y <- Y %>%
#    mutate(trait =  replace(trait, sample(row_number(),  
#           size = ceiling(0.3 * n()), replace = FALSE), NA) )

#run univariate mixed model with sommer 
##specifying random covariates of id and relationship
##specifying covariance structure by units
ans_gs <- mmer(trait~1,
               random=~vs(id,Gu=A),
               rcov=~units,
               data=Y) 

gebv.pb <- as.data.frame(ans_gs$U) #pull genomic estimated breeding values from mixed model
rownames(gebv.pb) <- factor(Y$id,levels = rownames(A)) #place id as rownames of gebv
gebv.pb <- subset(gebv.pb,rownames(gebv.pb)==simOutput[[cycle+1]]@id)
chosenParents <-tibble(id=subset(Y,Y$id==simOutput[[cycle+1]]@id), 
                       gebv=gebv.pb$trait) %>% 
  arrange(desc(gebv.pb$trait)) %>% 
  slice(1:5)
chosenParents <- simOutput[[cycle+1]][chosenParents$id$id]

offspringPop<-randCross(pop=chosenParents, 
                        nCrosses=5,nProgeny = 10)
# phenotype  new offspring
offspringPop<-setPheno(pop = offspringPop,h2=0.5,reps=2)
# add new offspring to simOutput list
simOutput[[cycle+2]]<-offspringPop


###################################################################
#THIRD CYCLE#
###################################################################

gen.mat.wrk <- pullSnpGeno(simOutput[[cycle+2]]) #pull genotypes from offspringPop object
x <- rbind(x,gen.mat.wrk)
A <- sommer::A.mat(x) #create A genomic relationship matrix
Y <- as.data.frame(cbind(simOutput[[cycle]]@id,simOutput[[cycle]]@pheno)) #create dataframe of traits and id for model
Y <- rbind(Y,cbind(simOutput[[cycle+1]]@id,simOutput[[cycle+1]]@pheno),cbind(simOutput[[cycle+2]]@id,simOutput[[cycle+2]]@pheno))

colnames(Y) <- c("id","trait") #rename columns something meaningful
Y$trait <- as.numeric(Y$trait) #make response (trait) numeric
#Y <- Y %>%
#    mutate(trait =  replace(trait, sample(row_number(),  
#           size = ceiling(0.3 * n()), replace = FALSE), NA) )

#run univariate mixed model with sommer 
##specifying random covariates of id and relationship
##specifying covariance structure by units
ans_gs <- mmer(trait~1,
               random=~vs(id,Gu=A),
               rcov=~units,
               data=Y) 

gebv.pb <- as.data.frame(ans_gs$U) #pull genomic estimated breeding values from mixed model
rownames(gebv.pb) <- factor(Y$id,levels = rownames(A)) #place id as rownames of gebv
gebv.pb <- subset(gebv.pb,rownames(gebv.pb)==simOutput[[cycle+2]]@id)
chosenParents <-tibble(id=subset(Y,Y$id==simOutput[[cycle+2]]@id), 
                       gebv=gebv.pb$trait) %>% 
  arrange(desc(gebv.pb$trait)) %>% 
  slice(1:5)
chosenParents <- simOutput[[cycle+2]][chosenParents$id$id]

offspringPop<-randCross(pop=chosenParents, 
                        nCrosses=5,nProgeny = 10)
# phenotype  new offspring
offspringPop<-setPheno(pop = offspringPop,h2=0.5,reps=2)
# add new offspring to simOutput list
simOutput[[cycle+2]]<-offspringPop
