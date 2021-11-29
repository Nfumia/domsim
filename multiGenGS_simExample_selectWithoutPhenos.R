## I'm going to add a few bells and whistles to the multigen GS.
## First, we are going to genotype _before_ phenotyping
## Second, let's add some control over the number of years of data
### that are used in the genomic prediction
### if you simulated 50 years, you could have a huge matrix
### Which might be a compute burden you do not want
library(AlphaSimR)
library(tidyverse)
library(sommer)
rm(list=ls())

###################
# simulation input parameters
###################
nCycles<-10 # loop length / number of cycles of selection
nParents<-5 # number of parents to select and cross
nCrosses<-5
nProgeny<-10
nYearsToUse<-5 # use up to and including this number of previous years for training / genomic prediction
###################

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

# Genotype, generate GEBV and select founders to make initial crosses
# We WILL NOT phenotype the progeny until after genotyping them and predicting
# their GEBV
## in reality, you might not phenotype everyone you genotype unless they get
## selected... will avoid that nuance
gen.mat <- pullSnpGeno(founders) #pull genotypes from offspringPop object
A <- sommer::A.mat(gen.mat) #create A genomic relationship matrix
#create dataframe of traits and id for model
Y <- tibble(id=founders@id,
            trait=as.numeric(founders@pheno)) %>%
  # make "id" a factor with levels = rownames(A)
  # will ensure prediction of all genotyped indivs
  dplyr::mutate(id=factor(id,levels=rownames(A)))
#run univariate mixed model with sommer
ans_gs <- mmer(trait~1,
               random=~vs(id,Gu=A),
               rcov=~units,
               data=Y)
# make a tibble of the GEBV from the mixed model
gebv <- tibble(id=names(ans_gs$U$`u:id`$trait),
               gebv=as.numeric(ans_gs$U$`u:id`$trait))
# Choose parents with the highest GEBV
chosenParents <- gebv %>%
  arrange(desc(gebv)) %>%
  slice(1:nParents) # created a param to set outside the loop, how many parents to choose
chosenParents <- founders[chosenParents$id] #subset genomic selected chosenParents from POP to be used in cross
# make crosses
offspringPop<-randCross(pop=chosenParents,
                        nCrosses=nCrosses,
                        nProgeny = nProgeny)
# add new offspring, no phenotypes, to simOutput list
# very simple container for each cycles sim output
simOutput<-list(founders,offspringPop)

# ENTER THE SELECTION CYCLE LOOP
#cycle<-1
for(cycle in 1:nCycles){
  cat(paste0(" C",cycle))

  # concatenate the list of sims
  # looks like tail() is a convenient way to
  # choose up to and including the last nYearsToUse
  # to be clear: after nYearsToUse cycles, the earliest cycles
  # will be dropped from the analysis
  trainingCyclesOutput<-purrr::reduce(simOutput %>% tail(nYearsToUse),`c`)

  gen.mat <- pullSnpGeno(trainingCyclesOutput) #pull genotypes from offspringPop object
  A <- sommer::A.mat(gen.mat) #create A genomic relationship matrix

  #create dataframe of traits and id for model
  Y <- tibble(id=trainingCyclesOutput@id,
              trait=as.numeric(trainingCyclesOutput@pheno)) %>%
    # make "id" a factor with levels = rownames(A)
    # will ensure prediction of all genotyped indivs
    # even the progeny not-yet-phenotyped
    dplyr::mutate(id=factor(id,levels=rownames(A)))

  #run univariate mixed model with sommer
  ans_gs <- mmer(trait~1,
                 random=~vs(id,Gu=A),
                 rcov=~units,
                 data=Y)

  # make a tibble of the GEBV from the mixed model

  gebv <- tibble(id=names(ans_gs$U$`u:id`$trait),
                 gebv=as.numeric(ans_gs$U$`u:id`$trait))

  # Choose parents with the highest GEBV
  chosenParents <- gebv %>%
    arrange(desc(gebv)) %>%
    slice(1:nParents) # created a param to set outside the loop, how many parents to choose
  chosenParents <- trainingCyclesOutput[chosenParents$id] #subset genomic selected chosenParents from POP to be used in cross
  # make crosses
  offspringPop<-randCross(pop=chosenParents,
                          nCrosses=nCrosses,
                          nProgeny = nProgeny)

  # phenotype new offspring from the PREVIOUS CYCLE
  ## NOTICE: we are phenotyping these progeny AFTER genotyping/predicting them
  ## Not phenotyping the `offspringPop` above
  simOutput[[length(simOutput)]]<-setPheno(pop = simOutput[[length(simOutput)]],
                                           h2=0.5,reps=2)
  # add new offspring to simOutput list
  simOutput[[length(simOutput)+1]]<-offspringPop
}

# Tidy up to output
tidySimOutput<-tibble(Cycle=0:(length(simOutput)-1),
                      Sims=simOutput) %>%
  mutate(meanG=map_dbl(Sims,~mean(.@gv)),
         varG=map_dbl(Sims,~var(.@gv)))
tidySimOutput

tidySimOutput %>% ggplot(.,aes(x=Cycle,y=meanG)) + geom_line()
