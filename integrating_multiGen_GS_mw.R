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


nCycles<-5 # loop length / number of cycles of selection
nParents<-5 # number of parents to select and cross
nCrosses<-5
nProgeny<-10

# very simple container for each cycles sim output
simOutput<-list(founders)

#cycle<-1
for(cycle in 1:nCycles){
  cat(paste0(" C",cycle))

  # concatenate the list of sims
  allCyclesOutput<-purrr::reduce(simOutput,`c`)


  gen.mat <- pullSnpGeno(allCyclesOutput) #pull genotypes from offspringPop object
  A <- sommer::A.mat(gen.mat) #create A genomic relationship matrix

  #create dataframe of traits and id for model
  Y <- tibble(id=allCyclesOutput@id,
              trait=as.numeric(allCyclesOutput@pheno)) %>%
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
  chosenParents <- allCyclesOutput[chosenParents$id] #subset genomic selected chosenParents from POP to be used in cross
  # make crosses
  offspringPop<-randCross(pop=chosenParents,
                          nCrosses=nCrosses,
                          nProgeny = nProgeny)
  # phenotype new offspring
  ## NOTICE: we are phenotyping these progeny BEFORE genotyping and predicting
  ## GS enables us to actually predict GEBV and select offspring without first phenotyping them

  offspringPop<-setPheno(pop = offspringPop,h2=0.5,reps=2)
  # add new offspring to simOutput list
  simOutput[[cycle+1]]<-offspringPop

}

# Tidy up to output
tidySimOutput<-tibble(Cycle=0:nCycles,
                      Sims=simOutput) %>%
  mutate(meanG=map_dbl(Sims,~mean(.@gv)),
         varG=map_dbl(Sims,~var(.@gv)))
tidySimOutput
