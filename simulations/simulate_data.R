####################################################################################################
## Package : CARD
## Version : 1.0.1
## Date    : 2022-5-7 09:10:08
## Title   : Main Simulation framework for CARD.
## Authors : Ying Ma
## Contacts: yingma@umich.edu
##           University of Michigan, Department of Biostatistics
####################################################################################################
## This code provide the simulation procedure for simulating a spatial transcriptomics dataset
## Detailed description about the simulation procejure can be found in the main manuscript and supplementary materials

## set the parameters, 
## randomly set the seed i.e.,
iseed <- 3 
imix <- 0
ntotal <- 10

## load single cell RNAseq dataset
## Here, we split the scRNA-seq data into two sets: 
## one set  (50% cells, denoted as split1) was used to simulate the spatial transcriptomics count 
## while the other set (50% cells, denoted as split2) was used to evaluate the performance of deconvolution methods.
## All data for simulation is stored in https://drive.google.com/drive/folders/1wRPxn1YI7f1oUw8eC42htXMjTUqyIT1g due to its large size

setwd("./data/")
load("./split.scRNAseq.forSimu.RData")
eset.sub.split1 = split$eset.sub.split1 ## simulate the data
eset.sub.split2 = split$eset.sub.split2 ## use as the reference for downstream deconvolution analysis
ct.varname = "cellType"
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
sample.varname = "sampleID"

## load the predefined layer label data and ST MOB location data
load("./pattern_gp_label.RData")
load("./sim_MOB_location.RData")

## functions that used to simulate the data

#### generate random numbers from Dirichelet distribution
generateMultiN <- function(pattern_gp_label,ipt,ntotal,mix,ct.select){
    library(MCMCpack)
    message(paste0("Generating cell type proportions for pattern",ipt))
    nobs = sum(pattern_gp_label == ipt)
    sample = matrix(0,nrow = nobs,ncol = length(ct.select))
    colnames(sample) = ct.select
    sampleNames = names(pattern_gp_label)[pattern_gp_label == ipt]
    prop = matrix(0,nrow = nobs,ncol = length(ct.select))
    colnames(prop) = ct.select
    ## sample total number of cell types in each layer: main cell type + colocalized cell types
    numCT = sample(1:length(ct.select),1) 
    if(ipt == 1){
      main = "Neurons" ### defined one dominant cell type in layer 1
      concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
    }else if (ipt == 2){ 
      main = "Astrocytes" ### defined one dominant cell type in layer 2
      concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
    }else if(ipt == 3){
      main = "Oligos" ### defined one dominant cell type in layer 3
      concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
    }
    propSample = rdirichlet(nrow(prop),concen)
    ct.select.sub.sample = sample(ct.select[ct.select != main],numCT-1)
    ct.select.sub = c(main,ct.select.sub.sample)
    Fix_Dirichlet = sample(1:nrow(sample),round(nobs * mix[1]))
    mix_Dirichlet = c(1:nrow(sample))[!(c(1:nrow(sample)) %in% Fix_Dirichlet)] 
    if(length(mix_Dirichlet) > 0){
    propSample[mix_Dirichlet,] = rdirichlet(length(mix_Dirichlet),rep(1,numCT))
    }
    print(ct.select.sub)
    if(length(Fix_Dirichlet) > 0){
    propSample[Fix_Dirichlet,] = t(sapply(Fix_Dirichlet,function(i){
      propSample[i,][order(propSample[i,],decreasing = T)] ### for non-noisy positions, order the proportion to assign the largest proportion to the dominant cell type
      }))
  }
    colnames(propSample) = ct.select.sub
    prop[,ct.select.sub] = propSample
    sample = round(ntotal * prop,digits = 0)
##### to avoid no numbers sample
    index = which(rowSums(sample) == 0)
    if(length(index) > 0){
    sample[index,] = t(sapply(index,function(i){rmultinom(1, ntotal, prob = prop[i,])}))
    }
####
  return(list(sample = sample)
}

generateSpatial_norep_fixedProp <- function(seed, eset.sub.split1, ct.varname, sample.varname,ct.select,sample.withRep = F,pattern_gp_label, ntotal,mix1,mix2,mix3){
  #### using split 1 to sample the single cell RNAseq data
  phenoData <- eset.sub.split1@phenoData@data
  # number of cell type of interest
  k <- length(unique(ct.select))
  message(paste('Using',k,'cell types to generate pseudo spatial dataset'))
  # select donors for each pseudo bulk sample.varname
  Sample_random = matrix(data = 0,ncol = length(ct.select), nrow = length(pattern_gp_label))
  rownames(Sample_random) = names(pattern_gp_label)
  colnames(Sample_random) = ct.select
  ##### Total number of cells in the subset
  ct.id <- droplevels(as.factor(eset.sub.split1@phenoData@data[,ct.varname]))
  library(Hmisc)
  ### random generate the number of cells for each spatial location in each layer
  pattern1 = generateMultiN(pattern_gp_label,1,ntotal,mix1,ct.select)
  pattern2 = generateMultiN(pattern_gp_label,2,ntotal,mix2,ct.select)
  pattern3 = generateMultiN(pattern_gp_label,3,ntotal,mix3,ct.select)
  
  Sample_random[pattern_gp_label == 1,] = pattern1$sample
  Sample_random[pattern_gp_label == 2,] = pattern2$sample
  Sample_random[pattern_gp_label == 3,] = pattern3$sample
  
  message(paste0("Generating pseudo spatial dataset for ",length(unique(pattern_gp_label))," patterns"))
  ##### calculate the number of cells we need to sample in each cell type
  set.seed(seed)
  temp.exprs <- exprs(eset.sub.split1)
  temp.nct <- Sample_random
  true.p = sweep(temp.nct,1,rowSums(temp.nct),"/")
  ##### use mcapply since pbmcapply will randomly assign the seed
  temp.pseudo = pbmclapply(1:nrow(temp.nct),function(isample){
      temp.sample = temp.nct[isample,]
      temp.sample.pseudo = sapply(ct.select,function(ict){
          temp.vec <- temp.exprs[,ct.id %in% ict]
          if(temp.nct[isample,ict] > ncol(temp.vec)){
            sample.withRep = T
          }
          temp.id <- sample(1:ncol(temp.vec), temp.nct[isample,ict], replace = sample.withRep)
          if(length(temp.id) > 1){
          rowSums(temp.vec[,temp.id])
          }else if(length(temp.id) == 1){
             temp.vec[,temp.id]
          }else if(length(temp.id) == 0){
             rep(0,nrow(temp.vec))
          }
      })
      rowSums(temp.sample.pseudo)
},mc.cores = 70,mc.set.seed = F)
  temp.pseudo = do.call("cbind",temp.pseudo)
  colnames(temp.pseudo) = rownames(Sample_random)
  rownames(temp.pseudo) = rownames(temp.exprs)
  return(list(pseudo.data = temp.pseudo, true.p = true.p, ntotal = ntotal,Sample_random = Sample_random))
}  

## simulate the dataset

library(pbmcapply)
library(SingleCellExperiment)
## we set the proportion of noisy locations to be 0, 0.2, 0.4, 0.6, corresponidngly imix = 0, 1, 2, 3
mix1 = mix2 = mix3 = c(1 - (0.2 * imix),0.2*imix)
set.seed(iseed)
spatial.pseudo = generateSpatial_norep_fixedProp(
seed = iseed,
eset.sub.split1 = eset.sub.split1,
ct.varname = ct.varname,
sample.varname = sample.varname,
ct.select = ct.select,
sample.withRep = F,
pattern_gp_label = pattern_gp_label,
ntotal = ntotal,
mix1 = mix1,
mix2 = mix2,
mix3 = mix3)

## check the simulated dataset, average number of cells in each spot
print(round(mean(rowSums(spatial.pseudo$Sample_random),0)))
10
## dimension of simulated count data
print(dim(spatial.pseudo$pseudo.data))
[1] 18263   260

## column means of simulated proportion across spatial locations
print(colMeans(spatial.pseudo$true.p))
Astrocytes    Neurons     Oligos   Vascular     Immune  Ependymal 
0.07065657 0.25750971 0.38714744 0.05901515 0.10253594 0.12313520 

## Due to the raondomness, for all the simulated data generated in Figure 2, you can find it in https://drive.google.com/drive/folders/1wRPxn1YI7f1oUw8eC42htXMjTUqyIT1g
## Then you can run CARD


