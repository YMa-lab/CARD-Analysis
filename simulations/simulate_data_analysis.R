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
## load single cell RNAseq dataset
## Here, we split the scRNA-seq data into two sets: 
## one set  (50% cells, denoted as split1) was used to simulate the spatial transcriptomics count 
## while the other set (50% cells, denoted as split2) was used to evaluate the performance of deconvolution methods.
## All data for simulation is stored in https://drive.google.com/drive/folders/1wRPxn1YI7f1oUw8eC42htXMjTUqyIT1g due to its large size

setwd("./data/")
load("./split.scRNAseq.forSimu.RData")
#eset.sub.split1 = split$eset.sub.split1 ## simulate the data
eset.sub.split2 = split$eset.sub.split2 ## use as the reference for downstream deconvolution analysis
ct.varname = "cellType"
## For ssimulation scenario I, we use the same set of cell types that generated the simulated count to do the deconvolution
## for other simulation scenarios, please follow the materials and methods and supplementary materials to modify the reference data
## i.e. in Scenario II, missing one cell type at a time, while in scenario III, add one more "Blood" cell type that is left in split2
## in Scenario IV, we merged two cell types into one cell type, so there are 15 combinations
## in scenario V, similar scRNA-seq data but sequenced from different platform was used for deconvolution. Details about the data is available 

## Here is an example of scenario I
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
sample.varname = "sampleID"
# load the simulated count data
load("./simulated_data_Figure2/sim.pseudo.MOB.n10.cellType6.Mixnoise0.repeat3.RData")
## since we use the first split to generate the simulation, for evaluating the robustness of CARD, we use the second split to do the devconovluton analysis
spatial_location = cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                   y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
rownames(spatial_location) = colnames(spatial.pseudo$pseudo.data)
library(CARD)
CARD_obj_S1 = createCARDObject(
    sc_count = exprs(eset.sub.split2),
    sc_meta = pData(eset.sub.split2),
    spatial_count = spatial.pseudo$pseudo.data,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = ct.select,
    sample.varname = "sampleID",
    minCountGene = 100,
    minCountSpot = 5) 
CARD_obj_S1 = CARD_deconvolution(CARD_object = CARD_obj_S1)
## The estimated proportion is stored in CARD_obj_S1@Proportion_CARD
