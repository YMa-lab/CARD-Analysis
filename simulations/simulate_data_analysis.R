####################################################################################################
## Package : CARD
## Version : 1.0.1
## Date    : 2022-06-21 23:10:08
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

## Here is an example of scenario 1
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


## Here is an example of scenario 2: missing one cell type, i.e. missing the neurons
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
ct.select = ct.select[ct.select != "Neurons"]
sample.varname = "sampleID"
# load the simulated count data
load("./simulated_data_Figure2/sim.pseudo.MOB.n10.cellType6.Mixnoise0.repeat3.RData")
## since we use the first split to generate the simulation, for evaluating the robustness of CARD, we use the second split to do the devconovluton analysis
spatial_location = cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                   y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
rownames(spatial_location) = colnames(spatial.pseudo$pseudo.data)
library(CARD)
CARD_obj_S2 = createCARDObject(
    sc_count = exprs(eset.sub.split2),
    sc_meta = pData(eset.sub.split2),
    spatial_count = spatial.pseudo$pseudo.data,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = ct.select,
    sample.varname = "sampleID",
    minCountGene = 100,
    minCountSpot = 5) 
CARD_obj_S2 = CARD_deconvolution(CARD_object = CARD_obj_S2)

## Scenario 3: more cell type, 
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
ct.select = c(ct.select,"Blood")
sample.varname = "sampleID"
# load the simulated count data
load("./simulated_data_Figure2/sim.pseudo.MOB.n10.cellType6.Mixnoise0.repeat3.RData")
## since we use the first split to generate the simulation, for evaluating the robustness of CARD, we use the second split to do the devconovluton analysis
spatial_location = cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                   y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
rownames(spatial_location) = colnames(spatial.pseudo$pseudo.data)
library(CARD)
CARD_obj_S3 = createCARDObject(
    sc_count = exprs(eset.sub.split2),
    sc_meta = pData(eset.sub.split2),
    spatial_count = spatial.pseudo$pseudo.data,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = ct.select,
    sample.varname = "sampleID",
    minCountGene = 100,
    minCountSpot = 5) 
CARD_obj_S3 = CARD_deconvolution(CARD_object = CARD_obj_S3)

## Scenario 4: miss-classified cell type, i.e. merged two cell types, here is an example of merged two cell types: mergeAstrocytes_Ependymal
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
Combinations = combn(ct.select,2)
#### for example, when we merge the cell types of Astrocytes Ependymal into one cell type
icom = 5
imerge = paste(Combinations[,icom],collapse = "_")
#### create the new single cell RNAseq with the merged cell types
eset.sub.split2.merged = eset.sub.split2
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
ct.select.merged = c(ct.select[!(ct.select %in% Combinations[,icom])],"Merged")
pData(eset.sub.split2.merged)[,ct.varname] <- as.character(pData(eset.sub.split2.merged)[,ct.varname])
pData(eset.sub.split2.merged)[,ct.varname][pData(eset.sub.split2.merged)[,ct.varname] %in% Combinations[,icom]] <- "Merged" 

load("./simulated_data_Figure2/sim.pseudo.MOB.n10.cellType6.Mixnoise0.repeat3.RData")
## since we use the first split to generate the simulation, for evaluating the robustness of CARD, we use the second split to do the devconovluton analysis
spatial_location = cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                   y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
rownames(spatial_location) = colnames(spatial.pseudo$pseudo.data)
library(CARD)

CARD_obj_S4 = createCARDObject(
    sc_count = exprs(eset.sub.split2.merged),
    sc_meta = pData(eset.sub.split2.merged),
    spatial_count = spatial.pseudo$pseudo.data,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = ct.select.merged,
    sample.varname = "sampleID",
    minCountGene = 100,
    minCountSpot = 5) 
CARD_obj_S4 = CARD_deconvolution(CARD_object = CARD_obj_S4)

## Scenario 5, a different scRNAseq reference sequenced from a different platform
## The similar scRNAseq reference dataset is available here: https://drive.google.com/drive/folders/18vSET90a3_cvcOwLdyc2beqlM6gJxZwX
sc_eset.scenario5 = readRDS("./GSE109447.sceset.scenario5.RDS")
## we used the following cell types that were highly similar to the cell types in the scRNAseq data that were used to simulate the dataset f
ct.select.similar = c("Astrocytes","Neuron","Oligodendrocytes","Endothelial","Microglia","Ependymal")
## the ct.select.similar matched with the ct.select.original, see supplementary figure 94 in the supplementary material of CARD
## ct.select.original = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
ct.varname = "cellType"
sample.varname = "sampleID"
sc_eset.scenario5 = sc_eset.scenario5[,pData(sc_eset.scenario5)[,ct.varname] %in% ct.select.similar]
load("./simulated_data_Figure2/sim.pseudo.MOB.n10.cellType6.Mixnoise0.repeat3.RData")
## since we use the first split to generate the simulation, for evaluating the robustness of CARD, we use the second split to do the devconovluton analysis
spatial_location = cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                   y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
rownames(spatial_location) = colnames(spatial.pseudo$pseudo.data)
library(CARD)

CARD_obj_S5 = createCARDObject(
    sc_count = exprs(sc_eset.scenario5),
    sc_meta = pData(sc_eset.scenario5),
    spatial_count = spatial.pseudo$pseudo.data,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = ct.select.similar,
    sample.varname = "sampleID",
    minCountGene = 100,
    minCountSpot = 5) 
CARD_obj_S5 = CARD_deconvolution(CARD_object = CARD_obj_S5)





