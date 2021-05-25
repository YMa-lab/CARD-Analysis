generateSpatial_norep_fixedProp <- function(seed,sc_eset, ct.varname, ct.select,sample.withRep = F,pattern, ntotal,con1,con2,con3,mix1,mix2,mix3,ncores){
  #### First QC on the scRNAseq data to simulate the spatial count data more reasonably
  counts = exprs(sc_eset)
  coldf = pData(sc_eset)
  nfeatures <- Matrix::colSums(x = counts > 0)
  counts <- counts[, which(x = nfeatures >= 300)]
  coldf <- coldf[which(x = nfeatures >= 300),]
  # filter genes on the number of cells expressing
  num.cells <- Matrix::rowSums(x = counts > 0)
  counts <- counts[which(x = num.cells >= 3), ]
  fdata = as.data.frame(rownames(counts))
  rownames(fdata) = rownames(counts)
  sc.eset <- Biobase::ExpressionSet(counts, phenoData=AnnotatedDataFrame(coldf),AnnotatedDataFrame(fdata))
  sc.eset.sub <- sc.eset[,sc.eset@phenoData@data[,ct.varname] %in% ct.select]
  # qc: remove non-zero genes
  sc.eset.sub <- sc.eset.sub[rowSums(exprs(sc.eset.sub)) > 0,]  
  phenoData <- sc.eset.sub@phenoData@data
  # number of cell type of interest
  k <- length(unique(ct.select))
  message(paste('Using',k,'cell types to generate pseudo spatial dataset'))
  # select donors for each pseudo bulk 
  Sample_random = matrix(data = 0,ncol = length(ct.select), nrow = length(pattern_gp_label))
  rownames(Sample_random) = names(pattern_gp_label)
  ##### Total number of cells in the subset
  ct.id <- droplevels(as.factor(sc.eset.sub@phenoData@data[,ct.varname]))
  ### suppose each spot contains 500 cells
  ##### Pattern I #### modify this later to be more flexible
  #pr1 = c(0.7,0.10,0.10,0.10)
  pattern1 = generateMultiN(pattern_gp_label,1,con1,ntotal,mix1)
  pattern2 = generateMultiN(pattern_gp_label,2,con2,ntotal,mix2)
  pattern3 = generateMultiN(pattern_gp_label,3,con3,ntotal,mix3)

  Sample_random[pattern_gp_label == 1,] = pattern1$sample
  Sample_random[pattern_gp_label == 2,] = pattern2$sample
  Sample_random[pattern_gp_label == 3,] = pattern3$sample
  
  pattern1Samples = list(multinormial_samples = pattern1$multinormial_samples,dirichlet_samples=pattern1$dirichlet_samples)
  pattern2Samples = list(multinormial_samples = pattern2$multinormial_samples,dirichlet_samples=pattern2$dirichlet_samples)
  pattern3Samples = list(multinormial_samples = pattern3$multinormial_samples,dirichlet_samples=pattern3$dirichlet_samples)
  multinormial_samples = c(pattern1$multinormial_samples,pattern2$multinormial_samples,pattern3$multinormial_samples)
  dirichlet_samples=c(pattern1$dirichlet_samples,pattern2$dirichlet_samples,pattern3$dirichlet_samples)

  ##### Pattern II
  #pr2 = c(0.10,0.7,0.10,0.10)
  
  colnames(Sample_random) = ct.select
  message(paste0("Generating pseudo spatial dataset for ",length(unique(pattern_gp_label))," patterns"))
  ##### calculate the number of cells we need to sample in each cell type
  temp.exprs <- exprs(sc.eset.sub)
  temp.nct <- Sample_random
  true.p = sweep(temp.nct,1,rowSums(temp.nct),"/")
  set.seed(seed)
  temp.pseudo = pbmclapply(1:nrow(temp.nct),function(isample){
      temp.sample = temp.nct[isample,]
      temp.sample.pseudo = sapply(ct.select,function(ict){
          temp.vec <- temp.exprs[,ct.id %in% ict]
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
},mc.cores = ncores,mc.set.seed = F)
  temp.pseudo = do.call("cbind",temp.pseudo)
  colnames(temp.pseudo) = rownames(Sample_random)
  rownames(temp.pseudo) = rownames(temp.exprs)
  con = data.frame(con1 = con1, con2 = con2, con3 = con3)
  mix = data.frame(mix1 = mix1, mix2 = mix2, mix3 = mix3)
  rownames(con) = colnames(true.p)
  return(list(pseudo.data = temp.pseudo, true.p = true.p, ntotal = ntotal, con = con,Sample_random = Sample_random,
    pattern1Samples = pattern1Samples,pattern2Samples=pattern2Samples,pattern3Samples = pattern3Samples,
    multinormial_samples = multinormial_samples,
    dirichlet_samples = dirichlet_samples))
}

#### generate multinormial random numbers
generateMultiN <- function(pattern_gp_label,ipt,con,ntotal,mix){
    library(MCMCpack)
    message(paste0("Generating cell type proportions for pattern",ipt))
    nobs = sum(pattern_gp_label == ipt)
    sample = matrix(NA,nrow = nobs,ncol = length(con))
    sampleNames = names(pattern_gp_label)[pattern_gp_label == ipt]
    mix_multinormial = sample(1:nrow(sample),round(nobs * mix[1]))
    mix_Dirichlet = c(1:nrow(sample))[!(c(1:nrow(sample)) %in% mix_multinormial)] 
    if(length(mix_multinormial) > 0){
    sample[mix_multinormial,] = t(rmultinom(length(mix_multinormial), ntotal, prob = con))
    }
    if(length(mix_Dirichlet) > 0){
    prop_dirichlet = rdirichlet(length(mix_Dirichlet),rep(1,length(con)))
    sample[mix_Dirichlet,] = t(apply(prop_dirichlet,1,function(x){
      rmultinom(1, ntotal, prob = x)
      }))
  }
    return(list(sample = sample,multinormial_samples = sampleNames[mix_multinormial],dirichlet_samples = sampleNames[mix_Dirichlet]))
}
