###############################################################
## expressolog score calculation - consumes a LOT of memory! ##
###############################################################

calc_spec <- function (res, np, nL, label = 'spec'){    
  
  ranks = 1:nL
  mini = sum(ranks[1:np])
  range = np*(nL - np)
  
  temp1 = t(apply(res, 1, function(x) rank(x, na.last = 'keep', ties.method = "average")))      
  spec1 <- (temp1 - mini)/range   
  temp2 = t(apply(t(res), 1, function(x) rank(x, na.last = 'keep', ties.method = "average")))      
  spec2 <- (temp2 - mini)/range
  
  # birectional spec scores and ranks                
  if(label=='spec'){
    specna = (!is.na(spec1)) + (!is.na(spec2))
    spec1 = as.data.frame(spec1)
    spec2 = as.data.frame(spec2)
    spec1[is.na(spec1)] = 0        
    spec2[is.na(spec2)] = 0        
    spec = (spec1 + spec2) / specna        
    return(spec)
  }else if(label=='rank'){
    tempna = (!is.na(temp1)) + (!is.na(temp2))
    temp1 = as.data.frame(temp1)
    temp2 = as.data.frame(temp2)
    temp1[is.na(temp1)] = 0        
    temp2[is.na(temp2)] = 0        
    rankmat = (temp1 + temp2) / tempna        
    return(rankmat)
  }else{   # rank-standardize
    tempna = (!is.na(temp1)) + (!is.na(temp2))
    temp1 = as.data.frame(temp1)
    temp2 = as.data.frame(temp2)
    temp1[is.na(temp1)] = 0        
    temp2[is.na(temp2)] = 0        
    rankmat = (temp1 + temp2) / tempna
    rankmat <- rankmat/max(rankmat, na.rm = T)
    return(rankmat)
  }
}  


# calculate specificity at each hierarchical level                    
calc_spec_list <- function (mlist, np, nL, startid, endid){    
  
  smat = 0
  snum = 0
  for(ii in startid:endid){
    currmat = as.data.frame(eval(as.symbol(mlist[ii])))        
    snum = snum + (!is.na(currmat))
    currmat[is.na(currmat)] = 0
    smat = smat + currmat
  }
  res = smat/snum  # resulting matrix of average ranks for each gene pair
  
  ranks = 1:nL
  mini = sum(ranks[1:np])
  range = np*(nL - np)
  
  temp1 = t(apply(res, 1, function(x) rank(x, na.last = 'keep', ties.method = "average")))      
  spec1 <- (temp1 - mini)/range   
  temp2 = t(apply(t(res), 1, function(x) rank(x, na.last = 'keep', ties.method = "average")))      
  spec2 <- (temp2 - mini)/range
  
  # birectional spec scores and ranks   
  specna = (!is.na(spec1)) + (!is.na(spec2))
  spec1 = as.data.frame(spec1)
  spec2 = as.data.frame(spec2)
  spec1[is.na(spec1)] = 0        
  spec2[is.na(spec2)] = 0        
  spec = (spec1 + spec2) / specna
  return(spec)   
}            

# calculate expressolog score for 14k genes in a pair of species across different
# cell type groups
get_expressolog_score <- function(sp1, sp2){
  
  # cell type groups over which we calculate expressolog scores for 14k orthologs
  glia_groups <- 1:7
  glia1 = 1:4
  glia2 = 5:7
  
  gaba_groups1 <- 8:23
  gaba11 = 1:3
  gaba12 = 4:9
  gaba13 = 10:16
  
  gaba_groups2 <- 24:38
  gaba21 = 1:5
  gaba22 = 6:11
  gaba23 = 12:15
  
  glut_groups1 <- 39:47
  glut11 = 1:4
  glut12 = 5:9
  
  glut_groups2 <- 48:57
  glut21 = 1:5
  glut22 = 6:10
  
  
  # load processed scRNA-seq - pseudo-bulk
  matlist = c('gliacorrspec', 'glutcorrspec', 'gabacorrspec', 'gliacorrspec1', 'gliacorrspec2', 'glutcorrspec1',
              'glutcorrspec2', 'gabacorrspec1', 'gabacorrspec2', 'glutcorrspec11', 'glutcorrspec12', 'glutcorrspec21', 'glutcorrspec22',
              'gabacorrspec11', 'gabacorrspec12', 'gabacorrspec13', 'gabacorrspec21', 'gabacorrspec22', 'gabacorrspec23')
  
  spe = c('human', 'chimp', 'gorilla', 'rhesus', 'marmoset')
  matlist1 = c('sch_cls', 'scc_cls', 'scg_cls', 'scr_cls', 'scm_cls')
  
  load(paste0('/data/suresh/allen/pseudobulk/', sp1, '_cluster_expression_new.Rdata'))
  load(paste0('/data/suresh/allen/pseudobulk/', sp2, '_cluster_expression_new.Rdata'))
  mat1 = eval(as.symbol(matlist1[match(sp1, spe)]))
  mat2 = eval(as.symbol(matlist1[match(sp2, spe)]))
  
  options(warn = -1)
  
  ## ------------------------------ overall ----------------------------- ##
  allcorr = cor(t(mat1), t(mat2), method = 'pearson')    
  allcorrspec = calc_spec(allcorr, 1, dim(allcorr)[1])   # specificity scores
  
  
  ## ------------------------- class-level ----------------------------- ##
  gliacorr = cor(t(mat1[,glia_groups]), t(mat2[,glia_groups]), method = 'pearson')
  glutcorr = cor(t(mat1[,c(glut_groups1,glut_groups2)]), t(mat2[,c(glut_groups1,glut_groups2)]), method = 'pearson')
  gabacorr = cor(t(mat1[,c(gaba_groups1,gaba_groups2)]), t(mat2[,c(gaba_groups1,gaba_groups2)]), method = 'pearson')
  
  # specificity scores
  gliacorrspec = calc_spec(gliacorr, 1, dim(gliacorr)[1], 'rank')
  glutcorrspec = calc_spec(glutcorr, 1, dim(glutcorr)[1], 'rank')
  gabacorrspec = calc_spec(gabacorr, 1, dim(gabacorr)[1], 'rank')
  
  classavgspec = calc_spec_list(matlist, 1, dim(allcorr)[1], 1, 3)
  
  
  ## ------------------------- subclass-level ----------------------------- ##
  gliacorr1 = cor(t(mat1[,glia_groups[glia1]]), t(mat2[,glia_groups[glia1]]), method = 'pearson')
  gliacorr2 = cor(t(mat1[,glia_groups[glia2]]), t(mat2[,glia_groups[glia2]]), method = 'pearson')
  glutcorr1 = cor(t(mat1[,glut_groups1]), t(mat2[,glut_groups1]), method = 'pearson')
  glutcorr2 = cor(t(mat1[,glut_groups2]), t(mat2[,glut_groups2]), method = 'pearson')
  gabacorr1 = cor(t(mat1[,gaba_groups1]), t(mat2[,gaba_groups1]), method = 'pearson')
  gabacorr2 = cor(t(mat1[,gaba_groups2]), t(mat2[,gaba_groups2]), method = 'pearson')
  
  # specificity scores
  gliacorrspec1 = calc_spec(gliacorr1, 1, dim(gliacorr1)[1], 'rank')
  gliacorrspec2 = calc_spec(gliacorr2, 1, dim(gliacorr2)[1], 'rank')
  glutcorrspec1 = calc_spec(glutcorr1, 1, dim(glutcorr1)[1], 'rank')
  glutcorrspec2 = calc_spec(glutcorr2, 1, dim(glutcorr2)[1], 'rank')
  gabacorrspec1 = calc_spec(gabacorr1, 1, dim(gabacorr1)[1], 'rank')
  gabacorrspec2 = calc_spec(gabacorr2, 1, dim(gabacorr2)[1], 'rank')
  
  subclassavgspec = calc_spec_list(matlist, 1, dim(allcorr)[1], 4, 9)
  
  
  
  ## ------------------------- cluster-level ----------------------------- ##
  ## glut ##
  glutcorr11 = cor(t(mat1[,glut_groups1[glut11]]), t(mat2[,glut_groups1[glut11]]), method = 'pearson')
  glutcorr12 = cor(t(mat1[,glut_groups1[glut12]]), t(mat2[,glut_groups1[glut12]]), method = 'pearson')
  glutcorr21 = cor(t(mat1[,glut_groups2[glut21]]), t(mat2[,glut_groups2[glut21]]), method = 'pearson')
  glutcorr22 = cor(t(mat1[,glut_groups2[glut22]]), t(mat2[,glut_groups2[glut22]]), method = 'pearson')
  
  ## gaba ##
  gabacorr11 = cor(t(mat1[,gaba_groups1[gaba11]]), t(mat2[,gaba_groups1[gaba11]]), method = 'pearson')
  gabacorr12 = cor(t(mat1[,gaba_groups1[gaba12]]), t(mat2[,gaba_groups1[gaba12]]), method = 'pearson')
  gabacorr13 = cor(t(mat1[,gaba_groups1[gaba13]]), t(mat2[,gaba_groups1[gaba13]]), method = 'pearson')
  gabacorr21 = cor(t(mat1[,gaba_groups2[gaba21]]), t(mat2[,gaba_groups2[gaba21]]), method = 'pearson')
  gabacorr22 = cor(t(mat1[,gaba_groups2[gaba22]]), t(mat2[,gaba_groups2[gaba22]]), method = 'pearson')
  gabacorr23 = cor(t(mat1[,gaba_groups2[gaba23]]), t(mat2[,gaba_groups2[gaba23]]), method = 'pearson')
  
  
  
  # specificity scores
  glutcorrspec11 = calc_spec(glutcorr11, 1, dim(glutcorr11)[1], 'rank')
  glutcorrspec12 = calc_spec(glutcorr12, 1, dim(glutcorr12)[1], 'rank')
  glutcorrspec21 = calc_spec(glutcorr21, 1, dim(glutcorr21)[1], 'rank')
  glutcorrspec22 = calc_spec(glutcorr22, 1, dim(glutcorr22)[1], 'rank')
  
  gabacorrspec11 = calc_spec(gabacorr11, 1, dim(gabacorr11)[1], 'rank')
  gabacorrspec12 = calc_spec(gabacorr12, 1, dim(gabacorr12)[1], 'rank')
  gabacorrspec13 = calc_spec(gabacorr13, 1, dim(gabacorr13)[1], 'rank')
  gabacorrspec21 = calc_spec(gabacorr21, 1, dim(gabacorr21)[1], 'rank')
  gabacorrspec22 = calc_spec(gabacorr22, 1, dim(gabacorr22)[1], 'rank')
  gabacorrspec23 = calc_spec(gabacorr23, 1, dim(gabacorr23)[1], 'rank')
  
  
  clusteravgspec = calc_spec_list(matlist, 1, dim(allcorr)[1], 10, 19)
  
  
  
  # ........................................................................ #
  # create data.frame, and save scores for 1-1 orthologs only
  dfsub1 = data.frame(rownames(mat1), rownames(mat2), diag(as.matrix(allcorrspec)), 
                      diag(as.matrix(classavgspec)), diag(as.matrix(subclassavgspec)), diag(as.matrix(clusteravgspec)))
  colnames(dfsub1) <- c(sp1, sp2, 'overall_spec', 'class_spec', 'subclass_spec', 'cluster_spec')
  rownames(dfsub1) <- NULL
  
  
  # ........................................................................ #
  # convert full rank matrices to spec scores
  
  gliacorrspec = calc_spec(gliacorr, 1, dim(gliacorr)[1])
  gliacorrspec1 = calc_spec(gliacorr1, 1, dim(gliacorr1)[1])
  gliacorrspec2 = calc_spec(gliacorr2, 1, dim(gliacorr2)[1])
  
  glutcorrspec = calc_spec(glutcorr, 1, dim(glutcorr)[1])
  glutcorrspec1 = calc_spec(glutcorr1, 1, dim(glutcorr1)[1])
  glutcorrspec2 = calc_spec(glutcorr2, 1, dim(glutcorr2)[1])
  
  gabacorrspec = calc_spec(gabacorr, 1, dim(gabacorr)[1])
  gabacorrspec1 = calc_spec(gabacorr1, 1, dim(gabacorr1)[1])
  gabacorrspec2 = calc_spec(gabacorr2, 1, dim(gabacorr2)[1])
  
  
  ### sub-dividing gaba and glut groups ###
  ## glut ##
  glutcorrspec11 = calc_spec(glutcorr11, 1, dim(glutcorr11)[1])
  glutcorrspec12 = calc_spec(glutcorr12, 1, dim(glutcorr12)[1])
  glutcorrspec21 = calc_spec(glutcorr21, 1, dim(glutcorr21)[1])
  glutcorrspec22 = calc_spec(glutcorr22, 1, dim(glutcorr22)[1])
  
  
  ## gaba ##
  gabacorrspec11 = calc_spec(gabacorr11, 1, dim(gabacorr11)[1])
  gabacorrspec12 = calc_spec(gabacorr12, 1, dim(gabacorr12)[1])
  gabacorrspec13 = calc_spec(gabacorr13, 1, dim(gabacorr13)[1])
  gabacorrspec21 = calc_spec(gabacorr21, 1, dim(gabacorr21)[1])
  gabacorrspec22 = calc_spec(gabacorr22, 1, dim(gabacorr22)[1])
  gabacorrspec23 = calc_spec(gabacorr23, 1, dim(gabacorr23)[1])
  
  
  # create data.frame, and save scores for 1-1 orthologs only
  dfsub2 = data.frame(rownames(mat1), rownames(mat2), diag(as.matrix(allcorrspec)), 
                      diag(as.matrix(gliacorrspec)), diag(as.matrix(glutcorrspec)), diag(as.matrix(gabacorrspec)),
                      diag(as.matrix(gliacorrspec1)), diag(as.matrix(gliacorrspec2)), diag(as.matrix(glutcorrspec1)),
                      diag(as.matrix(glutcorrspec2)), diag(as.matrix(gabacorrspec1)), diag(as.matrix(gabacorrspec2)),
                      diag(as.matrix(glutcorrspec11)), diag(as.matrix(glutcorrspec12)), diag(as.matrix(glutcorrspec21)),
                      diag(as.matrix(glutcorrspec22)), diag(as.matrix(gabacorrspec11)), diag(as.matrix(gabacorrspec12)), 
                      diag(as.matrix(gabacorrspec13)), diag(as.matrix(gabacorrspec21)), diag(as.matrix(gabacorrspec22)),
                      diag(as.matrix(gabacorrspec23)))
  colnames(dfsub2) <- c(sp1, sp2, 'overall_spec', 'glia_spec', 'glut_spec', 'gaba_spec',
                        'glia_spec1', 'glia_spec2', 'glut_spec1', 'glut_spec2',
                        'gaba_spec1', 'gaba_spec2', 'glut_spec11', 'glut_spec12',
                        'glut_spec21', 'glut_spec22', 'gaba_spec11', 'gaba_spec12',
                        'gaba_spec13', 'gaba_spec21', 'gaba_spec22', 'gaba_spec23')
  rownames(dfsub2) <- NULL
  
  
  return(cbind(dfsub1, dfsub2[,4:22]))
  
}