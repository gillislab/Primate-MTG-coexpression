###############################################################
## get coexpression neighborhood conservation for 14k genes  ## 
## across a pair of primates                                 ##
###############################################################

# this script can also be adapted to calculate coexpression neighborhood
# conservation of 1-1 orthologs across any pair of networks

get_top_pairs_mat <- function (mat.rank, n){ 
  f = mat.rank>dim(mat.rank)[1]-n
  mat.pairs = mat.rank*0
  mat.pairs[f] = 1
  return(mat.pairs)
}    

auroc_analytic2 <-function (mat, np, nL) { 
  ranks = 1:nL
  mini = sum(ranks[1:np])
  maxi = sum(ranks[(nL-np+1):nL])    
  range = maxi - mini    
  auroc <- (mat - mini)/range
  return(auroc)
}

roc_bg_predict <- function (bg_predict_res, species1, species2, np, nL){
  rocs = matrix(0, ncol = 2, nrow = dim(bg_predict_res[[1]])[2])
  colnames(rocs) = c(species1,species2)
  
  ranks = 1:nL
  mini = sum(ranks[1:np])
  maxi = sum(ranks[(nL-np+1):nL])    
  range = maxi - mini
  
  for(i in 1:2){
    temp = t(apply(bg_predict_res[[i]], 1, function(x) rank(x, ties.method = "average")))        
    rocs[,i] <- (diag(temp) - mini)/range       
  }
  return(rocs)
}

get_coexpression_conservation <- function(sp1, sp2){
  
  spe = c('human', 'chimp', 'gorilla', 'rhesus', 'marmoset')
  spe_bulk = c('human', 'chimp', 'gorilla', 'rhesusm', 'marmoset')
  path0 = '/data/suresh/allen/coexp_networks/aggregate_nets/'
  path1 = '/data/suresh/allen/coexp_networks/'
  
  # get list of 14131 1-1 orthologs across all 5 species
  common_genes = read.delim('species5_common_OrthoDb_1-1.csv', sep = ',')
  
  
  # subset to common genes
  load(paste0(path0, sp1, '_cls_aggregate_top4500_new.Rdata'))
  diag(aggNet) = 0
  coexp_species1 = aggNet    
  
  load(paste0(path0, sp2, '_cls_aggregate_top4500_new.Rdata'))
  diag(aggNet) = 0
  coexp_species2 = aggNet
  
  cids = intersect(match(rownames(coexp_species1), common_genes[,sp1]), match(rownames(coexp_species2), common_genes[,sp2]))
  genelist1 = common_genes[cids, sp1]
  genelist2 = common_genes[cids, sp2]
  coexp_species1 <- coexp_species1[genelist1, genelist1]
  coexp_species2 <- coexp_species2[genelist2, genelist2]
  
  # calculate ranks and top 10s
  rank_species1 = apply(coexp_species1, 1, rank, ties.method = 'average')
  rank_species2 = apply(coexp_species2, 1, rank, ties.method = 'average')
  
  mat1 = get_top_pairs_mat(t(rank_species1), 10)
  mat2 = get_top_pairs_mat(t(rank_species2), 10)
  
  gg_sp12 = mat1%*%rank_species2    # sp1 genes predict sp2 genes this well
  gg_sp21 = mat2%*%rank_species1    # sp2 genes predict sp1 genes this well
  
  res.mat1 = auroc_analytic2(gg_sp12, 10, dim(mat1)[2])
  res.mat2 = auroc_analytic2(gg_sp21, 10, dim(mat2)[2])
  res = vector("list", length = 2)
  res[[1]] = res.mat1
  res[[2]] = res.mat2
  
  
  # function conservation scores - "spec_aurocs" #
  fncons_aurocs = matrix(0, ncol = 2, nrow = dim(res[[1]])[2])
  colnames(fncons_aurocs) = c(sp2, sp1)
  fncons_aurocs[,1] = diag(res.mat1)
  fncons_aurocs[,2] = diag(res.mat2)
  
  spec_aurocs = roc_bg_predict(res, sp2, sp1, 1, dim(mat1)[1])
  
  geneList = cbind(rownames(coexp_species1), rownames(coexp_species2))
  colnames(geneList) = c(sp1, sp2)
  
  roc_scores = data.frame(geneList, rowMeans(spec_aurocs))
  colnames(roc_scores) <- c('sp1_gene', 'sp2_gene', 'coexp_cons')
  
  return(roc_scores)
}
