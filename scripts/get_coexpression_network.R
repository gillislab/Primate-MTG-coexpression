###############################################################
## create coexpression networks from gene expression dataset ##
###############################################################

# build SCE from species Rdata files
source('Rdata_2SCE.R')
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater)
library(matrixStats)

## create bulk samples from cluster-specific sc data ##
get_bulk_samples <- function(mat, number_of_samples){
  nsamples = round(dim(mat)[2]/number_of_samples)    # divide into samples of ~20 cells each
  bulk_samples = matrix(NA, nrow = dim(mat)[1], ncol = nsamples)
  
  for(ii in 1:(nsamples-1)){      
    ids = sample(dim(mat)[2], number_of_samples, replace = FALSE)
    bulk_samples[,ii] = rowSums(mat[,ids])
    mat <- mat[,-ids]        
  }
  bulk_samples[,nsamples] = rowSums(mat)    
  rownames(bulk_samples) = rownames(mat)
  return(bulk_samples)
}

## create coexpression network ## 
build_coexp_network <- function(net, method = "spearman", flag = "rank"){
  
  # Calculate correlation coefficients
  genes = rownames(net)
  net = net[rowSums(net) > 0, ]
  
  #     print(paste0("... ",dim(net)[1]," of ",length(genes), " genes have non-zero expression" ))
  
  net = cor(t(net), method = method)
  
  # Create network
  temp = net[upper.tri(net, diag = T)]
  if (flag == "abs"){
    temp = abs(temp)
  }
  
  #     print("...ranking")
  temp = rank(temp, ties.method = "average")  
  
  #     print("...reconstructing matrix")
  net[upper.tri(net, diag = T)] = temp
  
  net = t(net)
  net[upper.tri(net, diag = T)] = temp
  
  net = net / max(net, na.rm = T)
  med = median(net, na.rm = T)
  
  ind = setdiff(genes, rownames(net))   # which genes are missing?
  
  temp = matrix(med, length(ind), dim(net)[2])
  rownames(temp) = ind
  
  net = rbind(net, temp)
  temp = matrix(med, dim(net)[1], length(ind))
  colnames(temp) = ind
  net = cbind(net, temp)
  
  # reorder to original
  net = net[genes, genes]
  diag(net) = 1
  
  return(net)
}

# rank-standardize
get_rankstd <- function(aggNet){
  temp = aggNet[upper.tri(aggNet, diag = T)]
  temp = rank(temp, ties.method = "average")  
  aggNet[upper.tri(aggNet, diag = T)] = temp
  aggNet = t(aggNet)
  aggNet[upper.tri(aggNet, diag = T)] = temp
  aggNet = aggNet / max(aggNet, na.rm = T)
  return(aggNet)
}


# get coexpression network from single cell data
get_coexpression_network <- function(species_name){
  
  # get list of 14131 1-1 orthologs across all 5 species
  common_genes = read.delim('species5_common_OrthoDb_1-1.csv', sep = ',', header = T, stringsAsFactors = F)
  clstab = read.delim('cross_species_cluster_aggregated_new.txt', sep = ',', header = T, stringsAsFactors = F)
  
  # 57 common clusters
  cls_groups <- c('glia|9', 'glia|1', 'glia|12', 'glia|13', 'glia|7', 'glia|5', 'glia|15',
                  'it_types|10', 'it_types|4', 'it_types|9', 'it_types|1', 'it_types|7',
                  'it_types|16', 'it_types|13', 'it_types|12', 'it_types|19', 'it_types|17',
                  'l5et_l56np_l6ct_l6b|1', 'l5et_l56np_l6ct_l6b|2', 'l5et_l56np_l6ct_l6b|3', 'l5et_l56np_l6ct_l6b|4',
                  'l5et_l56np_l6ct_l6b|7', 'l5et_l56np_l6ct_l6b|8', 'l5et_l56np_l6ct_l6b|9', 'l5et_l56np_l6ct_l6b|10', 'l5et_l56np_l6ct_l6b|13',
                  'sst_sst_chodl_pvalb|10', 'sst_sst_chodl_pvalb|12', 'sst_sst_chodl_pvalb|1', 'sst_sst_chodl_pvalb|3', 'sst_sst_chodl_pvalb|6',           
                  'sst_sst_chodl_pvalb|32', 'sst_sst_chodl_pvalb|30', 'sst_sst_chodl_pvalb|15', 
                  'sst_sst_chodl_pvalb|18', 'sst_sst_chodl_pvalb|22', 'sst_sst_chodl_pvalb|20', 'sst_sst_chodl_pvalb|13',
                  'sst_sst_chodl_pvalb|25', 'sst_sst_chodl_pvalb|26', 'sst_sst_chodl_pvalb|29',
                  'lamp5_sncg_vip|3', 'lamp5_sncg_vip|4', 'lamp5_sncg_vip|2', 'lamp5_sncg_vip|5', 'lamp5_sncg_vip|6', 'lamp5_sncg_vip|7',
                  'lamp5_sncg_vip|11', 'lamp5_sncg_vip|9', 'lamp5_sncg_vip|10', 'lamp5_sncg_vip|15', 'lamp5_sncg_vip|27', 'lamp5_sncg_vip|26', 'lamp5_sncg_vip|25',
                  'lamp5_sncg_vip|22', 'lamp5_sncg_vip|24', 'lamp5_sncg_vip|17')
  
  # load processed scRNA-seq data
  spe = c('human', 'chimp', 'gorilla', 'rhesus', 'marmoset')

  sc1 = Rdata_2SCE(paste0('/data/suresh/allen/', species_name, '_data_new.Rdata'))
  sc1 <- sc1[,sc1$exclude_in_final=='no']
  assayNames(sc1) <- 'counts'
  cpm(sc1) <- calculateCPM(sc1)
  
  
  # get a list of highly variable genes for downstream analyses
  upper_limit = ceiling((dim(sc1)[2])/1e4)
  pb = txtProgressBar(min = 0, max = upper_limit, initial = 0)
  qmat = matrix(0, nrow = dim(sc1)[1])
  
  for(yy in 1:upper_limit){
    
    if(yy<upper_limit){
      startid = 1 + 10000*(yy-1)
      endid = 10000*yy
    }else{
      startid = endid + 1
      endid = dim(sc1)[2]
    }
    
    temp1 = as.matrix(cpm(sc1)[,startid:endid])
    qq = colQuantiles(temp1, probs = 0.8)
    tempqq = matrix(rep(qq, each = dim(sc1)[1]), nrow = dim(sc1)[1])
    
    # aggregate counts to see how many genes are expressed in no. of cells
    qmat = qmat + rowSums(temp1>tempqq)
    
    if(yy == 1){
      qnts = qq
    }else{
      qnts = c(qnts, qq)
    }
    
    setTxtProgressBar(pb, yy)
  }
  
  rankqmat = rank(-qmat)
  
  
  # threshold to top 4.5k genes expressed in at least 25k cells
  geneids = which(rankqmat<=4500)
  sc1 = sc1[geneids,]
  
  # build bulk samples, and coexpression networks from pseudo-bulk samples
  
  # can also be adapted to calculate networks for each donor/tech
  # and then aggregated
  pb = txtProgressBar(min = 0, max = length(cls_groups), initial = 0)
  agg_coexp_net = matrix(0, nrow = 4500, ncol = 4500)
  
  for(oo in 1:length(cls_groups)){
    currcls = cls_groups[oo]
    temp1 = cpm(sc1)[,sc1$cross_species_cluster==currcls]
    
    # get pseudo-bulk samples with 20 cells each
    bulk1 = get_bulk_samples(temp1, 20) 
    cor1 = build_coexp_network(bulk1, method = 'spearman')
    
    agg_coexp_net = agg_coexp_net + cor1
  }
  
  
  # rank-standardize
  temp = agg_coexp_net[upper.tri(agg_coexp_net, diag = T)]
  temp = rank(temp, ties.method = "average")  
  agg_coexp_net[upper.tri(agg_coexp_net, diag = T)] = temp
  agg_coexp_net = t(agg_coexp_net)
  agg_coexp_net[upper.tri(agg_coexp_net, diag = T)] = temp
  agg_coexp_net = agg_coexp_net / max(agg_coexp_net, na.rm = T)
  
  
  return(agg_coexp_net)
}
