###############################################################
## get cluster-level pseudo-bulk data from species-specific  ##
## snRNA-seq MTG datasets                                    ##
###############################################################

# build SCE from species Rdata files
library(MetaMarkers)
source('Rdata_2SCE.R')
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater)

# get pseudo-bulk datasets from CPM data
get_cpm_pseudo_bulk <- function(scdata, new_labels, nbds){
  sc_sub2 = matrix(NA, nrow = dim(scdata)[1], ncol = length(nbds))
  
  for (i in 1:length(nbds)){
    cols = which(new_labels==nbds[i])
    
    if(length(cols)>1){
      sc_sub2[,i] = rowSums(cpm(scdata)[,cols])/length(cols)
    }else if(length(cols)==1){            
      sc_sub2[,i] = cpm(scdata)[,cols]
    }else{
      sc_sub2[,i] = NA
    }
  }    
  rownames(sc_sub2) = rownames(scdata)
  colnames(sc_sub2) = nbds
  return(sc_sub2)
}

# get pseudo-bulk data from snRNA-seq MTG datasets for each primate
get_pseudobulk_data <- function(species_name){
  
  # 57 consensus cluster labels
  cls_groups <- c('glia|9', 'glia|1', 'glia|12', 'glia|13', 'glia|7', 'glia|5', 'glia|15',
                  'lamp5_sncg_vip|3', 'lamp5_sncg_vip|4', 'lamp5_sncg_vip|2', 'lamp5_sncg_vip|5', 'lamp5_sncg_vip|6', 'lamp5_sncg_vip|7',
                  'lamp5_sncg_vip|11', 'lamp5_sncg_vip|9', 'lamp5_sncg_vip|10', 'lamp5_sncg_vip|15', 'lamp5_sncg_vip|27', 'lamp5_sncg_vip|26',
                  'lamp5_sncg_vip|25', 'lamp5_sncg_vip|22', 'lamp5_sncg_vip|24', 'lamp5_sncg_vip|17',
                  'sst_sst_chodl_pvalb|10', 'sst_sst_chodl_pvalb|12', 'sst_sst_chodl_pvalb|1', 'sst_sst_chodl_pvalb|3', 'sst_sst_chodl_pvalb|6',           
                  'sst_sst_chodl_pvalb|32', 'sst_sst_chodl_pvalb|30', 'sst_sst_chodl_pvalb|15', 
                  'sst_sst_chodl_pvalb|18', 'sst_sst_chodl_pvalb|22', 'sst_sst_chodl_pvalb|20', 'sst_sst_chodl_pvalb|13',
                  'sst_sst_chodl_pvalb|25', 'sst_sst_chodl_pvalb|26', 'sst_sst_chodl_pvalb|29',
                  'l5et_l56np_l6ct_l6b|1', 'l5et_l56np_l6ct_l6b|2', 'l5et_l56np_l6ct_l6b|3', 'l5et_l56np_l6ct_l6b|4',
                  'l5et_l56np_l6ct_l6b|7', 'l5et_l56np_l6ct_l6b|8', 'l5et_l56np_l6ct_l6b|9', 'l5et_l56np_l6ct_l6b|10', 'l5et_l56np_l6ct_l6b|13',
                  'it_types|10', 'it_types|4', 'it_types|9', 'it_types|1', 'it_types|7',
                  'it_types|16', 'it_types|13', 'it_types|12', 'it_types|19', 'it_types|17')
  
  # get species name
  scdata = Rdata_2SCE(paste0('/data/suresh/allen/', species_name, '_data_new.Rdata'))
  scdata <- scdata[,scdata$exclude_in_final=='no']  # remove poorly annotated cells
  assay(scdata, "cpm") = convert_to_cpm(assay(scdata))
  
  # create pseudo-bulk matrices of 14k orthologs x 57 cell types for each species
  scdata_bulk = get_cpm_pseudo_bulk(scdata, scdata$cross_species_cluster, cls_groups)

  return(scdata_bulk)
}