###############################################################
## plot coexpression conservation scores for a gene across   ##
## 21 animals and yeast                                      ##
###############################################################

library(corrplot)
library(viridis)

plot_coexpression_conservation_heatmap <- function(genename){
  
  # get heatmap across 22 species a given gene #

  my_palette = viridis::magma(10)
  
  spe = c('human', 'chimp', 'rhesus-macaque', 'crab-eating-macaque',
          'mouse', 'rat', 'rabbit', 'boar', 'cow', 'dog', 'horse', 'goat',
          'sheep', 'chicken', 'zebrafish', 'salmon', 'trout',
          'fruitfly', 'roundworm', 'bee', 'silkworm', 'yeast')
  
  bulkspec = read.delim('coexp_cons_across_22_species.csv', sep = ',')
  
  
  # get table for heatmap 
  scores = matrix(NA, nrow = length(spe), ncol = length(spe))
  scores[lower.tri(scores, diag = F)] <- unlist(bulkspec[bulkspec$gene==genename,-1])
  scores = t(scores)
  scores[lower.tri(scores, diag = F)] <- unlist(bulkspec[bulkspec$gene==genename,-1])
  rownames(scores) = spe
  colnames(scores) = spe
  
  
  # plot
  corrplot(scores, method = 'color', is.corr = F, col = my_palette, #col.lim = c(0,1),
           tl.col = 'black', na.label = ' ', na.label.col = 'white')
  mtext(genename, at = 9, line = 0.4, cex = 1.2)

  return(scores)
}