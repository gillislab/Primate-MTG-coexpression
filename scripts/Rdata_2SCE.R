Rdata_2SCE <- function(filepath){    
  # given as Rdata
  load(filepath, userDat <- new.env())
  flag = c(0,0)
  # Are we given SEO?
  isSEO =  lapply(userDat,class) == 'SummarizedExperiment'
  isSCE =  lapply(userDat,class) == 'SingleCellExperiment'
  if (sum(isSEO)>0) { # a SEO object is found
    SEO = userDat[[names(which( isSEO ) )[1]  ]]
    # convert to singlecellexperiment
    SCE = as(SEO , 'SingleCellExperiment')
    if (ncol(colData(SCE))==0 ) { 
      flag = c(1,0)
    } else{
      flag = c(1,1)
    }
  }
  # are we given SCE?
  if (sum(isSCE) > 0 ){ # an SCE object is found
    SCE = userDat[[names(which( isSCE ) ) ]]
    if (ncol(colData(SCE))==0 ) { 
      # shinyalert('colData field was empty.', type='error',closeOnClickOutside=TRUE) 
      flag = c(1,0)
    } else{
      flag=c(1,1)
    }
  }   
  # Are we given a matrix and a dataframe (for cells)?
  if (sum(flag)==0 ){ # nothing found yet - no SCE or SEO
    ismat = names(userDat) == 'mat'
    iscol = names(userDat) == 'colData'
    if (sum(ismat)> 0 ) {   # mat is named mat as it should be
      mat = userDat[[ 'mat' ]]
      flag=c(1,0)
    } else {        # else, we have to guess based on which is a matrix
      ismat = lapply(userDat ,class)%in%c('matrix','dgCMatrix')
      names(ismat) = names(userDat)
      if(sum(ismat ) > 0 ){
        mat = userDat[[ names(which( ismat >0 ) )  ]]
        flag=c(1,0)
      } else{
        shinyalert('Unable to find an expression matrix in uploaded data. Be sure that it is a matrix, SummarizedExperiment, or SingleCellExperiment class.',type= 'error',closeOnClickOutside=TRUE)
      }
    }
    # look for colData
    if (sum(flag) ==1){     # matrix found
      if (sum(iscol)>0 ) {
        colDat = userDat[['colData']]
        flag =c(1,1)
      } else{     # guess using which is a dataframe of the right size
        iscol = lapply(userDat ,class)%in% c('DataFrame','data.frame','DFrame')
        nCells = ncol(mat) 
        names(iscol)=names(userDat)
        # which dataframe is the right size?
        if (sum(iscol) > 0 ){
          for (i in length(names(iscol)[iscol]) ){
            colDat =  userDat[[names(which(iscol)==T)[i] ]]
            if (nrow(colDat) ==nCells ){
              flag = c(1,1)
              break
            }
          }
        } else {
          shinyalert('Unable to find cell types in uploaded data. Be sure that it is a dataframe.',type= 'error',closeOnClickOutside=TRUE)
        }
      }
      # convert to SCE
      if (sum(flag)==2 ){
        SCE = SingleCellExperiment(
          mat, 
          colData = colDat
        )
      } else{shinyalert('Unable to find cell types!')}
    }
  }
  return(SCE)
}