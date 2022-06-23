
#' Performe massive Delta analysis 
#' 
#' @description
#'`multiDelta()` performs the analysis of Delta for multiple raters against a goldstandard.
#' 
#'
#'@param data a data.frame whose first column is the goldstandard.
#'@param which.measure character string indicating the measure of interest to retrieve from the analysis. Valid options are:
#' \code{"Delta"} (returns global Delta and SE), \code{"Agreement"} (returns measure of agreement of each category and SE), 
#' \code{"Conformity"} (returns measure of conformity of each category and SE), 
#' \code{"Predictivity"} (returns measure of predictivity of each category and SE), 
#' and \code{"Consistency"} (returns measure of consistency of each category and SE).
#'@param tol the desired tolerance applied to find the root of the unknown constant B, needed to estimate the model parameters.
#'@param mxits the maximum numer of iterations applied to find the root of the unknown constant B, needed to estimate the model parameters.
#' 
#'@return An object of class \code{"multiDelta"}, which is a list of as many elements as measures selected.
#'@details A print method is available for \code{"multiDelta"} objects. 
#'The results can be reported as plain tex (tex = F) or LaTeX formatted (tex = T). 
#' In the latter case, the table can be transposed (transpose = T).
#'@export
#'@examples 
#' # Create a data.frame for 1 goldstandards and 9 raters
#' dat = data.frame(replicate(10, sample(1:3, 120, replace = TRUE)))
#' 
#' # Compute de Delta model and return the Consistency a Conformity measures
#' mDelta = multiDelta(dat, which.measure = c('Consistency', 'Conformity'))
#' print(mDelta, tex = TRUE, transpose = TRUE)




multiDelta = function(data, which.measure = c('Delta','Agreement', 'Conformity', 'Predictivity', 'Consistency'), tol = 1e-7, mxits = 100){
  
  simpleCap <- function(x) {
    unname(sapply(x, function(x){
      s <- strsplit(x, " ")[[1]]
      paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
            sep="", collapse=" ")
    }))

  }
  
  which.measure = simpleCap(which.measure)
  
  if(!all(which.measure%in%c('Delta','Agreement', 'Conformity', 'Predictivity', 'Consistency'))){
    stop('Argument "which.measure" should be one of "Delta", "Agreement", "Conformity", "Predictivity", "Consistency"')
  }

  
  nraters = ncol(data)
  names = paste(colnames(data)[1], colnames(data)[-1], sep = '.')
  
  # Compute delta (all against column 1)
  res = lapply(2:nraters, function(i){delta(data[,c(1,i)], standard = T, rawdata = T, tol = tol, mxits = mxits)})
  
  # Retrieve all measures
  all.measures = lapply(res, '[[', 'all.measures')
  names(all.measures) = names
  
  

    
  # Give format to list object containing all measures
  measures = list()
  B = list()
  A = lapply(all.measures, '[', c('Agreement', 'Conformity', 'Predictivity', 'Consistency'))
  
  delta_measure = lapply(lapply(all.measures, '[[', 'Global'), '[', 'Delta')
  delta_se = lapply(lapply(all.measures, '[[', 'Global'), attr, 'SE_Delta')

  for(i in 1:length(delta_measure)){
    B$Delta[[i]] = list(Delta = list(Measure = delta_measure[[i]]$Delta, 
                      SE_samplingI = delta_se[[i]]$SE_samplingI, 
                      SE_samplingII = delta_se[[i]]$SE_samplingII))
    
    measures[[i]] = append(B$Delta[[i]], A[[i]])
  }

  
  names(measures) = names
  

  # Get selected measures
  selectedMeasures = lapply(which.measure, function(i){
    lapply(measures, '[[',i)
  })
  names(selectedMeasures) = which.measure
  class(selectedMeasures) = c('multiDelta', 'deltaMAN')
  
  return(selectedMeasures)
}

#' Print object of class multiDelta
#' 
#' multiDelta contains the analysis of Delta of multiples raters against a standard.
#' @param x an object of class multiDelta
#' @param digits an integer indicating the significant digits to be used
#' @param tex a logical value indicating whether to generate formatted LaTeX output or plain output
#' @param transpose if tex = logical, the output table can be transposed
#' @importFrom xtable xtable print.xtable
#' @noRd

print.multiDelta = function(x, digits, tex = F, transpose = F){
  
  if(missing(digits)){
    digits = 4
  }
  
  which.measure = names(x)
  
  result = lapply(which.measure, function(i){
    
    obj = x[[i]]
    classes =  if(i == 'Delta'){
      'Global'
    }else{
      names(obj[[1]]$Measure)
    }

    
    res = .tableMultiDelta(obj, classes = classes, digits = digits)
    if(tex){
      cat(paste0("\\section*{",i,"}"))
     res =  .plain2latex(res, classes = classes, transpose = transpose)
    }
    res
  })
  
  if(!tex){
    names(result) = which.measure
    print(result, row.names = F)
  }

}


.tableMultiDelta = function(obj, classes, digits = 4){
  estimates = lapply(lapply(obj, '[[',1), function(x).round(unname(x),digits))
  
  
  # SE for sampling type I
  SE_1 = lapply(lapply(obj, '[[',2), function(x)if(!is.null(x)).round(unname(x),(digits+1)))
  
  # SE for sampling type II
  SE_2 = lapply(lapply(obj, '[[',3), function(x)if(!is.null(x)).round(unname(x),(digits+1)))
  
  
  
  A = lapply(SE_1, function(x){sapply(x, function(x)ifelse(is.null(x), '', paste0(' (', x, ')* ')))})
  lst1 = lapply(1:length(obj), function(i){paste0(estimates[[i]], A[[i]])})
  df = do.call(cbind, lst1)
  colnames(df) = names(SE_1)
  
  df1 = data.frame(Class = classes, df);df1
  
  
  # MUESTREO TIPO II:
  SE_2[which(sapply(SE_2, is.null))] = ' - '
  A = lapply(SE_2, function(x){sapply(x, function(x)ifelse(x == ' - ', ' ', paste0(' (', x, ')**')))});A
  
  df = do.call(cbind, A);df
  rownames(df) = NULL
  
  df2 = data.frame(Class = ' ', df);df2
  if(all(df2 == ' ')){
    # return(df1)
    df = df1
  }else{
    # CREAR DATA.FRAME CON LAS MEDIDAS DE ACUERDO POR CLASE
    df = data.frame()
    for(i in 1:length(classes)){
      df = rbind(df, df1[i,],df2[i,])
    }
    rownames(df) = NULL
  }
  
  return(df)

}


.plain2latex = function(obj, classes, transpose = F){
  
  A = do.call(cbind, lapply(obj, .p2pm));A
  
  raters_names = colnames(obj[-1])
  
  if(transpose){
    
    r = c()
    se = c()
    for(i in 2:ncol(A)){
      if(i%%2 != 0){
        next
      }
      
      r = rbind(r, A[,i])
      se= rbind(se,A[,i+1])
    }
    
    if(ncol(r)==length(classes)*2){
      n = ncol(r)
      se1 = se[,seq(1, n, 2), drop = F]
      se2 = se[,seq(2, n, 2), drop = F]
      
      SEs = c()
      for(i in 1:nrow(se)){
        SEs = rbind(SEs, se1[i,],se2[i,])
      }
      
      B = rep(r, each = 2)
      
      r = matrix(B, nrow = nrow(SEs))
      r[seq(2,nrow(SEs),2),] = ' '
      
      se = SEs
      r[,seq(2, ncol(r),2)] = se
      r
      multiR = 2
      
      raters_names_LTX = paste0("\\multirow{",multiR,"}{*}{",raters_names,"}")
      rowNams = rbind(raters_names_LTX, 
                 rep(" ", length(raters_names_LTX)))

    }else{
      multiR = 1
      r = cbind(r, se)[, order(c(seq(ncol(r)), seq(ncol(se))))]
      
      rowNams = paste0("\\multirow{",multiR,"}{*}{",raters_names,"}")
    }
    
    
    r = cbind(as.vector(rowNams), r)
    
    namesLTX = paste0("& \\multicolumn{2}{c}{",classes,"}")
    namesLTX[length(namesLTX)] = paste(namesLTX[length(namesLTX)], '\\\\\n')
    
    namesLTX0 = " "
    colNams = c(namesLTX0, namesLTX)
    
    addtorow <- list()
    addtorow$pos <- as.list(c(rep(0, length(colNams))))
    addtorow$command = c(colNams)
    
    
    print(xtable(r),  
          add.to.row = addtorow, include.colnames = FALSE, include.rownames=FALSE,
          floating = F, sanitize.text.function = identity, scalebox='0.8',comment=F)
  }else{
    namesLTX = paste0("& \\multicolumn{2}{c}{",raters_names,"}")
    namesLTX[length(namesLTX)] = paste(namesLTX[length(namesLTX)], '\\\\\n')
    
    namesLTX0 = "Class"
    colNams = c(namesLTX0, namesLTX)
    
    addtorow <- list()
    addtorow$pos <- as.list(c(rep(0, length(colNams))))
    addtorow$command = c(colNams)
    
    
    print(xtable(A), add.to.row = addtorow, include.colnames = FALSE,include.rownames=FALSE,
          floating = F, sanitize.text.function = identity, scalebox='0.8',comment=F)
  }
  
}
