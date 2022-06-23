

.p2pm = function(x){
  a = gsub('\\(','ww\\(', x)
  b = gsub(')', ')', a)
  d = strsplit(b, 'ww')
  d[lengths(d) ==0] = ''
  do.call(rbind, d)
}

.round = function(x, digits){
  format(round(x, digits), nsmall = digits)
}


#'Print object of class deltaMAN.
#'
#'
#' \code{print} method for class \code{"deltaMAN"}. 
#' @param x an object of class deltaMAN.
#' @param ... optional arguments passed to print for other classes created in the deltaMAN package. 
#' Currently, it supports \code{"digits"} (the significant digits to be used) 
#' and \code{"tex"} (a logical value indicating whether to generate formatted LaTeX output or plain output). 
#' Also, optional argument \code{"transpose"} is available for objects of class \code{"multiDelta"}.
#' @return No return value. Results are printed on console.
#' @export
print.deltaMAN = function(x, ...){
  print(x, ...)
}

# Informe Delta -----
#'Print object of class delta.
#'
#' \code{print} method for class \code{"delta"}. 
#' @noRd
print.delta = function(obj){
  
  messages = attr(obj$Data$Input, 'messages')
  
  if(length(messages)>0){
    cat(paste(messages, collapse = '\n'), '\n')
  }
  
  # if(obj$standard){
  #   cat('There is a gold standard classification by rows.\n')
  # }else{
  #   cat('None of the raters is a gold standard.\n')
  # }
  # 
  # if(obj$fixedRows){
  #   cat('Row marginal frequencies have been set beforehand (type II sampling).\n')
  # }else{
  #   cat('Row marginal frequencies have not been set beforehand (type I sampling).\n')
  # }
  
  
  
  cat('\nDelta = ', obj$Delta$Estimates$Delta, ', SE = ', obj$Delta$SE$Delta, sep = '')
}



#' Print object of class defP
#' 
#' defP contains the definition of the problem to be solved in order to estimate Delta
#' @param obj an object of class delP
#' @param digits an integer indicating the significant digits to be used
#' @param tex a logical value indicating whether to generate formatted LaTeX output or plain output
#' @importFrom xtable xtable print.xtable
#' @noRd
print.defP <- function(obj, digits, tex = F){
  if(missing(digits)){
    digits = 4
  }
  
  X = obj$info
  x = c(X$h, X$B0, X$PType, X$yM, X$dyM, X$yP, X$dyP, X$ns, X$cr2x)
  x = .round(x, digits)
  names(x) = c('h', 'B0', 'PType', 'Y-', 'dY-', 'Y+', 'dY+', 'n-s', 'c+r-2x')
  
  if(tex ==F){
    
    cat(' B = ', .round(obj$B$root, digits), ' (',obj$B$iterations,' iterations)\n\n', sep = '')
    
    print(x, quote = F)
  }else{
    
    cat('B = ', .round(obj$B$root, digits), ' (',obj$B$iterations,' iterations)\n\n', sep = '')
    
    m = matrix(x, ncol = 9, dimnames = list(c(''), c('h', 'B0', 'PType', '$Y^-$', '$dY^-$', '$Y^+$', '$dY^+$', '$n-s$', '$c+r-2x$')))
    print(xtable(m),  sanitize.text.function=identity, comment=F, floating = F)
  }
  
  
}



#' Print object of class psetDelta
#' 
#' psetDelta contains the estimated parameters of the model computed by get_parameters().
#' @param obj an object of class psetDelta
#' @noRd

print.psetDelta <- function(obj){
  classes = names(obj$partial_delta)
  Delta = as.vector(unname(obj$partial_delta))
  Pi = as.vector(unname(obj$proportions))
  
  print(data.frame(Class = classes, Delta, Pi))
  
  cat('Overall Delta: ', obj$Delta)
}



#' Print object of class chisqTestDelta
#' 
#' chisqTestDelta computes the goodness of fit of Delta model
#' @param obj an object of class chisqTestDelta
#' @noRd

print.chisqTestDelta <- function(obj){
  cat('Observed frequencies: \n')
  print(obj$Observed)
  
  cat('\n Expected frequencies: \n')
  print(obj$Expected)
  
  cat('\n Chi-square Goodness of fit:')
  if(attr(obj, 'validity')==T){
    cat('\n','X-squared = ', obj$statistic, ', df = ', obj$df, ', p-value = ', obj$pval)
  }else{
    cat('\n',attr(obj, 'validity'))
  }
  
}




#' Print object of class cov.Delta
#' 
#' cov.Delta contains the covariance matrices of the partial deltas and proportions of responses at random.
#' @param obj an object of class delP
#' @param digits an integer indicating the significant digits to be used
#' @param tex a logical value indicating whether to generate formatted LaTeX output or plain output
#' @importFrom xtable xtable print.xtable
#' @noRd

print.cov.Delta <- function(obj, digits, tex = F){
  names(dimnames(obj$Cov_Delta)) = NULL
  names(dimnames(obj$cov_Mix)) = NULL
  names(dimnames(obj$cov_Pi)) = NULL
  
  if(missing(digits)){
    digits = 4
  }
  
  if(tex == F){
    cat('-- Cov(D[i], D[j]) -- \n')
    print(.round(obj$Cov_Delta, digits), quote = F)
    
    
    cat('\n-- Cov(D[i], Pi[j]) -- \n')
    print(.round(obj$cov_Mix, digits), quote = F)
    
    
    cat('\n-- Cov(Pi[i], Pi[j]) -- \n')
    print(.round(obj$cov_Pi, digits), quote = F)
  }else{
    
    k = ncol(obj$Cov_Delta)
    
    # cat('-- $Cov(D_i, D_j)$ -- \n\n')
    colNams = paste0("& \\multicolumn{",k,"}{c}{$Cov(\\hat\\Delta_i,\\hat\\Delta_j)$}\\\\\n")
    addtorow <- list()
    addtorow$pos <- as.list(c(rep(-1, length(colNams))))
    addtorow$command = c(colNams)
    # print(xtable(obj$Cov_Delta), floating = F,comment=F)
    print(xtable(obj$Cov_Delta, digits = digits), add.to.row = addtorow, 
          floating = F, sanitize.text.function = identity, scalebox='0.9',comment=F)
    
    
    
    # cat('\n\n-- $Cov(D_i, Pi_j)$ -- \n\n')
    colNams = paste0("& \\multicolumn{",k,"}{c}{$Cov(\\hat\\Delta_i,\\hat\\pi_j)$}\\\\\n")
    addtorow <- list()
    addtorow$pos <- as.list(c(rep(-1, length(colNams))))
    addtorow$command = c(colNams)
    print(xtable(obj$cov_Mix,digits = digits), add.to.row = addtorow, 
          floating = F, sanitize.text.function = identity, scalebox='0.9',comment=F)
    
    
    # cat('\n\n-- $Cov(Pi_i, Pi_j)$ -- \n\n')
    # print(xtable(obj$cov_Pi), floating = F,comment=F)
    
    colNams = paste0("& \\multicolumn{",k,"}{c}{$Cov(\\hat\\pi_i,\\hat\\pi_j)$}\\\\\n")
    addtorow <- list()
    addtorow$pos <- as.list(c(rep(-1, length(colNams))))
    addtorow$command = c(colNams)
    print(xtable(obj$cov_Pi,digits = digits), add.to.row = addtorow, 
          floating = F, sanitize.text.function = identity, scalebox='0.9',comment=F)
    
  }
  
  
}




#'Print all measures of agreement
#' 
#' Print object of class measuresDelta
#' 
#' measuresDelta contains all the measures of agreement of a given contingency table, computed by get_measures().
#' @param obj an object of class measuresDelta
#' @param digits an integer indicating the significant digits to be used
#' @param tex a logical value indicating whether to generate formatted LaTeX output or plain output
#' @importFrom xtable xtable print.xtable
#' @noRd

print.measuresDelta <- function(obj, digits, tex = F){
  if(missing(digits)){
    digits = 4
  }
  classes = names(obj$Agreement$Measure)
  
  case = attr(obj, 'DA_type')
  # Delta parameters
  params = obj$Global
  Delta = .round(as.vector(unname(params$partial_delta)), digits)
  Pi = .round(as.vector(unname(params$proportions)), digits)
  
  SE_delta = attr(params, 'SE_Delta')
  if(!all(sapply(SE_delta, is.null))){ 
    SE_delta = sapply(SE_delta, .round, digits = (digits+1))
    SE_1_delta = SE_delta[1]
    SE_2_delta = SE_delta[2]
  }else{
    SE_1_delta = NULL
    SE_2_delta = NULL
  }
  
  
  # obj2 = obj
  obj = obj[-1]
  measures = lapply(lapply(obj, '[[',1), function(x).round(unname(x),digits))
  
  
  # SE for sampling type I
  SE_1 = lapply(lapply(obj, '[[',2), function(x)if(!is.null(x)).round(unname(x),(digits+1)))
  
  # SE for sampling type II
  SE_2 = lapply(lapply(obj, '[[',3), function(x)if(!is.null(x)).round(unname(x),(digits+1)))
  
  
  # MUESTREO TIPO I:
  A = lapply(SE_1, function(x){sapply(x, function(x)ifelse(is.null(x), '', paste0(' (', x, ')* ')))})
  lst1 = lapply(1:4, function(i){paste0(measures[[i]], A[[i]])})
  df = do.call(cbind, lst1)
  colnames(df) = names(SE_1)
  df1 = data.frame(Class = classes, Delta, Pi, df);df1
  
  se1 = ifelse(is.null(SE_1_delta), '', paste0(' (',SE_1_delta,')* '))
  
  # MUESTREO TIPO II:
  SE_2[which(sapply(SE_2, is.null))] = ' - '
  A = lapply(SE_2, function(x){sapply(x, function(x)ifelse(x == ' - ', ' ', paste0(' (', x, ')**')))});A
  
  df = do.call(cbind, A);df
  rownames(df) = NULL
  df2 = data.frame(Class = classes, Delta = '', Pi = '', df)
  df2[1] = ''
  
  # CREAR DATA.FRAME CON LAS MEDIDAS DE ACUERDO POR CLASE
  df = data.frame()
  
  for(i in 1:length(classes)){
    df = rbind(df, df1[i,],df2[i,])
  }
  rownames(df) = NULL
  
  
  # CREAR DATA.FRAME CON DELTA GLOBAL
  se2 = ifelse(is.null(SE_2_delta), '', paste0('(',SE_2_delta,')**'))
  df_overall = data.frame('Overall Delta (SE):', .round(params$Delta, digits), se1)
  df_overall[2,] = c('','',se2)
  colnames(df_overall)=NULL
  
  
  
  
  if(tex){
    A = do.call(cbind, lapply(df, .p2pm));A
    n = ncol(A)
    measure_names = colnames(df[-c(1:3)])
    namesLTX = paste0("& \\multicolumn{2}{c}{",measure_names,"}")
    namesLTX[length(namesLTX)] = paste(namesLTX[length(namesLTX)], '\\\\\n')
    
    namesLTX0 = paste("& ", colnames(df)[1:3])
    namesLTX0[1] = "Class"
    colNams = c(namesLTX0, namesLTX)
    
    # se1 = ifelse(is.null(SE_1_delta), '', paste0('&$\\pm$',SE_1_delta,'*'))
    # se2 = ifelse(is.null(SE_2_delta), '', paste0('&$\\pm$',SE_2_delta,'**'))
    
    se1 = ifelse(is.null(SE_1_delta), '', paste0(' & (',SE_1_delta,')*'))
    se2 = ifelse(is.null(SE_2_delta), '', paste0(' & (',SE_2_delta,')**'))
    overall = c(paste0('\\hline \\multicolumn{3}{c}{Overall Delta:} &', .round(params$Delta, digits), se1, 
                       paste(rep('&',n-5),collapse = ''), '\\\\\n'),
                paste0(paste(rep('&',3),collapse = ''), se2, paste(rep('&',n-5),collapse = ''), '\\\\\n'))
    r = nrow(A)
    addtorow <- list()
    addtorow$pos <- as.list(c(rep(0, length(colNams)), rep(r, length(overall))))
    addtorow$command = c(colNams, overall)
    
    # addtorow$pos <- list(0, 0, 0,0,0,0,0,r,r)
    # addtorow$command = c(namesLTX0, namesLTX, overall)
    
    print(xtable(A), add.to.row = addtorow, include.colnames = FALSE,include.rownames=FALSE,
          floating = F, sanitize.text.function = identity, scalebox='0.95',comment=F)
    
  }else{
    print(df, row.names = F)
    cat(paste(rep('_', sum(nchar(df[1,]))+20), collapse = ''))
    
    print(df_overall, row.names = F)
  }
  
  
  
  if(is.null(case) || case != 'DA0'){
    if(tex ==F){
      cat(paste(rep('_', sum(nchar(df[1,]))+20), collapse = ''), sep = '')
    }
    
    cat('\n','* Valid under type-I sampling (columns and rows totals are at random)')
    if(tex){
      cat("\\newline")  
    }
    cat('\n** Valid under type-II sampling (rows totals are fixed beforehand)\n')
  }else{
    cat('At least one estimation falls on the boundary of the parametric space. \nFor this reason, the SE values have no been obtained.\n')
    
  }
}




#' Print object of class validMeasuresDelta
#' 
#' validMeasuresDelta contains the valid measures of agreement under selected model (sampling type and goldstandard), computed by validMeasures().
#' @param obj an object of class validMeasuresDelta
#' @param digits an integer indicating the significant digits to be used
#' @param tex a logical value indicating whether to generate formatted LaTeX output or plain output
#' @importFrom xtable xtable print.xtable
#' @noRd

print.validMeasuresDelta = function(obj, digits = 4, tex = F){
  
  classes = names(obj$Estimates$partial_delta)
  
  # Delta parameters
  params = obj$Estimates
  Delta = .round(as.vector(unname(params$partial_delta)), digits)
  Pi = .round(as.vector(unname(params$proportions)), digits)
  
  measures = obj$Estimates[-c(1:3)]
  SE = obj$SE[-1]
  n = length(measures)
  A = lapply(SE, function(x){sapply(x, function(x)ifelse(is.null(x), '', paste0(' (', .round(x,(digits+1)), ')')))})
  lst1 = lapply(1:n, function(i){paste0(.round(measures[[i]], digits), A[[i]])})
  
  df = do.call(cbind, lst1)
  colnames(df) = names(SE)
  df1 = data.frame(Class = classes, Delta, Pi, df)
  se = ifelse(is.null(obj$SE$Delta), '', paste0(' (',.round(obj$SE$Delta, (digits+1)),')'))
  

  
  if(tex){
    
    
    A = do.call(cbind, lapply(df1, .p2pm));A
    n = ncol(A)
    measure_names = colnames(df1[-c(1:3)])
    namesLTX = paste0("& \\multicolumn{2}{c}{",measure_names,"}")
    namesLTX[length(namesLTX)] = paste(namesLTX[length(namesLTX)], '\\\\\n')
    
    namesLTX0 = paste("& ", colnames(df1)[1:3])
    namesLTX0[1] = "Class"
    colNams = c(namesLTX0, namesLTX)
    # se = ifelse(is.null(obj$SE$Delta), '', paste0('&$\\pm$',.round(obj$SE$Delta, (digits+1))))
    se = ifelse(is.null(obj$SE$Delta), '', paste0(' & (',.round(obj$SE$Delta, (digits+1)),")"))
    
    overall = c(paste0('\\hline \\multicolumn{3}{c}{Overall Delta:} &', .round(params$Delta, digits), se, 
                       paste(rep('&',n-5),collapse = ''), '\\\\\n'))
    r = nrow(A)
    addtorow <- list()
    addtorow$pos <- as.list(c(rep(0, length(colNams)), rep(r, length(overall))))
    addtorow$command = c(colNams, overall)
    
    print(xtable(A), add.to.row = addtorow, include.colnames = FALSE,include.rownames=FALSE,
          floating = F, sanitize.text.function = identity, scalebox='0.95',comment=F)
    
    
    
  }else{
    print(df1, row.names=F)
    
    # cat('Overall Delta (SE): ', .round(params$Delta, digits), se,'\n', sep = '')
    df_overall = data.frame('Overall Delta (SE): ', .round(params$Delta, digits), se)
    colnames(df_overall) = NULL
    print(df_overall, row.names = F)
  }
  
  
  case = attr(obj, 'DA_type')
  if(!is.null(case) && case == 'DA0'){
    cat('At least one estimation falls on the boundary of the parametric space. \nFor this reason, the SE values have no been obtained.\n')
  }
}