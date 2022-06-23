
#' Compute the Delta coefficient
#' 
#' @description
#'`delta()` computes Delta coefficient, or proportion of agreements that are not due to chance, 
#' which is used to measure nominal agreement between two raters.
#' 
#'`delta()` and `Delta()` are synonyms.
#'
#'@references Andrés, A. M., & Marzo, P. F. (2004). Delta: A new measure of agreement between two raters. British journal of mathematical and statistical psychology, 57(1), 1-19.
#'@references Andrés, A. M., & Marzo, P. F. (2005). Chance-corrected measures of reliability and validity in KK tables. Statistical methods in medical research, 14(5), 473-492.
#'@param data either a contingency table or raw data.
#'@param standard a logical value indicating whether the observer on the rows of the contingency table (or first column of raw data) is a goldstandard (i.e., gives the correct responses).
#'@param fixedRows a logical value indicating whether the row marginals are fixed in advance (sampling type II) or not (sampling type I).
#'@param rawdata a logical value indicating whether the data is raw (TRUE) or a contingency table (FALSE). If not specified, the function will try to guess the data type.
#'@param tol the desired tolerance applied to find the root of the unknown constant B, needed to estimate the model parameters.
#'@param mxits the maximum numer of iterations applied to find the root of the unknown constant B, needed to estimate the model parameters.
#' 
#'@return An object of class \code{"delta"}, which is a list of 9 elements (or 10 if the dimension of the contingency table is 2x2). See details.
#'
#'@details The allowed input data type are (1) contingency tables (of class \code{"table"} or \code{"matrix"}) or (2) raw data (of class \code{"data.frame"}).
#' If the data is of type (1), the empty classes (if any) are removed. If the data is of class (2), the function checks the number of columns (n) and throws an
#' error if n < 2 or n > 3. If n = 2, a frequency table is computed. If n = 3, it is assumed that one column represents the row index (normally, it is expected to be the first one).
#' Once the row index column is identified, it is removed and a frequency table is computed using the remaining two columns. 
#' In all cases, the result is always a squared matrix, which will be used in the subsequent computation of the Delta coefficient.
#' The observer on the rows will be referred to as observer (or rater) R and the one on the columns will be referred to as observer (or rater) C.
#' 
#' The function returns a list of 9 elements (if the number of classes is >2):
#' * Delta: This is a list of 2 elements: 
#' (1) the estimates of the model: overall delta (Delta), partial delta for class j (partial_delta) and the distribution of responses made at random by observer C (proportions); 
#' and the agreement measurements: agreement, conformity, predictivity and consistency (only some of the are shown, depending on the model assumed)
#' (2) the standard error of the estimates under the model assumed (sampling type I or II).
#' * Kappa: The estimate and standard error of the Cohen's Kappa coffiecient.
#' * Data: Input data and analyzed data (may be the same)
#' * GOF: Goodness of Fit for the Delta model. The chi squares statistic is computed. If the performed test is significative (p-value < alpha), the model Delta is not suitable for the data.
#' * fixedRows: logical value that matches the "fixedRows" argument.
#' * standard: logical value that matches the "standard" argument.
#' * all.measures: This is a list including all the estimates and standard errors (disregarding the model assumed by the user).
#' * problem.parameters: This list contains information about the estimation of the auxiliary constant B, needed to estimate the model parameters.
#' * cov: This list contains 3 elements: the covariance matrix of the partial delta estimates; the covariance matrix of the proportion estimates; and the covariance matrix of the partial delta and proportion estimates.
#' 
#' If the number of classes is k = 2, another element is added to the aforementioned list, including the asymptotic analysis (asymptoticDelta).
#' @export
#' @seealso [summary.delta()] for the summary method created for objects of class delta, and [print.deltaMAN()] for the print method.
#' @examples 
#' # Create a 3x3 matrix
#' m = matrix(c(15, 5, 0, 4, 21, 1, 3, 4, 25), ncol = 3)
#' # Compute the Delta coefficient assuming the rater on the rows 
#' # is a goldstandard and type II sampling.
#' obj = delta(m, standard = TRUE, fixedRows = TRUE)
#' # Get the complete report
#' summary(obj, fullReport = TRUE)
#' 
#' # Create a 2x2 matrix
#' m = matrix(c(15, 7, 3, 21), ncol = 2)
#' # Compute the Delta coefficient assuming no one is a goldstandard 
#' # and type I sampling.
#' obj = delta(m, standard = FALSE, fixedRows = FALSE)
#' # Get the report
#' summary(obj, fullReport = FALSE)
#' 
delta = function(data, standard = FALSE, fixedRows = FALSE, rawdata = NULL, tol = 1e-7, mxits = 100){

  # obtener matriz de contingencia depurada
  M1 = getM1(data, is.raw = rawdata)
  
  # número de clases válidas
  K = length(attr(M1, 'valid_classes'))
  
  
  if(K<2){# fin del problema
    stop('Insuficient number of valid classes')
  }else if(K == 2){# análisis 2x2
    # Add an extra class
    M = getM2(M1)
    
    # Perform asymptotic Delta
    # 1. Original data
    A1 = asymptoticDelta(M1)
    # 2. Original data increased by +1
    A2 = asymptoticDelta(getM4(M1))
    
    # all measures of agreement obtained by the asymptotic analysis.
    DA_all = list(DA_A1 = A1, DA_A2 = A2)
    
    # valid measures of asymptotic Delta under choosen conditions
    DA = lapply(DA_all, validMeasures, standard = standard, fixedRows = fixedRows)
    
    for(i in 1:length(DA)){
      attr(DA[[i]], 'DA_type') = attr(DA_all[[i]], 'DA_type')
    }
  }else{ # análisis K>2

    # Determinar si algún marginal es nulo
    # Elementos de la matriz 
    me = matrixElements(M1)
    margins = c(me$ri, me$ci)
    
    if(any(margins==0)){
      # Add +0.5
      M = getM3(M1)
      class(M)
    }else{
      M = M1
    }
  }
  
  
  # obtener kappa
  kappa = Kappa(M1)$k
  se_kappa = Kappa(M1)$se
  
  
  # Definir problema: getB() devuelve B0, B, h, s, tipo... matriz analizada...
  defP = getB(M, tol = tol, mxits = mxits)
  defP
  
  B = defP$B$root
  
  # La matriz analizada puede variar dentro de getB() (+0.5)
  Mx = defP$analyzedData
  
  delta_params = get_parameters(Mx, defP)
  

  all_measures = get_measures(Mx, B, delta_params)

 
  xi = chisq.test.Delta(Mx, delta_params)
  
  # Medidas válidas bajo el modelo asumido
  DN = validMeasures(all_measures, standard = standard, fixedRows = fixedRows)
  
  Estimates = DN$Estimates
  SE = DN$SE
  
  results = list(Delta = DN,
                 Kappa = list(Estimate = kappa, SE = se_kappa),
                 Data = list(Input = M1, Analyzed = Mx),
                 GOF = xi, 
                 fixedRows = fixedRows,
                 standard = standard)
  
  
  # Devolver todas las medidas para ambos tipos de muestreo, parámetros del modelo y matrices de covarianzas
  # Mostrar en el informe solo si fullReport == T
  results$all.measures = all_measures
  results$problem.parameters = defP
  results$cov = cov.Delta(Mx)
  
  
  # DELTA ASINTOTICO
  if('M2' %in% class(Mx)){
    results$asymptoticDelta = list(validMeasures = DA, allMeasures = DA_all)
  }
  

  class(results) = c('delta', 'deltaMAN')
  return(results)
  
}

#' @rdname delta
#' @export
Delta <- delta


.title = function(txt, col = 32, tex = F, sectioning = 'subsubsection'){
  if(tex == F){
    cat(paste0("\033[0;", col, "m","\n--- ", txt, " ---\n","\033[0m","\n"))
  }else{
    cat(paste0("\\",sectioning,"*{",txt,"}"))
  }
}


.section = function(txt, col = 42, tex = F){
  if(tex == F){
    cat(paste0("\033[0;", col, "m","\n------------", txt, "------------", "\033[0m","\n"))
  }else{
    cat(paste0("\\subsection*{",txt,"}"))
  }
}

#'Summary of Delta model
#'
#' \code{summary} method for class \code{"delta"}. 
#' This functions creates a report with the information contained in an object of class \code{"delta"}.
#' @param object and object of class \code{"delta"}
#' @param fullReport a logical value indicating whether to generate an exhaustive report of the results (TRUE) or not (FALSE, default).
#' @param digits an integer value indicating the significant digits to be used.
#' @param tex a logical value indicating whether to generate formatted LaTeX output (TRUE) or plain output (FALSE, default).
#' @param ... further arguments passed to summary. No one else is currently available.
#' @return No return value. Results are printed on console.
#' @export
#' @importFrom xtable xtable print.xtable
summary.delta = function(object, fullReport = FALSE, digits = 4, tex = FALSE,...){
  
  .section('SELECTED MODEL CONDITIONS', tex = tex)
  if(object$standard){
    cat('\nThere is a gold standard classification by rows.\n')
    txt2 = '   The rater by rows is a gold standard.\n'
  }else{
    cat('\nNone of the raters is a gold standard.\n')
    txt2 = '   None of the raters is a gold standard.\n'
  }
  
  if(object$fixedRows){
    cat('Row marginal frequencies have been set beforehand (type II sampling).\n')
    txt1 = '   Measures (SE) under type-II sampling.\n'
  }else{
    cat('Row marginal frequencies have not been set beforehand (type I sampling).\n')
    txt1 = '   Measures (SE) under type-I sampling.\n'
  }
  
  .section('STANDARD ANALYSIS', tex = tex)

  .title('Input data', tex = tex)
  
  messages = attr(object$Data$Input, 'messages')
  
  if(length(messages)>0){
    if(tex){
      cat(paste0(messages, '. '), sep = '\\\\ ')
    }else{
      cat(paste(messages, collapse = '\n'), '\n') 
    }
  }
  cat('Valid classes detected:', length(attr(object$Data$Input, 'valid_classes')),'\n\n')
  
  if(tex){
    print(xtable(object$Data$Input), floating = F, comment=F)
  }else{
    print(object$Data$Input)
  }
  
  if(fullReport){

    if(!all(class(object$Data$Input) == class(object$Data$Analyzed))){
      if('M2' %in% class(object$Data$Analyzed)){
        cat('\nThe problem with K = 2 classes has infinite solutions.\nTo make the problem solvable:\n\t (1) an extra ficticious third class is added with X[3,3] = 1, and\n\t (2) all the data are increaded by +0.5.')
      }
      if('M3' %in% class(object$Data$Analyzed)){
        cat('\nAt least one estimator lies at the boundary of the parametric space. \nTo turn the problem solvable, data are increased by +0.5.')
      }
      
      .title('Analyzed data', tex = tex)
      
      if(tex){
        print(xtable(object$Data$Analyzed), floating = F, comment=F)
      }else{
        print(object$Data$Analyzed)
      }
      
    }
  }

  if(fullReport){
    .title('Expected frecuencies', tex = tex)
    
    if(tex){
      print(xtable(.round(object$GOF$Expected, digits)), floating = F, comment=F)
    }else{
      print(.round(object$GOF$Expected, digits), quote = F)
    }
  }
  
  .title('Summary: Goodness of fit, Kappa and Delta ', tex = tex)
  chi_val = .round(object$GOF$statistic, digits)
  chi_df = object$GOF$df
  chi_pval = .round(object$GOF$pval, digits)
  
  Kappa_val = .round(object$Kappa$Estimate, digits)
  Kappa_SE = .round(object$Kappa$SE, (digits+1))
  Delta_val = .round(object$Delta$Estimates$Delta, digits)
  Delta_SE = .round(object$Delta$SE$Delta, (digits+1))
  
  if(tex){
    if(attr(object$GOF, 'validity')==T){
      cat(paste0(
        "\n\\begin{itemize}
          \\item $\\chi^2 =",chi_val,"$, $df =  ", chi_df,"$, $p-value = ",chi_pval,"$
          \\item $\\kappa~(SE) = ",Kappa_val,"~(",Kappa_SE,")$
          \\item $\\Delta~(SE) = ",Delta_val,"~(",Delta_SE,")$
      \\end{itemize}\n"
      ))
    }else{
      x = gsub('>20%', '$>20\\\\%$', attr(object$GOF, 'validity'))
      
      cat(paste0(
        "\n\\begin{itemize}
          \\item $\\chi^2$ ",tolower(x),"
          \\item $\\kappa~(SE) = ",Kappa_val,"~(",Kappa_SE,")$
          \\item $\\Delta~(SE) = ",Delta_val,"~(",Delta_SE,")$
      \\end{itemize}\n"
      ))
    }


  }else{
    
    if(attr(object$GOF, 'validity')==T){
      cat(' X-squared = ', chi_val, ', df = ', chi_df, ', p-value = ', chi_pval, '\n')
    }else{
      cat(' X-squared',tolower(attr(object$GOF, 'validity')),'\n')
    }
    cat(' Kappa (SE) = ', Kappa_val, ' (', Kappa_SE, ')\n', sep= '')
    cat(' Delta (SE) = ', Delta_val, ' (', Delta_SE, ')\n', sep= '')
    
  }
  
  if(fullReport){
    .title('Determination of constants B0 and B', tex = tex)
    print(object$problem.parameters, digits = digits, tex = tex)
    
    .title('Variances and covariances', tex = tex)
    print(object$cov, digits = digits, tex = tex)
    
    .title('Parameters of the Model (Delta and Pi) and ALL the Measures (SE) of Concordance', tex = tex)
    print(object$all.measures, digits = digits, tex = tex)
  }
  
  .title('Parameters of the Model (Delta and Pi) and the Measures (SE) of Concordance under the selected conditions', tex = tex)
    # cat(txt1,txt2, '\n', sep = '')
    print(object$Delta, digits = digits, tex = tex)
    
  
  
  if('M2' %in% class(object$Data$Analyzed)){
    
    .section('ASYMPTOTIC ANALYSIS', tex = tex)
    
    validDA = object$asymptoticDelta$validMeasures
    allDA = object$asymptoticDelta$allMeasures

      .title('ASYMPTOTIC SOLUTIONS. I- Based on the original data', col = 31, tex = tex)
    
      if(fullReport){
        .title('Analyzed data',tex = tex, sectioning = 'subsubsection')
        
        if(tex){
          print(xtable(object$Data$Input), floating = F)
        }else{
          print(object$Data$Input)
        }
        
        .title('Parameters of the Model (Delta and Pi) and ALL the Measures (SE) of Concordance', tex = tex, sectioning = 'subsubsection')
        print(allDA$DA_A1, digits = digits, tex = tex)
      }
      
      .title('Parameters of the Model (Delta and Pi) and the Measures (SE) of Concordance under the selected conditions', tex = tex, sectioning = 'subsubsection')
      # cat(txt1,txt2, '\n', sep = '')
      print(validDA$DA_A1, digits = digits, tex = tex)
      

      
      
      .title('ASYMPTOTIC SOLUTIONS. II- Based on the original data incremented +1', col = 31,tex = tex)

      if(fullReport){
        .title('Analyzed data', tex = tex, sectioning = 'subsubsection')
        
        if(tex){
          print(xtable(object$Data$Input+1), floating = F)
        }else{
          print(object$Data$Input+1)
        }
        .title('Parameters of the Model (Delta and Pi) and ALL the Measures (SE) of Concordance', tex = tex, sectioning = 'subsubsection')
        print(allDA$DA_A2, digits = digits, tex = tex)
      }
      
      .title('Parameters of the Model (Delta and Pi) and the Measures (SE) of Concordance under the selected conditions', tex = tex, sectioning = 'subsubsection')
      # cat(txt1,txt2, '\n', sep = '')
      print(validDA$DA_A2, digits = digits, tex = tex)
      
    
  }
  cat('\n')
}


