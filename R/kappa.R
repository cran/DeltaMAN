#' Compute the Cohen's kappa coefficient
#' 
#' Kappa() computes de Cohen's kappa coefficient for nominal or ordinal data. 
#' If data is ordinal, weigthed kappa can be applied to allow disagreements to be weighted differently.
#' @param m a squared matrix of frequencies between two observers.
#' @param r an integer (0, 1, or 2) to create a matrix of weigths. See details.
#' @param alternative a string specifying the alternative hypothesis to construct the confidence interval: 
#' either "two.sided" (default), "greater" or "less".
#' @param conf.level confidence level of the interval.
#' @param partial a logical value indicating whether to evaluate the degree of agreement of each category by collapsing the contingency table.
#' 
#' @return A list of 3 elements containing the kappa statistic, the standard error and the confidence interval. 
#' If \code{"partial = TRUE"}, a data.frame containing 3 columns (the class, the unweighted partial kappa coefficient 
#' for each class and the standards error of each estimate) is added to the list.
#' @details 
#' The weighted kappa can be computed when data are ordinal and the argument \code{r} is eigher 1 or 2:
#' * if r = 0, unweighted kappa is computed (used for nominal variables)
#' * if r = 1, weighted kappa with linear formula is applied
#' * if r = 2, weighted kappa with quadratic formula is applied
#' 
#' @export
#' @importFrom stats qnorm
#' 
#' @examples 
#' # Create a 3x3 matrix
#' m = matrix(c(15, 5, 0, 4, 21, 1, 3, 4, 25), ncol = 3)
#' # Compute the Kapa coefficient for nominal data
#' Kappa(m, r = 0, partial = TRUE)
#' # Compute the Kapa coefficient for ordinal data, using linear formula
#' Kappa(m, r = 1)
Kappa = function(m, r = 0, alternative = c("two.sided", "less", "greater"), 
                 conf.level = 0.95, partial = FALSE){
  if(r<0 | r>2){
    stop('Argument "r" must be either 0, 1 or 2.')
  }
  if(missing(alternative)){
    alternative = 'two.sided'
  }
  
  value = .kappa(m, r = r)
  
  se = .SE_kappa(m, r = r)
  
  if(alternative == "less") {
    conf.int <- c(-Inf, value + qnorm(conf.level)*se )
  }
  else if (alternative == "greater") {
    
    conf.int <- c(value - qnorm(conf.level)*se, Inf)
  }
  else {
    alpha = (1-conf.level)/2
    z = abs(qnorm(alpha))
    conf.int = value + c(- z*se, z*se)
  }
  
  result = list(k = value, se = se, conf.int = conf.int)
  if(partial == TRUE){
    pkappa = .partial_kappa(m)
    result$partial.kappa = pkappa
  }
  return(result)

}

#' Compute the Cohen's kappa coefficient for partial agreement
#' 
#' partial_kappa() evaluates the degree of agreement by category by collapsing the 
#' contingency table for each one.
#' 
#' @param m a squared matrix of frequencies between two observers.
#' @noRd
#' @return A data.frame containing 3 columns: the class, the unweighted partial kappa coefficient for each class and the standards error of each estimate.
#' 
#' @examples 
#' # Create a 3x3 matrix
#' m = matrix(c(15, 5, 0, 4, 21, 1, 3, 4, 25), ncol = 3)
#' # Compute the partial kapa coefficient for each category
#' partial_kappa(m)
.partial_kappa = function(m){
  K = ncol(m)
  cm = list()
  
  partial_k = list()
  for(i in 1:K){
    cm[[i]] = .colapseMatrix(m, i)
    partial_k[[i]] = Kappa(cm[[i]])
  }
  v_names = dimnames(m)[[1]]
  if(is.null(v_names)){
    v_names = 1:K
  }
  names(cm) = v_names
  
  k = sapply(partial_k, '[[',1)
  se = sapply(partial_k, '[[',2)
  
  result = data.frame(Class = v_names, kappa = k, SE = se)
  return(result)
}

#' Coeficiente kappa
#' @noRd
.kappa = function(m, r = 0){

  K = ncol(m)
  w = .weightMatrix(K, r = r)
  
  # Frecuencias observadas
  Io = sum(w*m)/sum(m)
  
  # Frecuencias esperadas debidas al azar
  x = .expected(m)
  Ie = sum(x*w)
  
  k = (Io - Ie)/(1-Ie)
  return(k)
}

# Matriz de frecuencias esperadas
.expected = function(m){
  margin.rows = prop.table(marginSums(m,1))
  margin.cols = prop.table(marginSums(m,2))
  
  x = outer(margin.rows, margin.cols)
  
  return(x)
}

#' Matriz de pesos en coeficiente ponderado
#' @noRd
.weightMatrix = function(K, r = 0){

  w = diag(K)
  
  for(i in 1:K){
    for(j in 1:K){
      if(i == j){
        next
      }
      w[i,j] = 1 - (abs(i-j)/(K-1))^r
    }
  }
  return(w)
}


#' Error estandar de estimacion
#' @noRd
.SE_kappa = function(m, r = 0){

  pm = prop.table(m)
  
  k = .kappa(m, r)
  
  # marginal pi.
  pi. = rowSums(pm)
  
  # marginal p.j
  p.j = colSums(pm)
  
  # pii (digonal)
  pii = diag(pm)
  
  K = ncol(m)
  
  w = .weightMatrix(K, r = r)
  
  # marginal wi.
  wi. = as.vector(w %*% p.j)
  # marginal w.j
  w.j = as.vector(w %*%pi.)
  
  
  A = sum(pm*(w-outer(wi.,w.j, FUN = '+')*(1-k))^2)
  
  # Expected frequencies (weigthed)
  xi = .expected(m)*w
  Ie = sum(xi)
  B = (k - Ie*(1-k))^2
  
  se = sqrt((A-B)/(sum(m)*(1-Ie)^2))
  
  return(se)
}


#' Collapse classification matrix to compute partial kappa
#' @noRd
.colapseMatrix = function(m, class){
  i = class
  mc = matrix(c(m[i,i], sum(m[i,-i]), 
           sum(m[-i,i]), sum(m[-i,-i])), ncol = 2, byrow = T)
  dimnames(mc) = list(c('Yes', 'No'), c('Yes', 'No'))
  return(mc)
}
