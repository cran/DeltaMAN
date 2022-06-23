#' Get M2
#' (M1+una clase extra con x33=c3=r3=1)+0.5: solo cuando M1=2x2.
#' @noRd

getM2 = function(M1){
  # M2 = matrix(0,3,3)
  # M2[1:2,1:2] = M1
  # M2
  M2 = rbind(cbind(M1,0),0)
  M2[3,3] = 1
  
  # Atributos
  attribute_list = attributes(M1)
  attribute_list$dim = c(3,3)
  
  dimnams = c(attribute_list$dimnames[[1]],'extCls')
  attribute_list$dimnames[[1]] = dimnams
  attribute_list$dimnames[[2]] = dimnams
  attributes(M2) = attribute_list
  M2 = M2+0.5
  class(M2) <- c('table', 'M2')
  return(M2)
}

#' M3 = M1+0.5.
#' @noRd
getM3 = function(M1){
  M3 = M1+0.5
  class(M3) <- c('table', 'M3')
  return(M3)
}


#' M4 = M1+1: solo cuando M1=2x2.
#' @noRd
getM4 = function(M1){
  M4 = M1+1
  class(M4) <- c('table', 'M4')
  return(M4)
}
