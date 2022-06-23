# Auxiliary functions to compute delta
# These are all non-exported functions (hidden)


#' Get elements of matrix
#' @noRd
matrixElements = function(M){
  # k = attr(M, 'valid_classes')
  # ri = marginSums(M, 1)[k]
  # ci = marginSums(M, 2)[k]
  # n = sum(r_i)
  
  ri = marginSums(M,1)
  ci = marginSums(M,2)
  xii = diag(M)
  S = sum(diag(M))
  n = sum(M)
  
  return(list(ri = ri, ci = ci, xii = xii, S = S, n = n))
}




# Obtención de la constante desconocida B: getB() ----
#' Estimate the value B
#' 
#' getB() computes the unknown constant B, needed to estimate the parameters of the model.
#' @param M an object of class \code{"M1"}, \code{"M2"} or \code{"M3"}.
#' @param tol the desired tolerance applied to find the root of the unknown constant B.
#' @param mxits the maximum numer of iterations applied to find the root of the unknown constant B.
#' @return An object of class \code{"defP"}. 
#' This object is a list of 3: the root of B, the analyzed matrix and a list containing information related to the determination of B.
#' @noRd

getB = function(M, tol = 1e-7, mxits = 100){
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Definir funciones internas
  # Las funciones internas llevan un '.' delante del nombre
  .getB0 = function(M){
    
    b0 = (sqrt(ci-xii) + sqrt(ri-xii))^2
    
    B0 = max(b0)
    # B0_i = which(b0 == max(b0)) # todos los índices que proporcionan a B0
    B0_i = which.max(b0) # primer índice que proporciona a B0
    
    attr(B0, 'h') = B0_i
    return(B0)
  }
  
  
  # El único argumento es B, el resto de objetos se toman del scope de la función principal getB()
  .y = function(B){
    
    # prevenir la aparición de raíces negativas por falta de precisión
    # si arg < 0 y abs(arg)< tol, tomar arg = 0
    arg = (B+(ci-ri))^2 - 4*B*(ci-xii)
    
    # tol = 1e-5
    # if(all(arg>= 0 | (arg<0 & abs(arg)<tol))){
    #   arg = pmax(arg, 0)
    # }else{
    #   stop('Negative roots')
    # }
    
    arg = pmax(arg, 0)
    B = (K-2)*B + sum(s*sqrt(arg))
  
    return(B)
  }
  

  
  # derivada y'(B)
  .dydb = function(B0){
    # incremento
    inc = 0.01
    
    d = (.y(B0+inc)-.y(B0))/inc
    return(d)
  }
  
  

  # resolver ecuacion no lineal con newton raphson
  .newton_raphson = function(B0){
    # tol = 1e-7
    # mxits = 10000
    
    # valor inicial
    b0 = B0
    attributes(b0) = NULL
    # vector de resultados
    v = c()
    
    # iterar desde 1 hasta mxits
    for(i in 1:mxits){
      
      # Obtener aproximación de la derivada y'(B)
      dB = .dydb(b0)
      # Evaluar y() en b0
      yB = .y(b0)
      
      # Buscar raíz 
      b1 = b0-yB/dB
      # Guardar resultado de iteración i
      v[i] = b1
      
      # Comprobar si se cumple el criterio de parada
      dif = abs(b1-b0)
      if(dif< tol){
        break
      }
      # Si no se cumple el criterio de parada, ir a la iteración i+1, con b0 actualizado
      b0 = b1
    }
    # Si llegamos al número máx de iteraciones, devolver último resultado
    return(list(root = b1, iterations = length(v), tol = dif))
  }
  

  
  # Obtener tipo de problema para buscar la solución en y+(B) = 0 o y-(B) = 0
  .getType = function(B0){

    Y = n-S
    Z = ci[h] + ri[h] - 2*M[h,h]
    sgn = -1
    yB0 = .y(B0)
    if(yB0>0){
      tipo = 1
      DN = 'DN2'
      
      # solución y-(B0) = 0
      # s[h] = -1
      # B = newton_raphson(B0)
    }else if(yB0<0){
      if(Y>Z){
        tipo = 2
        DN = 'DN2'
        
        # solución y+(B0) = 0
        # s[h] = +1
        # B = newton_raphson(B0)
      }else{
        tipo = 3
        DN = 'DN1'
        # B = Inf
      }
    }else if(yB0==0){
      if(Y>Z){
        tipo = 4
        DN = 'DN0'
        # B = B0
      }else if(Y==Z & Y>0){
        tipo = 5
        DN = 'DN3'
        # B = NaN
      }else if(Y == Z & Y == 0){
        tipo = 6
        DN = 'DN0'
        # B = 0
      }
      
    }
    # return(list(tipo = tipo, DN = DN, B = B))
    return(tipo)
  }
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  # Elementos de M
  me = matrixElements(M)
  ri <- me$ri
  ci <- me$ci
  xii <- me$xii
  n <- me$n
  S <- me$S
  
  
  K <- length(attr(M, 'valid_classes'))
  if('M2'%in% class(M)){
    K = 3
  }
  
  # Obtener B0 (atributo 'h' es la categoría que proporcional el valor B0)
  B0 <- .getB0(M)
  h <- attr(B0, 'h')
  
  # Determinar tipo de problema
  s <- rep(-1, K)
  type = .getType(B0)
  
  if(type > 2){
    M = getM3(M)
    message('Data are increased by +0.5 to turn the problem solvable')
    result = getB(M)
    return(result)
    
  }
  
  # Obtener B
  yB0 <- .y(B0)
  if(yB0<0){
    s[h] <- 1
  }
  B = .newton_raphson(B0)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Información para el report
  s_save = s # vector s(i)
  
  attributes(B0) = NULL
  s[h] = -1
  yM = .y(B0)
  dyM = .dydb(B0)
  
  s[h] = 1
  yP = .y(B0)
  dyP = .dydb(B0)
  
  Y = n-S
  Z = ci[h] + ri[h] - 2*M[h,h]
  
  info = list(B0 = B0, h = h, s = s_save, PType = type, 
              yM = yM, dyM = dyM, yP = yP, dyP = dyP,
              'ns' = Y, 'cr2x' = Z)

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  result = list(B = B, analyzedData = M,info = info)
  # class(result) = 'defP'
  class(result) = c('defP', 'deltaMAN')
  return(result)
}



# Obtener parámetros -----
#' Obtain the parameters of the Delta model
#' 
#' Compute the global Delta, partial delta for each category and the distributions of responses made at random.
#' 
#' @param M an object of class \code{"M1"}, \code{"M2"} or \code{"M3"}.
#' @param defP An object of class \code{"defP"}. 
#' @return An object of class \code{"psetDelta"}. This is a list of 3 elements:
#' the global Delta coefficient, the partial delta coefficient for each category and the distribution of responses made at random by observer C (proportions)
#' @noRd
get_parameters = function(M, defP){
  .get_pset = function(){

    # proportions
    arg = pmax((B + (ci-ri))^2 - 4*B*(ci-xii),0)
    pi_i = (B + (ci-ri) + s*sqrt(arg))/(2*B)
    
    # Partial delta
    delta_i = (xii-ri*pi_i)/(ri*(1-pi_i))
    
    
    # Overall delta
    if('M2' %in% class(M)){
      k = attr(M, 'valid_classes')
      delta = sum(ri[k]*delta_i[k])/sum(ri[k])
    }else{
      delta = 1-(B/n)
    }
    
    return(list(Delta = delta, proportions = pi_i, partial_delta = delta_i))
  }
  
  me = matrixElements(M)
  ri <- me$ri
  ci <- me$ci
  xii <- me$xii
  n <- me$n
  S <- me$S
  
  
  B = defP$B$root
  h = defP$info$h
  s = defP$info$s # s(i), ya viene definido si s(h) es +1 o -1

  type = defP$info$PType
  
  pset = .get_pset()
  
  if(type == 3){
    partial_delta = xii/ri
    partial_delta[h]  = -Inf
    pset$partial_delta = partial_delta
    
    pset$Delta = -Inf
    
  }else if(type == 6){
    pset$Delta = 1
    k = length(pset$partial_delta)
    pset$partial_delta[1:k] =1
  }
  
  # class(pset) = 'psetDelta'
  class(pset) = c('psetDelta', 'deltaMAN')
  return(pset)
}




# Bondad de ajuste -------
#' chisq.test.Delta computes the goodness of fit of the Delta model for a given contingency table.
#' 
#' @param M a matrix of type M1, M2 or M3
#' @param psetDelta the set of parameters obtained as a result of function result of get_parameters(): 
#' i.e., overall delta, partial delta and distribution of responses at random (proportions). 
#' If not supplied, it is computed by running functions getB() and get_parameters().
#' 
#' @details 
#' If the parameter psetDelta is given, M must be the analyzed matrix returned by getB().
#' If it is not given, the parameters are computed and the matrix M is substituted by the analyzed matrix returned by getB().
#' @importFrom stats pchisq
#' @noRd
chisq.test.Delta <- function(M, psetDelta){
  

  if(missing(psetDelta)){
    suppressMessages({
      defP = getB(M)
    })
    M = defP$analyzedData
    psetDelta = get_parameters(M, defP)
  }
  
  
  
  K = dim(M)[1]
  
  k = attr(M, 'valid_classes')
  
  
  ri = marginSums(M, 1)
  
  delta_i = psetDelta$partial_delta
  pi = psetDelta$proportions
  
  # Expected frequencies for Delta
  E_ij = outer(ri*(1-delta_i), pi)
  diag(E_ij) = diag(M)
  
  E_ij = E_ij[k,k]
  M = M[k,k]

  # Get chi-square statistic
  xi = sum((M-E_ij)^2/E_ij)
  
  # Degrees of freedom
  df = (K-1)*(K-2)-1
  
  # p-value
  pvalor = pchisq(xi, df, lower.tail = F)
  
  result = list(statistic = xi, df = df, pval = pvalor, Observed = M, Expected = E_ij)
  
  # Check test validity
  # Cuando más de un 20% de las frecuencias esperadas por un modelo es menor a 5, o bien alguna es menor a 1, el test no es válido
  lowfreq = which(E_ij<5);lowfreq
  E_ij[lowfreq]
  propE_ij = prop.table(E_ij)*100;propE_ij
  if(sum(propE_ij[lowfreq])>20 | any(E_ij<1)){
    result$statistic = NA
    result$df = NA
    result$pval = NA
    if(any(E_ij<1)){
      validity = paste('Test not aplicable, min(Eij) = ', round(min(E_ij),4))
    }else{
      validity = paste('Test not aplicable, >20% of expected frequencies are lower than 5')
    }
  }else{
    validity = T
  }
  attr(result, 'validity') = validity
  class(result) = c('chisqTestDelta', 'deltaMAN')
  return(result)
}



# Covarianzas -------

#' cov.Delta computes the covariance matrices of (1) partial deltas, (2) proportions and (3) partial deltas and proportions,
#' where proportions refers to the distribution of responses made at random.
#' 
#' parameters: 
#' @param M a matrix of type M1, M2 or M3
#' @param B the value of the unknown constant B (returned by getB()). If not supplied, it is computed.
#' @param psetDelta the set of parameters obtained as a result of function result of get_parameters(): 
#' i.e., overall delta, partial delta and distribution of responses at random (proportions). If not supplied, it is computed by running function get_parameters().
#' 
#' @details 
#' If the parameters B and psetDelta are given, M must be the analyzed matrix returned by getB().
#' @noRd

cov.Delta <- function(M, B, psetDelta){
  
  if(missing(psetDelta) | missing(B)){
    suppressMessages({
      defP = getB(M)
    })
    B = defP$B$root
    M = defP$analyzedData
    psetDelta = get_parameters(M, defP)
  }

  
  r_i = marginSums(M, 1)
  xii = diag(M)
  
  
  delta_i = psetDelta$partial_delta
  pi_i = psetDelta$proportions
  
  v_i = (1-delta_i)/(1-pi_i)
  
  E_i = pi_i/(B-r_i*v_i)
  E = sum(E_i)
  
  cov_D_iD_j= -outer(v_i*E_i, v_i*E_i)/E+diag(v_i*(xii/r_i^2+(v_i*E_i)))
  
  cov_D_iPi_j = outer(v_i*E_i, E_i)/E - diag(v_i*E_i)
  
  cov_Pi_iPi_j = -outer(E_i, E_i)/E+diag(E_i) # no sale igual que en https://wpd.ugr.es/~bioest/delta.php
  
  result = list(Cov_Delta = cov_D_iD_j, cov_Mix = cov_D_iPi_j, cov_Pi = cov_Pi_iPi_j)
  
  # class(result) = 'cov.Delta'
  class(result) = c('cov.Delta', 'deltaMAN')
  return(result)
  
}




# Varianzas -----
# Depende del tipo de muestreo (I: nada prefijado y II: marginal de filas prefijado)

#' var.Delta computes the variance of the different delta measures (overall, agreement, conformity, predictivity and consistency)
#' under type-I and type-II samplings.
#' 
#' @param M a matrix of type M1, M2 or M3
#' @param B the value of the unknown constant B (returned by getB()). If not supplied, it is computed.
#' @param psetDelta the set of parameters obtained as a result of function result of get_parameters(): 
#' i.e., overall delta, partial delta and distribution of responses at random (proportions). If not supplied, it is computed by running function get_parameters().
#' @details 
#' If the parameters B and psetDelta are given, M must be the analyzed matrix returned by getB().
#' @noRd

var.Delta = function(M, B, psetDelta){

  
  if(missing(psetDelta) | missing(B)){
    suppressMessages({
      defP = getB(M)
    })
    B = defP$B$root
    M = defP$analyzedData
    psetDelta = get_parameters(M, defP)
  }
  
  
  
  # M = defP$analyzedData
  K = length(attr(M, 'valid_classes'))
  k = attr(M, 'valid_classes')
  r_i = marginSums(M, 1)[k]
  c_i = marginSums(M, 2)[k]
  x_ii = diag(M)[k]
  
  
  # delta_params = get_parameters(M, defP)
  delta_i = psetDelta$partial_delta[k]
  pi_i = psetDelta$proportions[k]
  delta = psetDelta$Delta
  
  v_i = (1-delta_i)/(1-pi_i)
  
  E_i = pi_i/(B-r_i*v_i)
  E = sum(E_i)
  n = sum(r_i)
  
  Vij = cov.Delta(M, B, psetDelta)$Cov_Delta
  Vii = diag(Vij)[k]
  
  ## Delta
  if('M2'%in%class(M)){
    Vij = Vij[k,k]
    # Muestreo Tipo I:
    V1_delta = (sum(outer(r_i,r_i)*Vij) + sum(r_i*delta_i^2) - n*delta^2)/n^2
    # Muestreo Tipo II:
    V2_delta = sum(outer(r_i,r_i)*Vij)/n^2
  }else{
    # Muestreo Tipo I:
    V1_delta = 1/n^2*(n-1/E-n*delta^2)
    # SE
    sqrt(V1_delta)
    
    # Muestreo Tipo II:
    V2_delta = 1/n^2*(n-1/E-sum(r_i*delta_i^2))
    # SE
    sqrt(V2_delta)
  }

  
  ## Conformidad: F_i = delta_i
  # Muestreo Tipo I == Muestreo Tipo II
  V_F = Vii
  sqrt(V_F)
  
  
  ## Predictividad: P_i (solo tipo I, no tiene sentido para tipo II)
  # Muestreo Tipo I:
  V1_P = (r_i/c_i)^2*(Vii + (c_i-r_i)/(c_i*r_i)*delta_i^2)
  # SE
  sqrt(V1_P)
  
  
  ## Acuerdo global: A_i
  # Muestreo Tipo I:
  V1_A = (r_i/n)^2*(Vii+(n-r_i)/(n*r_i)*delta_i^2)
  sqrt(V1_A)
  
  
  # Muestreo Tipo II:
  V2_A = (r_i/n)^2*Vii
  sqrt(V2_A) # no sale igual que en https://wpd.ugr.es/~bioest/delta.php
  
  
  ## Consistencia: S_i (solo tipo I, no tiene sentido para tipo II)
  # Muestreo Tipo I:
  V1_S = (2*r_i/(r_i+c_i))^2*(Vii + delta_i^2/(r_i+c_i)*(c_i/r_i - 2 + (2*x_ii)/(r_i+c_i)))
  # SE
  sqrt(V1_S)
  
  sampling1 = list(Delta = V1_delta, Agreement = V1_A, Conformity = V_F, Predictivity = V1_P, Consistency = V1_S)
  sampling2 = list(Delta = V2_delta, Agreement = V2_A, Conformity = V_F)
  
  result = list(Sampling_I = sampling1, Sampling_II = sampling2)
  
  # class(result) = 'var.Delta'
  class(result) = c('var.Delta', 'deltaMAN')
  return(result)
}

# Medidas para Delta Normal ----
#' Estimate measures of agreement based on Delta
#' 
#' get_measures() computes the different delta measures (overall, agreement, conformity, predictivity and consistency).
#' It returns their estimates and their standard errors (SE) under type-I and type-II samplings.
#' 
#' parameters: 
#' @param M a matrix of type M1, M2 or M3
#' @param B the value of the unknown constant B (returned by getB()). If not supplied, it is computed.
#' @param psetDelta the set of parameters obtained as a result of function result of get_parameters(): 
#' i.e., overall delta, partial delta and distribution of responses at random (proportions). If not supplied, it is computed by running function get_parameters().
#' @details 
#' If the parameters B and psetDelta are given, M must be the analyzed matrix returned by getB().
#' @return An object of class \code{"measuresDelta"}
#' @noRd

get_measures = function(M, B, psetDelta){
  
  if(missing(psetDelta) | missing(B)){
    suppressMessages({
      defP = getB(M)
    })
    B = defP$B$root
    M = defP$analyzedData
    psetDelta = get_parameters(M, defP)
  }
  
  k = attr(M, 'valid_classes')
  r_i = marginSums(M, 1)[k]
  c_i = marginSums(M, 2)[k]

  n = sum(r_i)
  
  D = psetDelta$Delta
  d_i = psetDelta$partial_delta[k]
  p_i = psetDelta$proportions[k]
  
  F_i = d_i
  P_i = r_i*d_i/c_i
  S_i = 2*r_i*d_i/(r_i+c_i)
  A_i = r_i*d_i/n
  
  # Varianzas MUESTREO TIPO I
  variance = var.Delta(M, B, psetDelta)
  var1 = variance$Sampling_I
  SE_1 = lapply(var1, function(x)sqrt(unname(x)))
  
  
  # Varianzas MUESTREO TIPO II
  var2 = variance$Sampling_II
  
  SE_2 = lapply(var2, function(x)sqrt(unname(x)))
  
  psetDelta$proportions = p_i
  psetDelta$partial_delta = d_i
  # delta_params$Delta = list(Delta = D, SE_samplingI = SE_1$Delta, SE_samplingII = SE_2$Delta)
  
  attr(psetDelta, 'SE_Delta') = list(SE_samplingI = SE_1$Delta, SE_samplingII = SE_2$Delta)
  # results = list(Agreement = A_i, Conformity = F_i, Predictivity = P_i, Consistency = S_i)
  results = list(Global = psetDelta, Agreement = list(Measure = A_i, SE_samplingI = SE_1$Agreement, SE_samplingII = SE_2$Agreement), 
       Conformity = list(Measure = F_i, SE_samplingI = SE_1$Conformity, SE_samplingII = SE_2$Conformity), 
       Predictivity = list(Measure = P_i, SE_samplingI = SE_1$Predictivity, SE_samplingII = NULL), 
       Consistency = list(Measure = S_i, SE_samplingI = SE_1$Consistency, SE_samplingII = NULL))
  
  # class(results) = 'measuresDelta'
  class(results) = c('measuresDelta', 'deltaMAN')
  # attr(results, 'delta') = psetDelta

  return(results)
}





# Medidas válidas para el modelo asumido ----
#' Estimate valid measures of agreement based on Delta under the model assumed.
#' 
#' get_measures() computes the valid delta measures (overall, agreement, conformity, predictivity and consistency)
#' for the model assumed (sampling type I or II, and presence/absence of goldstandard).
#' It returns the valid estimates and their standard errors (SE) under the model assumed.
#' 
#' parameters: 
#' @param measuresDelta an object of class \code{"measuresDelta"}. 
#' @param standard a logical value indicating whether the observer on the rows of the contingency table (or first column of raw data) is a goldstandard.
#' @param fixedRows a logical value indicating whether the row marginals are fixed in advance (sampling type II) or not (sampling type I).
#' @return An object of class \code{"validMeasuresDelta"}
#' @noRd

validMeasures = function(measuresDelta, standard, fixedRows){
  global = measuresDelta$Global
  SE_delta = attr(global, 'SE_Delta')
  
  if(standard == F){ # SIN GOLDSTANDARD
    
    if(fixedRows == F){# Muestreo TIPO I 
      # AGREEMENT Y CONSISTENCY
      measures = measuresDelta[c('Agreement', 'Consistency')]
      m_Estimate = lapply(measures, '[[', 'Measure')
      m_SE = lapply(measures, '[[', 'SE_samplingI')
      
      SE_D = SE_delta$SE_samplingI
      
    }else{# Muestreo TIPO II 
      # AGREEMENT 
      measures = measuresDelta[c('Agreement')]
      m_Estimate = lapply(measures, '[[', 'Measure')
      m_SE = lapply(measures, '[[', 'SE_samplingII')
      
      SE_D = SE_delta$SE_samplingII
    }
    
  }else{ # CON GOLDSTANDARD (en filas)
    if(fixedRows == F){# Muestreo TIPO I 
      # AGREEMENT, CONFORMITY Y PREDICTIVITY
      
      measures = measuresDelta[c('Agreement', 'Conformity', 'Predictivity')]
      m_Estimate = lapply(measures, '[[', 'Measure')
      m_SE = lapply(measures, '[[', 'SE_samplingI')
      
      SE_D = SE_delta$SE_samplingI
      
    }else{# Muestreo TIPO II 
      # AGREEMENT Y CONFORMITY
      measures = measuresDelta[c('Agreement', 'Conformity')]
      m_Estimate = lapply(measures, '[[', 'Measure')
      m_SE = lapply(measures, '[[', 'SE_samplingII')
      
      SE_D = SE_delta$SE_samplingII
    }
  }
  
  
  Estimates = unlist(list(global, m_Estimate), recursive = F)
  SE = unlist(list(Delta = list(SE_D), m_SE), recursive = F)
  
  result = list(Estimates = Estimates, SE = SE)
  
  # class(result) = 'validMeasuresDelta'
  class(result) = c('validMeasuresDelta', 'deltaMAN')
  return(result)
}





