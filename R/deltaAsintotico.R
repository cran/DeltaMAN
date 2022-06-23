
# Delta Asintótico------
# Solo para matrices 2x2

# Medidas para Delta asintótico ----
#' Estimate measures of agreement based on the asymptotic solution of Delta.
#' 
#' asymptoticDelta() is used for 2x2 matrices. 
#' It applies an aymptotic solution to compute the different delta measures (overall, agreement, conformity, predictivity and consistency).
#' It returns their estimates and their standard errors (SE) under type-I and type-II samplings.
#' 
#' parameters: 
#' @param M a matrix of class M1 or M4
#' @return An object of class \code{"measuresDelta"}
#' @details A matrix of class M1 is an squared matrix (2x2) containing the original frequencies. 
#' A matrix of class M4 is an squared matrix (2x2) containing the original frequencies increased by +1. 
#' @noRd

asymptoticDelta = function(M){
  
  .DA_params_SE = function(x){
    n = sum(x)
    xii = diag(x)
    
    r_i = marginSums(x, 1)
    c_i = marginSums(x, 2)
    offdiag = c(x[2,1], x[1,2])
    
    # PARÁMETROS DEL MODELO
    pi_i = sqrt(offdiag)/sum(sqrt(offdiag))
    
    
    A_i = (xii-sqrt(prod(offdiag)))/n
    
    delta_i = n*A_i/r_i
    F_i = delta_i
    
    P_i = n*A_i/c_i 
    P_i
    
    S_i = (2*n*A_i)/(r_i+c_i)
    
    Delta = (sum(xii) - 2*sqrt(prod(offdiag)))/n
    
    
    # VARIANZAS (solo en el caso DA1)
    if(case == 'DA1'){
      Vii = 1/r_i^2*(xii*(1-delta_i) + 1/4*(sum(offdiag) - prod(n,offdiag)/prod(r_i)))
      
      
      ## Conformidad: F_i = delta_i
      # Muestreo Tipo I == Muestreo Tipo II
      
      V2_F = Vii
      V1_F = 1/r_i^2*(xii*(1-delta_i) + 1/4*sum(offdiag))
      
      ## Predictividad: P_i (solo tipo I, no tiene sentido para tipo II)
      # Muestreo Tipo I:
      V1_P = 1/c_i^2*(xii*(1-P_i) + 1/4*sum(offdiag))
      
      
      ## Acuerdo global: A_i
      # Muestreo Tipo I:
      V1_A = 1/n^2*(xii + 1/4*sum(offdiag) - n*A_i^2)
      
      
      # Muestreo Tipo II:
      V2_A = (r_i/n)^2*Vii
      
      
      ## Consistencia: S_i (solo tipo I, no tiene sentido para tipo II)
      # Muestreo Tipo I:
      V1_S = n*(1-Delta)/(r_i+c_i)^2 * (2 - n*(1-Delta)*sum(offdiag)/(r_i+c_i)^2)
      
      ## Delta
      # Muestreo Tipo I:
      V1_delta = 1/n*(1-Delta)*(1+Delta)
      # Muestreo Tipo II:
      V2_delta = (1-Delta)/n * sum(xii/r_i)
      
      
      sampling1 = list(Delta = V1_delta, Agreement = V1_A, Conformity = V1_F, Predictivity = V1_P, Consistency = V1_S)
      sampling2 = list(Delta = V2_delta, Agreement = V2_A, Conformity = V2_F)
      
      SE_1 = lapply(sampling1, sqrt)
      SE_2 = lapply(sampling2, sqrt)
      
    }else{
      SE_1 = list(Delta = NULL, Agreement = NULL, Conformity = NULL, Predictivity = NULL, Consistency = NULL)
      SE_2 = list(Delta = NULL, Agreement = NULL, Conformity = NULL)
    }
    
    
    delta_params = list(Delta = Delta, proportions = pi_i, partial_delta = delta_i)
    
    attr(delta_params, 'SE_Delta') = list(SE_samplingI = SE_1$Delta, SE_samplingII = SE_2$Delta)
    
    class(delta_params) = 'psetDelta'
    
    results = list(Global = delta_params,
                   Agreement = list(Measure = A_i, SE_samplingI = SE_1$Agreement, SE_samplingII = SE_2$Agreement), 
                   Conformity = list(Measure = F_i, SE_samplingI = SE_1$Conformity, SE_samplingII = SE_2$Conformity), 
                   Predictivity = list(Measure = P_i, SE_samplingI = SE_1$Predictivity, SE_samplingII = NULL), 
                   Consistency = list(Measure = S_i, SE_samplingI = SE_1$Consistency, SE_samplingII = NULL))
    
    # class(results) = 'measuresDelta'
    class(results) = c('measuresDelta', 'deltaMAN')
    return(results)
  }
  
  
  
  r_i = marginSums(M, 1)
  c_i = marginSums(M, 2)
  offdiag = c(M[2,1], M[1,2])
  
  # CASO DA0 o DA1
  if(all(offdiag == 0) | any(r_i == 0) | any(c_i == 0)){
    # CASO DA0
    case = 'DA0'
  }else{
    # CASO DA1
    case = 'DA1'
  }
  
  # A1 = .DA_params_SE(M)
  # 
  # if(case != 'DA0'){
  #   # Obtener matriz M4 (+1)
  #   M4 = getM4(M)
  #   A2 = .DA_params_SE(M4)
  #   results = list(DA_A1 = A1, DA_A2 = A2)
  #   attr(results, 'DA_type') = 'DA1'
  # 
  # }else{
  #   results = list(DA_A1 = A1, DA_A2 = NULL)
  #   attr(results, 'DA_type') = 'DA0'
  # }

  
  results = .DA_params_SE(M)
  
  if(case != 'DA0'){
    attr(results, 'DA_type') = 'DA1'
    
  }else{
    attr(results, 'DA_type') = 'DA0'
  }
  
 return(results)

}
  
  
  
  
  
  
  
  
  
  
  
  
  