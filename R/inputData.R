

#' Get contingency table
#' 
#' @param x raw data or frequency table
#' @param is.raw the user can give information about the data type. 
#'         If not defined, the function will try to guess the type
#'         
#' Returns object of class 'M1'
#' @noRd

getM1 <- function(x, is.raw = NULL){
  if(!is.null(is.raw)){
    if(is.raw == T){ # datos en bruto
      if(is.matrix(x)){ # si el objeto es matrix, cambiar a data.frame
        x = data.frame(x)
      }
      ct = table.ct(x)
    }else{ # datos en tabla de frecuencias
      if(is.matrix(x)){# si el objeto es matrix, cambiar a table
        x = as.table(x)
      }
      ct = checkTable(x)
      
    }
    return(ct)
  }
  
  # GUESS TYPE OF X
  # comprobar dimensiones del objeto de entrada
  # si el número de filas y el de columnas es distinto, tenemos datos en bruto

  squared = length(unique(dim(x))) == 1
  messages = NULL
  # tenemos datos en bruto si:
  #   * la matriz no es cuadrada (squared == F), 0
  #   * (caso extremo) es cuadrada 2x2 y los elementos no son numéricos (categorías cualitativas no codificadas a números)
  # Además, x no es de tipo 'table' (clase correspondiente a tabla de frecuencias)
  # si tenemos datos en bruto, comprobar numero de columnas:
  if((squared == F | !is.numeric(unlist(x))) & !is.table(x)){
    # message('Assuming raw data')
    messages = c(messages, 'Assuming raw data')
    nc = ncol(x)
    if(nc < 2){# no puede haber menos de 2 columnas. Arrojar error
      stop('Number of columns <2')
      
    }else if (nc >3){# no puede haber más de 3 columnas. Arrojar error
      stop('Number of columns >3')
      
    }else{ # si hay 2 o 3 columnas, se calcula la tabla de frecuencias con la función table.ct()
      # que comprueba qué columna eliminar (en caso de nc == 3)
      
      # si el objeto es matrix, transformar a data.frame
      if(is.matrix(x)){
        x = data.frame(x)
      }
      
      ct = table.ct(x, messages = messages)
    }
  }else{
    # si la tabla es cuadrada, puede ser una tabla de contingencia o no (2 observadores clasifican 2 objetos)
    # Casos indistinguibles donde se asume que es una tabla de contingencia
    #   * si la dimensión es 2x2 y los elementos numéricos
    #   * si la dimensión es 3x3 (con primera columna indicando número de filas) y los elementos numéricos
    
    if(is.matrix(x)){
      x = as.table(x)
    }
    
    n = max(dim(x))
    if((n==2 | n == 3) & is.numeric(unlist(x))){
      # message(paste0('Assuming ', n, 'x',n, ' contingency table'))
      messages = c(messages, paste0('Assuming ', n, 'x',n, ' contingency table'))
    }
    ct = checkTable(x, messages = messages)
  }
  
  return(ct)
}





# AUXILIARY FUNCTIONS ----------
#' Compute frequency table from raw data.
#' If ncol(x)==3, one column is dropped.
#' After that, both columns are converted to factors 
#' with the same levels to ensure that we get a squared matrix
#' 
#' Returns object of class 'M1'
#' @noRd

table.ct <- function(x, empty_class = NULL, messages = NULL){
  nr = nrow(x)
  nc = ncol(x)
  if(nc == 3){ # si hay 3 columnas, hay que eliminar una (posiblemente índice de filas)
    # se asume que la que tenga un mayor número de valores distintos es la columna sobrante
    # en caso de empate, es la primera

    nval = sapply(x, function(x){
      length(unique(x)) # número de valores distintos de cada columna
    })
    

    # el indice de fila debe ser único (tantos como numero de filas)
    pos_max_nval = which(nval == nr) # posición de las columnas con tantos valores como numero de filas
    
    # Si ninguna columna tiene tantos valores como filas, arrojar error
    if(length(pos_max_nval)==0){
      stop("Only 2 raters are allowed and ID column must have unique values for each row.")
    }
    
    # Si solo hay una, eliminarla y mandar mensaje
    if(length(pos_max_nval) == 1){
      # message(paste('Column', names(pos_max_nval), 'is assumed to be an ID column and will be removed'))
      messages = c(messages, paste('Column', names(pos_max_nval), 'is assumed to be an ID column and will be removed'))
      rm = pos_max_nval
    }
    
    # Si hay empates (más de una columna con tantos valores como filas), lanzar aviso
    if(length(pos_max_nval)>1){
      # Si la primera columna está incluida, elegirla como ID
      if(1%in%pos_max_nval){
        # warning(paste0('First column, named ', names(pos_max_nval)[1], ', is assumed to be an ID column and will be removed'))
        messages = c(messages, paste0('First column, named ', names(pos_max_nval)[1], ', is assumed to be an ID column and will be removed'))
        rm = 1
      }else{# si no está la primera, elegir la última como ID
        # warning(paste0('Last column, named ', names(which(pos_max_nval==3)), ', is assumed to be an ID column and will be removed'))
        messages = c(messages, paste0('Last column, named ', names(which(pos_max_nval==3)), ', is assumed to be an ID column and will be removed'))
        rm = 3
      }
      
    }
    
    x = x[,-rm] # eliminar columna

  }
  
  # niveles (categorías): nos aseguramos de obtener una matriz cuadrada
  lev = unique(sort(unlist(sapply(x, unique))))
  
  x1 = as.data.frame(lapply(x, factor, levels = lev))
  x1 = table(x1)
  class(x1) <- c('table', 'M1')
  
  # Add attributes related to valid and empty classes
  classes = sort(c(rownames(x1),  as.character(empty_class)))
  attr(x1, "classes") = classes
  attr(x1, "valid_classes") = rownames(x1)
  attr(x1, "empty_classes") = as.character(empty_class)
  attr(x1, "messages") = as.character(messages)
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Change name of raters (if no name was given by user) 
  if(all(names(dimnames(x1)) == c("Var1", "Var2"))){
    names(dimnames(x1)) = c('R', 'C')
  }
  
  return(x1)
}


#' If object is of class 'table', make sure it is a squarred matrix
#' @noRd
checkTable = function(x, messages = NULL){
  if(!is.table(x)){
    stop("Object is not of class 'table'")
  }
  
  # Transformar table a data.frame
  # Devuelve una columna 'Freq'
  x1 = data.frame(x)
  
  # Almacenar clases vacías para informe
  ri = marginSums(x, 1)
  ci = marginSums(x, 2)
  checkClass = ri==0 & ci ==0
  if(any(checkClass)){
    empty_class = names(ri[checkClass])
    # message("The following empty classes have been detected: ", paste(empty_class, collapse = ', '))
    messages = c(messages, paste0("The following empty classes have been detected: ", paste(empty_class, collapse = ', ')))
    
  }else{
    empty_class = NULL
  }
  
  # Obtener datos en 'bruto', repitiendo las filas según Freq
  x2 = as.data.frame(lapply(x1, rep, x1$Freq))
  x2 = x2[,-3]
  
  # Obtener tabla de frecuencias. 
  # table.ct() convierte las columnas a factores con el mismo número de levels
  # (categorías en las que se clasifican los objetos). El resultado siempre es 
  # una matriz cuadrada
  x2 = table.ct(x2, empty_class = empty_class, messages = messages)

  # Permitir decimales
  if(all(dim(x2)==dim(x))){
    cte = unique(x-x2)
    if(length(cte) == 1 & cte !=0){
      x2 = x2+cte
      
    }
  }
  return(x2)
}

