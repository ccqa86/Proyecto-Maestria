crear_poblacion <- function(n_poblacion, n_variables, nmax,miu,escenario,des, limite_inf = NULL,
                            limite_sup = NULL, verbose = TRUE) {
  
  
  n_poblacion = n_poblacion
  n_variables = n_variables
  nmax        = nmax
  limite_inf  = limite_inf
  limite_sup  = limite_sup
  verbose     = verbose %in% c(2)
  
  
  
  # Esta función crea una matriz en la que, cada fila, está formada por una
  # combinación de valores numéricos aleatorios. El rango de posibles valores
  # para cada variable puede estar acotado.
  #
  # ARGUMENTOS
  # ============================================================================
  # n_poblacion: número total de individuos de la población.
  # n_variables: longitud de los individuos.
  # limite_inf:  vector con el límite inferior de cada variable. Si solo se
  #              quiere imponer límites a algunas variables, emplear NA para
  #              las que no se quiere acotar.
  # limite_sup:  vector con el límite superior de cada variable. Si solo se
  #              quiere imponer límites a algunas variables, emplear NA para
  #              las que no se quieren acotar.
  # verbose:     mostrar información del proceso por pantalla.
  #   
  # RETORNO
  # ============================================================================
  # Una matriz de tamaño n_poblacion x n_variables que representa una población.
  
  # COMPROBACIONES
  # ----------------------------------------------------------------------------
  if (!is.null(limite_inf) & (length(limite_inf) != n_variables)) {
    stop(paste(
      "limite_inf debe tener un valor por cada variable.",
      "Si para alguna variable no se quiere límite, emplear NA.",
      "Ejemplo: lim_sup = c(10, NA, 10)"
    ))
  } else if (!is.null(limite_sup) & length(limite_sup) != n_variables) {
    stop(paste(
      "limite_sup debe tener un valor por cada variable.",
      "Si para alguna variable no se quiere límite, emplear NA.",
      "Ejemplo: lim_sup = c(10, NA, 10)"
    ))
  } else if (is.null(limite_sup) | is.null(limite_inf)) {
    warning(paste(
      "Es altamente recomendable indicar los límites dentro de los",
      "cuales debe buscarse la solución de cada variable.",
      "Por defecto se emplea [-10^3, 10^3]."
    ))
  } else if (any(any(is.na(limite_sup)), any(is.na(limite_inf)))) {
    warning(paste(
      "Los límites empleados por defecto cuando no se han definido son:",
      " [-10^3, 10^3]."
    ))
    cat("\n")
  }
  
  # Si no se especifica limite_inf, el valor mínimo que pueden tomar las variables
  # es -10^3.
  if (is.null(limite_inf)) {
    limite_inf <- rep(x = -10^3, times = n_variables)
  }
  
  # Si no se especifica limite_sup, el valor máximo que pueden tomar las variables
  # es 10^3.
  if (is.null(limite_sup)) {
    limite_sup <- rep(x = 10^3, times = n_variables)
  }
  
  # Si los límites no son nulos, se reemplazan aquellas posiciones NA por el valor
  # por defecto -10^3 y 10^3
  if (!is.null(limite_inf)) {
    limite_inf[is.na(limite_inf)] <- -10^3
  }
  
  if (!is.null(limite_sup)) {
    limite_sup[is.na(limite_sup)] <- 10^3
  }
  
  # CREAR POBLACIÓN
  # ----------------------------------------------------------------------------
  # Matriz donde almacenar los individuos generados.
  poblacion <- matrix(data = NA, nrow = n_poblacion, ncol = n_variables)

    # Bucle para crear cada individuo.
  for (i in 1:n_poblacion) {
   
    # Se crea un vector de NA que representa el individuo.
    individuo <- rep(NA, times = n_variables)
    if (escenario==1){

        # Para cada posición, se genera un valor aleatorio dentro del rango permitido
        # para cada variable.
        individuo[1] <- runif(1, min=miu, max=miu+4*des)
        individuo[2] <- runif(1, min=miu, max=individuo[1])
        individuo[3] <- runif(1, min=miu, max=individuo[2])
        individuo[4] <- runif(1, min=miu, max=individuo[3])
        individuo[5] <- runif(1, min=0.1, max=20)
        individuo[6] <- runif(1, min=0.1, max=individuo[5])
        individuo[7] <- as.integer(runif(1, min=1, max=nmax-1))
        individuo[8] <- as.integer(runif(1, min=individuo[7]+1, max=nmax+1))
        individuo[9] <- sample(limite_inf[9]:limite_sup[9], 1, replace = FALSE)
        
    }else if (escenario==2){
    # Escenario 2 para la carta VSSIWL

        # Para cada posición, se genera un valor aleatorio dentro del rango permitido
        # para cada variable.
        individuo[1] <- runif(1, min=miu, max=miu+4*des)
        individuo[2] <- runif(1, min=miu, max=individuo[1])
        individuo[3] <- individuo[1]
        individuo[4] <- runif(1, min=miu, max=individuo[3])
        individuo[5] <- runif(1, min=0.1, max=20)
        individuo[6] <- runif(1, min=0.1, max=individuo[5])
        individuo[7] <- as.integer(runif(1, min=1, max=nmax-1))
        individuo[8] <- as.integer(runif(1, min=individuo[7]+1, max=nmax+1))
        individuo[9] <- sample(limite_inf[9]:limite_sup[9], 1, replace = FALSE)
        
    }else if (escenario==3){
      #Escenario 3 para la carta VSSI
      
        # Para cada posición, se genera un valor aleatorio dentro del rango permitido
        # para cada variable.
        individuo[1] <- runif(1, min=miu, max=miu+4*des)
        individuo[2] <- runif(1, min=miu, max=individuo[1])
        individuo[3] <- individuo[1]
        individuo[4] <- individuo[2]
        individuo[5] <- runif(1, min=0.1, max=20)
        individuo[6] <- runif(1, min=0.1, max=individuo[5])
        individuo[7] <- as.integer(runif(1, min=1, max=nmax-1))
        individuo[8] <- as.integer(runif(1, min=individuo[7]+1, max=nmax+1))
        individuo[9] <- sample(limite_inf[9]:limite_sup[9], 1, replace = FALSE)
    }else if (escenario==4){
      #Escenario 4 para la carta VSI

        # Para cada posición, se genera un valor aleatorio dentro del rango permitido
        # para cada variable.
        individuo[1] <- runif(1, min=miu, max=miu+4*des)
        individuo[2] <- runif(1, min=miu, max=individuo[1])
        individuo[3] <- individuo[1]
        individuo[4] <- individuo[2]
        individuo[5] <- runif(1, min=0.1, max=20)
        individuo[6] <- runif(1, min=0.1, max=individuo[5])
        individuo[7] <- as.integer(runif(1, min=1, max=nmax-1))
        individuo[8] <- individuo[7]
        individuo[9] <- sample(limite_inf[9]:limite_sup[9], 1, replace = FALSE)
    }else if (escenario==5){
      #Escenario 5 para la carta VSS
        
        # Para cada posición, se genera un valor aleatorio dentro del rango permitido
        # para cada variable.
        individuo[1] <- runif(1, min=miu, max=miu+4*des)
        individuo[2] <- runif(1, min=miu, max=individuo[1])
        individuo[3] <- individuo[1]
        individuo[4] <- individuo[2]
        individuo[5] <- runif(1, min=0.1, max=20)
        individuo[6] <- individuo[5]
        individuo[7] <- as.integer(runif(1, min=1, max=nmax-1))
        individuo[8] <- as.integer(runif(1, min=individuo[7]+1, max=nmax+1))
        individuo[9] <- sample(limite_inf[9]:limite_sup[9], 1, replace = FALSE)
    }else{
      #Escenario 6 para la carta X tradicional
      
      # Para cada posición, se genera un valor aleatorio dentro del rango permitido
      # para cada variable.
      individuo[1] <- runif(1, min=miu, max=miu+4*des)
      individuo[2] <- runif(1, min=miu, max=individuo[1])
      individuo[3] <- individuo[1]
      individuo[4] <- individuo[2]
      individuo[5] <- runif(1, min=0.1, max=20)
      individuo[6] <- individuo[5]
      individuo[7] <- as.integer(runif(1, min=1, max=nmax-1))
      individuo[8] <- individuo[7]
      individuo[9] <- sample(limite_inf[9]:limite_sup[9], 1, replace = FALSE)
    }

    # Se añade el nuevo individuo a la población.
    poblacion[i, ] <- individuo
  }
  # INFORMACIÓN ALMACENADA EN LOS ATRIBUTOS
  # ----------------------------------------------------------------------------
  attr(poblacion, 'fecha_creacion')    <- Sys.time()
  attr(poblacion, 'numero_individuos') <- n_poblacion
  attr(poblacion, "class") <- c("matrix", "poblacion")
  
  if (verbose) {
    cat("Población inicial creada", "\n")
    cat("------------------------", "\n")
    cat("Fecha creación:", as.character(Sys.time()), "\n")
    cat("Número de individuos =", n_poblacion, "\n")
    cat("Límites inferiores de cada variable =", paste(limite_inf, collapse = ", "), "\n")
    cat("Límites superiores de cada variable =", paste(limite_sup, collapse = ", "), "\n")
    cat("\n")
  }
  
  return(poblacion)
}
