library("loadeR")
library('transformeR')
library('visualizeR')
library('parallel')


#' Download and Load Fire Data
#'
#' Downloads a `.nc` file with monthly burned area data for Spain,
#' stores it in a temporary directory, and loads it as a `grid` object.
#'
#' @return A `grid` object containing the fire data.
#' @export
func.downloadData <- function(){
    temp_dir <- tempdir()
    temp_nc <- file.path(temp_dir, "fireData.nc")
    
    url <- "https://zenodo.org/records/14644902/files/firespain025_WGS84_v0.nc?download=1"
    
    download.file(url, temp_nc, mode = "wb")
    
    grid <- loadGridData(temp_nc, var = 'ba')

    return (grid)
}

#' Extract Coordinates with Non-zero Data
#'
#' Iterates over all grid coordinates and returns those that contain non-zero data values.
#'
#' @param grid A `grid` object containing fire data.
#'
#' @return A list with two vectors: \code{x} and \code{y} containing coordinates with valid data.
#' @export
func.coordenadasConDatos <- function(grid){
    coordX <- c()
    coordY <- c()
    for (x in grid$xyCoords$x){ 
        for (y in grid$xyCoords$y){ 
            subGrid = subsetGrid(grid = grid,  lonLim  = x, latLim  = y) 
            if (any(0 != subGrid$Data)){
                    coordX <- c(coordX, x)
                    coordY <- c(coordY, y)   
            }
        }
    }
    return(list('x' = coordX, 'y' = coordY))
}

#' Calculate Monthly Statistics at a Coordinate
#'
#' Calculates the monthly average (or other statistic) of burned area for a given location.
#'
#' @param grid A `grid` object with fire data.
#' @param x Longitude coordinate.
#' @param y Latitude coordinate.
#' @param func Aggregation function to apply (default is \code{mean}).
#'
#' @return A numeric vector of length 12, with one value per month.
#' @export
func.mediasMensuales <- function(grid, x, y, func = mean){
        results <- c()
        for (season in 1:12){
            new_grid <- subsetGrid(grid = grid,season = season, lonLim  = x, latLim  = y)
            var_name <- paste("season", season, sep = "_")
            assign(var_name,new_grid)
            results <- c(results, func(get(var_name)$Data))
        }
        return(results)
}

#' Parallel Wrapper to Apply a Function to Monthly Subsets
#'
#' Auxiliary function to apply \code{func.applyDataFrame} on a subset of coordinates for a given month.
#'
#' @param i Index of the month (1 to 12).
#'
#' @return A vector of results for the specified month's coordinates.
#' @export
parallel_apply <- function(i) {
  inicio <- filas_df_params_mes_inicios[i]
  final <- filas_df_params_mes_finales[i]
  parametros_mes <- df_params[inicio:final, ]
  resultado_mes <- apply(parametros_mes, 1, func.applyDataFrame)
  return(resultado_mes)
}

#' Apply Function to a Subgrid Based on Parameter Row
#'
#' Takes a row with coordinates and month, extracts the corresponding subgrid,
#' and applies a function (default: \code{mean}) to its data.
#'
#' @param row A row containing fields \code{Mes}, \code{coord_x}, and \code{coord_y}.
#' @param func A function to apply to the subgrid values (default is \code{mean}).
#'
#' @return A numeric result of applying the function to the selected grid cell.
#' @export
func.applyDataFrame <- function (row, func = mean){
    subgrid <- subsetGrid(grid, season = row["Mes"], lonLim = row["coord_x"], latLim  = row["coord_y"])
    serie <- func(subgrid$Data)
    return (serie)
}

#' Get Index Ranges for Monthly Batches
#'
#' Returns the start and end row indices for each month, based on grid combinations of longitude and latitude.
#'
#' @param grid A `grid` object with fire data.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{filas_df_params_mes_inicios}: Start row indices for each month.
#'   \item \code{filas_df_params_mes_finales}: End row indices for each month.
#'   \item \code{df_params}: A data frame with all (lon, lat, month) combinations.
#' }
#' @export
func.iniciosFinales_mes<- function (grid){
    lon = grid$xyCoords$x
    lat = grid$xyCoords$y
    df_params <- expand.grid( "coord_x" = lon, "coord_y" = lat, "Mes" = c(1:12))
    filas_df_params_mes <- (seq(0, nrow(df_params), length(lon) * length(lat)))
    
    filas_df_params_mes_finales = c()
    for (i in filas_df_params_mes[2:length(filas_df_params_mes)]){
        filas_df_params_mes_finales = c(filas_df_params_mes_finales, i)
    }
    filas_df_params_mes_inicios = c(1)
    for (i in seq(1,12)){
            lista_finales = filas_df_params_mes_finales[2:length(filas_df_params_mes_finales)-1]
            filas_df_params_mes_inicios = c(filas_df_params_mes_inicios, lista_finales[i]+1)    
    }
    filas_df_params_mes_inicios = na.omit(filas_df_params_mes_inicios)
    return (list(filas_df_params_mes_inicios,filas_df_params_mes_finales,df_params))
}


#' Convert Grid Data to Data Frame of Time Series
#'
#' Aggregates fire data from a spatial grid into a data frame where each row corresponds 
#' to a spatial coordinate and columns represent monthly aggregated values.
#' Parallel computation is used to process each month's data efficiently.
#'
#' @param grid A `grid` object containing fire data.
#' @param func A function to apply to the grid values at each coordinate and month 
#'   (default is \code{mean}).
#'
#' @return A data frame with columns \code{coord_x}, \code{coord_y}, and one column for each 
#'   month from 1 to 12 containing the aggregated values.
#'
#' @details This function uses parallel computing to accelerate the extraction of 
#'   monthly statistics. It automatically detects the number of available cores 
#'   and applies the aggregation function (e.g., \code{mean}) to each grid cell.
#'
#' @export
func.toDataFrame <- function(grid, func = mean){
    iniciosFinalesMes <- func.iniciosFinales_mes(grid)
    filas_df_params_mes_inicios <- iniciosFinalesMes[[1]]
    filas_df_params_mes_finales <- iniciosFinalesMes[[2]]
    lon = grid$xyCoords$x
    lat = grid$xyCoords$y
    df_params <- expand.grid( "coord_x" = lon, "coord_y" = lat, "Mes" = c(1:12))
    df.seriesTemporales <-  expand.grid( "coord_x" = lon, "coord_y" = lat)

    num_cores <- detectCores()
    cl <- makeCluster(num_cores)
    clusterExport(cl, c("filas_df_params_mes_inicios", "filas_df_params_mes_finales", "df_params", "func.applyDataFrame","subsetGrid","grid"), envir = environment())
    resultados <- parLapply(cl, 1:12, parallel_apply)

    stopCluster(cl)
    df.seriesTemporales <- cbind(df.seriesTemporales, resultados)
    colnames(df.seriesTemporales) <- c('coord_x','coord_y', '1','2','3','4','5','6','7','8','9','10','11','12')
    df.seriesTemporales <- df.seriesTemporales[order(df.seriesTemporales$'coord_x'),]
    
    return(df.seriesTemporales)
    
}



#####FIRE SEASON########

#' Identify Fire Season Months
#'
#' This function determines the months of the fire season from a monthly series of average burned area.
#' It uses a cumulative threshold approach to identify significant months.
#'
#' @param sereTemporalMedias A numeric vector of length 12 representing average burned area for each month.
#' @param umbral A numeric threshold between 0 and 1 (default is 0.8) to determine the fire season.
#'
#' @return A vector of months considered as fire season, or NA/NaN if not applicable.
#' @export
func.fireSeason <- function(sereTemporalMedias, umbral = 0.8){
  if (any(is.nan(sereTemporalMedias))) {
    return(NaN)
  }

  # Si hay al menos un valor distinto de 0
  if (any(sereTemporalMedias != 0)) {
    
    proporcionAreaQuemada <- ifelse(sereTemporalMedias != 0, 
                                     sereTemporalMedias / sum(sereTemporalMedias), 
                                     0)
    vector_acumulado <- cumsum(proporcionAreaQuemada)
    vector_acumulado_ordenado <- sort(vector_acumulado)
    
    # Intentar sucesivamente con distintos umbrales
    for (umbral_local in c(umbral, 0.9, 0.99, 1)) {
      vector_filtrado <- vector_acumulado_ordenado[
        vector_acumulado_ordenado > 0 & 
        vector_acumulado_ordenado <= umbral_local
      ]
      
      meses_vector <- unique(match(vector_filtrado, vector_acumulado))
      
      if (length(meses_vector) != 0) {
        return(meses_vector)
      }
    }
    
    # Si no se encuentra ningún mes, devolver NA
    return(NA)
    
  } else {
    return(NA)
  }
}

#' Reconfigure Fire Season for Bimodal Series
#'
#' Detects and separates bimodal fire seasons based on top values not being consecutive.
#'
#' @param serie A numeric vector of length 12 representing monthly values.
#'
#' @return A sorted vector of months representing the split fire seasons.
#' @export
reconfigurarFireSeason <- function(serie){ #Para bimodales de segundo filtro
    #Esto lo hacemos para dividir la fire season en dos
    #Como previamente ya sabemos que estas series tienen fire season y cumplen con el criterio del 80% ahora lo reordenamos para que se divida en dos
    maximo1 = which(serie == sort(serie, decreasing = TRUE)[1])
    maximo2 = which(serie == sort(serie, decreasing = TRUE)[2])
    maximos = unname(sort(c(maximo1, maximo2)))
    if (identical(maximos,seq(min(maximos), max(maximos)))){
        maximo3 = which(serie == sort(serie, decreasing = TRUE)[3])
        maximos = c(maximos, maximo3)
        maximos = unname(sort(maximos))
        if (identical(maximos,seq(min(maximos), max(maximos)))){
            maximo4 = which(serie == sort(serie, decreasing = TRUE)[4])
            maximos = c(maximos, maximo4)
            maximos = unname(sort(maximos))
            if (identical(maximos,seq(min(maximos), max(maximos)))){
                maximo5 = which(serie == sort(serie, decreasing = TRUE)[5])
                maximos = c(maximos, maximo5)
                maximos = unname(sort(maximos))
                if (identical(maximos,seq(min(maximos), max(maximos)))){
                    maximo6 = which(serie == sort(serie, decreasing = TRUE)[6])
                    maximos = c(maximos, maximo6)
                    maximos = unname(sort(maximos))
                    if (identical(maximos,seq(min(maximos), max(maximos)))){
                        maximo7 = which(serie == sort(serie, decreasing = TRUE)[7])
                        maximos = c(maximos, maximo7)
                        maximos = unname(sort(maximos))
                        if (identical(maximos,seq(min(maximos), max(maximos)))){
                            maximo8 = which(serie == sort(serie, decreasing = TRUE)[8])
                            maximos = c(maximos, maximo8)
                            maximos = unname(sort(maximos))
                             if (identical(maximos,seq(min(maximos), max(maximos)))){
                                    maximo9 = which(serie == sort(serie, decreasing = TRUE)[9])
                                    maximos = c(maximos, maximo9)
                                    maximos = unname(sort(maximos))
                                    if (identical(maximos,seq(min(maximos), max(maximos)))){
                                        maximo10 = which(serie == sort(serie, decreasing = TRUE)[10])
                                        maximos = c(maximos, maximo10)
                                        maximos = unname(sort(maximos))
                                        if (identical(maximos,seq(min(maximos), max(maximos)))){
                                            maximo11 = which(serie == sort(serie, decreasing = TRUE)[11])
                                            maximos = c(maximos, maximo11)
                                            maximos = unname(sort(maximos))
                                            if (identical(maximos,seq(min(maximos), max(maximos)))){
                                                maximo12 = which(serie == sort(serie, decreasing = TRUE)[12])
                                                maximos = c(maximos, maximo12)
                                                maximos = unname(sort(maximos))
                                                }else{
                                                    return (sort(maximos))
                                                }
                                            }else{
                                                return (sort(maximos))
                                        }
                                        }else{
                                           return (sort(maximos))
                                    }
                                }else{
                                    return (sort(maximos))
                                 }
                            }else{
                                return (sort(maximos))
                            }
                        }else{
                            return (sort(maximos))
                        }
                    }else{
                       return (sort(maximos))
                }
            }else{
                return (sort(maximos))
            }
        }else{
            return (sort(maximos))
        }
    }else{
        return (sort(maximos))
    }
}

#' First Filter for Bimodality Detection
#'
#' Checks whether the fire season months are non-consecutive.
#'
#' @param fireSeason_serie A list containing a vector of fire season months.
#'
#' @return TRUE if bimodal, FALSE otherwise.
#' @export
isBimodal_filtro1 <- function(fireSeason_serie) {
  # Primero chequeo NaN y NA
  if (is.na(fireSeason_serie) || (is.atomic(fireSeason_serie) && length(fireSeason_serie) == 1 && is.nan(fireSeason_serie))) {
    return(FALSE)
  } else {
    fireSeason_serie <- fireSeason_serie[[1]]
    vector_2 <- seq(min(fireSeason_serie), max(fireSeason_serie))
    if (identical(fireSeason_serie, vector_2)) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}

#' Second Filter for Bimodality Detection Based on Series Shape
#'
#' Detects bimodality from a time series by analyzing turning points and peak similarity.
#'
#' @param serie A numeric vector of monthly values.
#' @param umbral_entre_maximos Threshold for relative difference between peaks (default 0.2).
#' @param umbral_entre_incrementos Threshold for relative amplitude between sign changes (default 0.2).
#'
#' @return TRUE if series is considered bimodal, FALSE otherwise.
#' @export
isBimodal_filtro2 <- function(serie, umbral_entre_maximos = 0.2, umbral_entre_incrementos = 0.2) {
  # Chequeo NaN y NA en el vector completo
  if (any(is.na(serie)) || any(is.nan(serie))) {
    return(FALSE)
  }
  
  if (all(serie == 0)) {  # Evitamos vectores todo ceros
    return(FALSE)
  } else {
    vector_signos <- sign(diff(serie))
    cambio_signo <- c()
    maximos <- c()
    for (i in 1:(length(vector_signos)-1)) {
      if (vector_signos[i] != 0 & vector_signos[i + 1] != 0) {
        if (vector_signos[i] != vector_signos[i + 1]) {
          cambio_signo <- c(cambio_signo, i+1)
        }
        if (vector_signos[i] > vector_signos[i + 1]) {
          maximos <- c(maximos, i + 1)
        }
      }
    }
    
    if (length(cambio_signo) == 4) {
      max1 <- serie[maximos[1]]
      max2 <- serie[maximos[2]]
      semejanza_maximos <- (max2 - max1) / max2
      incrementos <- c()
      if (abs(semejanza_maximos) < umbral_entre_maximos) {
        for (i in 1:(length(cambio_signo)-1)) {
          incremento <- (serie[cambio_signo[i+1]] - serie[cambio_signo[i]]) / sum(serie)
          incrementos <- c(incrementos, incremento)
        }
      }
      # Verificamos que incrementos no esté vacío y no tenga NA/NaN antes de usar abs()
      if (length(incrementos) > 0) {
        incrementos <- incrementos[!is.na(incrementos) & !is.nan(incrementos)]
        incrementos_significativos <- incrementos[abs(incrementos) >= umbral_entre_incrementos]
        if (length(incrementos_significativos) == length(incrementos)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      } else {
        return(FALSE)
      }
      
    } else if (length(cambio_signo) == 3 & length(maximos) == 1) {
      max1 <- serie[maximos[1]]
      max2 <- serie[length(serie)]
      incrementos <- c()
      if (which(serie == max1) != 1) {
        semejanza_maximos <- (max2 - max1) / max2
        if (abs(semejanza_maximos) < umbral_entre_maximos) {
          for (i in 1:(length(cambio_signo)-1)) {
            incremento <- (serie[cambio_signo[i+1]] - serie[cambio_signo[i]]) / sum(serie)
            incrementos <- c(incrementos, incremento)
          }
        }
        # Mismo chequeo aquí
        if (length(incrementos) > 0) {
          incrementos <- incrementos[!is.na(incrementos) & !is.nan(incrementos)]
          incrementos_significativos <- incrementos[abs(incrementos) >= umbral_entre_incrementos]
          if (length(incrementos_significativos) == length(incrementos)) {
            return(TRUE)
          } else {
            return(FALSE)
          }
        } else {
          return(FALSE)
        }
      } else {
        return(FALSE)
      }
      
    } else {
      return(FALSE)
    }
  }
}

#' Apply Bimodality Filters to Fire Season Data
#'
#' Applies both bimodality filters to each row of the input data and reconfigures fire seasons accordingly.
#'
#' @param df.seriesTemporales_conCoords A data frame with coordinates and monthly time series (columns 3:14).
#'
#' @param df.fireSeason A data frame with fire season results.
#'
#' @return Updated data frame with fire season and bimodality annotations.
#' @export
func.bimodalidad <- function(df.seriesTemporales_conCoords, df.fireSeason){
    #Pasamos el primer filtro
    bimodales_1 <- apply(df.fireSeason, 1, isBimodal_filtro1)
    coord_x = df.seriesTemporales_conCoords$coord_x
    coord_y = df.seriesTemporales_conCoords$coord_y
    df_bimodales_1 = data.frame(coord_x, coord_y, bimodales_1)
    #Segundo filtro
    bimodales_2 <- c()
    for (i in 1:nrow(df.seriesTemporales_conCoords)){
                bimodales_2 <- c(bimodales_2, isBimodal_filtro2(unlist(df.seriesTemporales_conCoords[i,3:14])))
        }
    df_bimodales = cbind(df_bimodales_1,bimodales_2)
    df_bimodales$Bimodal <- df_bimodales$bimodales_1 | df_bimodales$bimodales_2
    df.fireSeason <- cbind(coord_x, coord_y, df.fireSeason, 'Bimodal'=df_bimodales$Bimodal,  'bimodales_1'=df_bimodales$bimodales_1 , 'bimodales_2'= df_bimodales$bimodales_2)
    #Reconfiguramos fire seasons nuevas
    for (i in 1:nrow(df.fireSeason)){
        if (df.fireSeason$bimodales_1[i] == FALSE & df.fireSeason$bimodales_2[i] == TRUE){
            x = df.fireSeason$coord_x[i]
            y = df.fireSeason$coord_y[i]
            serie = unlist(df.seriesTemporales_conCoords[df.seriesTemporales_conCoords$coord_x == x &  df.seriesTemporales_conCoords$coord_y == y, ][,3:14])
            df.fireSeason$FireSeason[i] = list(reconfigurarFireSeason(serie))
    
        }
    }
    df.fireSeason <- subset(df.fireSeason, select = -c(bimodales_1, bimodales_2))
    return (df.fireSeason)
    }

#' Compute Angular Positions (Sigmas) for Fire Seasons
#'
#' This function computes the angular positions (in radians) for a specified number of seasons.
#'
#' @param numberSeasons Integer. The number of seasons (e.g., 12 for monthly).
#'
#' @return A numeric vector of angular positions.
#' @export
sigma_m <- function(numberSeasons){
    sigmas = c()
    for (m in 1:numberSeasons){
        sigmas = c(sigmas, 2*pi * (m-1)/numberSeasons)
    }
    return (sigmas)
}

#' Characterize Fire Season Using Seasonal Concentration and Timing
#'
#' Calculates the seasonal concentration (C) and seasonal timing (P) of fire occurrence.
#'
#' @param sereTemporalMedias Numeric vector of average monthly values.
#' @param numberSeasons Integer. Default is 12.
#'
#' @return A list containing 'C' (concentration) and 'P' (timing) or NA.
#' @export
func.caracterizacion_fireSeason <- function(sereTemporalMedias, numberSeasons = 12) {
  # Comprobar NA o NaN en el vector
  if (any(is.na(sereTemporalMedias)) || any(is.nan(sereTemporalMedias))) {
    return(NA)
  }
  
  if (any(0 != sereTemporalMedias)) {
    sigmas = sigma_m(numberSeasons = numberSeasons)
    # mediasMensuales es una lista de vectores en donde cada vector son 12 medias mensuales
    x <- sereTemporalMedias
    L_x_vector <- c()
    L_y_vector <- c()
    for (m in 1:12) {
      L_x_vector <- c(L_x_vector, x[m] * cos(sigmas[m]))
      L_y_vector <- c(L_y_vector, x[m] * sin(sigmas[m]))
    }
    L_x = sum(L_x_vector)
    L_y = sum(L_y_vector)
    # seasonal concentration
    C = (sqrt(L_x^2 + L_y^2)) / sum(x)
    # seasonal timing 
    P = atan(L_x / L_y)  
    ####Modulo de la seasonal timing
    # P = (P + 2*pi) %% (2*pi) 
    return(list('C' = C, 'P' = P)) 
  } else {
    return(NA)
  }
}

#' Convert Fire Season Timing Phase to Month Index
#'
#' Transforms the computed phase into a month-like index using trigonometric quadrant checks.
#'
#' @param sereTemporalMedias Numeric vector of average monthly values.
#' @param numberSeasons Integer. Default is 12.
#'
#' @return Numeric value indicating the estimated month.
#' @export
func.phase2meses <- function(sereTemporalMedias,numberSeasons = 12){
    if (any(0 != sereTemporalMedias)){
        sigmas = sigma_m(numberSeasons = numberSeasons)
        #mediasMensuales es una lista de vectores en donde cada vector son 12 medias mensuales
        x <- sereTemporalMedias
        L_x_vector <- c()
        L_y_vector <- c()
        for (m in 1:12){
            L_x_vector <- c(L_x_vector, x[m] * cos(sigmas[m]))
            L_y_vector <- c(L_y_vector, x[m] * sin(sigmas[m]))
        }
        L_x = sum(L_x_vector)
        L_y = sum(L_y_vector)
        P = atan(L_x / L_y)  
        
        if (L_x > 0 & L_y > 0){
            m = 0.5*pi - P
        }else if (L_x < 0 & L_y > 0){
            m = 0.5*pi + P
        }else if (L_x < 0 & L_y < 0){
            m = 1.5*pi - P
        }else if (L_x > 0 & L_y < 0){
            m = 1.5*pi + P
        }
        return(m) 
        
    }else{
        return (NA)
    }
}

#' Extract Coordinates from Data Frame Row Names
#'
#' Parses the row names of a data frame to extract coordinate pairs.
#'
#' @param df A data frame with row names in the format 'x_y'.
#'
#' @return A list with numeric vectors 'x' and 'y'.
#' @export
getCoordsFromDataFrame <- function(df){
        names = rownames(df)
        vector_x <- c()
        vector_y <- c()
        for (name in names){
            x = strsplit(name,split = '_')[[1]][1]
            y = strsplit(name,split = '_')[[1]][2]
            vector_x <- c(vector_x, as.numeric(x))
            vector_y <- c(vector_y, as.numeric(y))
        
       
    }
     return(list('x' = vector_x, 'y' = vector_y))
}

#' Standardize Fire Season Data Frame Columns
#'
#' Cleans and standardizes the columns in a fire season data frame.
#'
#' @param df A data frame with fire season data.
#'
#' @return A cleaned data frame with standardized numeric columns.
#' @export
procesar_fireSeason_dataframe <- function(df) {
  
  # Estandarizar la columna FireSeason
  df$FireSeason <- lapply(df$FireSeason, function(x) {
    if (is.null(x) || all(is.na(x)) || all(is.nan(x))) {
      return(NaN)
    } else {
      return(x)
    }
  })
    # Estandarizar la columna SeasonalConcentration
    df$SeasonalConcentration <- sapply(df$SeasonalConcentration, function(x) {
      if (is.null(x) || all(is.na(x)) || all(is.nan(x))) {
        return(NaN)
      } else {
        return(as.numeric(x))
      }
    })
    # Estandarizar la columna SeasonalTiming
  df$SeasonalTiming <- sapply(df$SeasonalTiming, function(x) {
    if (is.null(x) || all(is.na(x)) || all(is.nan(x))) {
      return(NaN)
    } else {
      return(x)
    }
  })
  
  # Extraer MainFirstSeason y MainLastSeason
  df$MainFirstSeason <- sapply(df$FireSeason, function(x) {
      return(x[1])
    })
  
  df$MainLastSeason <- sapply(df$FireSeason, function(x) {
      return(x[length(x)])
    })
  
  # Convertir Bimodal: TRUE → 1, FALSE → 0, NA → NaN
  df$Bimodal <- sapply(df$Bimodal, function(x) {
    if (is.na(x)) {
      return(NaN)
    } else if (x) {
      return(1)
    } else {
      return(0)
    }
  })
  
  return(df)
}


#' Generate Fire Season Data Frame from Grid
#'
#' Converts a spatial grid to a fire season data frame with temporal and seasonal metrics.
#'
#' @param grid A climate grid object containing fire data.
#'
#' @return A data frame containing fire season statistics.
#' @export
func.FromGridToFireSeasonDF <- function (grid){
    ##Hacemos la media de todos los eneros, de todos los febreros,... de todos los meses para cada gridBox:
    df.seriesTemporales_conCoords <- func.toDataFrame(grid = grid,func = mean)
    
    df.seriesTemporales <- df.seriesTemporales_conCoords[,3:14]
    #####Calculamos la Fire Season usando las series temporales (las medias de cada mes)
    df.fireSeason <- data.frame(t(data.frame(t(apply(df.seriesTemporales, 1, func.fireSeason)))))
    names(df.fireSeason)[ncol(df.fireSeason)] <- 'FireSeason'
    rownames(df.fireSeason) <- NULL
    df.fireSeason <- func.bimodalidad(df.seriesTemporales_conCoords, df.fireSeason)
    ##Calculamos la Seasonal Concentration y el Seasonal Timing y los incluimos en el data frame de la fire season
    vector_c <- c()
    vector_p <- c()
    for (i in 1:nrow(df.seriesTemporales)){
    carFS <- func.caracterizacion_fireSeason(unlist(df.seriesTemporales[i,]))
    if (is.na(carFS[1])){
        vector_c <- c(vector_c, 0)
        vector_p <- c(vector_p, 0)
    }else{
        C = carFS$C
        P = carFS$P
        vector_c <- c(vector_c, C)
        vector_p <- c(vector_p, P)  
    }
    }
    
    df.fireSeason <- cbind(df.fireSeason,vector_c)
    names(df.fireSeason)[ncol(df.fireSeason)] <- 'SeasonalConcentration'
    df.fireSeason <- cbind(df.fireSeason,vector_p)
    names(df.fireSeason)[ncol(df.fireSeason)] <- 'SeasonalTiming'
    df.fireSeason$SeasonalTiming[is.na(df.fireSeason$FireSeason)] <- NA
    df.fireSeason$SeasonalConcentration[is.na(df.fireSeason$FireSeason)] <- NA
    df.fireSeason$Bimodal[is.na(df.fireSeason$FireSeason)] <- NA
    df.fireSeason = procesar_fireSeason_dataframe(df.fireSeason)
    return (df.fireSeason)
}


#' Convert Quantity Vector to Grid Object
#'
#' Maps a vector of quantities onto a grid object for visualization or export.
#'
#' @param quantity Numeric vector of values to map.
#' @param what Description of the quantity (e.g., 'mean').
#' @param ref.grid A reference grid to match dimensions.
#' @param backperm Optional permutation index.
#'
#' @return A grid object with mapped data.
#' @export
quantity2clim <- function(quantity, what, ref.grid, backperm = NULL) {
  if(!is.null(backperm)){quantity <- quantity[backperm]}
  mat <- matrix(quantity, nrow = 1)  
  ref.grid$Data <- mat2Dto3Darray(mat, x = ref.grid$xyCoords$x , y = ref.grid$xyCoords$y)
  attr(ref.grid$Data, "climatology:fun") <- what
  return(ref.grid)
}


#' Analyze clustering structure of fire season data
#'
#' Performs exploratory clustering analysis on fire season features.
#' It filters invalid observations, scales the selected variables, 
#' applies k-means for multiple values of k, and plots both the elbow method 
#' and a hierarchical clustering dendrogram to help assess the number of clusters.
#'
#' @param df A data frame containing fire season features (e.g., `SeasonalTiming`, `SeasonalConcentration`, etc.).
#' @param max_k Maximum number of clusters to try in the elbow method (default is 10).
#' @param conCoords Logical. Whether to include spatial coordinates (`coord_x`, `coord_y`) in the clustering (default is FALSE).
#'
#' @return No value is returned. The function produces two plots:
#'   - Elbow method for k-means clustering
#'   - Hierarchical clustering dendrogram
#'
#' @export
analisis_clustering_fireSeason <- function(df, max_k = 10, conCoords = FALSE) {
    
    # 1. Filtrar filas con FireSeason no NaN
    is_valid_fireseason <- sapply(df$FireSeason, function(x) {
    !(is.numeric(x) && length(x) == 1 && is.nan(x))
    })
    
    df_valid <- df[is_valid_fireseason, ]
    
    # 2. Seleccionar variables relevantes
    if (conCoords) {
    vars_to_cluster <- df_valid[, !(names(df_valid) %in% c("FireSeason", "Cluster"))]
    } else {
    vars_to_cluster <- df_valid[, !(names(df_valid) %in% c("FireSeason", "Cluster", "coord_x", "coord_y"))]
    }
    
    # 3. Asegurar que los valores sean numéricos válidos
    vars_to_cluster$SeasonalConcentration <- sapply(df_valid$SeasonalConcentration, function(x) ifelse(is.null(x) || is.nan(x), NA, x))
    vars_to_cluster$SeasonalTiming <- sapply(df_valid$SeasonalTiming, function(x) ifelse(is.null(x) || is.nan(x), NA, x))
    
    # 4. Eliminar filas con NA
    complete_rows <- complete.cases(vars_to_cluster)
    df_clean <- df_valid[complete_rows, ]
    vars_clean <- vars_to_cluster[complete_rows, ]
    
    # 5. Escalar las variables
    vars_scaled <- scale(vars_clean)
                                   
    # 7. Graficar el método del codo
    inercias <- numeric(max_k)
    for (k in 1:max_k) {
    inercias[k] <- kmeans(vars_scaled, centers = k, nstart = 10)$tot.withinss
    }
    
    par(cex.axis = 1.5)
    plot(1:max_k, inercias, type = "b", pch = 19, frame = FALSE,
    xlab = "Número de Clusters (k)",
    ylab = "Inercia total (tot.withinss)",
    main = "Método del Codo")
    # 7. Clustering jerárquico
    d <- dist(vars_scaled)
    hc <- hclust(d)
    
    # 8. Dendrograma
    par(cex.axis = 1.5)
    plot(hc, main = paste("Dendrograma jerárquico"),
    xlab = "", sub = "", cex = 0.6)

  
}

#' Perform k-means clustering on fire season data
#'
#' Applies k-means clustering to fire season features, excluding missing or invalid rows,
#' and optionally including spatial coordinates. Returns the original data frame with an 
#' additional column `Cluster` indicating the cluster assignment.
#'
#' @param df A data frame containing fire season data and variables such as `SeasonalConcentration` and `SeasonalTiming`.
#' @param k Number of clusters to form.
#' @param conCoords Logical. Whether to include spatial coordinates (`coord_x`, `coord_y`) in the clustering (default is FALSE).
#'
#' @return A data frame identical to `df` with an additional `Cluster` column indicating the assigned cluster for each observation.
#'
#' @export
clustering_fireSeason <- function(df, k, conCoords = FALSE) {
  
    #Identificar qué filas tienen FireSeason como NaN
    is_nan_fireseason <- sapply(df$FireSeason, function(x) {
    is.numeric(x) && length(x) == 1 && is.nan(x)
    })
    
    # Inicializar columna de Cluster
    df$Cluster <- NA
    
    # Asignar Cluster 0 a las filas sin FireSeason
    df$Cluster[is_nan_fireseason] <- NA
    
    #Filtrar las filas con FireSeason válida
    df_valid <- df[!is_nan_fireseason, ]
    
    #Seleccionar solo columnas  relevantes
    if (conCoords == TRUE){
        vars_to_cluster <- df_valid[, !(names(df_valid) %in% c("FireSeason",'Cluster'))]
    }else{
        vars_to_cluster <- df_valid[, !(names(df_valid) %in% c("FireSeason",'Cluster','coord_x','coord_y'))]
    }
    
    # Convertir listas a vectores numéricos (SeasonalTiming, SeasonalConcentration)
    vars_to_cluster$SeasonalConcentration <- sapply(df_valid$SeasonalConcentration, function(x) ifelse(is.null(x) || is.nan(x), NA, x))
    vars_to_cluster$SeasonalTiming <- sapply(df_valid$SeasonalTiming, function(x) ifelse(is.null(x) || is.nan(x), NA, x))
    
    # Eliminar filas con NAs en las variables seleccionadas
    complete_cases <- complete.cases(vars_to_cluster)
    df_valid_clean <- df_valid[complete_cases, ]
    vars_to_cluster_clean <- vars_to_cluster[complete_cases, ]
    
    #Aplicar K-means
    kmeans_result <- kmeans(scale(vars_to_cluster_clean), centers = k)
    
    #Asignar los clusters a las observaciones originales
    df_valid_clean$Cluster <- kmeans_result$cluster
    
    #Combinar los resultados
    df_final <- df
    df_final[rownames(df_valid_clean), "Cluster"] <- df_valid_clean$Cluster
    
    return(df_final)
}

