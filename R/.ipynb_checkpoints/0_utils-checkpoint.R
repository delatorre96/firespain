

#' @title remove_parenthesis
#' @description Elimina los paréntesis de un nombre
#' @param name El nombre a limpiar
#' @return El nombre sin paréntesis
#' @author juaco
#' @date 2024-04-05
#' @examples remove_parenthesis("Noáin (Valle de Elorz)/Noain (Elortzibar)") 

remove_parenthesis <- function(name) {
    gsub("\\(|\\)", "", name)
}


#' @title clean_name
#' @description Limpia un nombre de determinantes y preposiciones
#' @param name El nombre a limpiar
#' @return El nombre limpio
#' @author juaco
#' @date 2024-03-31
#' @examples #' clean_name("O Grove")
#' clean_name("BÁRDENAS, LAS")

clean_name <- function(name, trail.accents = TRUE) {
    # Lista de determinantes y preposiciones a eliminar
    stopwords <- c("el", "la", "los", "las",
                   "els", "les", "s'", "d'",  "l'", # cat
                    "o", "a", "os", "as", # gal
                   "s", "sta", "sant") # abreviatura de San o Santa
    # Divide el nombre en palabras
    words <- strsplit(name, "\\s|,")[[1]]
    # Filtra las palabras que no están en la lista de stopwords
    cleaned_words <- words[!tolower(words) %in% tolower(stopwords)]
    # Une las palabras limpias nuevamente en un nombre
    cleaned_name <- paste(cleaned_words, collapse = " ")
    # Quitar las p**** tildes y otros signos diacríticos
    if (trail.accents) {
        cleaned_name <- chartr("ÁÉÍÓÚÀÈÌÒÙÜ", "AEIOUAEIOUU", toupper(cleaned_name))
    }
    # Reemplazar las barras por operador lógico OR
    cleaned_name <- gsub("/", "|", cleaned_name)
    # Lo mismo para el signo "-", por ejemplo en San Sebastian-Donostia
    cleaned_name <- gsub("-", "|", cleaned_name)
    # Se eliminan espacios en blanco que puedan haber quedado al principio o al final
    cleaned_name <- gsub("^\\s+|\\s+$", "", cleaned_name)
    # Se eliminan los paréntesis
    cleaned_name <- remove_parenthesis(cleaned_name)
    
    return(cleaned_name)
}



#' @title get UTM zone
#' @date 2024-04-01
#' @description Devuelve la zona UTM correspondiente a una provincia
#' @param provincia El nombre de la provincia, tal y como se codifica
#' en la base de datos de fuegos original.
#' @return La zona UTM correspondiente
#' @examples get_zone("JAEN")
#' @author juaco
#' @details ## Hay 15 registros con anotación defectuosa del huso:
#' Los corregimos manualmente a partir de la provincia con esta función:
#'  (ver \url{https://images.app.goo.gl/hXfDsk6uXFYcd9RY8}

get_zone <- function(provincia) {
    switch(provincia, 
           "JAEN" = 30,
           "ALMERIA" = 30,
           "LA RIOJA" = 30,
           "ASTURIAS" = 30,
           "ALICANTE" = 30,
           "BADAJOZ" = 29,  ## coincide extremo occidental
           "ILLES BALEARS" = 31,
           "VALENCIA" = 30,
           "BURGOS" = 30)
}

#' @title write_muni_check_file
#' @description Escribe un archivo CSV con los municipios localizados en la
#'  construcción del raster y los identificados en la base de datos de incendios
#' indicando la distancia entre cadenas de caracteres según el método especificado.
#' @param logfile El archivo de log con los municipios localizados
#' @param outfile El archivo de salida
#' @param method El método de comparación de cadenas de caracteres.
#' Ver \code{help(stringdist::stringdist)} para más detalles.
#' @return NULL. La fución escribe el fichero de salida indicado
#'  en \code{outfile}, en ormato csv.
#'  @author juaco
#'  @date 2024-04-05

write_muni_check_file <- function(logfile, outfile, method = "osa") {
    method <- match.arg(method,
                        choices = c("osa", "lv", "dl", "hamming",
                                    "lcs", "qgram", "cosine",
                                    "jaccard", "jw", "soundex"))
    a <- read.delim(logfile, header = FALSE,
                    stringsAsFactors = FALSE, sep = ";")
    # Se retienen los campos de interés
    b <- a[,c(2,4)]
    a <- NULL; gc()
    # Se eliminan los duplicados
    bmat <- unique(b)
    b <- NULL; gc()
    dist <- sapply(1:nrow(bmat), function(i) {
        stringdist::stringdist(bmat[i,1], bmat[i,2])
    })
    ind.sort <- sort(dist, decreasing = TRUE, index.return = TRUE)
    # rm.ind <- which(ind.sort$x == 0)
    bmat.sorted <- bmat[ind.sort$ix,]
    bmat <- NULL; gc()
    bmat.sorted <- cbind.data.frame(bmat.sorted, "dist" = dist[ind.sort$ix])
    write.csv(bmat.sorted, file = outfile, row.names = FALSE)
    message("Saved", outfile)
}
