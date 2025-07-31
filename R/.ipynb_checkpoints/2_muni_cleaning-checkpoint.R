library(magrittr)
setwd("/lustre/gmeteo/WORK/juaco/PTI-clima")
source("fireDatabase/R/0_utils.R")
load("fireDatabase/consultas/burned_area_curated_1983-2021_WGS84.RData", verbose = TRUE)


## SE MODIFICA EL LOG PARA INSERTAR TABULADOR ANTES DE LOS MUNICIPIOS
a <- read.delim("fireDatabase/raster_outputs/fdb_muni_not_found.log",
                header = FALSE, stringsAsFactors = FALSE)
mun_not_found <- (a$V2)
mun_not_found <- gsub("\\s+en .*", "", mun_not_found, ignore.case = FALSE) 

unique(mun_not_found) 
write.csv(unique(mun_not_found), file = "fireDatabase/raster_outputs/mun_not_found.csv", row.names = FALSE)


str(fiV2str(fires)

findMuniByName <- function(name) {
  return(fires[grepl(name, fires$Municipio_Nombre, ignore.case = TRUE),])
}

findMuniByName("URROTZ")

findMuniByName("MONDARIZ")
findMuniByName("SANTIBAÑEZ DEL VALLE")

?t.test

grep("EGÜÉS", fires$Municipio_Nombre, value = TRUE, ignore.case = TRUE)

clean_name("TORRE DEL CAMPO")

oldName <- "LLOMBAY"
newName <- "LLOMBAI"

replaceMuniName <- function(oldName, newName) {
    ind <- grepl(oldName, fires$Municipio_Nombre, ignore.case = TRUE)
    fires[ind,]$Municipio_Nombre <- newName
    message(paste("Replaced", oldName, "by", newName, "in",
                  nrow(fires[ind,]), "records"))
    invisible(fires)
}


oldName <- "LLOMBAY"
newName <- "LLOMBAI"
replaceMuniName(oldName, newName)
findMuniByName(oldName)

save(fires, file = "fireDatabase/consultas/burned_area_curated_1983-2021_WGS84.RData")

## Matrix database name / shapefile name
c("LLOMBAY", "LLOMBAI",
#  "PRESES", "LES PRESES",
  "BENLLOCH", "BENLLOC",
#  "NEVES", "AS NEVES",
## "CAÑIZA" ## este mete CAÑIZARES en el grep

clean_name("ROBLEDA-CERVANTES")

findMuniByName("ROBLEDA-CERVANTES")

findMuniByName("FOZ-CALANDA")
clean_name("FOZ-CALANDA")

## ESTOS DATOS SE PIERDEN
findMuniByName("")
findMuniByName("PORTUGAL")


%>% str()
grep("NEVES", clean_name("NEVES, AS"), value = TRUE)

## A partir del log inicial, vamos a ordenar el municipio de la base de datos 
## y el identificado en el shapefile por distancias netre cadenas de caracteres,
## usando el método de Levehnstein. Así, vemos al principio cuáles son
## sospechosos de estar mal identificados, e ignoramos las distancias de cero
## (perfect matching entre ambas fuentes de datos)

library(stringdist)

a <- read.delim("fireDatabase/raster_outputs/fdb_muni_found_025.log",
                header = FALSE, stringsAsFactors = FALSE, sep = ";")
# Se retienen los campos de interés
b <- a[,c(2,4)]
# Se eliminan los duplicados
bmat <- unique(b)

dist <- sapply(1:nrow(bmat), function(i) {
    stringdist::stringdist(bmat[i,1], bmat[i,2])
})

ind.sort <- sort(dist, decreasing = TRUE, index.return = TRUE)
# rm.ind <- which(ind.sort$x == 0)
bmat.sorted <- bmat[ind.sort$ix,]
df <- cbind.data.frame(bmat.sorted, "dist" = dist[ind.sort$ix])
write.csv(df, file = "fireDatabase/raster_outputs/mun_found_check_stringdists.csv", row.names = FALSE)


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

write_muni_check_file("fireDatabase/raster_outputs/fdb_muni_found_025_new.log",
                      "fireDatabase/raster_outputs/mun_found_check_stringdists.csv")

## adicion de municpios cvon guines en el diccionario

errs <- read.delim("fireDatabase/raster_outputs/ERRORES_MUNIS", sep = ";", header = FALSE) 
problemas <- errs[,4]
p <- sapply(problemas, "clean_name")
tmp <- tempfile()
write.csv(p, file = tmp, row.names = FALSE)
tmp
