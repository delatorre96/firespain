
## /////////////////////////////////////////////////////////////////////////////
## SpatialPolygonsDataFrame de centroides de municipios ------------------------
## /////////////////////////////////////////////////////////////////////////////

# # Nota: centroides calculado en QGIS por simplicidad
# ## load shapefile
# library(rgdal)
# munis <- readOGR("/lustre/gmeteo/WORK/juaco/PTI-clima/fireDatabase/carto/myQGIS/muni_centroids.shp")
# str(munis)
# head(munis)
# munis@data[,"NAMEUNIT"]
# class(munis)
# munis <- munis[,-c(1,2,3,4,5,7,8,9)]

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save(munis, file = "/lustre/gmeteo/WORK/juaco/PTI-clima/fireDatabase/carto/municipios_SPDF.RData")
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# # Idem para los centroides en UTM -ETRS89
# 
# library(rgdal)
# munis <- readOGR("/lustre/gmeteo/WORK/juaco/PTI-clima/fireDatabase/carto/myQGIS/muni_centroids_UTM_ETRS89.shp")
# str(munis)
# head(munis)
# munis@data[,"NAMEUNIT"]
# class(munis)
# munis <- munis[,-c(1,2,3,4,5,7,8,9)]
# 
# ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save(munis, file = "/lustre/gmeteo/WORK/juaco/PTI-clima/fireDatabase/carto/municipios_ETRS89_SPDF.RData")
# ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


## /////////////////////////////////////////////////////////////////////////////
## Creación de un raster regular 0.25 de referencia ----------------------------
## /////////////////////////////////////////////////////////////////////////////
# library(sp)
# library(raster)
# library(magrittr)
# setwd("/lustre/gmeteo/WORK/juaco/PTI-clima")
# ## Extent de la capa municipios como referencia
# load("/lustre/gmeteo/WORK/juaco/PTI-clima/fireDatabase/carto/municipios_SPDF.RData",
#      verbose = TRUE)
# ext <- extent(munis)
# 
#  # Crear un objeto raster con resolución diferente en x (longitud) y y (latitud)
#  res_x <- 0.25 # 0.5 #1.0 # Resolución en grados para la longitud
#  res_y <- 0.25 # 0.5 #1.0 # Resolución en grados para la latitud
#  # Asignamos valor cero por defecto 
#  raster_base <- raster(ext, res = c(res_x, res_y), vals = 0)
# 
# ## Mascara de Tierra
#  library(rgdal)
#  muni_shp <- readOGR("/lustre/gmeteo/WORK/juaco/PTI-clima/fireDatabase/carto/lineas_limite/SHP_ETRS89/recintos_municipales_inspire_peninbal_etrs89/recintos_municipales_inspire_peninbal_etrs89.shp")
#  pts <- as(raster_base, "SpatialPoints")
#  proj4string(pts) <- CRS(proj4string(muni_shp))
#  ov <- over(pts, muni_shp)
#  str(ov)
#  na.mask <- which(is.na(ov$NAMEUNIT))
#  values(raster_base)[na.mask] <- NA
# 
# # Comprobar que la máscara es correcta:
#  pdf("raster_base.pdf", width = 8, height = 8)
#  plot(raster_base) %>% print()
#  dev.off()
# 
# # #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save(raster_base,
#      file = paste0("/lustre/gmeteo/WORK/juaco/PTI-clima/fireDatabase/carto/raster_base_",
#      res_x, "deg_landmasked.RData"))
# # #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## /////////////////////////////////////////////////////////////////////////////
## Ajustes en la base de datos de fuegos (fechas etc.) y limpieza de datos------
## /////////////////////////////////////////////////////////////////////////////

# library(magrittr)
# setwd("/lustre/gmeteo/WORK/juaco/PTI-clima")
# bdfile <- "fireDatabase/consultas/consulta_basica_2.txt"
#  
# fires <- read.table(file = bdfile,
#                     header = TRUE, sep = ";", quote = "\"", dec = ",",
#                     stringsAsFactors = FALSE, fileEncoding = "latin1")
# # str(fires)
#  
# ## Manejo de fechas
# fires$start_date <- as.Date(fires$deteccion, format = "%d/%m/%Y %H:%M:%S")
# fires$end_date <- as.Date(fires$extinguido, format = "%d/%m/%Y %H:%M:%S")
# 
# ## duracion del evento:
# fires$duracion <- fires$end_date - fires$start_date
# 
# # fires$duracion %>%  table()
# ## Hay un buen numero de fuegos cuya duracion no es creible. Se trata de errores
# ## en la anotación de las fechas de extinción.
# ## Ejemplos:
# # 
# # which(fires$duracion < 0)
# # fires[which(fires$duracion < 0),] ## fuegos que terminan en una fecha imposible
# # 
# # which(fires$duracion > 30)
# # fires[which(fires$duracion > 30),] ## fuegos que duran más de 30 días
# 
# ## Eliminar filas innecesarias
# keep.rows <- which(fires$start_date >= as.Date("1983-01-01"))
# fires <- fires[keep.rows,]
# # str(fires)
# 
# ## Se elimina Canarias (problemas de proyecciones y poco representativo a 0.25) --------
# fires$Provincia_Nombre %>% unique()
# can.ind <- which(fires$Provincia_Nombre == "LAS PALMAS" | fires$Provincia_Nombre == "S.C. TENERIFE")
# # fires[can.ind,"Provincia_Nombre"] 
# # str(fires)
# fires <- fires[-can.ind,]
# 
# 
# ## Ajuste coordenadas geográficas ----------------------------------------------
# library(sp)
# library(rgdal)
# 
# # https://gis.stackexchange.com/a/331057
# # For data in Spain, the string would likely be:
# # "+proj=utm +zone=29 +ellps=GRS80 +units=m +no_defs "
# # or
# # "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
# 
# # En nuestro caso, la zona UTM está indicada en la columna fires$huso
# unique(fires$huso)
# fires$huso %>% table()
# ## Corrección de husos mal anotados --------------------------------------------
# ## Hay registros con anotación defectuosa del huso:
# ## Los corregimos manualmente a partir de la provincia:
# ## https://images.app.goo.gl/hXfDsk6uXFYcd9RY8
# 
# huso.defectuoso <- which(fires$huso < 29 | fires$huso > 31)
# source("fireDatabase/R/0_utils.R")
# for (i in 1:length(huso.defectuoso)) {
#     ind <- huso.defectuoso[i]
#     fires$huso[ind] <- get_zone(fires$Provincia_Nombre[ind])
# }
# 
# # if (length(which(fires$huso < 27 | fires$huso > 31)) > 0) {
# #     stop("Aun hay husos defectuosos")
# # }
# 
# ## Transformación de proyectadas ETRS89 a greográficas WGS84 -------------------
# 
# ## Inicializa campos lon-lat
# fires$lon <- fires$lat <- NA
# str(fires)
# # Comprobacion, solo queda ETRS89 y husos 29, 30 y 31
# fires$Descripcion %>%  table() 
# fires$huso %>% table()
# 
# ## Tarda mucho punto a punto. Dividimos en chunks por zonas UTM.
# zn <- c(29,30,31)
# for (z in zn) {
#     message("Transforming coords Zone UTM ", z)
#     xy.rows <- which(!is.na(fires$x) & !is.na(fires$y) & fires$huso == z) 
#     xymat.UTM <- fires[xy.rows, c("x","y")] %>% as.matrix()
#     xysp <- SpatialPoints(coords = xymat.UTM,
#                           proj4string = CRS(paste0("+proj=utm +zone=",
#                                                    as.character(z),
#                                                    " +ellps=GRS80 +units=m +no_defs")))
#     xysp <- spTransform(xysp, CRS("+proj=longlat +datum=WGS84"))
#     fires$lon[xy.rows] <- xysp@coords[,1]
#     fires$lat[xy.rows] <- xysp@coords[,2]
# }
# 
# ## Eliminar columnas innecesarias
# str(fires)
# rm.ind <- c(1:6,10,11)
# fires <- fires[,-rm.ind]
# str(fires)

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 # save(fires,
 #      file = "fireDatabase/consultas/burned_area_curated_1983-2021_WGS84.RData")
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## /////////////////////////////////////////////////////////////////////////////
## PROGRAMA PRINCIPAL ----------------------------------------------------------
## /////////////////////////////////////////////////////////////////////////////

rm(list=ls())
gc()

library(magrittr)
library(raster)
library(sp)
library(log4r)

# setwd("/lustre/gmeteo/WORK/juaco/PTI-clima")
source("fireDatabase/R/0_utils.R")

load("fireDatabase/consultas/burned_area_curated_1983-2021_WGS84.RData",
     verbose = TRUE)

# Generar un rango de fechas continuo 1983-01-01 a 2021-12-31
# fecha_minima <- min(fires$start_date)
fecha_minima <- as.Date("1983-01-01")
# fecha_maxima <- as.Date("1983-03-31")
fecha_maxima <- as.Date("2021-12-31") # max(fires$start_date)
rango_fechas <- seq(from = fecha_minima, to = fecha_maxima, by = "days")

## rango_fechas es un vector continuo diario desde el inicio hasta el fin de la
## serie de datos. Se trata de iterar sobre el día a día para ir extrayendo las
## áreas quemadas diariamente


## El objeto munis contiene los centroides de todos los municipios 
## de España peninsular y Baleares
load("fireDatabase/carto/municipios_SPDF.RData", verbose = TRUE)
ref_muni_names <- sapply(as.character(munis$"NAMEUNIT"), "clean_name") 

## La tabla voc.munis contiene un vocabulario de municipios que mapea
## nombres conflictivos entre la base de datos de fuegos
## y el shapefile de municipios
voc.munis <- read.csv("fireDatabase/municipios_vocabulary.csv",
                      header = TRUE, sep = ",",
                      stringsAsFactors = FALSE, na.strings = "")[,1:2] 
# tail(voc.munis)
# str(voc.munis)

## El objeto raster_base contiene la rejilla de referencia,
## con todas las celdas de tierra con valor inicial cero por defecto
load("fireDatabase/carto/raster_base_0.25deg_landmasked.RData",
     verbose = TRUE)

## Preparacion del cubo de datos vacío:
arr <- array(data = NA, dim = c(length(rango_fechas),
                                nrow(raster_base),
                                ncol(raster_base)))
attr(arr, "dimensions") <- c("time", "lat", "lon")

## Log especial para detectar problemas de nomeclatura de municipios
logfile1 <- 'fireDatabase/raster_outputs/test_fdb_muni_not_found_025.log'
log_notfound <- create.logger(logfile = logfile1, level = "INFO")

logfile2 <- 'fireDatabase/raster_outputs/test_fdb_muni_found_025_new.log'
log_munis <- create.logger(logfile = logfile2, level = "INFO")

## Iterar sobre las fechas

raster.list <- lapply(1:length(rango_fechas), function(i) {
    # i=22
    # i=1 ## Ejemplo de día ausente en base de datos (se asignan 0's a todo)
    # i=2 ## Ejemplo de día con un fuego en la base de datos
    # i=5 ## Ejemplo de día con varios fuegos en la base de datos
    # i=13570 ## Ejemplo de día con varios fuegos geolocalizados
    ## grep("Llombay", fires$Municipio_Nombre, ignore.case = TRUE)[1]
    ## Ejemplo de municipio con nombre conflictivo
    ## i=6
    
    r <- raster_base
    
    # Aqui he tenido que convertir las fechas a character string
    ## para asegurar el matching (usando Dates no coincidían -??-)
    fecha <- rango_fechas[i] %>% as.character()
    loc.ind <- which(as.character(fires[["start_date"]]) == fecha)

    # fires[loc.ind,]    
    
    if (length(loc.ind) > 0) {
        
        for (j in 1:length(loc.ind)) {
            
            # j=1
            ## Localización por coordenadas exactas
            if (!is.na(fires$lon[loc.ind[j]]) & !is.na(fires$lat[loc.ind][j])) {
                co <- c(fires$lon[loc.ind[j]], fires$lat[loc.ind[j]])
                celda <- raster::cellFromXY(r, co)
                valor_punto <- fires$area_total[loc.ind[j]]
                ## sum con na.rm=TRUE para evitar NA's en Menorca, Ibiza
                ## y quizás en otros puntos costeros,
                ## inicialmente enmascarados
                values(r)[celda] <- sum(values(r)[celda],
                                        valor_punto, na.rm = TRUE)
             
            }
            
            ## Localización por coordenada del municipio  
            else {
                # j=5 con i=6 Pilla LLombay (ejemplo de uso de vocabulario)
                mun <- fires$Municipio_Nombre[loc.ind[j]] %>% clean_name()
                
                if(mun == "INDETERMINADO") {
                    info(log_notfound,
                         paste("Municipio NO DETERMINADO en", fecha))
                }
                
                ## Municipios mapeados a vocabulario:
                ind.voc <- match(mun, voc.munis$cleaned_name)
                
                ## Aqui tenemos cuidado con los parentesis dentro de las regexp
                if (!is.na(ind.voc)) {
                    mun <- voc.munis$muni_shp_name[ind.voc] %>%
                        remove_parenthesis()
                    loc.ind2 <- grep(paste0("^", mun, "$"),
                                     remove_parenthesis(names(ref_muni_names)))
                } 
                
                ## Municipios no mapeados (encontramos con pattern matching
                ## usando clean_name) 
                else { 
                    loc.ind2 <- grep(paste0("^", mun, "$"), ref_muni_names)
                    
                }
                
                ## Si no se encuentra ese municipio se anota en el log
                if (length(loc.ind2) == 0) {
                    info(log_notfound,
                         paste("Municipio no encontrado:", mun, "en", fecha))
                }
                
                else {
                    
                    ## Se asume el primer match (cuidado, revisar el log!)
                    ## se ha creado la funcion 'write_muni_check_file' para
                    ## revisar a posteriori
                    
                    loc.ind2 <- loc.ind2[1] 
                    
                    info(log_munis,
                         paste0("Tomado ;", ref_muni_names[loc.ind2],
                                "; como ;",
                                fires[["Municipio_Nombre"]][loc.ind[j]],
                                "; el ", fecha))
                    
                    co <- coordinates(munis)[loc.ind2,]
                    celda <- raster::cellFromXY(r, co)
                    # fires[loc.ind[j],]
                    valor_punto <- fires$area_total[loc.ind[j]]
                    
                    ## sum con na.rm=TRUE para evitar NA's en Menorca, Ibiza
                    ## y quizás en otros puntos costeros,
                    ## inicialmente enmascarados
                    values(r)[celda] <- sum(values(r)[celda],
                                            valor_punto, na.rm = TRUE)
                
                }
            }
        }
    }
    
    # Fecha no encontrada en base de datos (dia sin fuegos)
    else {
        message(paste("No fires on", fecha))
    }
    return(r)
    # arr[i,,] <- r
})


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save(raster.list,
#      file = "fireDatabase/raster_outputs/r025/raster_list_025deg.RData")
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


load("fireDatabase/raster_outputs/r025/raster_list_025deg.RData")

raster_stack <- raster::stack(raster.list)

message("[",
        Sys.time(),
        "] 3D Raster Stack successfully created\nWriting netCDF file...")


raster.list <- NULL
gc()


library(ncdf4)


# output filename
filename <- "fireDatabase/raster_outputs/r025/fire_025deg.nc"


# Longitude and Latitude data
xvals <- unique(values(init(raster_stack, "x")))
yvals <- unique(values(init(raster_stack, "y")))
nx <- length(xvals)
ny <- length(yvals)
lon <- ncdim_def("longitude", "degrees_east", xvals)
lat <- ncdim_def("latitude", "degrees_north", yvals)


# Missing value to use
mv <- -9999

# Time component
ntimes <- dim(raster_stack)[3]
time <- ncdim_def(name = "time", 
                  units = "days since 1983-01-01T00:00:00Z", 
                  vals = 1:ntimes, 
                  unlim = TRUE,
                  longname = "Day")

# Define the precipitation variables
var_ba <- ncvar_def(name = "ba",
                    units = "ha",
                    dim = list(lon, lat, time),
                    longname = "Total_daily_burned_area",
                    missval = mv,
                    compression = 9)

# Add the variables to the file
ncout <- nc_create(filename,
                   list(var_ba), ## Esta lista puede contener varias variables
                   force_v4 = TRUE)
message(paste("The file has", ncout$ndim,"dimensions"))

# add some global attributes
ncatt_put(ncout, 0,
          "Title",
          "Daily burned area 0.25 regular gridded dataset for Spain")
ncatt_put(ncout, 0,
          "Source",
          "Estadística General de Incendios Forestales (EGIF), MITECO")
ncatt_put(ncout, 0,
          "References",
          "https://www.miteco.gob.es/es/biodiversidad/temas/incendios-forestales/estadisticas-datos.html")
ncatt_put(ncout, 0,
          "Creator",
          "https://orcid.org/0000-0001-6219-4312")
ncatt_put(ncout, 0,
          "Created on", date())

# Place ba in the file
# need to loop through the layers to get them 
# to match to correct time index

for (i in 1:nlayers(raster_stack)) { 
    #message("Processing layer ", i, " of ", nlayers(prec))
    ncvar_put(nc = ncout, 
              varid = "ba", 
              vals = values(raster_stack[[i]]), 
              start = c(1, 1, i), 
              count = c(-1, -1, 1))
}
# Close the netcdf file when finished adding variables
nc_close(ncout)

message("[",
        Sys.time(),
        "] SUCCESS!\nNetCDF file written to ", filename)


q("no")

# # #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save(arr, file = "fireDatabase/raster_outputs/test_fires_cube_025deg.RData")
# # #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# # #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save(arr, file = "fireDatabase/raster_outputs/r025/fires_cube_025deg.RData")
# # #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



 # pdf("r_prueba.pdf", width = 8, height = 8)
 # plot(r) %>% print()
 # dev.off()

# mat <- arr[5,,]
# pdf("mat_prueba.pdf", width = 8, height = 8)
# raster::plot(r) %>% print()
# dev.off()
# 
# 
# ?transformeR::array3Dto2Dmat


