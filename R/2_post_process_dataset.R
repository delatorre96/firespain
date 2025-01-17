rm(list=ls())
gc()

library(magrittr)
library(raster)

setwd("/lustre/gmeteo/WORK/juaco/PTI-clima")

load("fireDatabase/raster_outputs/test_fires_cube_025deg.RData", verbose = TRUE)

library(transformeR)
arrmat <- array3Dto2Dmat(arr)
str(arrmat)


x <- coordinates(raster_base)[,1] %>%
    unique() %>%
    sort(decreasing = FALSE)
y <- coordinates(raster_base)[,2] %>%
    unique() %>%
    sort(decreasing = FALSE)

arr3d <- mat2Dto3Darray(arrmat, x, y)

str(arr3d)
arr3d[5,,] 
# load("fireDatabase/raster_outputs/r025/fires_cube_025deg.RData", verbose = TRUE)
# clim <- apply(arr, MAR = c(2,3), FUN = "mean", na.rm=TRUE)
dim(arr)
str(arr)
# attr(arr, "dimensions") <- c("time", "lon", "lat")

pdf("test.pdf", width = 8, height = 8)
a <- arr[1,,] %>% t() 
a[,34:1] %>% image()
dev.off()

## Inicializa el objeto
ba <- list()

## Variable

Variable <- list()
Variable$varName <- "ba"
Variable$level <- "surface"
attr(Variable, "units") <- "ha"
attr(Variable, "longname") <- "Burned area"
attr(Variable, "daily_agg_cellfun") <- "sum"
attr(Variable, "verification_time") <- "DD"
attr(Variable, "description") <- "Total burned area in ha"
ba$Variable <- Variable

## Data cube
# ba$Data <- arr
ba$Data <- arr3d

## Spatial coordinates
load("fireDatabase/carto/raster_base_0.25deg_landmasked.RData", verbose = TRUE)

xyCoords <- list()
xyCoords$x <- coordinates(raster_base)[,1] %>%
    unique() %>%
    sort(decreasing = FALSE)
xyCoords$y <- coordinates(raster_base)[,2] %>%
    unique() %>%
    sort(decreasing = FALSE)

attr(xyCoords, "resX") <- res(raster_base)[1]
attr(xyCoords, "resY") <- res(raster_base)[2]
ba$xyCoords <- xyCoords
attr(ba$xyCoords, "proj4string") <- "+proj=longlat +datum=WGS84"
    
## Dates
fecha_minima <- as.Date("1983-01-01")
fecha_maxima <- as.Date("1983-03-31")
# fecha_maxima <- as.Date("2021-12-31")
rango_fechas <- seq(from = fecha_minima, to = fecha_maxima, by = "days")
Dates <- list()
Dates$start <- rango_fechas
Dates$end <- rango_fechas + 1
ba$Dates <- Dates

## Atributos globales
attr(ba, "source") <- "Estadística General de Incendios Forestales (EGIF), MITECO"
attr(ba, "description") <- "Base de datos de incendios forestales en España"
attr(ba, "creator") <- "J. Bedia <https://orcid.org/0000-0001-6219-4312>"

str(ba)


library(transformeR)
library(visualizeR)

pdf("fireDatabase/ba_025deg_test.pdf", width = 10, height = 10)
climatology(ba) %>% spatialPlot(backdrop.theme = "coastline", rev.colors = TRUE) %>% print()
dev.off()


data("NCEP_Iberia_ta850")
str(NCEP_Iberia_ta850)


str(ba)

## Atributos globales 
attr(ba, "source") <- "Estadística General de Incendios Forestales (EGIF)"
    

object.size(raster_stack) %>% print(units = "Mb")


writeRaster(raster_stack, "fireDatabase/rstack_test.nc", 
            format = "CDF", overwrite = TRUE, zname = "time")

library(loadeR)
di <- dataInventory("fireDatabase/rstack_test.nc")





library(magrittr)
library(loadeR)
library(transformeR)
library(visualizeR)

ds <- "fireDatabase/raster_outputs/r025/fire_025deg_v0.1.nc"
di <- dataInventory(ds)
str(di)

ba <- loadGridData(ds, var = "ba") 
ba.an <- aggregateGrid(ba, aggr.y = list(FUN = "sum", na.rm = TRUE))

spatialPlot(climatology(ba.an), backdrop.theme = "countries", rev.colors = TRUE,
            main = "Mean annual Total Burned Area (ha) 1983-2001",
            at = seq(0,2000,50), set.max = 2000) 

mysum <- function(x) {
    if (all(is.na(x))) {
        NA
    }
    else {
        x <- na.omit(x)
        which(x > 0) %>% length()    
    }
}

count <- function(x) {
    if (all(is.na(x))) {
        NA
    }
    else {
        x <- na.omit(x)
        which(x > 0) %>% length()    
    }
}

count <- function(x) {
        x <- na.omit(x)
        which(x > 0) %>% length()    
    }
}

ba.count <- climatology(ba, clim.fun = list(FUN = "count"))

spatialPlot(ba.count, backdrop.theme = "countries", 
            color.theme = "YlOrRd",
            set.max = 3500,
            at = seq(0,3500,100),
            main = "Total number of fire events 1983-2021") %>% print()











