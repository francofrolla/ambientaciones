
## 5. MULTIVARIATE SITE CLASSIFICATION ---------------------------------------------------

# After all the variables have been processed from step 1 up to the interpolation with the same prediction grid (step 5),
# the different data sets obtained should be concatenated using the cbind function.
# Below the R code was deactivated because a new database that has been concatenated is used.

# Pred <- cbind(PredECa30[,1:3], PredECa90[,3], PredElev[,3],PredSd[,3])
# names(Pred)[3]<-paste("ECa30")
# names(Pred)[4]<-paste("ECa90")
# names(Pred)[5]<-paste("Elev")
# names(Pred)[6]<-paste("Sd")

require(raster)
require(tcltk)
require(rgeos)

ruta <-tk_choose.dir()
setwd(ruta)
dir(ruta)

bar <- function() {print('ATENCION: Hubo un error que detuvo la ejecucion del programa')}


ambientes<-multispatipca() 

mapa_vectorizado<-tryCatch(vectorizacion(ambientes,3,0.5,poligono,5),  error = function(e){ bar()} )


vectorizacion<-function(ambientes,n_ambientes,superficie_minima,poligono,celda){
  
  #Controlo si se ejecuto multispatipca
  
  if(!exists("ambientes")){stop("Primero debe ejecutar multispatipca")}
  
  #Controlo valor logico de ventana de suavizado
  
  if(celda != 3 && celda != 5 && celda != 9 && celda != 13){
    stop("No se pueden hacer ese suavizado, valores posibles 3,5,9,13")
  }
  
  
  if(n_ambientes == 2){
    rasterambientes<- raster(ambientes[[1]],layer=1)
  }
  if(n_ambientes == 3){
    rasterambientes<- raster(ambientes[[1]],layer=2)
  }
  if(n_ambientes == 4){
    rasterambientes<- raster(ambientes[[1]],layer=3)
  }
  if(n_ambientes == 5){
    rasterambientes<- raster(ambientes[[1]],layer=4)
  }
  if(n_ambientes == 6){
    rasterambientes<- raster(ambientes[[1]],layer=5)
  }
  if(n_ambientes == 7){
    rasterambientes<- raster(ambientes[[1]],layer=6)
  }
  if(n_ambientes != 2 && n_ambientes != 3 && n_ambientes != 4 && n_ambientes != 5 && n_ambientes != 6 && n_ambientes != 7){
        stop("No se pueden hacer ese cantidad de ambientes, valores posibles 2,3,4,5,6,7")
  }
  
  #Existe la superficie minima
  
  if(!exists("superficie_minima")){superficie_minima<-0.5}
   
  ambientesvectorizados<-rasterToPolygons(rasterambientes,na.rm=TRUE,dissolve=TRUE)
  ambientesvectorizados<-disaggregate (ambientesvectorizados)
  plot(ambientesvectorizados)
  
  
  plot(ambientesvectorizados["MC3"],col=viridis(4),border="NA")
  
  
  #Remuevo aquella areas menores a 1/2 ha
  ambientesvectorizados@data$area<-area(ambientesvectorizados)/10000
  ambientesvectorizados<-ambientesvectorizados[ambientesvectorizados@data$area>superficie_minima,]
  proj4string(ambientesvectorizados) <- crs(poligono)
  

  #Convierto a raster
 
  r.polys <- rasterize(ambientesvectorizados[,1], ambientes[[2]],field= names(ambientesvectorizados@data)[1])
  
  #Remuestreo la capa para ver donde nos quedamos sin datos.
  
  #Grilla para muestreas capas
  #Genero grilla para hacer los clusters
  
  
  grid <- makegrid(buffer(poligono,50), cellsize = 5)
  coordinates(grid)<-c("x1","x2")
  proj4string(grid) <- crs(poligono)
  recortegrilla <- raster::crop(grid,buffer(poligono,50))
  
  data<-try(extract(r.polys,recortegrilla),silent=TRUE)
  
  #Remuevo outliers
  datos<-data.frame(data,recortegrilla@coords[,1],recortegrilla@coords[,2])
  colnames(datos)<-c("Ambientes","x","y")
  #datos<-na.omit(datos)
  
  #teniendo los datos a una escala de muestreo inferior a la original, rasterizo, suavizo
  #y por ultimo vectorizo.
  
  amb_raster <- SpatialPixelsDataFrame(points = datos[c("x", "y")], data =datos)
  amb_raster<- raster(amb_raster,layer=1)
  
  
  #Genero la ventana de suavizado
  ventana<-matrix(1,nrow=celda,ncol=celda)
  
  mapa<-focal(amb_raster[[1]],w=ventana,fun=modal,na.rm=TRUE)
  
  #Fuerzo rellenar nulls hasta ventana maxima de 13x13
  ventana<-matrix(1,nrow=13,ncol=13)
  mapa<-focal(mapa,w=ventana,fun=modal,na.rm=TRUE,NAonly=TRUE)
  
  
  #mapa<-crop(mapa,poligono)
  #mapa<-mask(mapa,poligono)
  

  mapa_vectorizado<-rasterToPolygons(mapa,na.rm=TRUE,dissolve=TRUE)
  
  
  gIsValid(mapa_vectorizado, reason = T)
  mapa_vectorizado <- gBuffer(mapa_vectorizado, width=0, byid = T)
  print("generando buffer para solucionar problemas de geometria")
  gIsValid(mapa_vectorizado, reason = T)
  
  crs(mapa_vectorizado)<-crs(poligono)
  
  mapa_vectorizado<-raster::crop(mapa_vectorizado,buffer(poligono,20))
  
  spplot(mapa_vectorizado,col.regions=viridis(100))
  
  
  return(mapa_vectorizado)
  }
  
writeOGR(mapa_vectorizado, layer = "mapa_vectorizado_4ambientes_2018", dsn=getwd(), driver="ESRI Shapefile",overwrite_layer=TRUE)

multispatipca<-function(){
  #Lllamar a todas las capas que se van a usar....
  require(raster)
  require(rgdal)
  require(ade4)
  require(spdep)
  require(adespatial)
  require(e1071)
  
  setwd(tk_choose.dir())
  
  cea30<-raster("2018.3.tif")
  cea90<-raster("2018.4.tif")
  dem<-raster("2018.5.tif")
  ndvi<-raster("2018.6.tif")

  #radarvh<-raster("RadarVV.tif")
  #twi <-raster("TWI.tif")
  #conv20<-raster("convergencia 20m.tif")
  #conv10<-raster("convergencia 10m.tif")
  #twimdear<-raster("twi mdear.tif")
  #profvalle<-raster("profvalle.tif")
  listacapas<-list(cea30,cea90,dem,ndvi)
  
  #Reproyecto a 3857
  
  for (i in (1:5)){
    newproj <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null
    +wktext +no_def"
    
    listacapas[[i]]<- projectRaster(listacapas[[i]],crs = newproj)
    listacapas[[i]][listacapas[[i]] <= 0]<-NA
  }
  
  
  
  #listacapas<-list(cea30,cea90,dem,ndvi,radarvv,radarvh,twi)
  #listacapas<-list(ndvi,radarvv,radarvh)
  
  newproj <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null
    +wktext +no_def"
  poligono<-readOGR("C:/Users/francofrolla/Downloads/LaMadrir/LaMadrir/Establecimiento_1/03_DatosdeCampo/Lote1/Lote1_Perimetro.shp")
  
  poligono<- spTransform(poligono,crs(newproj))
  
  #FUNCION DE CONTROL, controlar el SRC de la primera capas con las restantes.

  
  for (i in 1:length(listacapas)){
    crs_control<-projection(listacapas[[1]])
    if(crs_control != projection(poligono)){print("el poligono y la capa raster no son iguales");stop("CRS: diferente crs para poligonos y raster, por favor reproyecte las capas a un mismo sistema en metros")}
    if(crs_control == "+proj=longlat +datum=WGS84 +no_defs"){print("una de las esta en grados");stop("CRS: Error en proyeccion capa en grados, por favor reproyecte las capas a un mismo sistema en metros")}
    if(crs_control != projection(listacapas[[i]])){print("una de las capas tiene diferente SRC");stop("CRS: Error en proyeccion en una de las capas, por favor reproyecte las capas a un mismo sistema en metros")}
     
  }
  
  #tomo la primer capa y la resampleo
  for (i in 2:length(listacapas)){
    listacapas[[i]]<-resample(listacapas[[i]],listacapas[[1]])
  }
  
  #Grilla para muestreas capas
  #Genero grilla para hacer los clusters
  grid <- makegrid(buffer(poligono,50), cellsize = 10)
  coordinates(grid)<-c("x1","x2")
  proj4string(grid) <- crs(poligono)
  recortegrilla <- raster::crop(grid,buffer(poligono,50))
  
  
  
  #Extraigo los valores de las capas y genero data.frame
  
  nombres<-c("")
  matrizdatos<-seq(1,length(recortegrilla),1)
  
  for (i in 1:length(listacapas)){
    print(i)
    data<-extract(listacapas[[i]],recortegrilla)
    matrizdatos<-cbind(matrizdatos,data)
    nombre<-listacapas[[i]]@data@names
    nombres<-c(nombres,nombre)
    
  }
  
  
  datos<-data.frame(matrizdatos,recortegrilla@coords[,1],recortegrilla@coords[,2])
  colnames(datos)<-c(nombres,"x","y")
  
  datos<-na.omit(datos)
  
  #controlo que halla datos
  if(nrow(datos)<1){stop("No se generaron datos en la grilla")}
  
  #multispati-pca
  #datos sin x e y
  datospca<-datos[,1:(ncol(datos)-2)]
  datospca<-datospca[,2:ncol(datospca)]

  pca <- dudi.pca(datospca, center=T,scannf = FALSE,  nf = 5)
  
  
  cord_1 <- coordinates(datos[,(ncol(datos)-1):(ncol(datos))])
  gri_1 <- dnearneigh(cord_1,0,25)
  lw_1 <- nb2listw(gri_1, style = "W")
  
  ms <- multispati(pca, lw_1, scannf = F, nfposi = 5)
  s.arrow(ms$c1,xax = 1, yax = 2, clabel = 1)
  
  # Extraction of spatial principal components
  sPC <- ms$li[,1:2]
  PredMA <- cbind(datos,sPC) 
  PredMA
  
  #
  parazonificar<-data.frame(PredMA$CS1,PredMA$CS2)
  #  Fuzzy k-means cluster analysis
  MC_2<-cmeans(parazonificar,2,100,method="cmeans",m=1.3)
  MC_3<-cmeans(parazonificar,3,100,method="cmeans",m=1.3)
  MC_4<-cmeans(parazonificar,4,100,method="cmeans",m=1.3)
  MC_5<-cmeans(parazonificar,7,100,method="cmeans",m=1.3)
  MC_6<-cmeans(parazonificar,7,100,method="cmeans",m=1.3)
  MC_7<-cmeans(parazonificar,7,100,method="cmeans",m=1.3)
  
  # Mapas con clases delimitadas
  MC_22 <-as.data.frame(MC_2$cluster)
  MC_33 <-as.data.frame(MC_3$cluster)
  MC_44 <-as.data.frame(MC_4$cluster)
  MC_55 <-as.data.frame(MC_5$cluster)
  MC_66 <-as.data.frame(MC_6$cluster)
  MC_77 <-as.data.frame(MC_7$cluster)
  
  
  baseMC <- cbind(datos$x,datos$y,MC_22,MC_33,MC_44,MC_55,MC_66,MC_77)
  colnames(baseMC)<-c("x","y","MC2","MC3","MC4","MC5","MC6","MC7")
  
  coordinates(baseMC) <- ~x+y
  gridded(baseMC) <- T
  
  library(viridis)
  
  spplot(baseMC["MC2"],col.regions=viridis(100),colorkey = F)
  
  spplot(baseMC["MC3"],col.regions=viridis(100),colorkey = F)
  
  spplot(baseMC["MC4"],col.regions=viridis(100),colorkey = F)
  
  spplot(baseMC["MC5"],col.regions=viridis(100),colorkey = F)
  
  spplot(baseMC["MC6"],col.regions=viridis(100),colorkey = F)
  
  
  spplot(baseMC["MC7"],col.regions=viridis(100),colorkey = F)
  
  
  return(list(baseMC,listacapas[[1]],poligono))
  
  
}

#Controlo el crs de la

