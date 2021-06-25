
#Script para cargar KML a Geemao y GEE
def imprimir(output):
    with output:
     print("hola")
    
def ingresar_poligono(nombrelote,ruta):

  ###Para ingresar un KML
  import geopandas as gpd
  import fiona
  # Enable fiona driver
  gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
  # Read file
  df = gpd.read_file(ruta, driver='KML')
  import os
  #os.chdir("/content/lotes")
  direccion = os.getcwd()
  rutafinal = direccion+"/salida1.geojson"
  df.to_file(rutafinal, driver='GeoJSON')

  f=open(rutafinal, "r")
  contents =f.read()
  geojson1 = contents.replace(", 0.0", "")
  #df.plot()
  df1 = df.to_crs("EPSG:3857")
  df1['dissolvefield'] = 1
  df1 = df1.dissolve(by='dissolvefield')
  df1 = df1["geometry"]
  return geojson1,df1,nombrelote

#script para hacer composiciones mensuales
def busqueda_imagenes(output,lote,año_inicio,año_fin,mes_inicio,mes_fin):
  #los valores pord defecto son 0.3,0.75 y 100
  import config
  import ee
  import json
  from pandas.io.json import json_normalize
  import pandas as pd
  import numpy as np
  import fiona
  import seaborn as sns
  import geopandas as gpd
  import matplotlib.pyplot as plt
  from IPython.display import HTML, display
  import time
  from pyproj import CRS

  
  #d = json.loads(geojson)
  #geometria = (d['features'][0]["geometry"]["coordinates"])
  #lote = ee.Geometry.Polygon(geometria[0])

    
  loteentero = lote
  lote = lote.buffer(-35)

  coleccionfiltrada = ee.ImageCollection('COPERNICUS/S2').filterBounds(lote).filter(ee.Filter.calendarRange(año_inicio,año_fin,'year'))
  
  lista = coleccionfiltrada.toList(coleccionfiltrada.size())
  imagen = ee.Image(lista.get(0))
  lista = lista.add(imagen)

  def detectar_duplicador(imagen):
        esduplicado = ee.String("")
        numero = lista.indexOf(imagen)
        imagen1 = ee.Image(lista.get(numero.add(1)))
        #Compare the image(0) in the ImageCollection with the image(1) in the List
        fecha1 = imagen.date().format("Y-M-d")
        fecha2 = imagen1.date().format("Y-M-d")
        estado = ee.Algorithms.IsEqual(fecha1,fecha2)
        esduplicado = ee.String(ee.Algorithms.If(estado,"duplicado","no duplicado"));
        imagen = imagen.set("duplicado", esduplicado)
        return imagen

  coleccionfiltrada = coleccionfiltrada.map(lambda image: detectar_duplicador(image))
  coleccionfiltrada = coleccionfiltrada.filter(ee.Filter.eq('duplicado', "no duplicado"))
  with output:
   print("Filtrado de imagenes con la misma fecha")



  coleccionfiltrada = coleccionfiltrada.map(lambda img: img.clip(loteentero))

  def agregar_nubes(image): 
      meanDict = image.reduceRegion(
      reducer= ee.Reducer.anyNonZero(),
      geometry= lote,
      scale= 10,
      )
      image = image.set("mascara",meanDict.get("QA60"))
      return image
  
  with output:
   print("Pasando por filtro de nubes")
  coleccionfiltrada = coleccionfiltrada.map(lambda image: agregar_nubes(image));
  coleccionfiltrada = coleccionfiltrada.filterMetadata('mascara', 'equals', 0);
  
  with output:
   print("NDVI para coleccion")

  def ndvi(img):
    ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI')
    img = img.addBands(ndvi)
    return img

  coleccionfiltrada = coleccionfiltrada.map(lambda img: ndvi(img))
    
  with output:
   print("Genero un mosaico mensual")
 
  #Genero un mosaico mensual
  #Mosaico mensual para NDVI. 
  # Modificado de https://gis.stackexchange.com/questions/258344/reduce-image-collection-to-get-annual-monthly-sum-precipitation
  #Se modificca para coincidir con los meses
 
    

  def sacar_meses(imagen):
    fecha = imagen.date().format("MM")
    imagen = imagen.set("meses",fecha)
    return imagen
   
  coleccion2017 = coleccionfiltrada.filter(ee.Filter.calendarRange(2017, 2017, 'year'))
  coleccion2018 = coleccionfiltrada.filter(ee.Filter.calendarRange(2018, 2018, 'year'))
  coleccion2019 = coleccionfiltrada.filter(ee.Filter.calendarRange(2019, 2019, 'year'))
  coleccion2020 = coleccionfiltrada.filter(ee.Filter.calendarRange(2020, 2020, 'year'))

  coleccion2017 = coleccion2017.map(lambda imagen:sacar_meses(imagen))
  coleccion2018 = coleccion2018.map(lambda imagen:sacar_meses(imagen))
  coleccion2019 = coleccion2019.map(lambda imagen:sacar_meses(imagen))
  coleccion2020 = coleccion2020.map(lambda imagen:sacar_meses(imagen))

  meses2017 = coleccion2017.aggregate_array("meses").distinct()   
  meses2018 = coleccion2018.aggregate_array("meses").distinct()   
  meses2019 = coleccion2019.aggregate_array("meses").distinct()   
  meses2020 = coleccion2020.aggregate_array("meses").distinct()   
    
  def convertir_string(valor):
    return ee.Number.parse(valor)

  meses2017 = meses2017.map(lambda valor:convertir_string(valor))
  meses2018 = meses2018.map(lambda valor:convertir_string(valor))
  meses2019 = meses2019.map(lambda valor:convertir_string(valor))
  meses2020 = meses2020.map(lambda valor:convertir_string(valor))

    
    
    

  def mosaico_mensual_2017(m):
      coleccion = coleccion2017.filter(ee.Filter.calendarRange(m, m, 'month')).select(['NDVI','B4', 'B3', 'B2'])
      fecha = coleccion.first().date().format("YYYY-MM")
      año = coleccion.first().date().format("YYYY")
      imagen = coleccion.mean()
      imagen = imagen.set('month', m)  
      #fecha = imagen.select("NDVI").date().format("YYYY-MM")
      imagen = imagen.set("system:time_start", fecha) 
      imagen = imagen.set("year", año) 
      return imagen
    
  def mosaico_mensual_2018(m):
      coleccion = coleccion2018.filter(ee.Filter.calendarRange(m, m, 'month')).select(['NDVI','B4', 'B3', 'B2'])
      fecha = coleccion.first().date().format("YYYY-MM")
      año = coleccion.first().date().format("YYYY")
      imagen = coleccion.mean()
      imagen = imagen.set('month', m)  
      #fecha = imagen.select("NDVI").date().format("YYYY-MM")
      imagen = imagen.set("system:time_start", fecha) 
      imagen = imagen.set("year", año) 
      return imagen

  def mosaico_mensual_2019(m):
      coleccion = coleccion2019.filter(ee.Filter.calendarRange(m, m, 'month')).select(['NDVI','B4', 'B3', 'B2'])
      fecha = coleccion.first().date().format("YYYY-MM")
      año = coleccion.first().date().format("YYYY")
      imagen = coleccion.mean()
      imagen = imagen.set('month', m)  
      #fecha = imagen.select("NDVI").date().format("YYYY-MM")
      imagen = imagen.set("system:time_start", fecha) 
      imagen = imagen.set("year", año) 
      return imagen
    
  def mosaico_mensual_2020(m):
      coleccion = coleccion2020.filter(ee.Filter.calendarRange(m, m, 'month')).select(['NDVI','B4', 'B3', 'B2'])
      fecha = coleccion.first().date().format("YYYY-MM")
      año = coleccion.first().date().format("YYYY")
      imagen = coleccion.mean()
      imagen = imagen.set('month', m)  
      #fecha = imagen.select("NDVI").date().format("YYYY-MM")
      imagen = imagen.set("system:time_start", fecha) 
      imagen = imagen.set("year", año) 
      return imagen

  
  coleccionfiltrada2017 = ee.ImageCollection.fromImages(meses2017.map(lambda m: mosaico_mensual_2017(m)))
  coleccionfiltrada2018 = ee.ImageCollection.fromImages(meses2018.map(lambda m: mosaico_mensual_2018(m)))
  coleccionfiltrada2019 = ee.ImageCollection.fromImages(meses2019.map(lambda m: mosaico_mensual_2019(m)))
  coleccionfiltrada2020 = ee.ImageCollection.fromImages(meses2020.map(lambda m: mosaico_mensual_2020(m)))
                                                     
  coleccionfiltrada = coleccionfiltrada2017.merge(coleccionfiltrada2018)         
  coleccionfiltrada = coleccionfiltrada.merge(coleccionfiltrada2019)
  coleccionfiltrada = coleccionfiltrada.merge(coleccionfiltrada2020)
  
  with output:
   print("Calculo de estadisticas")

  def ndvi_medio(image):
      image1 = image.select("NDVI").rename("NDVI_medio")
      reduced = image1.reduceRegion(geometry=lote, reducer=ee.Reducer.mean(), scale=10)
      image = image.set(reduced)
      return image
  def ndvi_sd(image):
      image1 = image.select("NDVI").rename("NDVI_sd")
      reduced = image1.reduceRegion(geometry=lote, reducer=ee.Reducer.stdDev(), scale=10)
      image = image.set(reduced)
      return image
  def normalidad(image):
      image1 = image.select("NDVI").rename("normalidad")
      reduced = ee.Number(ee.Number(image1.get("NDVI_medio")).divide(ee.Number(image1.get("NDVI_mediana"))))
      image = image.set("normalidad",reduced)
      return image
  def ndvi_cv(image):
      sd = ee.Number(image.get("NDVI_sd"))
      medio = ee.Number(image.get("NDVI_medio"))
      cv = sd.divide(medio).multiply(100)
      image = image.set("cv",cv)
      return image
  

  coleccionfiltrada = coleccionfiltrada.map(lambda imagen: ndvi_medio(imagen))
  coleccionfiltrada = coleccionfiltrada.map(lambda imagen: ndvi_sd(imagen))
  coleccionfiltrada = coleccionfiltrada.map(lambda imagen: ndvi_cv(imagen))
  
  with output:
   print("Armando diccionario y descargando datos del servidor")

  fechas = coleccionfiltrada.aggregate_array("system:time_start")
  años = coleccionfiltrada.aggregate_array("year")

  NDVI_medio = coleccionfiltrada.aggregate_array("NDVI_medio")
  NDVI_sd = coleccionfiltrada.aggregate_array("NDVI_sd")
  meses = coleccionfiltrada.aggregate_array("month")
    
  test_dict = ee.Dictionary.fromLists(['system:time_start', 'NDVI_medio','NDVI_sd','años','mes'], [fechas, NDVI_medio,NDVI_sd,años,meses])
  featureCollection = ee.FeatureCollection([ee.Feature(None, test_dict)])
  
    
 
    
    
  link = featureCollection.getDownloadURL(filetype="CSV", selectors=None, filename=None)
  #Generamos el diccionario con los valores.
  import csv, urllib.request
  response = urllib.request.urlopen(link)
  lines = [l.decode('utf-8') for l in response.readlines()]
  reader = csv.DictReader(lines)

  with output:
   print("Covirtiendo a dataframe")

  data = list(reader)
  # Converting string to list
  fechas_cliente = data[0]["system:time_start"].strip('][').split(', ')
  ndvi_cliente = data[0]["NDVI_medio"].strip('][').split(', ')
  sd_ndvi_cliente = data[0]["NDVI_sd"].strip('][').split(', ')
  años_cliente = data[0]["años"].strip('][').split(', ')
  meses_cliente = data[0]["mes"].strip('][').split(', ')
  
  fechas_cliente = map(str, fechas_cliente)
  fechas = list(fechas_cliente)
  
  ndvi_cliente = map(float, ndvi_cliente)
  sd_ndvi_cliente = map(float, sd_ndvi_cliente)
  años_cliente = map(float, años_cliente)
  meses_cliente = map(float, meses_cliente)

  #Genero dataframe con los datos
  df = pd.DataFrame(list(zip(fechas, ndvi_cliente,sd_ndvi_cliente,años_cliente,meses_cliente)),columns =['Fecha', 'NDVI_medio',"sd_ndvi","años","meses"])
    
    
  config.coleccion = coleccionfiltrada
  config.lote = loteentero
  config.df = df

  with output:
   print("Exito!")
  return df, coleccionfiltrada,loteentero


#Graficamos lo que vemos
def graficar_series(output1,datos_lote):
    from matplotlib import pyplot as plt
    import numpy as np
    import pandas as pd

    #Si algun mes no tiene datos lo reemplazo or NULL
    for a in range(2017,2021):
     datos = datos_lote[datos_lote['años'] == a]
     for i in range(1,12):
      condicion = i in list(datos["meses"])
      if condicion is False:
       fecha =  (str(a)+"-"+str(i))
       df = pd.DataFrame([[fecha, None,None,a,i]],columns =['Fecha', 'NDVI_medio',"sd_ndvi","años","meses"])
       datos_lote = datos_lote.append(df)


    datos2017 = datos_lote[datos_lote['años'] == 2017].sort_values(by=['meses'])
    datos2018 = datos_lote[datos_lote['años'] == 2018].sort_values(by=['meses'])
    datos2019 = datos_lote[datos_lote['años'] == 2019].sort_values(by=['meses'])
    datos2020 = datos_lote[datos_lote['años'] == 2020].sort_values(by=['meses'])
    fechas= pd.Series([1,2,3,4,5,6,7,8,9,10,11,12])




    with output1:
        #fig, axs = plt.subplots(2,1)
        fig, ax = plt.subplots(1,1,figsize=(15,10))

        #2017
        transparencia=0.4
        colorlinea = '#CC4F1B'
        edgecolor= '#CC4F1B'
        facecolor= '#CC4F1B'
        ndvi =  datos2017["NDVI_medio"]
        sdmas =  datos2017["NDVI_medio"]+datos2017["sd_ndvi"]
        sdmenos =  datos2017["NDVI_medio"]-datos2017["sd_ndvi"]

        ax.plot(fechas, ndvi,'k-',color=colorlinea,label=2017)
        ax.fill_between(fechas,sdmenos,sdmas, alpha=transparencia, edgecolor=edgecolor, facecolor=facecolor,
            linewidth=4, antialiased=True)

        #2018
        transparencia=0.4
        colorlinea = '#1B2ACC'
        edgecolor= '#1B2ACC'
        facecolor= '#1B2ACC'
        ndvi =  datos2018["NDVI_medio"]
        sdmas =  datos2018["NDVI_medio"]+datos2018["sd_ndvi"]
        sdmenos =  datos2018["NDVI_medio"]-datos2018["sd_ndvi"]

        ax.plot(fechas, ndvi,'k-',color=colorlinea,label=2018)
        ax.fill_between(fechas,sdmenos,sdmas, alpha=transparencia, edgecolor=edgecolor, facecolor=facecolor,
            linewidth=4, antialiased=True)

        #2019
        transparencia=0.4
        colorlinea = '#3F7F4C'
        edgecolor= '#3F7F4C'
        facecolor= '#3F7F4C'
        ndvi =  datos2019["NDVI_medio"]
        sdmas =  datos2019["NDVI_medio"]+datos2019["sd_ndvi"]
        sdmenos =  datos2019["NDVI_medio"]-datos2019["sd_ndvi"]

        ax.plot(fechas, ndvi,'k-',color=colorlinea,label=2019)
        ax.fill_between(fechas,sdmenos,sdmas, alpha=transparencia, edgecolor=edgecolor, facecolor=facecolor,
            linewidth=4, antialiased=True)

        #2020
        transparencia=0.4
        colorlinea = '#F6FF33'
        edgecolor= '#F6FF33'
        facecolor= '#F6FF33'
        ndvi =  datos2020["NDVI_medio"]
        sdmas =  datos2020["NDVI_medio"]+datos2020["sd_ndvi"]
        sdmenos =  datos2020["NDVI_medio"]-datos2020["sd_ndvi"]

        ax.plot(fechas, ndvi,'k-',color=colorlinea,label=2020)
        ax.fill_between(fechas,sdmenos,sdmas, alpha=transparencia, edgecolor=edgecolor, facecolor=facecolor,
            linewidth=4, antialiased=True)

        ax.set_title('Evolucion de NDVI en KML cargado')
        #ax.set(xlabel="MES",ylabel="NDVI",fontsize=18)
        ax.set_ylabel('NDVI', fontsize = 20.0) # Y label
        ax.set_xlabel('MES', fontsize = 20) # X label
        ax.legend(prop=dict(size=18))

        ax.tick_params(labelsize=15)
    return
