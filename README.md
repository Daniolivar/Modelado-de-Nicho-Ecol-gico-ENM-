# Modelado de Nicho Ecologico ENM
Modelado de Nicho Ecológico (ENM) del fitopatógeno Colletotrichum gloeosporioides en las Américas usando R y ENMeval. Implementación de un enfoque de fondo restringido por hospederos (Aguacate, Mango, Papaya, Fresa) para evitar sesgos geográficos y mejorar la precisión biológica.

🎯 Objetivo
Identificar zonas de alto riesgo climático para la antracnosis en las Américas, restringiendo el análisis a las áreas agrícolas donde sus hospederos principales están presentes, para evitar sesgos ecológicos triviales.

🛠️ Metodología 
Datos de Presencia: Descarga y limpieza de registros de GBIF para C. gloeosporioides (n = 1,404). Se implementó una limpieza de coordenadas personalizada para evitar conflictos de paquetes.

Definición del Área M (Fondo):

Se descartó el uso de un fondo global o continental.

Se construyó una máscara de "Hospederos Disponibles" basada en la distribución de Aguacate (Persea americana), Mango (Mangifera indica), Fresa (Fragaria) y Papaya (Carica papaya).

Variables Climáticas: Selección basada en biología y baja colinealidad:

bio10: Temperatura media del trimestre más cálido.

bio12: Precipitación anual.

bio15: Estacionalidad de la precipitación.

Modelado: Se utilizó el paquete ENMeval 2.0 con el algoritmo maxnet (Maxent sin Java). Se evaluaron 15 configuraciones de modelo (feature classes L, Q, H, LQ, LQH y regularización 1-3).

📊 Resultados Principales
Mejor Modelo: Configuración LQH con rm = 1.

Desempeño: AUC de Validación = 0.827. La diferencia entre AUC de entrenamiento y validación fue mínima (0.006), indicando ausencia de sobreajuste.

Test de Nulidad: El modelo es significativamente mejor que el azar (p < 0.01).





## 📦 1. Configuración del Entorno

Para garantizar la reproducibilidad del análisis, este flujo de trabajo utiliza un conjunto específico de librerías de R para la descarga de datos, el manejo espacial y el modelado ecológico.

Se utilizan paquetes como `maxnet` para evitar dependencias externas complejas (como Java) y `CoordinateCleaner` para asegurar la calidad de los datos biológicos.

```r
install.packages(c(
  "rgbif",              # Descarga de ocurrencias de GBIF
  "CoordinateCleaner",  # Limpieza automatizada de coordenadas
  "geodata",            # Descarga de clima (WorldClim)
  "terra",              # Manejo de datos raster y vectoriales
  "sf",                 # Operaciones espaciales simples
  "dplyr",              # Manipulación de datos
  "stringr",            # Manejo de cadenas de texto
  "ENMeval",            # Calibración y evaluación de modelos
  "maxnet"              # Algoritmo Maxent (versión R pura)
))



# --- 1. CARGAR LIBRERÍAS ----
library(rgbif)
library(CoordinateCleaner) # ¡Clave para la limpieza!
library(geodata)
library(terra)
library(sf)
library(dplyr)
library(stringr) #
```

1.1 Definición de Función de Descarga Robusta
Para asegurar la reproducibilidad y evitar errores comunes en la descarga de datos (como desconexiones o inconsistencias en las columnas entre continentes), definimos la función personalizada get_occurrences_america_V7.

Esta función implementa:

Descarga segura por continente (Norte y Sur América).

Manejo de errores (evita fallos si GBIF retorna NULL).

Unión inteligente de datos usando bind_rows.

Limpieza automatizada con CoordinateCleaner.

Conversión forzada a data.frame base para evitar conflictos con terra::extract en pasos posteriores.

Fragmento de código

# 2. Definir Función de Descarga Robusta (V7.0)
get_occurrences_america_V7 <- function(species_name, limit = 10000, clean_data = TRUE) {
  
  print(paste("Descargando:", species_name, "(Norteamérica)..."))
  data_na <- occ_search(
    scientificName = species_name,
    continent = "north_america",
    hasCoordinate = TRUE,
    limit = limit
  )$data
  
  print(paste("Descargando:", species_name, "(Sudamérica)..."))
  data_sa <- occ_search(
    scientificName = species_name,
    continent = "south_america",
    hasCoordinate = TRUE,
    limit = limit
  )$data
  
  if (is.null(data_na)) data_na <- data.frame()
  if (is.null(data_sa)) data_sa <- data.frame()

  print("Combinando datos...")
  df <- dplyr::bind_rows(data_na, data_sa)
  
  if (nrow(df) == 0) {
    print("No se encontraron registros.")
    return(data.frame(decimalLongitude=numeric(), decimalLatitude=numeric()))
  }

  target_cols <- c("decimalLongitude", "decimalLatitude", "species")
  available_cols <- intersect(target_cols, names(df))
  
  df <- df %>%
    # ¡AQUÍ ESTÁ LA PRIMERA CORRECCIÓN!
    dplyr::select(all_of(available_cols)) %>%
    dplyr::filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

  if (nrow(df) == 0) {
    print("No se encontraron registros tras el filtro inicial.")
    return(data.frame(decimalLongitude=numeric(), decimalLatitude=numeric()))
  }
  
  print(paste("Total de registros brutos:", nrow(df)))
  
  if (clean_data) {
    print("Limpiando coordenadas...")
    if (!"species" %in% names(df)) {
      df$species <- species_name
    }
    df_clean <- clean_coordinates(
      x = df,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "species",
      tests = c("centroids", "equal", "gbif", "institutions", "zeros")
    ) %>% filter(.summary == TRUE)
    
    print(paste("Total de registros limpios:", nrow(df_clean)))
    
    # --- ¡LA CORRECCIÓN CLAVE QUE ARREGLA TODO! ---
    # Forzamos la conversión a data.frame y usamos dplyr::select
    df_final <- as.data.frame(df_clean) %>%
      dplyr::select(decimalLongitude, decimalLatitude)
    # ---------------------------------------------
      
  } else {
    print("Omitiendo CoordinateCleaner.")
    df_final <- df %>%
      dplyr::select(decimalLongitude, decimalLatitude)
  }
  return(df_final)
}
1.2 Adquisición de Datos del Patógeno


