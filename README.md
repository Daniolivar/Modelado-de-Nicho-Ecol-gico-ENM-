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







## 🛠️ 2. Función de Descarga y Limpieza Robusta
Para solucionar los problemas de desconexión y formatos de objetos incompatibles (el error spatialvalid), definimos una función personalizada get_occurrences_america_V7. Esta función:

Descarga datos por continente para evitar tiempos de espera.

Maneja respuestas vacías (NULL) sin romper el script.

Aplica la limpieza de coordenadas y fuerza el resultado a un data.frame simple compatible con ENMeval.


```r
## --- 1. FUNCIÓN DE DESCARGA (V7.0 - CORREGIDA) ---
# (Esta función maneja 'NULL', 'bind_rows', Y el error 'spatialvalid')

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
```
##🦠 3. Adquisición de Datos del Patógeno
Descargamos los registros de Colletotrichum gloeosporioides utilizando la función robusta, asegurando que las presencias (occs) estén limpias y listas para el modelado.

```r
## --- 2. PREPARAR PRESENCIAS (Occs) - ¡CON LIMPIEZA! ---
print("Descargando patógeno (CON CoordinateCleaner)...")
# Usamos V7 con clean_data = TRUE
occs <- get_occurrences_america_V7("Colletotrichum gloeosporioides", clean_data = TRUE)
```
##🥑 4. Adquisición de Datos de Hospederos
Para restringir el espacio de fondo ("M") a zonas biológicamente relevantes, descargamos los registros de los principales hospederos productivos: Aguacate, Mango, Fresa y Papaya. Estos puntos se combinan en un solo objeto vectorial.

```r
## --- 3. PREPARAR HOSPEDEROS (M) - ¡CON LIMPIEZA! ---
print("Descargando hospederos (CON CoordinateCleaner)...")
occ_aguacate <- get_occurrences_america_V7("Persea americana", clean_data = TRUE)
occ_mango <- get_occurrences_america_V7("Mangifera indica", clean_data = TRUE)
occ_fresa <- get_occurrences_america_V7("Fragaria", clean_data = TRUE)
occ_papaya <- get_occurrences_america_V7("Carica papaya", clean_data = TRUE)

puntos_hospederos_df <- dplyr::bind_rows(occ_aguacate, occ_mango, occ_fresa, occ_papaya) %>%
  dplyr::filter(!is.na(decimalLongitude), !is.na(decimalLatitude))
puntos_hospederos_vect <- terra::vect(puntos_hospederos_df, 
                                      geom = c("decimalLongitude", "decimalLatitude"), 
                                      crs = "EPSG:4326")
print("¡Hospederos listos!")

```

##🌦️ 5. Preparación de Variables Climáticas
Descargamos los datos de WorldClim 2.1 y realizamos una selección de variables a priori para evitar la multicolinealidad. Se seleccionaron bio10, bio12 y bio15 por su relevancia fisiológica para el desarrollo fúngico (calor y humedad).

```r
## --- 4. PREPARAR CAPAS AMBIENTALES (env) - NUEVAS VARIABLES ---
print("Preparando capas ambientales...")
bio_global <- geodata::worldclim_global(var = "bio", res = 10, path = ".")
names(bio_global) <- paste0("bio", 1:19)
americas_extent <- ext(-170, -30, -55, 75)
env_americas <- terra::crop(bio_global, americas_extent)

# CAMBIO: Ahora seleccionas bio10, bio12, bio15
vars_seleccionadas <- c("bio10", "bio12", "bio15")
env <- env_americas[[vars_seleccionadas]]

print("¡Capas 'env' listas!")
print(env)

```


##🗺️ 6. Construcción del Fondo (Background M)
Generamos la máscara de fondo creando un buffer de 100 km alrededor de los cultivos. Luego, muestreamos 10,000 puntos aleatorios (bg_points) exclusivamente dentro de esta zona, evitando sesgos por comparar con climas extremos no agrícolas.


```r
## --- 5. CREAR MÁSCARA "M" y FONDO "BG" ---
print("Creando 'M' y 'BG'...")
m_hospederos <- terra::buffer(puntos_hospederos_vect, width = 100000) %>%
  terra::rasterize(env[[1]], background = NA) %>%
  terra::subst(from = 1, to = 1)

bg_points <- terra::spatSample(x = m_hospederos, size = 10000, method = "random", as.df = TRUE, xy = TRUE) %>%
  dplyr::rename(decimalLongitude = x, decimalLatitude = y)
print("¡'M' y 'BG' listos!")
```



##⚙️ 7. Ejecución del Modelo ENMeval
Ejecutamos la evaluación de modelos utilizando maxnet. Se prueban múltiples configuraciones de complejidad (Lineal, Cuadrática, Hinge) y regularización para encontrar el modelo óptimo, utilizando validación cruzada (k-fold) para evitar el sobreajuste.

```r
## --- 6. CORRER ENMEVAL ---
cat("\n=== VERIFICACIÓN PRE-VUELO ===\n")
stopifnot(is.data.frame(occs), is.data.frame(bg_points), ncol(occs)==2)

cat("\n=== EJECUTANDO ENMevaluate ===\n")
eval_results <- ENMeval::ENMevaluate(
  occs = occs,
  env = env,
  bg = bg_points,
  algorithm = 'maxnet',
  partitions = 'randomkfold',
  partition.settings = list(kfolds = 5),
  tune.args = list(
    fc = c("L", "Q", "H", "LQ", "LQH"),
    rm = c(1, 2, 3)
  ),
  quiet = FALSE
)

cat("\n=== ¡EVALUACIÓN COMPLETADA! ===\n")
print(eval_results)
```


##📊 8. Selección del Mejor Modelo y Predicción
Analizamos la tabla de resultados para seleccionar el modelo con el menor AICc (Criterio de Información de Akaike). Generamos el mapa de idoneidad final y lo exportamos como un archivo raster GeoTIFF.

```r
## --- VER RESULTADOS Y PREDECIR ---
results_table <- eval_results@results

# Identificar el mejor modelo (menor AICc)
best_idx <- which.min(results_table$AICc)
cat("\n=== MEJOR MODELO (menor AICc) ===\n")
print(results_table[best_idx, ])

# Extraer el mejor modelo
best_model <- eval_results@models[[best_idx]]

# Hacer predicción CON manejo de NAs
prediction <- predict(env, best_model, type = "cloglog", na.rm = TRUE)

# Visualizar
library(terra)
colores <- colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100)
plot(prediction, 
     main = "Idoneidad de hábitat - C. gloeosporioides",
     col = colores)
points(occs$decimalLongitude, occs$decimalLatitude, pch = 20, cex = 0.3, col = "blue")

# Guardar el raster
terra::writeRaster(prediction, "mapa_idoneidad_colletotrichum.tif", overwrite = TRUE)
cat("✓ Mapa guardado!\n")

```

<img width="831" height="586" alt="image" src="https://github.com/user-attachments/assets/4131c0cd-c1f8-45b9-823e-3fb41b9fed10" />



