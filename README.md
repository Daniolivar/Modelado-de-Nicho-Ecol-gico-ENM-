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
## 🦠 3. Adquisición de Datos del Patógeno
Descargamos los registros de Colletotrichum gloeosporioides utilizando la función robusta, asegurando que las presencias (occs) estén limpias y listas para el modelado.

```r
## --- 2. PREPARAR PRESENCIAS (Occs) - ¡CON LIMPIEZA! ---
print("Descargando patógeno (CON CoordinateCleaner)...")
# Usamos V7 con clean_data = TRUE
occs <- get_occurrences_america_V7("Colletotrichum gloeosporioides", clean_data = TRUE)
```
## 🥑 4. Adquisición de Datos de Hospederos

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

## 🌦️ 5. Preparación de Variables Climáticas

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
```r
## ============================================================================
## PANELES B y C: ANÁLISIS DE MULTICOLINEALIDAD (15 VARIABLES)
## Este análisis justifica la selección de solo 3 variables
## Excluye: bio8, bio9, bio18, bio19
## ============================================================================

library(ggplot2)
library(terra)
library(reshape2)
library(dplyr)

cat(">>> Analizando multicolinealidad de 15 variables bioclimáticas...\n")

# Cargar las 19 capas bioclimáticas
bio_global <- geodata::worldclim_global(var = "bio", res = 10, path = ".")
names(bio_global) <- paste0("bio", 1:19)
americas_extent <- ext(-170, -30, -55, 75)
env_19 <- terra::crop(bio_global, americas_extent)

# SELECCIONAR SOLO 15 VARIABLES (excluir bio8, bio9, bio18, bio19)
vars_to_use <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                 "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", 
                 "bio16", "bio17")

env_15 <- env_19[[vars_to_use]]

cat("Variables incluidas:", paste(vars_to_use, collapse = ", "), "\n")
cat("Variables excluidas: bio8, bio9, bio18, bio19\n")

# Extraer valores (muestra de 5000 puntos)
env_15_values <- as.data.frame(env_15, xy = FALSE, na.rm = TRUE)
if (nrow(env_15_values) > 5000) {
  set.seed(123)
  env_15_values <- env_15_values[sample(1:nrow(env_15_values), 5000), ]
}

# Calcular matriz de correlación
cor_matrix_15 <- cor(env_15_values, use = "complete.obs")

## --- PANEL B: MDS de 15 variables ---
cat("\n>>> Panel B: MDS de 15 variables...\n")

mds_result <- cmdscale(1 - abs(cor_matrix_15), k = 2)
mds_df <- as.data.frame(mds_result)
colnames(mds_df) <- c("MDS1", "MDS2")
mds_df$variable <- rownames(cor_matrix_15)

# Resaltar las 3 variables seleccionadas (bio10, bio12, bio15)
mds_df$selected <- ifelse(mds_df$variable %in% c("bio10", "bio12", "bio15"), 
                           "Selected", "Not selected")

panel_b <- ggplot(mds_df, aes(x = MDS1, y = MDS2, label = variable, 
                              color = selected, size = selected)) +
  geom_point(alpha = 0.7) +
  geom_text(size = 3, fontface = "bold", nudge_y = 0.02) +
  scale_color_manual(values = c("Selected" = "#d73027", 
                                "Not selected" = "gray60")) +
  scale_size_manual(values = c("Selected" = 5, "Not selected" = 3)) +
  labs(
    title = "(b) Variable Correlation Structure (MDS)",
    subtitle = "15 bioclimatic variables - selected variables (bio10, bio12, bio15) in red",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 9),
    legend.position = "none",
    panel.grid = element_line(color = "gray90")
  )

print(panel_b)
ggsave("Panel_B_MDS_15vars.png", panel_b, 
       width = 7, height = 6, dpi = 300, bg = "white")

## --- PANEL C: Heatmap de correlación (15 variables) ---
cat("\n>>> Panel C: Heatmap de correlación (15 variables)...\n")

cor_melted <- melt(cor_matrix_15)
colnames(cor_melted) <- c("Var1", "Var2", "value")

# Marcar pares con correlación alta (|r| > 0.7)
cor_melted$high_cor <- abs(cor_melted$value) > 0.7 & cor_melted$value != 1

# Marcar las 3 variables seleccionadas
cor_melted$selected_var <- (cor_melted$Var1 %in% c("bio10", "bio12", "bio15")) | 
                           (cor_melted$Var2 %in% c("bio10", "bio12", "bio15"))

panel_c <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  # Agregar borde rojo a las variables seleccionadas
  geom_tile(data = subset(cor_melted, selected_var), 
            color = "#d73027", linewidth = 0.8, fill = NA) +
  scale_fill_gradient2(
    low = "#4575b4", 
    mid = "white", 
    high = "#d73027",
    midpoint = 0, 
    limit = c(-1, 1),
    name = "Pearson\nr"
  ) +
  labs(
    title = "(c) Correlation Matrix (15 variables)",
    subtitle = "Variables with |r| > 0.7 are highly correlated | Selected vars outlined in red",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "gray40", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed()

print(panel_c)
ggsave("Panel_C_Heatmap_15vars.png", panel_c, 
       width = 8, height = 7, dpi = 300, bg = "white")

cat("\n✓ Paneles B y C guardados (justificación de selección de variables)\n")
cat("✓ Variables analizadas: 15 (excluidas bio8, bio9, bio18, bio19)\n")
cat("✓ Variables seleccionadas para modelo: bio10, bio12, bio15\n")

```

<img width="1335" height="584" alt="image" src="https://github.com/user-attachments/assets/526c9ebc-2351-4359-8219-aa309c9cae16" />


<img width="764" height="567" alt="image" src="https://github.com/user-attachments/assets/e88b47c1-64fe-47b0-ba38-c20711cfb727" />

<h3>🔍 Interpretación de la Selección de Variables</h3>

<p>
  Como se evidencia en el <b>MDS (Panel B)</b>, las variables seleccionadas (<span style="color: #d73027;"><b>puntos rojos</b></span>) están maximizadas en distancia, lo que indica que capturan dimensiones climáticas distintas (Temperatura, Precipitación y Estacionalidad).
</p>

<p>
  La <b>Matriz de Correlación (Panel C)</b> confirma que, mientras existe alta redundancia entre las variables térmicas (bloque bio1-bio11), las variables seleccionadas (<b>bio10, bio12, bio15</b>) presentan baja correlación cruzada (<i>|r| < 0.7</i>), reduciendo el riesgo de inflación de la varianza en el modelo.
</p>





## 🗺️ 6. Construcción del Fondo (Background M)

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



## ⚙️ 7. Ejecución del Modelo ENMeval

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


## 📊 8. Selección del Mejor Modelo y Predicción

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

<h3>🏆 Resumen del Mejor Modelo Seleccionado</h3>

<table border="0">
  <tr>
    <td width="40%"><b>Feature Class (fc):</b></td>
    <td><code>LQH</code> (Linear, Quadratic, Hinge)</td>
  </tr>
  <tr>
    <td><b>Regularization (rm):</b></td>
    <td><code>1</code></td>
  </tr>
  <tr>
    <td><b>AICc:</b></td>
    <td><code>9622.828</code></td>
  </tr>
  <tr>
    <td><b>AUC (Validación):</b></td>
    <td><code>0.8271</code> 🟢 <i>(Alto desempeño)</i></td>
  </tr>
  <tr>
    <td><b>AUC (Entrenamiento):</b></td>
    <td><code>0.8336</code></td>
  </tr>
  <tr>
    <td><b>Diferencia (Overfitting):</b></td>
    <td><code>0.006</code> <i>(Mínimo, modelo robusto)</i></td>
  </tr>
</table>

<p>✅ </b> Mapa de idoneidad generado y guardado exitosamente.</p>

<img width="796" height="578" alt="image" src="https://github.com/user-attachments/assets/18a3511e-ea5f-4c92-8a85-6acf552c5041" />


<h3>🌍 Interpretación Biogeográfica del Modelo</h3>

<p>
  La proyección espacial muestra una fuerte correspondencia entre la idoneidad climática predicha 
  (<span style="color: darkred;"><b>zonas rojas</b></span>) y las ocurrencias conocidas del patógeno 
  (<span style="color: blue;"><b>puntos azules</b></span>), validando visualmente el alto AUC (0.827).
</p>

<ul>
  <li><b>Patrón Latitudinal:</b> Se observa una clara restricción a zonas <b>tropicales y subtropicales</b>. El riesgo disminuye drásticamente en latitudes altas (>40°N y >40°S), lo cual es consistente con la limitación térmica del hongo (reflejada en la variable <i>bio10</i>).</li>
  
  <li><b>Hotspots de Riesgo Fitosanitario:</b> El modelo identifica zonas de máxima idoneidad (Clase > 0.8) en regiones clave para la producción de los hospederos estudiados:
    <ul>
      <li><b>Norteamérica:</b> Sureste de EE.UU. (Florida) y la costa del Golfo de México.</li>
      <li><b>Región Andina:</b> Zonas cafeteras y frutícolas de Colombia, Ecuador y Perú.</li>
      <li><b>Sudamérica:</b> El Cerrado y la Mata Atlántica en Brasil, zonas críticas para la producción de papaya y mango.</li>
    </ul>
  </li>

  <li><b>Validación del Fondo Restringido:</b> A diferencia de los modelos globales genéricos, este enfoque delimita correctamente la ausencia de riesgo en zonas áridas (ej. desiertos del norte de México/Sur de EE.UU.) donde, aunque la temperatura podría ser adecuada, la falta de precipitación (<i>bio12</i>, <i>bio15</i>) y la ausencia de hospederos impiden el establecimiento de la enfermedad.</li>
</ul>


<h3>🌱¿Qué controla al patógeno?</h3>

<p>
  El análisis de contribución de variables revela que la distribución de <i>C. gloeosporioides</i> en las Américas 
  está gobernada principalmente por la <b>energía térmica</b>, seguida de la disponibilidad hídrica.
</p>

<div align="center">
  <img src="https://github.com/user-attachments/assets/a89f2ab8-cdff-4cd6-b063-10502d428c1e" width="85%" alt="Gráfico de Importancia de Variables">
  <p><em>Fig 2. Importancia relativa de las variables bioclimáticas en el modelo Maxent.</em></p>
</div>

<p><b>Interpretación Biológica:</b></p>
<ul>
  <li>
    🔥 <b>Bio10 (Temperatura Media del Trimestre más Cálido):</b> Es la variable dominante (<i>Importance > 16</i>). 
    Esto confirma que el patógeno es dependediente de la temperatura. Su desarrollo está restringido principalmente por la falta de calor suficiente durante la temporada de crecimiento.
  </li>
  <li>
    💧 <b>Bio12 (Precipitación Anual):</b> Juega un papel secundario pero crítico (<i>Importance ~ 11.5</i>). 
    El hongo requiere niveles basales de humedad para la esporulación, excluyéndolo de zonas áridas aunque sean cálidas.
  </li>
  <li>
    📉 <b>Bio15 (Estacionalidad):</b> Tiene la menor influencia, indicando que el patógeno tolera cierta variabilidad en los patrones de lluvia siempre que se cumplan los requisitos de temperatura y humedad relativa.
  </li>
</ul>

## 📊 9. Evaluación de Calibración y Fiabilidad (Boyce & ECE)

Este bloque de código calcula métricas avanzadas para asegurar que el modelo no esté sobreajustado y genera el **Panel de Diagnóstico** (Panel D) con curvas de calibración y histogramas de discriminación.

<details>
<summary>📐 <strong>Click aquí para ver el código de Validación (R)</strong></summary>

```r
## ============================================================================
## PANEL D: CALIBRACIÓN DEL MODELO
## (Estilo publicación - con scatter plot, histogramas y métricas)
## ============================================================================

library(ggplot2); library(cowplot); library(terra); library(dplyr)

cat(">>> Generando Panel D: Calibración del Modelo...\n")

# --- 1. EXTRAER MEJOR MODELO Y HACER PREDICCIONES ---
best_idx <- which.min(eval_results@results$AICc)
best_model <- eval_results@models[[best_idx]]

# ... (Resto del código de cálculo de Boyce y ECE) ...

# --- 2. CALCULAR MÉTRICAS DE CALIBRACIÓN ---
# Índice de Boyce Continuo (simplificado)
calc_boyce <- function(pred_presence, pred_background, n_bins = 10) {
  # ... lógica de la función ...
  boyce <- cor(1:n_bins, ratio, method = "spearman")
  return(boyce)
}

# ... (Código de gráficos ggplot y cowplot) ...

ggsave("Panel_D_Model_Calibration.png", panel_d_final,
       width = 7, height = 6, dpi = 300, bg = "white")
```


<h3>📈 Validación Avanzada y Calibración</h3>

<p>
  Además del AUC, se evaluó la robustez del modelo mediante el <b>Índice de Boyce</b> (independiente del umbral) 
  y métricas de error de calibración (ECE/MCE) utilizando <i>n = 1,404</i> presencias de validación.
</p>

<div align="center">
  <table border="0">
    <tr>
      <td width="50%" align="center">
        <b>Discriminación (Presencias vs Fondo)</b><br>
        <img src="https://github.com/user-attachments/assets/a56..." width="100%" alt="Histogramas">
      </td>
      <td width="50%" align="center">
        <b>Curva de Calibración</b><br>
        <img src="<img width="633" height="557" alt="image" src="https://github.com/user-attachments/assets/85824cb7-4595-4016-a966-90e06c5beb35" />
" width="100%" alt="Boyce Index">
      </td>
    </tr>
  </table>
  
  <p>
    <code>Boyce Index: 0.976</code> 🌟 | 
    <code>ECE: 0.082</code> (Error Promedio Bajo) | 
    <code>MCE: 0.202</code> (Error Máximo)
  </p>
</div>

<h4>📌 Interpretación de Resultados:</h4>
<ul>
  <li>
    <b>Boyce Index (0.976):</b> Cercano al 1.0 teórico. Confirma que el modelo ordena la idoneidad casi perfectamente respecto a las presencias reales.
  </li>
  <li>
    <b>Calibración (ECE 0.08):</b> El <i>Expected Calibration Error</i> indica que, en promedio, la probabilidad predicha por el modelo solo se desvía un <b>8.1%</b> de la realidad observada.
  </li>
  <li>
    <b>Discriminación:</b> Los histogramas muestran una separación clara: el modelo asigna valores altos (>0.75) a la mayoría de las presencias (rojo), mientras mantiene el fondo (azul) en valores bajos.
  </li>
</ul>
