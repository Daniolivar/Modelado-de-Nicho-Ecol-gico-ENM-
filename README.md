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



📦 1. Configuración del Entorno
Para reproducir este análisis, es necesario instalar un ecosistema de paquetes de R especializados en bioinformática y análisis espacial. Este bloque asegura que todas las dependencias estén presentes.

Las librerías clave incluyen:

rgbif & geodata: Para la descarga automatizada de ocurrencias biológicas y capas climáticas.

CoordinateCleaner: Para la limpieza automatizada de errores geográficos comunes.

terra & sf: Para el manejo de datos raster y vectoriales (la base del análisis espacial).

ENMeval & maxnet: Para la calibración rigurosa del modelo y la ejecución del algoritmo Maxent (sin necesidad de Java).


```{r}

install.packages(c(
  "rgbif", 
  "CoordinateCleaner", 
  "geodata", 
  "terra", 
  "sf", 
  "dplyr", 
  "stringr", 
  "ENMeval", 
  "maxnet"
))

```


```{r}
# --- 1. CARGAR LIBRERÍAS ----
library(rgbif)
library(CoordinateCleaner) # ¡Clave para la limpieza!
library(geodata)
library(terra)
library(sf)
library(dplyr)
library(stringr) #
```

