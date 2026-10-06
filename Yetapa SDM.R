install.packages(c("sf", "dplyr", "readr", "ggplot2", "rgbif"))

source("01_yetapa_GBIF_AOI_M.R")

library(sf)

M <- st_read(
  "outputs/AOI_yetapa_M.gpkg",
  layer = "M_yetapa",
  quiet = TRUE
)

# Export as a shapefile
st_write(
  M,
  "outputs/AOI_yetapa_M.shp",
  driver = "ESRI Shapefile",
  delete_layer = TRUE,
  quiet = TRUE
)

# Zip file creation
archivos_shp <- list.files(
  "outputs",
  pattern = "^AOI_yetapa_M\\.(shp|shx|dbf|prj|cpg)$",
  full.names = TRUE
)

zip(
  "outputs/AOI_yetapa_M_shapefile.zip",
  files = archivos_shp,
  flags = "-j"
)

#----#

library(terra)

setwd(
  "C:/Users/Lucho/GeoAI Projects/SDM project 1/YETAPÁ COLLAR"
)

# ------------------------------------------------------------
# 1. LOCALIZAR Y LEER EL GEOTIFF
# ------------------------------------------------------------

archivos <- list.files(
  path = "data/predictors",
  pattern = "^yetapa_predictores_actuales_1km\\.(tif|tiff)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(archivos) != 1) {
  stop(
    "Se esperaba un único archivo .tif o .tiff llamado ",
    "yetapa_predictores_actuales_1km en data/predictors. ",
    "Se encontraron: ", length(archivos)
  )
}

predictores <- rast(archivos[1])

print(predictores)

if (nlyr(predictores) != 6) {
  stop("Se esperaban seis bandas. Revisar el archivo exportado.")
}

# Orden de las bandas según el script de exportación de GEE
names(predictores) <- c(
  "pastizal",
  "arbolado",
  "humedal_herbaceo",
  "bio1_temp_media",
  "bio12_precipitacion",
  "bio15_estacionalidad_precipitacion"
)

# ------------------------------------------------------------
# 2. REVISAR RANGOS
# ------------------------------------------------------------

rangos <- global(
  predictores,
  fun = c("min", "max"),
  na.rm = TRUE
)

print(rangos)

# ------------------------------------------------------------
# 3. COMPROBAR COBERTURA EN LAS OCURRENCIAS
# ------------------------------------------------------------

archivo_ocurrencias <- "outputs/yetapa_ocurrencias_limpias.csv"

if (!file.exists(archivo_ocurrencias)) {
  stop("No se encontró: ", archivo_ocurrencias)
}

ocurrencias <- read.csv(archivo_ocurrencias)

if (!all(c("longitude", "latitude") %in% names(ocurrencias))) {
  stop("El CSV debe contener longitude y latitude.")
}

puntos <- vect(
  ocurrencias,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

puntos <- project(puntos, crs(predictores))

valores <- extract(
  predictores,
  puntos,
  ID = FALSE
)

cat(
  "\nOcurrencias con los seis predictores:",
  sum(complete.cases(valores)),
  "de", nrow(ocurrencias), "\n"
)

cat(
  "Ocurrencias sin alguno de los predictores:",
  sum(!complete.cases(valores)), "\n"
)

# ------------------------------------------------------------
# 4. VISUALIZAR LAS COBERTURAS
# ------------------------------------------------------------

plot(predictores[[1:3]])

dir.create("outputs", showWarnings = FALSE)

archivos <- list.files(
  "data/predictors",
  pattern = "^yetapa_predictores_actuales_1km\\.(tif|tiff)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(archivos) != 1) {
  stop("Se esperaba un único GeoTIFF de predictores.")
}

predictores <- rast(archivos[1])

if (nlyr(predictores) != 6) {
  stop("El GeoTIFF debe tener seis bandas.")
}

bandas <- c(
  "pastizal",
  "arbolado",
  "humedal_herbaceo",
  "bio1_temp_media",
  "bio12_precipitacion",
  "bio15_estacionalidad_precipitacion"
)

names(predictores) <- bandas

M <- vect(
  "outputs/AOI_yetapa_M.gpkg",
  layer = "M_yetapa"
)

ocurrencias <- read.csv(
  "outputs/yetapa_ocurrencias_limpias.csv",
  colClasses = c(timestamp = "character")
)

stopifnot(
  all(c("longitude", "latitude", "timestamp") %in%
        names(ocurrencias))
)

# Identificador para mantener trazabilidad con el CSV de entrada
ocurrencias$registro_id <- seq_len(nrow(ocurrencias))

puntos <- vect(
  ocurrencias,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

# ------------------------------------------------------------
# 2. CREAR GRILLA DE 1.000 × 1.000 M
# ------------------------------------------------------------

crs_modelo <- "EPSG:32721"

M_utm <- project(M, crs_modelo)
puntos_utm <- project(puntos, crs_modelo)

# Límites ajustados a múltiplos de 1.000 m
limites <- ext(M_utm)

grilla <- rast(
  xmin = floor(xmin(limites) / 1000) * 1000,
  xmax = ceiling(xmax(limites) / 1000) * 1000,
  ymin = floor(ymin(limites) / 1000) * 1000,
  ymax = ceiling(ymax(limites) / 1000) * 1000,
  resolution = 1000,
  crs = crs_modelo
)

# ------------------------------------------------------------
# 3. REPROYECTAR Y GUARDAR LOS SEIS PREDICTORES
# ------------------------------------------------------------

pred_utm <- project(
  predictores,
  grilla,
  method = "bilinear"
)

# Conservar celdas cuyo centro está dentro de M
pred_1km <- mask(
  pred_utm,
  M_utm,
  touches = FALSE
)

names(pred_1km) <- bandas

writeRaster(
  pred_1km,
  "outputs/yetapa_predictores_UTM21S_1km.tif",
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = "COMPRESS=LZW"
  )
)

# ------------------------------------------------------------
# 4. ASIGNAR CELDA Y EXTRAER PREDICTORES
# ------------------------------------------------------------

ocurrencias$celda_id <- cellFromXY(
  pred_1km,
  crds(puntos_utm)
)

valores <- extract(
  pred_1km,
  puntos_utm,
  ID = FALSE
)

tabla <- cbind(ocurrencias, valores)

validas <- !is.na(tabla$celda_id) &
  complete.cases(tabla[, bandas])

# Guardar registros excluidos para revisarlos, si los hubiera
write.csv(
  tabla[!validas, ],
  "outputs/yetapa_presencias_sin_predictores_1km.csv",
  row.names = FALSE,
  na = ""
)

tabla_valida <- tabla[validas, ]

if (nrow(tabla_valida) == 0) {
  stop("No quedaron presencias con predictores completos.")
}

# ------------------------------------------------------------
# 5. SELECCIONAR UNA PRESENCIA POR CELDA
# ------------------------------------------------------------

set.seed(42)

# Aleatorizar el orden y conservar la primera de cada celda
orden <- sample.int(nrow(tabla_valida))
seleccion <- tabla_valida[orden, ]

presencias_1km <- seleccion[
  !duplicated(seleccion$celda_id),
]

presencias_1km <- presencias_1km[
  order(presencias_1km$celda_id),
]

presencias_1km$presencia <- 1L
rownames(presencias_1km) <- NULL

# Coordenadas originales, no centros de píxel
presencias_vector <- vect(
  presencias_1km,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

# ------------------------------------------------------------
# 6. GUARDAR RESULTADOS
# ------------------------------------------------------------

# CSV mínimo de coordenadas y fecha
write.csv(
  presencias_1km[, c("longitude", "latitude", "timestamp")],
  "outputs/yetapa_presencias_1km.csv",
  row.names = FALSE,
  na = ""
)

# Tabla completa para modelar y auditar la selección
write.csv(
  presencias_1km,
  "outputs/yetapa_presencias_1km_con_predictores.csv",
  row.names = FALSE,
  na = ""
)

writeVector(
  presencias_vector,
  "outputs/yetapa_presencias_1km.gpkg",
  overwrite = TRUE
)

# ------------------------------------------------------------
# 7. CONTROLES
# ------------------------------------------------------------

control <- data.frame(
  ocurrencias_entrada = nrow(ocurrencias),
  sin_predictores_o_celda = sum(!validas),
  redundantes_por_celda =
    nrow(tabla_valida) - nrow(presencias_1km),
  presencias_1km = nrow(presencias_1km),
  celdas_duplicadas =
    sum(duplicated(presencias_1km$celda_id))
)

print(control)

cat("\nResolución en metros:\n")
print(res(pred_1km))

cat("\nRangos después de reproyectar:\n")
print(global(pred_1km, c("min", "max"), na.rm = TRUE))

plot(
  pred_1km[["pastizal"]],
  main = "Pastizal y presencias seleccionadas a 1 km"
)

lines(M_utm, col = "black")

points(
  project(presencias_vector, crs_modelo),
  col = "red",
  pch = 20,
  cex = 0.5
)

#---------------------------------------------------------------------------#
# ------------------------------------------------------------
# 1. CARGAR LOS RESULTADOS DEL PASO ANTERIOR
# ------------------------------------------------------------

pred <- rast("outputs/yetapa_predictores_UTM21S_1km.tif")

bandas <- c(
  "pastizal",
  "arbolado",
  "humedal_herbaceo",
  "bio1_temp_media",
  "bio12_precipitacion",
  "bio15_estacionalidad_precipitacion"
)

stopifnot(all(bandas %in% names(pred)))
pred <- pred[[bandas]]

M <- vect("outputs/AOI_yetapa_M.gpkg", layer = "M_yetapa")
M <- project(M, crs(pred))

originales <- read.csv(
  "outputs/yetapa_ocurrencias_limpias.csv",
  colClasses = c(timestamp = "character")
)

presencias <- read.csv(
  "outputs/yetapa_presencias_1km_con_predictores.csv",
  colClasses = c(timestamp = "character")
)

# ------------------------------------------------------------
# 2. IDENTIFICAR TODAS LAS CELDAS OCUPADAS
# ------------------------------------------------------------

p_originales <- vect(
  originales,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

p_originales <- project(p_originales, crs(pred))

celdas_ocupadas <- cellFromXY(pred, crds(p_originales))
celdas_ocupadas <- sort(unique(
  celdas_ocupadas[!is.na(celdas_ocupadas)]
))

# Verificar que los IDs de presencia corresponden a esta grilla
p_presencias <- vect(
  presencias,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

p_presencias <- project(p_presencias, crs(pred))

ids_actuales <- cellFromXY(pred, crds(p_presencias))

stopifnot(
  !anyNA(ids_actuales),
  all(ids_actuales == presencias$celda_id),
  !anyDuplicated(ids_actuales)
)

# ------------------------------------------------------------
# 3. CREAR MÁSCARA DE CELDAS CANDIDATAS
# ------------------------------------------------------------

# 1 = seis predictores disponibles; NA = no elegible
n_disponibles <- app(!is.na(pred), fun = "sum")

mascara <- ifel(n_disponibles == length(bandas), 1, NA)
mascara <- mask(mascara, M, touches = FALSE)

# Excluir celdas con cualquiera de las ocurrencias originales
mascara[celdas_ocupadas] <- NA
names(mascara) <- "candidata_pseudoausencia"

writeRaster(
  mascara,
  "outputs/yetapa_mascara_pseudoausencias_1km.tif",
  overwrite = TRUE,
  wopt = list(datatype = "INT1U", gdal = "COMPRESS=LZW")
)

# ------------------------------------------------------------
# 4. MUESTREAR SIN REEMPLAZO
# ------------------------------------------------------------

candidatas <- which(
  values(mascara, mat = FALSE) == 1
)

semilla <- 42L
n_pseudo <- 2L * nrow(presencias)

if (length(candidatas) < n_pseudo) {
  stop("No hay suficientes celdas elegibles para el muestreo.")
}

set.seed(semilla)

celdas_pseudo <- candidatas[
  sample.int(length(candidatas), n_pseudo, replace = FALSE)
]

xy <- xyFromCell(pred, celdas_pseudo)

pseudo_utm <- vect(
  data.frame(x = xy[, 1], y = xy[, 2]),
  geom = c("x", "y"),
  crs = crs(pred)
)

lonlat <- crds(project(pseudo_utm, "EPSG:4326"))

valores_pseudo <- extract(pred, pseudo_utm, ID = FALSE)

pseudoausencias <- cbind(
  data.frame(
    muestra_id = paste0("PA_", seq_len(n_pseudo)),
    celda_id = celdas_pseudo,
    longitude = lonlat[, 1],
    latitude = lonlat[, 2],
    x_utm = xy[, 1],
    y_utm = xy[, 2],
    presencia = 0L,
    tipo = "pseudoausencia"
  ),
  valores_pseudo
)

# ------------------------------------------------------------
# 5. CONSTRUIR TABLA PARA EL MODELO
# ------------------------------------------------------------

xy_pres <- crds(p_presencias)

# Extraer nuevamente desde el raster definitivo
valores_pres <- extract(pred, p_presencias, ID = FALSE)

tabla_presencias <- cbind(
  data.frame(
    muestra_id = paste0("P_", presencias$registro_id),
    celda_id = ids_actuales,
    longitude = presencias$longitude,
    latitude = presencias$latitude,
    x_utm = xy_pres[, 1],
    y_utm = xy_pres[, 2],
    presencia = 1L,
    tipo = "presencia"
  ),
  valores_pres
)

datos_sdm <- rbind(tabla_presencias, pseudoausencias)

# ------------------------------------------------------------
# 6. VALIDAR ANTES DE GUARDAR
# ------------------------------------------------------------

control_pseudo <- data.frame(
  presencias = nrow(tabla_presencias),
  pseudoausencias = nrow(pseudoausencias),
  celdas_candidatas = length(candidatas),
  pseudo_en_celdas_ocupadas = sum(
    celdas_pseudo %in% celdas_ocupadas
  ),
  celdas_duplicadas = sum(duplicated(datos_sdm$celda_id)),
  muestras_sin_predictores = sum(
    !complete.cases(datos_sdm[, bandas])
  )
)

stopifnot(
  control_pseudo$pseudo_en_celdas_ocupadas == 0,
  control_pseudo$celdas_duplicadas == 0,
  control_pseudo$muestras_sin_predictores == 0
)

# ------------------------------------------------------------
# 7. GUARDAR TABLAS, PUNTOS Y METADATOS
# ------------------------------------------------------------

write.csv(
  pseudoausencias,
  "outputs/yetapa_pseudoausencias_1km.csv",
  row.names = FALSE
)

write.csv(
  datos_sdm,
  "outputs/yetapa_datos_sdm.csv",
  row.names = FALSE
)

write.csv(
  control_pseudo,
  "outputs/yetapa_control_pseudoausencias.csv",
  row.names = FALSE
)

writeVector(
  vect(
    pseudoausencias,
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  ),
  "outputs/yetapa_pseudoausencias_1km.gpkg",
  overwrite = TRUE
)

saveRDS(
  list(
    semilla = semilla,
    proporcion_pseudo_por_presencia = 2L,
    buffer_exclusion_m = 0,
    celdas_ocupadas = celdas_ocupadas,
    celdas_seleccionadas = celdas_pseudo,
    predictores = bandas,
    resolucion = res(pred),
    crs = crs(pred)
  ),
  "outputs/yetapa_configuracion_pseudoausencias.rds"
)

capture.output(
  sessionInfo(),
  file = "outputs/yetapa_sessionInfo_paso4.txt"
)

# ------------------------------------------------------------
# 8. MOSTRAR RESULTADOS
# ------------------------------------------------------------

print(control_pseudo)

plot(
  mascara,
  col = "grey90",
  legend = FALSE,
  main = "Presencias y pseudoausencias"
)

lines(M, col = "grey40")
points(pseudo_utm, col = "steelblue", pch = 20, cex = 0.4)
points(p_presencias, col = "red", pch = 20, cex = 0.5)

legend(
  "topright",
  legend = c("Presencias", "Pseudoausencias"),
  col = c("red", "steelblue"),
  pch = 20,
  bty = "n"
)

#----------------#

# ------------------------------------------------------------
# PASO 5A — AUTOCORRELACIÓN DE LOS PREDICTORES
# ------------------------------------------------------------

if (!requireNamespace("blockCV", quietly = TRUE)) {
  install.packages("blockCV")
}

library(terra)
library(blockCV)


dir.create("outputs", showWarnings = FALSE)

pred <- rast(
  "outputs/yetapa_predictores_UTM21S_1km.tif"
)

stopifnot(
  nlyr(pred) == 6,
  all(res(pred) == 1000),
  !is.lonlat(pred)
)

# Semilla fija para reproducir el muestreo
set.seed(42)

autocor <- blockCV::cv_spatial_autocor(
  r = pred,
  num_sample = 1000,
  plot = FALSE
)

# Resultados por predictor y resumen
print(autocor)

cat("\nTabla de rangos por predictor:\n")
print(autocor$range_table)

cat("\nMediana de los rangos, en km:\n")
print(autocor$range / 1000)

# Guardar el diagnóstico completo
saveRDS(
  autocor,
  "outputs/yetapa_autocorrelacion_predictores.rds"
)

write.csv(
  autocor$range_table,
  "outputs/yetapa_rangos_autocorrelacion.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "outputs/yetapa_sessionInfo_paso5a.txt"
)

#------#

# ------------------------------------------------------------
# PASO 5B — COMPROBAR TAMAÑOS DE BLOQUE VIABLES
# ------------------------------------------------------------

datos <- read.csv("outputs/yetapa_datos_sdm.csv")

stopifnot(
  all(c("x_utm", "y_utm", "presencia") %in% names(datos)),
  all(complete.cases(datos[, c("x_utm", "y_utm", "presencia")])),
  all(datos$presencia %in% c(0, 1))
)

evaluar_bloques <- function(tamano_km) {
  
  lado <- tamano_km * 1000
  
  # Grilla fija en UTM, con origen en (0, 0).
  # Todas las muestras de una celda de 1 km quedan juntas.
  bloque_id <- paste(
    floor(datos$x_utm / lado),
    floor(datos$y_utm / lado),
    sep = "_"
  )
  
  conteos <- aggregate(
    cbind(
      presencias = datos$presencia,
      pseudoausencias = 1L - datos$presencia
    ),
    by = list(bloque_id = bloque_id),
    FUN = sum
  )
  
  # Guardar composición de cada bloque para trazabilidad
  write.csv(
    conteos,
    paste0("outputs/yetapa_bloques_", tamano_km, "km_conteos.csv"),
    row.names = FALSE
  )
  
  data.frame(
    lado_km = tamano_km,
    bloques_con_muestras = nrow(conteos),
    bloques_con_presencias = sum(conteos$presencias > 0),
    bloques_con_ambas_clases = sum(
      conteos$presencias > 0 & conteos$pseudoausencias > 0
    ),
    max_presencias_en_un_bloque = max(conteos$presencias),
    porcentaje_presencias_en_bloque_mayor = round(
      100 * max(conteos$presencias) / sum(conteos$presencias),
      1
    )
  )
}

comparacion <- do.call(
  rbind,
  lapply(c(75, 160, 215), evaluar_bloques)
)

print(comparacion)

write.csv(
  comparacion,
  "outputs/yetapa_comparacion_tamanos_bloque.csv",
  row.names = FALSE
)


#----#

# ------------------------------------------------------------
# PASO 5C — DIVISIÓN ESPACIAL ENTRENAMIENTO / EVALUACIÓN
# ------------------------------------------------------------

datos <- read.csv("outputs/yetapa_datos_sdm.csv")

lado_m <- 75000
semilla <- 42L

# 1. Asignar bloques con el mismo origen del diagnóstico
datos$bloque_id <- paste(
  floor(datos$x_utm / lado_m),
  floor(datos$y_utm / lado_m),
  sep = "_"
)

bloques <- aggregate(
  cbind(
    presencias = datos$presencia,
    pseudoausencias = 1L - datos$presencia
  ),
  by = list(bloque_id = datos$bloque_id),
  FUN = sum
)

bloques <- bloques[order(bloques$bloque_id), ]

# 2. Buscar una partición equilibrada por clases
# Mantener aproximadamente 30 % de bloques en evaluación,
# tanto entre bloques con presencias como entre los restantes.

con_presencias <- bloques$bloque_id[bloques$presencias > 0]
sin_presencias <- bloques$bloque_id[bloques$presencias == 0]

n_test_con <- round(length(con_presencias) * 0.30)
n_test_sin <- round(length(sin_presencias) * 0.30)

n_candidatas <- 2000L
objetivo <- 0.30

set.seed(semilla)

mejor_desbalance <- Inf
bloques_test <- NULL

busqueda <- data.frame(
  iteracion = seq_len(n_candidatas),
  fraccion_presencias_test = NA_real_,
  fraccion_pseudoausencias_test = NA_real_,
  desbalance = NA_real_
)

for (i in seq_len(n_candidatas)) {
  
  candidatos <- c(
    con_presencias[
      sample.int(length(con_presencias), n_test_con)
    ],
    sin_presencias[
      sample.int(length(sin_presencias), n_test_sin)
    ]
  )
  
  seleccionados <- bloques$bloque_id %in% candidatos
  
  fraccion_pres <- sum(bloques$presencias[seleccionados]) /
    sum(bloques$presencias)
  
  fraccion_pseudo <- sum(bloques$pseudoausencias[seleccionados]) /
    sum(bloques$pseudoausencias)
  
  # Ambas clases reciben el mismo peso
  desbalance <- (fraccion_pres - objetivo)^2 +
    (fraccion_pseudo - objetivo)^2
  
  busqueda[i, 2:4] <- c(
    fraccion_pres, fraccion_pseudo, desbalance
  )
  
  if (desbalance < mejor_desbalance) {
    mejor_desbalance <- desbalance
    bloques_test <- candidatos
    mejor_iteracion <- i
  }
}

bloques$conjunto <- ifelse(
  bloques$bloque_id %in% bloques_test,
  "evaluacion",
  "entrenamiento"
)

datos$conjunto <- bloques$conjunto[
  match(datos$bloque_id, bloques$bloque_id)
]

train <- datos[datos$conjunto == "entrenamiento", ]
test  <- datos[datos$conjunto == "evaluacion", ]

# Registrar cómo se eligió la partición
saveRDS(
  list(
    semilla = semilla,
    n_candidatas = n_candidatas,
    objetivo_por_clase = objetivo,
    criterio = "Suma de desviaciones cuadradas respecto a 0.30",
    mejor_iteracion = mejor_iteracion,
    bloques_evaluacion = bloques_test,
    resultados = busqueda
  ),
  "outputs/yetapa_busqueda_particion_equilibrada.rds"
)

# 3. Controles de separación e integridad
bloques_compartidos <- intersect(
  unique(train$bloque_id),
  unique(test$bloque_id)
)

stopifnot(
  length(bloques_compartidos) == 0,
  !anyNA(datos$conjunto),
  nrow(train) + nrow(test) == nrow(datos),
  all(c(0, 1) %in% train$presencia),
  all(c(0, 1) %in% test$presencia)
)

resumen <- do.call(
  rbind,
  lapply(c("entrenamiento", "evaluacion"), function(grupo) {
    
    d <- datos[datos$conjunto == grupo, ]
    
    data.frame(
      conjunto = grupo,
      bloques = length(unique(d$bloque_id)),
      bloques_con_presencias =
        length(unique(d$bloque_id[d$presencia == 1])),
      presencias = sum(d$presencia == 1),
      pseudoausencias = sum(d$presencia == 0),
      total = nrow(d),
      porcentaje_muestras = round(100 * nrow(d) / nrow(datos), 1)
    )
  })
)

# 4. Guardar la partición para todas las etapas siguientes
write.csv(
  datos,
  "outputs/yetapa_datos_sdm_bloques.csv",
  row.names = FALSE
)

write.csv(
  train,
  "outputs/yetapa_entrenamiento.csv",
  row.names = FALSE
)

write.csv(
  test,
  "outputs/yetapa_evaluacion.csv",
  row.names = FALSE
)

write.csv(
  bloques,
  "outputs/yetapa_asignacion_bloques.csv",
  row.names = FALSE
)

write.csv(
  resumen,
  "outputs/yetapa_resumen_particion.csv",
  row.names = FALSE
)

saveRDS(
  list(
    lado_m = lado_m,
    origen_utm = c(0, 0),
    crs = "EPSG:32721",
    semilla = semilla,
    fraccion_evaluacion_objetivo = 0.30,
    metodo = "2000 asignaciones estratificadas; equilibrio de clases hacia 30%",
    asignacion = bloques
  ),
  "outputs/yetapa_configuracion_particion.rds"
)

# 5. Mostrar resultados
print(resumen)
cat("\nBloques compartidos:", length(bloques_compartidos), "\n")

colores <- c(
  entrenamiento = "steelblue",
  evaluacion = "darkorange"
)

plot(
  datos$x_utm,
  datos$y_utm,
  col = colores[datos$conjunto],
  pch = ifelse(datos$presencia == 1, 17, 1),
  cex = 0.6,
  asp = 1,
  xlab = "Este UTM (m)",
  ylab = "Norte UTM (m)",
  main = "Partición espacial — bloques de 75 km"
)

abline(
  v = seq(
    floor(min(datos$x_utm) / lado_m) * lado_m,
    ceiling(max(datos$x_utm) / lado_m) * lado_m,
    by = lado_m
  ),
  h = seq(
    floor(min(datos$y_utm) / lado_m) * lado_m,
    ceiling(max(datos$y_utm) / lado_m) * lado_m,
    by = lado_m
  ),
  col = "grey85"
)

legend(
  "topright",
  legend = c("Entrenamiento", "Evaluación"),
  col = colores,
  pch = 16,
  bty = "n"
)

#----#


