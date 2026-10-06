# ============================================================
# PASO 6 — RANDOM FOREST Y EVALUACIÓN ESPACIAL
# ============================================================

paquetes <- c("ranger", "PRROC")
for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(ranger)
library(PRROC)

setwd(
  "C:/Users/Lucho/GeoAI Projects/SDM project 1/YETAPÁ COLLAR"
)

dir.create("outputs/modelo", recursive = TRUE, showWarnings = FALSE)
salida <- "outputs/modelo"

bandas <- c(
  "pastizal",
  "arbolado",
  "humedal_herbaceo",
  "bio1_temp_media",
  "bio12_precipitacion",
  "bio15_estacionalidad_precipitacion"
)

train <- read.csv("outputs/yetapa_entrenamiento.csv")

stopifnot(
  all(complete.cases(train[, bandas])),
  all(train$presencia %in% c(0, 1)),
  !anyDuplicated(train$celda_id)
)

# Factor: imprescindible para clasificación probabilística
train$respuesta <- factor(train$presencia, levels = c(0, 1))

# ------------------------------------------------------------
# 1. VALIDACIÓN INTERNA: TRES GRUPOS DE BLOQUES
# ------------------------------------------------------------

bloques <- aggregate(
  cbind(
    presencias = train$presencia,
    pseudoausencias = 1L - train$presencia
  ),
  by = list(bloque_id = train$bloque_id),
  FUN = sum
)

bloques <- bloques[order(bloques$bloque_id), ]

k <- 3L
set.seed(2026)

# Buscar grupos razonablemente equilibrados sin dividir bloques.
# Solo usamos los conteos del entrenamiento, nunca el test.
mejor <- Inf
mejor_fold <- NULL

for (i in seq_len(1000)) {
  
  fold <- integer(nrow(bloques))
  
  for (tiene_presencia in c(TRUE, FALSE)) {
    ids <- which((bloques$presencias > 0) == tiene_presencia)
    fold[ids] <- sample(rep(seq_len(k), length.out = length(ids)))
  }
  
  totales <- sapply(seq_len(k), function(f) {
    c(
      sum(bloques$presencias[fold == f]),
      sum(bloques$pseudoausencias[fold == f])
    )
  })
  
  if (any(totales == 0)) next
  
  proporciones <- sweep(totales, 1, rowSums(totales), "/")
  desbalance <- sum((proporciones - 1 / k)^2)
  
  if (desbalance < mejor) {
    mejor <- desbalance
    mejor_fold <- fold
  }
}

if (is.null(mejor_fold)) {
  stop("No se pudieron crear tres grupos con ambas clases.")
}

bloques$fold <- mejor_fold
train$fold <- bloques$fold[match(train$bloque_id, bloques$bloque_id)]

write.csv(
  bloques,
  file.path(salida, "bloques_validacion_interna.csv"),
  row.names = FALSE
)

cat("\nMuestras por grupo interno y clase:\n")
print(table(fold = train$fold, presencia = train$presencia))

# ------------------------------------------------------------
# 2. FUNCIONES DE AJUSTE Y MÉTRICAS
# ------------------------------------------------------------

ajustar_rf <- function(datos, mtry, nodo, semilla,
                       importancia = "none") {
  
  ranger::ranger(
    x = datos[, bandas, drop = FALSE],
    y = datos$respuesta,
    probability = TRUE,
    num.trees = 500,
    mtry = mtry,
    min.node.size = nodo,
    importance = importancia,
    num.threads = 2,
    seed = semilla
  )
}

calcular_auc <- function(y, p) {
  
  stopifnot(
    all(c(0, 1) %in% y),
    all(is.finite(p))
  )
  
  positivos <- p[y == 1]
  negativos <- p[y == 0]
  
  c(
    AUC_ROC = PRROC::roc.curve(
      scores.class0 = positivos,
      scores.class1 = negativos
    )$auc,
    AUC_PR = PRROC::pr.curve(
      scores.class0 = positivos,
      scores.class1 = negativos
    )$auc.integral
  )
}

# ------------------------------------------------------------
# 3. COMPARAR SEIS CONFIGURACIONES
# ------------------------------------------------------------

config <- expand.grid(
  mtry = c(2L, 3L, 4L),
  min.node.size = c(5L, 15L)
)

config$config_id <- seq_len(nrow(config))

# Cada muestra recibe una predicción de un modelo que
# no utilizó su bloque durante el ajuste.
pred_oof <- matrix(
  NA_real_,
  nrow = nrow(train),
  ncol = nrow(config)
)

resultados_fold <- list()
contador <- 1L

for (j in seq_len(nrow(config))) {
  
  for (f in seq_len(k)) {
    
    idx_validacion <- train$fold == f
    
    modelo <- ajustar_rf(
      datos = train[!idx_validacion, ],
      mtry = config$mtry[j],
      nodo = config$min.node.size[j],
      semilla = 1000L + f
    )
    
    prob <- predict(
      modelo,
      data = train[idx_validacion, bandas, drop = FALSE]
    )$predictions[, "1"]
    
    pred_oof[idx_validacion, j] <- prob
    
    auc <- calcular_auc(train$presencia[idx_validacion], prob)
    
    resultados_fold[[contador]] <- data.frame(
      config_id = j,
      fold = f,
      AUC_ROC = unname(auc["AUC_ROC"]),
      AUC_PR = unname(auc["AUC_PR"])
    )
    
    contador <- contador + 1L
  }
  
  cat("Configuración", j, "de", nrow(config), "terminada.\n")
}

resultados_fold <- do.call(rbind, resultados_fold)

medias <- aggregate(
  cbind(AUC_ROC, AUC_PR) ~ config_id,
  data = resultados_fold,
  FUN = mean
)

comparacion_rf <- merge(config, medias, by = "config_id")
comparacion_rf <- comparacion_rf[
  order(-comparacion_rf$AUC_PR, comparacion_rf$config_id),
]

elegida <- comparacion_rf[1, ]
j_mejor <- elegida$config_id

print(comparacion_rf)

write.csv(
  comparacion_rf,
  file.path(salida, "comparacion_configuraciones.csv"),
  row.names = FALSE
)

write.csv(
  resultados_fold,
  file.path(salida, "metricas_por_fold.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. UMBRAL A PARTIR DE VALIDACIÓN INTERNA
# ------------------------------------------------------------

p_oof <- pred_oof[, j_mejor]
stopifnot(!anyNA(p_oof))

# Evaluar umbrales de 0 a 1, en pasos de 0.001
umbrales <- seq(0, 1, by = 0.001)

sensibilidad <- vapply(umbrales, function(u) {
  mean(p_oof[train$presencia == 1] >= u)
}, numeric(1))

especificidad <- vapply(umbrales, function(u) {
  mean(p_oof[train$presencia == 0] < u)
}, numeric(1))

tabla_umbrales <- data.frame(
  umbral = umbrales,
  sensibilidad = sensibilidad,
  especificidad = especificidad,
  TSS = sensibilidad + especificidad - 1
)

# En caso de empate, elegir el umbral menor:
# favorece sensibilidad.
umbral <- tabla_umbrales$umbral[
  which.max(tabla_umbrales$TSS)
]

cat("\nUmbral elegido con validación interna:", umbral, "\n")

train$prediccion_oof <- p_oof

write.csv(
  train,
  file.path(salida, "predicciones_validacion_interna.csv"),
  row.names = FALSE
)

write.csv(
  tabla_umbrales,
  file.path(salida, "seleccion_umbral.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. MODELO FINAL AJUSTADO SOLO CON ENTRENAMIENTO
# ------------------------------------------------------------

modelo_final <- ajustar_rf(
  datos = train,
  mtry = elegida$mtry,
  nodo = elegida$min.node.size,
  semilla = 42L,
  importancia = "permutation"
)

saveRDS(
  list(
    modelo = modelo_final,
    predictores = bandas,
    umbral = umbral,
    configuracion = elegida,
    semilla_modelo = 42L,
    semilla_folds = 2026L,
    n_arboles = 500L,
    criterio_seleccion = "Mayor AUC-PR media en 3 folds espaciales",
    criterio_umbral = "Maximizar TSS con predicciones internas OOF",
    bloques_internos = bloques
  ),
  file.path(salida, "yetapa_random_forest.rds")
)

# ------------------------------------------------------------
# 6. EVALUACIÓN FINAL — DATOS RESERVADOS
# ------------------------------------------------------------

test <- read.csv("outputs/yetapa_evaluacion.csv")

stopifnot(
  length(intersect(train$bloque_id, test$bloque_id)) == 0,
  length(intersect(train$celda_id, test$celda_id)) == 0,
  all(complete.cases(test[, bandas])),
  all(c(0, 1) %in% test$presencia)
)

test$suitability <- predict(
  modelo_final,
  data = test[, bandas, drop = FALSE]
)$predictions[, "1"]

test$prediccion_binaria <- as.integer(test$suitability >= umbral)

TP <- sum(test$presencia == 1 & test$prediccion_binaria == 1)
FN <- sum(test$presencia == 1 & test$prediccion_binaria == 0)
TN <- sum(test$presencia == 0 & test$prediccion_binaria == 0)
FP <- sum(test$presencia == 0 & test$prediccion_binaria == 1)

auc_test <- calcular_auc(test$presencia, test$suitability)

metricas <- data.frame(
  AUC_ROC = unname(auc_test["AUC_ROC"]),
  AUC_PR = unname(auc_test["AUC_PR"]),
  referencia_PR = mean(test$presencia),
  umbral = umbral,
  sensibilidad = TP / (TP + FN),
  especificidad = TN / (TN + FP),
  precision = if ((TP + FP) > 0) TP / (TP + FP) else NA_real_,
  TSS = TP / (TP + FN) + TN / (TN + FP) - 1,
  TP = TP, FN = FN, TN = TN, FP = FP
)

write.csv(
  metricas,
  file.path(salida, "metricas_evaluacion_final.csv"),
  row.names = FALSE
)

write.csv(
  test,
  file.path(salida, "predicciones_evaluacion_final.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. CURVAS PARA EL REPORTE
# ------------------------------------------------------------

roc_test <- PRROC::roc.curve(
  scores.class0 = test$suitability[test$presencia == 1],
  scores.class1 = test$suitability[test$presencia == 0],
  curve = TRUE
)

pr_test <- PRROC::pr.curve(
  scores.class0 = test$suitability[test$presencia == 1],
  scores.class1 = test$suitability[test$presencia == 0],
  curve = TRUE
)

pdf(
  file.path(salida, "curvas_ROC_PR_evaluacion.pdf"),
  width = 10,
  height = 5
)

par(mfrow = c(1, 2))

plot(
  roc_test$curve[, 1], roc_test$curve[, 2],
  type = "l", lwd = 2, col = "steelblue",
  xlim = c(0, 1), ylim = c(0, 1),
  xlab = "1 - Specificity", ylab = "Sensitivity",
  main = sprintf("ROC — AUC = %.3f", roc_test$auc)
)
abline(0, 1, lty = 2, col = "grey50")

plot(
  pr_test$curve[, 1], pr_test$curve[, 2],
  type = "l", lwd = 2, col = "darkgreen",
  xlim = c(0, 1), ylim = c(0, 1),
  xlab = "Recall / Sensitibity", ylab = "Precision",
  main = sprintf("PR — AUC = %.3f", pr_test$auc.integral)
)
abline(h = mean(test$presencia), lty = 2, col = "grey50")

dev.off()

saveRDS(
  list(ROC = roc_test, PR = pr_test),
  file.path(salida, "curvas_evaluacion.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(salida, "sessionInfo.txt")
)

cat("\nCONFIGURACIÓN ELEGIDA:\n")
print(elegida)

cat("\nMÉTRICAS DE EVALUACIÓN FINAL:\n")
print(metricas)

#----#

# ============================================================
# PASO 7 — MAPAS DE SUITABILITY Y HÁBITAT BINARIO
# ============================================================

library(terra)
library(ranger)

setwd(
  "C:/Users/Lucho/GeoAI Projects/SDM project 1/YETAPÁ COLLAR"
)

carpeta_mapas <- "outputs/mapas"
dir.create(carpeta_mapas, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. CARGAR MODELO Y PREDICTORES
# ------------------------------------------------------------

ajuste <- readRDS("outputs/modelo/yetapa_random_forest.rds")

modelo <- ajuste$modelo
bandas <- ajuste$predictores
umbral <- ajuste$umbral

pred <- rast("outputs/yetapa_predictores_UTM21S_1km.tif")

stopifnot(
  all(bandas %in% names(pred)),
  all(res(pred) == 1000),
  length(umbral) == 1,
  is.finite(umbral),
  umbral >= 0,
  umbral <= 1
)

# Orden exacto de los predictores usados por el modelo
pred <- pred[[bandas]]

M <- vect("outputs/AOI_yetapa_M.gpkg", layer = "M_yetapa")
M <- project(M, crs(pred))

# Respetar el dominio del modelo
pred <- mask(pred, M, touches = FALSE)

cat("Umbral utilizado:", umbral, "\n")

# ------------------------------------------------------------
# 2. FUNCIÓN DE PREDICCIÓN PARA RANGER
# ------------------------------------------------------------

predecir_suitability <- function(model, data, ...) {
  
  resultado <- predict(
    model,
    data = as.data.frame(data),
    num.threads = 2
  )
  
  data.frame(
    suitability = resultado$predictions[, "1"]
  )
}

# ------------------------------------------------------------
# 3. MAPA CONTINUO
# ------------------------------------------------------------

archivo_continuo <- file.path(
  carpeta_mapas,
  "yetapa_suitability_actual_1km.tif"
)

suitability <- terra::predict(
  object = pred,
  model = modelo,
  fun = predecir_suitability,
  na.rm = TRUE,
  cores = 1,
  filename = archivo_continuo,
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = "COMPRESS=LZW"
  )
)

# Leer el raster guardado: será la referencia para el binario
suitability <- rast(archivo_continuo)
names(suitability) <- "suitability"

# ------------------------------------------------------------
# 4. MAPA BINARIO
# ------------------------------------------------------------

# 1 = hábitat adecuado según el modelo y umbral
# 0 = hábitat no adecuado según el modelo y umbral
# NA = fuera del dominio o sin predictores

binario <- ifel(suitability >= umbral, 1, 0)
names(binario) <- "habitat_adecuado"

archivo_binario <- file.path(
  carpeta_mapas,
  "yetapa_habitat_binario_actual_1km.tif"
)

writeRaster(
  binario,
  archivo_binario,
  overwrite = TRUE,
  wopt = list(
    datatype = "INT1U",
    NAflag = 255,
    gdal = "COMPRESS=LZW"
  )
)

# ------------------------------------------------------------
# 5. CONTROLES
# ------------------------------------------------------------

stopifnot(compareGeom(suitability, binario))

rangos <- global(suitability, c("min", "max"), na.rm = TRUE)
print(rangos)

stopifnot(
  is.finite(rangos$min[1]),
  is.finite(rangos$max[1]),
  rangos$min[1] >= 0,
  rangos$max[1] <= 1
)

cat("\nNúmero de celdas por clase:\n")
print(freq(binario))

# Superficie estimada por clase usando el área de las celdas
area_celdas <- cellSize(binario, unit = "km", mask = TRUE)

superficie <- zonal(
  area_celdas,
  binario,
  fun = "sum",
  na.rm = TRUE
)

names(superficie) <- c("clase", "superficie_km2")

superficie$descripcion <- ifelse(
  superficie$clase == 1,
  "Habitat adecuado",
  "Habitat no adecuado"
)

superficie$porcentaje_area_modelada <- round(
  100 * superficie$superficie_km2 /
    sum(superficie$superficie_km2),
  2
)

write.csv(
  superficie,
  file.path(carpeta_mapas, "superficie_por_clase.csv"),
  row.names = FALSE
)

print(superficie)

# ------------------------------------------------------------
# 6. VISUALIZACIÓN Y PNG DE CONTROL
# ------------------------------------------------------------

dibujar_mapas <- function() {
  
  anterior <- par(mfrow = c(1, 2))
  on.exit(par(anterior))
  
  plot(
    suitability,
    col = hcl.colors(100, "YlGnBu"),
    range = c(0, 1),
    main = "Yetapá de collar — Suitability",
    axes = TRUE
  )
  lines(M, col = "grey30", lwd = 0.7)
  
  plot(
    binario,
    breaks = c(-0.5, 0.5, 1.5),
    col = c("grey90", "#238B45"),
    legend = FALSE,
    main = paste0("Hábitat binario — umbral ", umbral),
    axes = TRUE
  )
  lines(M, col = "grey30", lwd = 0.7)
  
  legend(
    "bottomleft",
    legend = c("0: No adecuado", "1: Adecuado"),
    fill = c("grey90", "#238B45"),
    bty = "n",
    cex = 0.8
  )
}

png(
  file.path(carpeta_mapas, "yetapa_dos_mapas_control.png"),
  width = 2400,
  height = 1400,
  res = 200
)

dibujar_mapas()
dev.off()

dibujar_mapas()

# ------------------------------------------------------------
# 7. DOCUMENTAR LAS SALIDAS
# ------------------------------------------------------------

writeLines(
  c(
    "Especie: Alectrurus risora",
    "Modelo: outputs/modelo/yetapa_random_forest.rds",
    "Ajustado con entrenamiento; sin reentrenar con evaluacion.",
    paste("Umbral:", umbral),
    "CRS: EPSG:32721; celdas de 1000 x 1000 m.",
    "Continuo: puntuacion de suitability entre 0 y 1.",
    "Binario: 1 adecuado; 0 no adecuado; NoData fuera del dominio.",
    "La clasificacion no confirma presencia ni ausencia real.",
    "Paisaje: ESA WorldCover 2021; clima: WorldClim V1 BIO.",
    "Sin mascara especifica de agua permanente.",
    "Superficies calculadas sobre celdas con prediccion valida."
  ),
  file.path(carpeta_mapas, "metadatos_mapas.txt")
)

cat("\nMapas guardados en:", carpeta_mapas, "\n")
