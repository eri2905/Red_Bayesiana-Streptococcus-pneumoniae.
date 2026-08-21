## =====================================================================
## Paqueterías
## =====================================================================

library(readxl)
library(bnlearn)
library(openxlsx)
library(Rgraphviz)
library(graph)
library(igraph)
library(readxl)
library(corrplot)
library(ordinal)
library(igraph)
library(ggraph)
library(tidygraph)


# ---------------------------------------------------------------------
# 1. Carga de datos
# ---------------------------------------------------------------------


archivo_entrada <- "C:\\Users\\mier_\\OneDrive\\Escritorio\\BGC_DATA.xlsx"

datos <- readxl::read_excel(archivo_entrada, sheet = "gcfs")
raw <- as.data.frame(datos)

rownames(raw) <- raw$assembly_ID
raw$assembly_ID <- NULL

# ---------------------------------------------------------------------
# 2. Variable de resistencia -> ordinal
# ---------------------------------------------------------------------
# Orden clínico CLSI: susceptible < intermediate < resistant
raw$Updated_phenotype_CLSI <- factor(
  raw$Updated_phenotype_CLSI,
  levels  = c("susceptible", "intermediate", "resistant"),
  ordered = TRUE
)

table(raw$Updated_phenotype_CLSI)  # revisar balance de clases

# ---------------------------------------------------------------------
# 3. Identificar los grupos de variables
# ---------------------------------------------------------------------
gcf_vars <- grep("^GCF", names(raw), value = TRUE)
aro_vars <- grep("^ARO", names(raw), value = TRUE)

length(gcf_vars)  # 61
length(aro_vars)  # 19

# ---------------------------------------------------------------------
# 4. Detectar variables sin varianza
# ---------------------------------------------------------------------
# bnlearn exige al menos 2 niveles por nodo; una columna constante
# rompe el aprendizaje de estructura y hay que excluirla.
constantes <- names(raw)[sapply(raw, function(x) length(unique(x)) == 1)]
constantes

datos <- raw[, !(names(raw) %in% constantes)]
gcf_vars <- gcf_vars[gcf_vars %in% names(datos)]
aro_vars <- aro_vars[aro_vars %in% names(datos)]

# ---------------------------------------------------------------------
# 5. Binarizar GCF/ARO a presencia/ausencia
# ---------------------------------------------------------------------

datos[gcf_vars] <- lapply(datos[gcf_vars], function(x) factor(as.integer(x > 0)))
datos[aro_vars] <- lapply(datos[aro_vars], function(x) factor(as.integer(x > 0)))

str(datos[, c(1:3, which(names(datos) == "Updated_phenotype_CLSI"))])
dim(datos)  # 300 cepas x (61 GCF + 18 ARO + 1 resistencia) = 300 x 80

# ---------------------------------------------------------------------
# 6. ver prevalencia mínima: qué variables puede sostener un bootstrap
# ---------------------------------------------------------------------

prevalencia <- sapply(c(gcf_vars, aro_vars), function(v) min(table(datos[[v]])))

umbral_min <- 10  # prob. de perder el nivel minoritario en una remuestra ~ exp(-10) ~ 4.5e-5

genotipo_incluidas  <- names(prevalencia)[prevalencia >= umbral_min]
genotipo_excluidas  <- names(prevalencia)[prevalencia <  umbral_min]

length(genotipo_incluidas)  # variables que entran a la red bootstrap
length(genotipo_excluidas)  # variables demasiado raras para estabilizarse
genotipo_excluidas          # revisa la lista -- nada se pierde, queda aquí guardada



## =====================================================================
## YA Aprendizaje de estructura
## =====================================================================

set.seed(2026)


datos_red <- as.data.frame(datos[, c(genotipo_incluidas, "Updated_phenotype_CLSI")])

# ---------------------------------------------------------------------
# 1. Blacklist: arcos prohibidos
# ---------------------------------------------------------------------
# El fenotipo no puede causar el genotipo: prohibimos cualquier arco
# que salga de Updated_phenotype_CLSI hacia mutaciones o BGCs.
bl <- data.frame(
  from = "Updated_phenotype_CLSI",
  to   = genotipo_incluidas
)

# ---------------------------------------------------------------------
# Verificación previa op.
# ---------------------------------------------------------------------
niveles_ok <- sapply(datos_red, function(x) nlevels(droplevels(x)))
niveles_ok[niveles_ok < 2]                 # debe imprimir 'named integer(0)'
stopifnot(all(niveles_ok >= 2))
sapply(datos_red, table)[1:3]              # inspección rápida de balance de clases

# ---------------------------------------------------------------------
# 2. Aprendizaje con bootstrap (selección por estabilidad)
# ---------------------------------------------------------------------


str_boot <- boot.strength(
  data = datos_red,
  R    = 300,
  algorithm = "hc",
  algorithm.args = list(
    blacklist = bl,
    score     = "bic",
    maxp      = 5
  )
)

# arcos más fuertes en general
head(str_boot[order(-str_boot$strength), ], 20)

# ---------------------------------------------------------------------
# 3. Red consenso
# ---------------------------------------------------------------------

red_consenso <- averaged.network(str_boot)
narcs(red_consenso)

# Resistencia 
str_boot[(str_boot$from == "Updated_phenotype_CLSI" |
          str_boot$to   == "Updated_phenotype_CLSI") &
         str_boot$strength > 0.5, ]

# ---------------------------------------------------------------------
# 4. Visualización
# ---------------------------------------------------------------------
strength.plot(red_consenso, str_boot, shape = "ellipse")


## =====================================================================
## Diagnostico de resultado
## =====================================================================

# Ver que el b. distinga bien
hist(str_boot$strength, breaks = 30,
     main = "Distribución de fuerzas de arco (bootstrap)")
summary(str_boot$strength)

# candidatos debiles
cands_resistencia <- str_boot[str_boot$to == "Updated_phenotype_CLSI", ]
cands_resistencia[order(-cands_resistencia$strength), ]

# ver sensibilidad al score
str_boot_bde <- boot.strength(
  data = datos_red, R = 100, algorithm = "hc",
  algorithm.args = list(blacklist = bl, score = "bde", iss = 1, maxp = 5)
)
str_boot_bde[str_boot_bde$to == "Updated_phenotype_CLSI" &
             str_boot_bde$strength > 0.3, ]


## =====================================================================
## Confirmación con regresión ordinal
## =====================================================================



modelo_confirmatorio <- clm(
  Updated_phenotype_CLSI ~ ARO3000186 + ARO3000375 + ARO3003041,
  data = datos, link = "logit"
)
summary(modelo_confirmatorio)
exp(coef(modelo_confirmatorio))  # odds ratios

# ---------------------------------------------------------------------
# Por qué ARO3000375 pierde significancia: confundida con tetM
# ---------------------------------------------------------------------
# 49 de 53 cepas con ARO3000375 (92%) también tienen tetM -- su señal
# univariada iba montada sobre el efecto de tetM. Esto confirma la
# red BIC original, que ya la conectaba como hija de ARO3000186 en
# vez de directo a la resistencia.
table(datos$ARO3000375, datos$ARO3000186)

# Modelo final sin la variable redundante
modelo_final <- clm(
  Updated_phenotype_CLSI ~ ARO3000186 + ARO3003041,
  data = datos, link = "logit"
)
anova(modelo_confirmatorio, modelo_final)  # ¿se justifica quitar ARO3000375?

# Supuesto de odds proporcionales (una sola pendiente entre los dos
# cortes) -- verificar por el efecto tan grande en tetM
nominal_test(modelo_final)

# IC de verosimilitud perfilada (más confiables que Wald dado el OR
# tan extremo de tetM)
confint(modelo_final)


## =====================================================================
## Odds proporcionales no se sostienen
## =====================================================================
# Dejamos que ambos tengan pendiente propia por corte:

modelo_nominal <- clm(
  Updated_phenotype_CLSI ~ 1,
  nominal = ~ ARO3000186 + ARO3003041,
  data = datos, link = "logit"
)
summary(modelo_nominal)

# Comparación de ajuste (menor AIC = mejor)
AIC(modelo_final, modelo_nominal)


## =====================================================================
## Grafo final
## =====================================================================

# Subgrafo enfocado en la resistencia
vecinos <- c("Updated_phenotype_CLSI", bnlearn::mb(red_consenso, "Updated_phenotype_CLSI"))
vecinos <- unique(c(vecinos, unlist(lapply(vecinos, bnlearn::parents, x = red_consenso))))
sub_red <- bnlearn::subgraph(red_consenso, vecinos)
strength.plot(sub_red, str_boot, shape = "ellipse")

# --- Opción publicación: red completa, coloreada por tipo de nodo,
#     con los arcos hacia la resistencia resaltados en rojo ---
# install.packages(c("igraph", "ggraph", "tidygraph"))


nodos <- data.frame(name = bnlearn::nodes(red_consenso))
nodos$tipo <- ifelse(nodos$name == "Updated_phenotype_CLSI", "Resistencia",
              ifelse(grepl("^GCF", nodos$name), "BGC", "Mutación (ARO)"))

arcos <- as.data.frame(bnlearn::arcs(red_consenso))
arcos <- merge(arcos, str_boot, by = c("from", "to"))
arcos$hallazgo <- arcos$to == "Updated_phenotype_CLSI" | arcos$from == "Updated_phenotype_CLSI"
arcos$tipo_arco <- "bootstrap (bic)"

# El manto de Markov de la resistencia en esta red es SOLO {ARO3000186}

arco_extra <- data.frame(
  from = "ARO3003041", to = "Updated_phenotype_CLSI",
  strength = 1, hallazgo = TRUE, tipo_arco = "regresión (no bootstrap)"
)
arcos_full <- dplyr::bind_rows(arcos, arco_extra)

g <- tbl_graph(nodes = nodos, edges = arcos_full, directed = TRUE)

set.seed(2026)
p <- ggraph(g, layout = "fr") +
  geom_edge_link(
    aes(width = strength, color = hallazgo, linetype = tipo_arco),
    arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
    end_cap = circle(3, "mm"), alpha = 0.75
  ) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  scale_edge_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey65"), guide = "none") +
  scale_edge_linetype_manual(values = c("bootstrap (bic)" = "solid",
                                         "regresión (no bootstrap)" = "22"), name = NULL) +
  geom_node_point(aes(color = tipo), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("BGC" = "#4C72B0", "Mutación (ARO)" = "#DD8452",
                                 "Resistencia" = "#C44E52"), name = NULL) +
  theme_void() +
  labs(title = "Red bayesiana consenso — resistencia a tetraciclina",
       caption = "Punteado: arco confirmado por regresión (clm), ausente en el bootstrap bic")

print(p)
ggsave("red_bayesiana_final.pdf", p, width = 10, height = 8)
