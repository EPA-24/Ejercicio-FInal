#paquetes necesarios

library(readxl)
library(ggplot2) 
library(sf) 
library(terra)
library(viridis) 
library(dplyr) 
library(ggspatial) 
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)

######################################################################################

### Mapa de muertes por cáncer en México

# Cargar base de datos de muertes por cancer a nivel estatal en México, tomada de: Estadísticas de Defunciones Registradas (EDR) 2023 - INEGI
bd1 <- read_xlsx("muertes_cancer.xlsx")

# Carga archivo vectorial de entidades de México
ent_mex <- st_read("entidades_mex.shp") 

# Realizar un left_join para fusionar los datos, asignar al archivo espacial
bd_join1 <- left_join(ent_mex, bd1, join_by("NOMGEO"=="NOM_ENT"))
bd_join1

#Mapa de muertes por cancer en mexcio

g1 <- ggplot(data = bd_join1) +
  geom_sf(aes(fill=std_MXC))+
  scale_fill_viridis(discrete = FALSE, option = "D") +
  labs(title = "Mortalidad de cáncer en México por estado",
       fill = "Muertes c/100 habitantes") +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.75, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering)+
  theme_minimal()
g1

##Volcano plot de genes diferencialmente expresados
# Cargar base de datos
bd2 <- read_xlsx("Base_lirio_completa.xlsx")
bd2

# Crear columna de selección diferencial
bd2 <- bd2 %>%
  mutate(Expression = case_when(
    Fold_change > 2 ~ "Sobreexpresado",
    Fold_change < -2 ~ "Subexpresado",
    TRUE ~ "ND"
  ))

# Se crea volcano plot
g2 <- ggplot(bd2, aes(x = Fold_change, y = -log10(p_value))) +
  geom_point(aes(color = Expression), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Sobreexpresado" = "red",
                                "Subexpresado" = "blue",
                                "ND" = "gray")) +
  geom_text_repel(data = bd2 %>%
                    slice_min(order_by = p_value, n = 10),
                  aes(label = GENE_ID),
                  size = 2.5,
                  box.padding = 0.2,
                  point.padding = 0.2,
                  max.overlaps = 20,) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 12)+
  labs(
    x = "Fold Change",
    y = "-Log10(p-value)",
    color = "Estado",
    title = "Expresión diferencial en células Hct-116 tratadas con Liriodenina",
    subtitle = "tratadas vs no tratadas"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top")+
  coord_cartesian(xlim = c(-50, 50), ylim = c(0, 10))

g2


####Barplot de genes sobre y sub expresados

g3 <- bd2%>%
  filter(Expression != "ND") %>%
  ggplot(aes(x = Expression, fill = Expression)) +
  geom_bar(width = 0.5, color = "black") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Sobreexpresado" = "red",
                               "Subexpresado" = "blue")) +
  ylab("Número de genes diferencialmente expresados") +
  ggtitle("Genes diferencialmente expresados tras tratamiento con Liriodenina") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
g3
#### Donut plot de DEGS y genes involucradps en ías de señalización de progresión tumoral en CCR.

# Filtrar base de datos con genes degs unucamente
bdDEG <- bd2 %>% filter(Expression != "ND")

# Cagar bases de datos de genes involuvcrados en vías de señalización
WNT_g <- read_xlsx("Wnt_genes.xlsx")
PI3K_g <- read_xlsx("PI3K_genes.xlsx")
MAPK_g <- read_xlsx("MAPK_genes.xlsx")
NOTCH_g <-read_xlsx("NOTCH_genes.xlsx") 
TGFB_g <- read_xlsx("TGFB_genes.xlsx")

# Crear conjuntos de genes en comun
WNT <- inner_join(bdDEG, WNT_g, by = "GENE_ID")
PI3K <- inner_join(bdDEG, PI3K_g, by = "GENE_ID")
MAPK <- inner_join(bdDEG, MAPK_g, by = "GENE_ID")
NOTCH <- inner_join(bdDEG, NOTCH_g, by = "GENE_ID")
TGFB <- inner_join(bdDEG, TGFB_g, by="GENE_ID")

# Crear un resumen con el número de genes por vía
datos_vias <- data.frame(
  Via_de_señalización = c("WNT", "PI3K", "MAPK", "NOTCH", "TGFB"),
  Genes = c(nrow(WNT), nrow(PI3K), nrow(MAPK), nrow(NOTCH), nrow(TGFB))
)
#Calcular el porcentaje de cada categoría
datos_vias <- datos_vias %>%
  mutate(Porcentaje = round(Genes / sum(Genes) * 100, 1),
         etiqueta = paste0(Via_de_señalización, "\n", Porcentaje, "%")
  )
# Crear el donut plot
g4 <- ggplot(datos_vias, aes(x = 2, y = Genes, fill = Via_de_señalización)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) + 
  theme_void() +    
  scale_fill_manual(values = c("WNT" = "blue",   
                               "PI3K" = "darkgreen", 
                               "MAPK" = "red",
                               "NOTCH"= "gold",
                               "TGFB"="orange")) +
  geom_text(aes(label = etiqueta),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5
            , fontface = "bold") +
  ggtitle("Porcentaje de involucrados en vías de señalización") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) 

g4

#### Heatmap de las vías mas representadas
PI3Khm<- PI3K %>%
  select(GENE_ID, B_1, B_2, B_3, LI_1, LI_2, LI_3)

# Preparar matriz de expresión 
PI3Kmat <- as.matrix(PI3Khm[,-1])  
rownames(PI3Kmat) <- PI3K$GENE_ID
PI3Kmat <- t(scale(t(PI3Kmat)))

#  Crear un vector para indicar tratamiento

col_anno <- c("Basal","Basal","Basal","Lirio","Lirio","Lirio")

# Definir paleta de colores
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Crear heatmap
g5 <-Heatmap(
  PI3Kmat,
  name = "Expresión",
  col = col_fun,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(
    Condición = col_anno,
    col = list(Condición = c("Basal" = "blue", "Lirio" = "red"))
  ),  # ← aquí cierra la anotación
  column_title = "Expresión diferencial en la vía PI3K/AKT",  
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

g5
