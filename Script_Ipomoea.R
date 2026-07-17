#'#################################################################
#' Script para análise da abundância de Ipomoea 
#' em diferentes geofácies
#' 
#' Perguntas
#i)	Há diferença na abundância e probabilidade de ocorrência Ipomoea cavalcantei conforme as diferentes geofácies?
#ii)A matriz circundante afeta a abundância de diferentes geo-ambientes? 
#' 
#'              ABRAÃO B. LEITE
#'#################################################################
library(bbmle)
library(brglm2)
library(DHARMa)
library(dplyr)
library(effects)
library(emmeans)
library(forcats)
library(FSA)
library(geobr)
library(ggeffects)
library(ggplot2)
library(ggrepel)
library(glmmTMB)
library(lmtest)
library(MASS)
library(multcomp)
library(MuMIn)
library(multcompView)
library(openxlsx)
library(pscl)
library(sandwich)
library(sf)
library(sjPlot)
library(terra)
library(tidyverse)
library(tidyr)
################################################################################
# SCRIPT DE ANÁLISE POPULACIONAL E ESPACIAL DE IPOMOEA EM CARAJÁS
# Módulo 1: Diferenças Locais Gerais (Geofácies vs. Platôs) com Zero-Inflação
# Módulo 2: Efeito do Microhabitat Circundante (Buffer Vetorial de 3m)
################################################################################


# Configuração global obrigatória para o pacote MuMIn (dredge) não falhar


# 1. CARREGAMENTO DOS DADOS ORIGINAIS
SnShape <- st_read("Data/PRJ112_DadosIpomoea_Densidades_v00.shp")
geof <- st_read("Data/Serra_Norte_2024_wLoc.shp")
IpomoeaData <- openxlsx::read.xlsx("PRJ152_DadosIpomoea_Densidades_v05.xlsx")%>%
mutate(ID_PARCELA= trimws(as.character(ID_PARCELA)),
       Geofacie= trimws(as.character(Geofacie)),
       Geofacie2= trimws(as.character(Geofacie2)),
       PLATO2= trimws(as.character(PLATO2)))

# 1. Ajustar o nome e o tipo da coluna no Shapefile para garantir o casamento dos dados
SnShape_df <- SnShape %>%
  st_drop_geometry() %>%
  rename(ID_PARCELA = ID_PARCELA) %>% # Garante o nome em maiúsculo se estiver diferente
  mutate(ID_PARCELA = trimws(as.character(ID_PARCELA)))
table(IpomoeaData $Geofacie, IpomoeaData $PLATO2)
# 2. Garantir que a coluna no Excel também seja texto
IpomoeaData <- IpomoeaData %>%
  mutate(ID_PARCELA = as.character(ID_PARCELA))

# 3. Fazer o join por ID_PARCELA, padronizar e filtrar
IpomoeaAbund <- SnShape_df %>%
  left_join(IpomoeaData, by = "ID_PARCELA") %>%
  
  # Padroniza as grafias usando a coluna do Excel (Geofacie2)
  mutate(Geofacie = ifelse(Geofacie2 %in% c("vegetação rupestre aberta", 
                                            "Vegetação Rupestre Aberta", 
                                            "Vegetação Rupestre Arbustiva"), 
                           "Vegetação Rupestre", Geofacie2))%>%
   # Remove o Campo Brejoso
  filter(Geofacie != "Campo Brejoso")%>%
mutate(Geofacie = case_match(Geofacie,
                             "Campo Graminoso"    ~ "Grassland",
                             "Lajedo"             ~ "Rocky outcrop",
                             "Mata Baixa"         ~ "Low forest",
                             "Mata baixa"         ~ "Low forest",
                             "Vegetação Rupestre" ~ "Rupestrian vegetation",
                             .default =Geofacie)) %>%
  mutate(Geofacie = droplevels(as.factor(Geofacie)),
         PLATO2 = as.factor(PLATO2))

# 

# Encontra o valor máximo para definir o último corte
max_val <- max(IpomoeaAbund$N_Ind_Tota, na.rm = TRUE)

# Sequências de 10 em 10 a partir do 0 
cortes_sequencia <- seq(0, ceiling(max_val / 10) * 10, by = 10)

# Une o -1 (para isolar o zero) com a sequência de 10 em 10
todos_cortes <- c(-1, cortes_sequencia)

# Rotulagem de forma automatizada (ex: "0", "1-10", "11-20"...)
rotulos <- c("0", paste0(cortes_sequencia[-length(cortes_sequencia)] + 1, "-", cortes_sequencia[-1]))

IpomoeaAbund$Faixa_10 <- cut(IpomoeaAbund$N_Ind_Tota, breaks = todos_cortes, labels = rotulos, 
                             include.lowest = TRUE)
nrow(IpomoeaAbund)
# Gerar e salvar o gráfico de distribuição de abundância por parcelas
ggplot(IpomoeaAbund, aes(x = Faixa_10, fill = Faixa_10 == "0")) +
  geom_bar(width = 1.0, color = "white", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#e63946", "FALSE" = "#4a90e2")) +
  scale_y_continuous(breaks = seq(0, 1000, by = 50), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Abundance", y = "Number of plots") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 9)) +
  stat_count(geom = "text", aes(label = after_stat(count)), 
             vjust = -0.5, fontface = "bold", size = 3.5)

ggsave("PlotParc.png", width = 8, height = 7, dpi = 300)

##########################################################################################################
#                                MODELAGEM ABUNDÂNCIA
###############################################################################################################

options(na.action = "na.fail") 
modAbundZI <- glmmTMB(N_Ind_Total ~ Geofacie*PLATO2, 
                      data = IpomoeaAbund,
                      ziformula = ~ Geofacie, 
                      family = nbinom2)
#ATENÇÃO: Mensagen de aviso:
#In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
 # Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
#Esse aviso é problemático, pois há geofácie sem amostragem, e isso afeta a convergência.
# Filtra removendo o Campo Brejoso
# --- 2. SELEÇÃO AUTOMATIZADA DE MODELOS VIA AIC ---
dredge_resultado <- dredge(modAbundZI, rank = "AIC")
print(dredge_resultado)

# --- 3. EXTRAÇÃO E CHECAGEM DO MELHOR MODELO (VENCEDOR AIC) ---
# Usamos colchetes duplos [[1]] para extrair o objeto glmmTMB de dentro da lista
best_model <- get.models(dredge_resultado, subset = 1)[[1]]

# Exibir o sumário estatístico do modelo vencedor
summary(best_model)

####################################################################################
#                          PÓS TESTE
##################################################################################
#Comparação por Geofácie
postCondGeof<- emmeans(best_model, specs = pairwise ~ PLATO2 | Geofacie, component = "cond", type = "response")
postCondGeof$emmeans
postCondGeof$contrasts
CondLettersGeof<- cld(postCondGeof, Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>%
  # Remove os valores que o modelo não conseguiu estimar (nonEst) para não quebrar o gráfico
  filter(!is.na(response))
CondLettersGeof<-CondLettersGeof%>%
  mutate(Geofacie = as.factor(Geofacie),
         .group = trimws(.group))
ggplot(CondLettersGeof, aes(x = PLATO2, y = response, color = PLATO2)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5) +
  geom_text(aes(y = asymp.UCL * 1.15, label = .group), size = 4.5, fontface = "bold", color = "black") +
  facet_wrap(~Geofacie, scales = "free_y") + 
  scale_color_brewer(palette = "Set1") +
  labs(x = "Geofácie",y = "Abundância de Indivíduos (Média ± IC 95%)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "none") # Esconde a legenda pois os nomes já estão no eixo X
ggsave("PlotGeo.tiff", width = 7, height = 6, dpi = 300, device = "tiff")
##################################################################################
#                      CONDITIONAL PLOT
##################################################################################
#Comparação por PLATO
postCondPlato <- emmeans(best_model, specs = pairwise ~Geofacie|PLATO2, component = "cond", type = "response")
postCondPlato$emmeans
postCondPlato$contrasts
CondLettersPlato<- cld(postCondPlato, Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>%
  # Remove os valores que o modelo não conseguiu estimar (nonEst) para não quebrar o gráfico
  filter(!is.na(response))
# Garante que o R trate a Geofacie como um fator e limpa as letrinhas
CondLettersPlato<- CondLettersPlato %>%
  mutate(Geofacie = as.factor(Geofacie),
         .group = trimws(.group))
ggplot(CondLettersPlato, aes(x = Geofacie, y = response, color = Geofacie)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5) +
  geom_text(aes(y = asymp.UCL * 1.15, label = .group), size = 4.5, fontface = "bold", color = "black") +
  facet_wrap(~ PLATO2, scales = "free_y") + 
  scale_color_brewer(palette = "Set1") +
  labs(x = "Geofácie",y = "Abundância de Indivíduos (Média ± IC 95%)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "none") # Esconde a legenda pois os nomes já estão no eixo X
ggsave("PlotPlato.tiff", width = 7, height = 6, dpi = 300, device = "tiff")

############################################################################################
#                       ZI PLOT
############################################################################################
#ComparaçãoZi
postZi <- emmeans(best_model, specs = pairwise ~ Geofacie, component = "zi", type = "response")
postZi$emmeans
postZi$contrasts
ZiLetters_Ordenado <- ZiLetters %>%
  mutate(prob_pct = response * 100,
         lcl_pct = asymp.LCL * 100,
         ucl_pct = asymp.UCL * 100) %>%
  mutate(Geofacie = as.factor(Geofacie),
         .group = trimws(.group)) %>%
  arrange(desc(prob_pct)) %>%
  mutate(Geofacie = factor(Geofacie, levels = unique(Geofacie)))
ggplot(ZiLetters_Ordenado, aes(x = Geofacie, y = prob_pct, color = Geofacie)) +
  geom_errorbar(aes(ymin = lcl_pct, ymax = ucl_pct), width = 0.15, linewidth = 1) +
  geom_point(size = 4.5) +
  scale_color_brewer(palette = "Set1") +
  geom_text(aes(y = ucl_pct + 4, label = .group), size = 5, fontface = "bold", color = "black") +
  scale_y_continuous(limits = c(0, 110)) +
  labs(x = "Geo-facies",y = "Structural Zero Probability (%)",) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(hjust = 1),
        legend.position = "none")
ggsave("ZiPlot.tiff", width = 7, height = 6, dpi = 300, device = "tiff")


##############################################################################
#                  MODELAGEM DA MATRIZ CIRCUNDANTE (5 METROS)
###############################################################################

SnShape <- st_read("Data/PRJ112_DadosIpomoea_Densidades_v00.shp")
geof    <- st_read("Data/Serra_Norte_2024_wLoc.shp")
IpomoeaData <- openxlsx::read.xlsx("PRJ152_DadosIpomoea_Densidades_v05.xlsx") %>%
  mutate(across(where(is.character), trimws))
SnShape_Unico <- SnShape %>% 
  mutate(ID_PARCELA = trimws(as.character(ID_PARCELA))) %>%
  distinct(ID_PARCELA, .keep_all = TRUE)
# Transforma os pontos únicos para o CRS métrico do mapa de fisionomias
pontos_utm <- st_transform(SnShape_Unico, st_crs(geof))

# Criando a zona de amortecimento (buffer) de 5 metros
buffers_parcelas <- st_buffer(pontos_utm, dist = 5)

# Interseção Vetorial (Recorta o mapa fisionômico no entorno de 5m)
intersecao <- st_intersection(buffers_parcelas, geof)
intersecao$area_pedaco <- as.numeric(st_area(intersecao))
df_proporcoes <- intersecao %>%
  st_drop_geometry() %>% 
  mutate(GeoF_2024 = dplyr::recode(GeoF_2024,
                                   "Campo Graminoso"              = "Grassland",
                                   "Lajedo"                       = "Rocky_outcrop",
                                   "Mata baixa"                   = "Low_forest",
                                   "Mata Baixa"                   = "Low_forest",
                                   "Vegetação Rupestre Arbustiva" = "Rupestrian_vegetation",
                                   "Vegetação Rupestre Aberta"    = "Rupestrian_vegetation",
                                   "vegetação rupestre aberta"    = "Rupestrian_vegetation")) %>%
  # Exclui o Campo Brejoso da paisagem circundante
  filter(GeoF_2024 != "Campo Brejoso") %>%
  group_by(ID_PARCELA, GeoF_2024) %>% 
  summarise(Area_Total_Classe = sum(area_pedaco), .groups = "drop") %>%
  group_by(ID_PARCELA) %>%
  mutate(Proporcao_Entorno = Area_Total_Classe / sum(Area_Total_Classe)) %>%
  ungroup()

# Pivotagem para formato Largo
df_matriz_larga <- df_proporcoes %>%
  dplyr::select(ID_PARCELA, GeoF_2024, Proporcao_Entorno) %>%
  pivot_wider(names_from = GeoF_2024, values_from = Proporcao_Entorno, values_fill = 0)
 
df_analise_restrito <- IpomoeaData %>%
  # 4.1 Padroniza as grafias locais do Excel usando a coluna Geofacie2
  mutate(Geofacie_Tratada = ifelse(Geofacie2 %in% c("vegetação rupestre aberta", 
                                                    "Vegetação Rupestre Aberta", 
                                                    "Vegetação Rupestre Arbustiva"), 
                                   "Vegetação Rupestre", Geofacie2)) %>%
  
  # 4.2 Remove o Campo Brejoso local baseado nos dados reais do Excel
  filter(Geofacie_Tratada != "Campo Brejoso") %>%
  
  # 4.3 Traduz os nomes locais para o inglês para casar com as matrizes
  mutate(Geofacie_Local = case_match(Geofacie_Tratada,
                                     "Campo Graminoso"    ~ "Grassland",
                                     "Lajedo"             ~ "Rocky outcrop",
                                     "Mata Baixa"         ~ "Low forest",
                                     "Mata baixa"         ~ "Low forest",
                                     "Vegetação Rupestre" ~ "Rupestrian vegetation",
                                     .default = Geofacie_Tratada
  )) %>%
  
  # 4.4 Junta as proporções de paisagem calculadas no entorno de 5m
  left_join(df_matriz_larga, by = "ID_PARCELA") %>%
  
  # 4.5 Seleciona e renomeia as colunas finais exatas do Excel para o modelo
  dplyr::select(ID_PARCELA, PLATO2, Geofacie_Local, N_Ind_Total, 
                Grassland, Rocky_outcrop, Rupestrian_vegetation) %>%
  dplyr::rename(PLATO = PLATO2,N_Ind_Total=N_Ind_Total)

# Preenche com zero fisionomias contínuas vazias e converte em fatores
df_analise_restrito[is.na(df_analise_restrito)] <- 0
df_analise_restrito <- df_analise_restrito %>%
  mutate(Geofacie_Local = as.factor(Geofacie_Local),
         PLATO = as.factor(PLATO))

# --- Diagnóstico visual do esforço amostral (Baseado puramente no Excel) ---
cat("\n--- Distribuição das Parcelas Locais (Referência Excel) ---\n")
print(table(df_analise_restrito$Geofacie_Local))

# 1. Filtra, consolida e transforma a proporção em porcentagem real
tabela_porcentagem <- df_proporcoes %>%
  left_join(df_analise_restrito %>% dplyr::select(ID_PARCELA, PLATO), by = "ID_PARCELA") %>%
  filter(!is.na(PLATO)) %>%
  
  mutate(GeoF_2024 = case_match(GeoF_2024,
                                "Rocky_outcrop"         ~ "Rocky outcrop",
                                "Rupestrian_vegetation" ~ "Rupestrian vegetation",
                                "Low_forest"            ~ "Low forest",
                                "Grassland"             ~ "Grassland",
                                .default = GeoF_2024
  )) %>%
  
  filter(GeoF_2024 %in% c("Grassland", "Low forest", "Rocky outcrop", "Rupestrian vegetation")) %>%
  
  group_by(PLATO, GeoF_2024) %>%
  summarise(
    Total_Area_m2 = sum(Area_Total_Classe, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  
  # MULTIPLICA POR 100 PARA GERAR A PORCENTAGEM
  mutate(Percentage_Val = (Total_Area_m2 / sum(Total_Area_m2)) * 100) %>%
  ungroup() %>%
  
  # Formata os números e adiciona o símbolo de % colado ao valor
  mutate(
    Total_Area_m2 = round(Total_Area_m2, 2),
    `Percentage (%)` = paste0(round(Percentage_Val, 1), "%")
  ) %>%
  
  # Seleciona e renomeia para o layout oficial em inglês do artigo
  dplyr::select(PLATO, GeoF_2024, Total_Area_m2, `Percentage (%)`) %>%
  dplyr::rename(
    Plateau = PLATO,
    `Geo-facies (Matrix)` = GeoF_2024,
    `Total Area (m²)` = Total_Area_m2
  ) %>%
  arrange(Plateau, `Geo-facies (Matrix)`)

# 2. Exibe a nova tabela com porcentagem no console
print(as.data.frame(tabela_porcentagem))

#Parcelas na Mata Baixa
LowForest<- glmmTMB(N_Ind_Total~ Rocky_outcrop + Grassland + Rupestrian_vegetation, 
                    ziformula = ~ Rocky_outcrop + Grassland + Rupestrian_vegetation, 
                    data = df_analise_restrito,family = nbinom2)
#“Essa fisionomia altera a abundância de forma estatisticamente relevante?”
cat("\n--- PARCELAS LOW FOREST ---\n")
car::Anova(LowForest, type = "III")

df_analise_restrito <- df_analise_restrito %>%
  mutate(Low_forest = 1 - (Rocky_outcrop + Grassland + Rupestrian_vegetation))

#Parcelas em campo graminoso
Grassland<- glmmTMB(N_Ind_Total~ Rocky_outcrop + Rupestrian_vegetation + Low_forest, 
                    ziformula = ~ Rocky_outcrop + Rupestrian_vegetation + Low_forest, 
                    data = df_analise_restrito, family = nbinom2)
cat("\n--- PARCELAS GRASSLAND ---\n")
car::Anova(Grassland, type = "III")

# Parcelas no Lajedo
RockyOut<- glmmTMB(N_Ind_Total ~ Grassland + Rupestrian_vegetation + Low_forest, 
                   ziformula = ~ Grassland + Rupestrian_vegetation + Low_forest, 
                   data = df_analise_restrito, family = nbinom2)
cat("\n--- PARCELAS ROCKY OUTCROP ---\n")
car::Anova(RockyOut, type = "III")
#Fazer gráfico

# Parcelas na Vegetação rupestre
Rupestrian<- glmmTMB(N_Ind_Total ~ Rocky_outcrop + Grassland + Low_forest, 
                             ziformula = ~ Rocky_outcrop + Grassland + Low_forest, 
                             data = df_analise_restrito, family = nbinom2)
cat("\n--- PARCELAS RUPESTRIAN VEGETATION ---\n")
car::Anova(Rupestrian, type = "III")
#Fazer gráfico

################################################################################
#                   PLOTAGEM CORRIGIDA (SEM ERROS DE STRING)
################################################################################
library(ggplot2)
library(dplyr)

# 1. PREDIÇÃO 1 e 2: Efeito das Matrizes Rochosas (Modelo LowForest)
grade_rochas <- expand.grid(
  Proporcao = seq(0, 1, length.out = 100),
  Fisionomia = c("Rocky outcrop matrix effect on Low forest plot", 
                 "Rupestrian vegetation matrix effect on Low forest plot")
)

# CORRIGIDO: Adicionada a palavra 'matrix' no segundo ifelse para match perfeito
grade_rochas_mod <- grade_rochas %>%
  mutate(
    Rocky_outcrop = ifelse(Fisionomia == "Rocky outcrop matrix effect on Low forest plot", Proporcao, 0),
    Rupestrian_vegetation = ifelse(Fisionomia == "Rupestrian vegetation matrix effect on Low forest plot", Proporcao, 0),
    Grassland = 0
  )

pred_rochas <- predict(LowForest, newdata = grade_rochas_mod, type = "link", se.fit = TRUE)
grade_rochas$Fit <- exp(pred_rochas$fit)
grade_rochas$LCL <- exp(pred_rochas$fit - 1.96 * pred_rochas$se.fit)
grade_rochas$UCL <- exp(pred_rochas$fit + 1.96 * pred_rochas$se.fit)

# ==============================================================================
# 2. PREDIÇÃO 3: Efeito do ganho de Floresta no Lajedo (Modelo RockyOut)
# ==============================================================================
grade_f_laj <- expand.grid(
  Proporcao = seq(0, 1, length.out = 100),
  Fisionomia = "Low forest matrix effect on Rocky outcrop plot"
)

grade_f_laj_mod <- grade_f_laj %>%
  mutate(Low_forest = Proporcao, Grassland = 0, Rupestrian_vegetation = 0)

pred_f_laj <- predict(RockyOut, newdata = grade_f_laj_mod, type = "link", se.fit = TRUE)
grade_f_laj$Fit <- exp(pred_f_laj$fit)
grade_f_laj$LCL <- exp(pred_f_laj$fit - 1.96 * pred_f_laj$se.fit)
grade_f_laj$UCL <- exp(pred_f_laj$fit + 1.96 * pred_f_laj$se.fit)

# ==============================================================================
# 3. PREDIÇÃO 4: Efeito do ganho de Floresta na Veg. Rupestre (Modelo Rupestrian)
# ==============================================================================
grade_f_rup <- expand.grid(
  Proporcao = seq(0, 1, length.out = 100),
  Fisionomia = "Low forest matrix effect on Rupestrian vegetation plot"
)

grade_f_rup_mod <- grade_f_rup %>%
  mutate(Low_forest = Proporcao, Grassland = 0, Rocky_outcrop = 0)

pred_f_rup <- predict(Rupestrian, newdata = grade_f_rup_mod, type = "link", se.fit = TRUE)
grade_f_rup$Fit <- exp(pred_f_rup$fit)
grade_f_rup$LCL <- exp(pred_f_rup$fit - 1.96 * pred_f_rup$se.fit)
grade_f_rup$UCL <- exp(pred_f_rup$fit + 1.96 * pred_f_rup$se.fit)

# ==============================================================================
# 4. UNIFICAÇÃO NO DATAFRAME LONGO DO PAINEL
# ==============================================================================
grade_painel_final <- rbind(grade_rochas, grade_f_laj, grade_f_rup)

grade_painel_final$Fisionomia <- factor(grade_painel_final$Fisionomia, 
                                        levels = c("Rocky outcrop matrix effect on Low forest plot", 
                                                   "Rupestrian vegetation matrix effect on Low forest plot",
                                                   "Low forest matrix effect on Rocky outcrop plot", 
                                                   "Low forest matrix effect on Rupestrian vegetation plot"))

# ==============================================================================
# 5. PLOTAGEM DO GRÁFICO FINAL
# ==============================================================================
Plot_Landscape_4_Panels <- ggplot(grade_painel_final, aes(x = Proporcao * 100, y = Fit)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.15, fill = "darkblue", color = NA) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  facet_wrap(~ Fisionomia, ncol = 2, scales = "free_y") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 101)) +
  labs(x = "Matrix Cover Proportion within 5m Radius (%)", 
       y = "Predicted Ipomoea abundance (Individuals / Plot)",
       title = "Surrounding matrix continuous Effects on Ipomoea abundance",
       subtitle = "All significant slopes extracted from the baseline-alternated models (Anova Type III, p < 0.05)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold", size = 8.5))

print(Plot_Landscape_4_Panels)

ggsave("Matrix.tiff", 
       plot = Plot_Landscape_4_Panels, width = 9, height = 7, dpi = 300, device = "tiff")
