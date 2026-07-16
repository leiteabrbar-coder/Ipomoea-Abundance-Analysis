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

##################################################################################
#                            MODELAGEM MATRIZ CIRCUNDANTE
##############################################################################################

# --- 1. UNIÃO E PREPARAÇÃO DOS DADOS ESPACIAIS ---
# Une a planilha de abundância externa ao Shapefile geográfico original das parcelas
dados_completos <- left_join(SnShape, IpomoeaData, by = "ID_PARCELA")

# Corrige o CRS dos pontos para UTM Zone 22S (Métrico), igualando ao mapa geológico/fisionômico 'geof'
pontos_utm <- st_transform(dados_completos, st_crs(geof))

# --- 2. ANÁLISE ESPACIAL DE MICRO-ENTORNO (SIG NO R) ---
# Criando a zona de amortecimento (buffer) com base no raio de 5 metros
buffers_parcelas <- st_buffer(pontos_utm, dist = 5)

# Interseção Geométrica Vetorial pura (Corta o mapa 'geof' nas bordas circulares de cada buffer de 5m)
intersecao <- st_intersection(buffers_parcelas, geof)

# Calcula a área exata em metros quadrados de cada fragmento recortado
intersecao$area_pedaco <- as.numeric(st_area(intersecao))

# --- 3. CÁLCULO DAS PROPORÇÕES DA MATRIZ EM 5M (AGRUPADO) ---
# Padroniza e agrupa as fisionomias rupestres e mata baixa vindas do mapa 'geof' antes de calcular as áreas
df_proporcoes <- intersecao %>%
  st_drop_geometry() %>% # Remove dados espaciais brutos para acelerar a manipulação
  mutate(GeoF_2024 = recode(GeoF_2024,
                            "Campo Graminoso"="Grassland",
                            "Lajedo"= "Rocky outcrop",
                            "Mata baixa" = "Low forest",
                            "Vegetação Rupestre Arbustiva" = "Rupestrian vegetation",
                            "Vegetação Rupestre Aberta"    = "Rupestrian vegetation",
                            "vegetação rupestre aberta"    = "Rupestrian vegetation")) %>%
  group_by(ID_PARCELA, GeoF_2024) %>% 
  summarise(Area_Total_Classe = sum(area_pedaco), .groups = "drop") %>%
  group_by(ID_PARCELA) %>%
  mutate(Proporcao_Entorno = Area_Total_Classe / sum(Area_Total_Classe)) %>%
  ungroup()

# Pivotagem para o formato "Largo" (Cada geofácie do entorno vira uma variável contínua de 0 a 1)
df_matriz_larga <- df_proporcoes %>%
  dplyr::select(ID_PARCELA, GeoF_2024, Proporcao_Entorno) %>%
  pivot_wider(names_from = GeoF_2024, values_from = Proporcao_Entorno, values_fill = 0)

# --- 4. ESTRUTURAÇÃO DA TABELA DE MODELAGEM FINAL ---
# Cruza as proporções calculadas no entorno de 5m com os dados populacionais locais
df_analise_final <- pontos_utm %>%
  st_drop_geometry() %>% # Remove a geometria para blindar e não travar o glmmTMB
  dplyr::select(ID_PARCELA, PLATO.x, Geofacie.x, N_Ind_Tota) %>% 
  dplyr::rename(PLATO = PLATO.x, Geofacie_Local = Geofacie.x, Abundancia = N_Ind_Tota) %>%
  left_join(df_matriz_larga, by = "ID_PARCELA")

# Preenche com zero caso alguma parcela não tenha interceptado nenhuma geofácie mapeada no raio de 5m
df_analise_final[is.na(df_analise_final)] <- 0

# --- 5. LIMPEZA DE ESCOPO E FILTRAGEM BIÓTICA ---
# 1. Filtramos as linhas aceitando também a string em minúsculo que foi identificada
df_analise_restrito <- df_analise_final %>%
  dplyr::filter(Geofacie_Local %in% c("Lajedo", 
                                      "Campo Brejoso", 
                                      "Vegetação Rupestre Arbustiva", 
                                      "Vegetação Rupestre Aberta", 
                                      "vegetação rupestre aberta", 
                                      "Mata baixa", 
                                      "Campo Graminoso"))

# 2. Consolidação Metodológica: Corrigindo as maiúsculas/minúsculas e agrupando as rupestres locais
df_analise_restrito <- df_analise_restrito %>%
  mutate(Geofacie_Local = recode(Geofacie_Local, 
                                 "Mata baixa"="Mata Baixa",
                                 "Vegetação Rupestre Arbustiva" = "Vegetação Rupestre",
                                 "Vegetação Rupestre Aberta"    = "Vegetação Rupestre",
                                 "vegetação rupestre aberta"    = "Vegetação Rupestre"))

# Verificação diagnóstica
cat("\n--- Distribuição Corrigida das Parcelas ---\n")
print(table(df_analise_restrito$Geofacie_Local))

# --- 6. MODELAGEM ESTATÍSTICA DO MICRO-ENTORNO (glmmTMB) ---
options(na.action = "na.fail")

# Nota: Omitimos `Mata Baixa` da fórmula abaixo para servir de base/referência espacial,
# evitando multicolinearidade exata (soma das proporções do círculo sempre igual a 1).
modMatrizZI <- glmmTMB(
  Abundancia ~ Lajedo + `Campo Graminoso` + `Vegetação Rupestre` + PLATO, 
  ziformula = ~ Lajedo + `Campo Graminoso` + `Vegetação Rupestre` + PLATO, 
  data = df_analise_restrito, 
  family = nbinom2
)

# --- 7. SELEÇÃO AUTOMATIZADA VIA CRITÉRIO DE INFORMAÇÃO (AIC) ---
dredge_matriz <- dredge(modMatrizZI, rank = "AIC")
cat("\n--- Ranking de Modelos Candidatos (AIC) ---\n")
print(dredge_matriz)

# Extração cirúrgica do melhor modelo da matriz de 5m usando colchetes duplos [[]]
best_model_matriz <- get.models(dredge_matriz, subset = 1)[[1]]

# --- 8. EXIBIÇÃO DO VEREDITO ESTATÍSTICO DO MELHOR MODELO ---
cat("\n--- Sumário Estatístico do Melhor Modelo Selecionado (Entorno 5m) ---\n")
summary(best_model_matriz)


# ==============================================================================
# MÓDULO 2.1: GRÁFICO MULTIVARIÁVEL COM CORES E SÍMBOLOS DISTINTOS (5M)
# ==============================================================================

# 1. CRIAR GRADE DE PREDIÇÃO PARA LAJEDO
grade_lajedo <- expand.grid(
  Lajedo = seq(0, 1, length.out = 100),
  `Campo Graminoso` = 0, 
  `Vegetação Rupestre` = 0,
  PLATO = unique(df_analise_restrito$PLATO)
)
pred_lajedo <- predict(best_model_matriz, newdata = grade_lajedo, type = "link", se.fit = TRUE)
grade_lajedo$Abundancia_Prevista <- exp(pred_lajedo$fit)
grade_lajedo$IC_Inferior        <- exp(pred_lajedo$fit - 1.96 * pred_lajedo$se.fit)
grade_lajedo$IC_Superior        <- exp(pred_lajedo$fit + 1.96 * pred_lajedo$se.fit)
grade_lajedo$Tipo_Matriz         <- "Lajedo (Rocky Outcrop)"
grade_lajedo$Proporcao           <- grade_lajedo$Lajedo

# 2. CRIAR GRADE DE PREDIÇÃO PARA CAMPO GRAMINOSO
grade_campo <- expand.grid(
  Lajedo = 0, 
  `Campo Graminoso` = seq(0, 1, length.out = 100),
  `Vegetação Rupestre` = 0,
  PLATO = unique(df_analise_restrito$PLATO)
)
pred_campo <- predict(best_model_matriz, newdata = grade_campo, type = "link", se.fit = TRUE)
grade_campo$Abundancia_Prevista <- exp(pred_campo$fit)
grade_campo$IC_Inferior        <- exp(pred_campo$fit - 1.96 * pred_campo$se.fit)
grade_campo$IC_Superior        <- exp(pred_campo$fit + 1.96 * pred_campo$se.fit)
grade_campo$Tipo_Matriz         <- "Campo Graminoso (Grassland)"
grade_campo$Proporcao           <- grade_campo$`Campo Graminoso`

# 3. CRIAR GRADE DE PREDIÇÃO PARA VEGETAÇÃO RUPESTRE
grade_rupestre <- expand.grid(
  Lajedo = 0,
  `Campo Graminoso` = 0,
  `Vegetação Rupestre` = seq(0, 1, length.out = 100),
  PLATO = unique(df_analise_restrito$PLATO)
)
pred_rupestre <- predict(best_model_matriz, newdata = grade_rupestre, type = "link", se.fit = TRUE)
grade_rupestre$Abundancia_Prevista <- exp(pred_rupestre$fit)
grade_rupestre$IC_Inferior        <- exp(pred_rupestre$fit - 1.96 * pred_rupestre$se.fit)
grade_rupestre$IC_Superior        <- exp(pred_rupestre$fit + 1.96 * pred_rupestre$se.fit)
grade_rupestre$Tipo_Matriz         <- "Vegetação Rupestre"
grade_rupestre$Proporcao           <- grade_rupestre$`Vegetação Rupestre`

# 4. UNIFICAR AS PREDIÇÕES EM UM ÚNICO BANCO LONGO
grade_total <- rbind(
  grade_lajedo   %>% dplyr::select(PLATO, Tipo_Matriz, Proporcao, Abundancia_Prevista, IC_Inferior, IC_Superior),
  grade_campo    %>% dplyr::select(PLATO, Tipo_Matriz, Proporcao, Abundancia_Prevista, IC_Inferior, IC_Superior),
  grade_rupestre %>% dplyr::select(PLATO, Tipo_Matriz, Proporcao, Abundancia_Prevista, IC_Inferior, IC_Superior)
)

# 5. ORGANIZAR OS DADOS REAIS DAS PARCELAS NO FORMATO LONGO
df_pontos_longo <- df_analise_restrito %>%
  dplyr::select(ID_PARCELA, PLATO, Abundancia, Lajedo, `Campo Graminoso`, `Vegetação Rupestre`) %>%
  pivot_longer(cols = c(Lajedo, `Campo Graminoso`, `Vegetação Rupestre`),
               names_to = "Tipo_Matriz", values_to = "Proporcao") %>%
  mutate(Tipo_Matriz = recode(Tipo_Matriz, 
                              "Lajedo" = "Lajedo (Rocky Outcrop)",
                              "Campo Graminoso" = "Campo Graminoso (Grassland)"))

# 6. GERAR E EXIBIR O GRÁFICO FINAL VIA GGPLOT2
grafico_final <- ggplot() +
  geom_ribbon(data = grade_total, aes(x = Proporcao, ymin = IC_Inferior, ymax = IC_Superior, fill = Tipo_Matriz), alpha = 0.15) +
  geom_line(data = grade_total, aes(x = Proporcao, y = Abundancia_Prevista, color = Tipo_Matriz), size = 1.2) +
  geom_point(data = df_pontos_longo, aes(x = Proporcao, y = Abundancia, color = Tipo_Matriz, shape = PLATO), alpha = 0.5, size = 2) +
  facet_wrap(~PLATO) +
  labs(x = "Proporção da Geofácie no Entorno (5m)", 
       y = "Abundância de Indivíduos", 
       title = "Efeito da Matriz Circundante na Abundância por Platô",
       color = "Tipo de Geofácie", fill = "Tipo de Geofácie", shape = "Platô") +
  theme_classic() +
  theme(legend.position = "bottom")

print(grafico_final)
