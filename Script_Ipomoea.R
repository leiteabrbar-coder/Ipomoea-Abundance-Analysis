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
library(multcompView)
library(openxlsx)
library(pscl)
library(sandwich)
library(sf)
library(sjPlot)
library(terra)
library(tidyr)
#==========================================================================================================================
IpomoeaData<-openxlsx::read.xlsx("IpomoeaAbund.xlsx")
IpomoeaData$Ocorrencia%>%table()
IpomoeaData$Geofaceis <- toupper(trimws(IpomoeaData$Geofaceis))
unique(IpomoeaData$Geofaceis)
#patch_data<-readRDS("patches_data.rds")
SnShape<-st_read("GeoAmbFac_2022_img21_Amplo_v03.shp")
SnShape<-SnShape%>%filter(!(GeoF_2021%in% c("Mata Alta","Floresta Ombrófila","Lagoa Temporária","Samambaial","Solo Exposto",
                                            "Buritizal","Pasto","Lagoa Permanente",
                          "Estruturas Relativas à Mineração")))%>%mutate(GeoF_2021= toupper(trimws(GeoF_2021)))
SnShape2<-SnShape%>%group_by(GeoF_2021) %>%summarize(area_total_km2 = sum(as.numeric(st_area(geometry)),
                                              na.rm = TRUE) / 1000000) %>%st_drop_geometry()
IpomoeaData2<-IpomoeaData %>%left_join(SnShape2, by = c("Geofaceis" = "GeoF_2021"))
IpomoeaData2<-IpomoeaData2%>%rename(AreaTot=area_total_km2)

#==========================================================================================================
#1)A geoface determina a abundância?
#Campo graminoso é um filtro para espécie, quando ocorre há uma alta densidade
#Pois ocorre nos gaps das gramineas
#Lajedo e cangas é abundante porém menos densa, seria devido ao tamanho da área?
#==========================================================================================================
IpomoeaAbund<-IpomoeaData2 %>% 
  group_by(Geofaceis,Parc,Local,X,Y) %>%
  summarize(AbundTot=sum(Ocorrencia, na.rm = TRUE), .groups = "drop")
# Adicionando a fórmula de zero-inflação (ziformula = ~1 assume que a inflação é constante)
modAbundZI<- glmmTMB(AbundTot ~Geofaceis +Local,data = IpomoeaAbund,
                      ziformula = ~Geofaceis,family = nbinom2)
# Certifique-se de que o df_plot existe
df_plot <- as.data.frame(emmeans(modAbundZI, ~ Geofaceis, type = "response"))

ggplot() +
  # 1. Pontos individuais (Dados Brutos)
  geom_jitter(data = IpomoeaAbund,aes(x = Geofaceis, y = AbundTot+1, color = Geofaceis), 
              width = 0.2, alpha = 0.2, size = 1) +
  
  # 2. Barras de Erro (Médias do Modelo)
  geom_errorbar(data = df_plot, 
                aes(x = Geofaceis, ymin = response - SE, ymax = response + SE, color = Geofaceis), 
                width = 0.1, size = 1) +
  
  # 3. Ponto da Média (Diamante)
  geom_point(data = df_plot, 
             aes(x = Geofaceis, y = response, color = Geofaceis), 
             shape = 18, size = 5) +
  geom_line() +
  scale_y_log10()+

  # 4. Ajustes de Eixo e Estética
  scale_y_continuous(trans = "pseudo_log")+
  #scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
  labs(x = "Geofaces", y = "Abundância (Média ± SE)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=13),legend.position = "none",
        axis.text.y= element_text(size=13),axis.title = element_text(size=17))
ggsave("AbundGeo.jpeg",plot = last_plot(),width = 15,height = 10, units = "in",dpi = 1000) 
modAbundZI$fit$convergence # Deve ser 0
summary(modAbundZI)
simAbund_ZI <- simulateResiduals(modAbundZI)
testZeroInflation(simAbund_ZI)
medias_geofaceis<-emmeans(modAbundZI, ~ Geofaceis)
comppar<-pairs(medias_geofaceis, adjust = "tukey")
print(comppar)
cld(medias_geofaceis, Letters = letters)
dfletras <- as.data.frame(cld(medias_geofaceis, Letters = letters))
colnames(dfletras) <- trimws(colnames(dfletras))
resumoabund <- IpomoeaAbund %>%
  group_by(Geofaceis) %>%
  summarise(Total = sum(AbundTot))
resumoabund <- IpomoeaAbund %>%group_by(Geofaceis) %>%  summarise(Total = sum(AbundTot))
ggplot(medias_geofaceis, aes(x = Geofaceis, y = Total)) +
  geom_col(fill = "lightgray") + 
  geom_text(aes(label = Total),vjust = -0.5,size = 5,fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # Dá espaço para o texto no topo
  labs(x = "Geofácies", y = "Abundância Total") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
medias_geofacies <- IpomoeaAbund %>%
  group_by(Geofaceis) %>%
  summarise(Abund_Media = mean(AbundTot, na.rm = TRUE))

# Visualizar o resultado
print(medias_geofacies)
ggplot(IpomoeaAbund, aes(x = Geofaceis, y =AbundTot, fill= Geofaceis)) + 
  geom_boxplot( outlier.shape = 16) + 
  scale_y_continuous(limits = c(0, 50),breaks = seq(0, 50, by = 5),
                     expand = expansion(mult = c(0, 0.05)))+
  labs(x = "Geofácies", y = "Abundância") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
ggsave("AbundGeo.jpeg",plot = last_plot(),width = 15,height = 10, units = "in",dpi = 1000) 
pred_zi <- ggpredict(modAbundZI, terms = "Geofaceis", type = "zi_prob")
pred_zi$predicted <- 1 - pred_zi$predicted
pred_zi$conf.low <- 1 - pred_zi$conf.low   # Nota: inverter limites se necessário
pred_zi$conf.high <- 1 - pred_zi$conf.high
ggplot(pred_zi, aes(x = x, y = predicted)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_point(size = 3, color = "darkgreen") +
  labs(title = "Probabilidade de Ocorrência de Ipomoea por Geoface",
    x = "Geoface",y = "Probabilidade de Presença (0 a 1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), )

#====================================================================================================================
#2)A perda de matriz circundante da geoface com maior abundancia diminuirá essa abundância?
#Lajedo foi a principal matriz circundante
#====================================================================================================================
# Criar dataframe de abundância original
IpomoeaAbund<- IpomoeaData2 %>% 
  group_by(Geofaceis, Parc, AreaTot, Local, X, Y) %>%
  summarize(AbundTot = sum(Ocorrencia, na.rm = TRUE), .groups = "drop")

# Converter para objeto espacial
parcelas <- st_as_sf(IpomoeaAbund, coords = c("X", "Y"), crs = 31982, remove = FALSE)

# Carregar e ajustar Shapefile
matriz_shp <- st_read("GeoAmbFac_2022_img21_Amplo_v03.shp") %>%
  st_transform(st_crs(parcelas)) %>%
  mutate(GeoF_2021 = toupper(trimws(GeoF_2021)))
buffer_10m <- st_buffer(parcelas, dist = 10)
intersecao <- st_intersection(buffer_10m, matriz_shp)
# Calcular porcentagem de cada fisionomia na vizinhança de 10m
matriz_larga <- intersecao %>%
  mutate(area_pedaco = st_area(geometry)) %>%
  st_drop_geometry() %>%
  group_by(X, Y, GeoF_2021) %>% 
  summarise(area_total_fisia = sum(area_pedaco), .groups = "drop_last") %>%
  mutate(porcentagem = as.numeric((area_total_fisia / sum(area_total_fisia)) * 100)) %>%
  ungroup() %>%
  dplyr::select(X, Y, GeoF_2021, porcentagem) %>%
  pivot_wider(names_from = GeoF_2021, values_from = porcentagem, values_fill = 0)
colnames(matriz_larga) <- make.names(colnames(matriz_larga))
# Unindo a matriz calculada (10m) com os dados originais de abundância
Ipom_Completo <- IpomoeaAbund %>%
  left_join(matriz_larga, by = c("X", "Y")) %>%
  # Renomeando as colunas principais para nomes simples
  rename(MatrizLaj  = LAJEDO,
         MatrizVrup = VEGETAÇÃO.RUPESTRE.ABERTA,
         MatrizVArb = VEGETAÇÃO.RUPESTRE.ARBUSTIVA)
# O coeficiente (Estimate) indicará qual tem maior efeito na abundância.
# Hipótese 1: O Lajedo (MatrizLaj) é o principal driver
mod_lajedo <- glmmTMB(AbundTot ~ MatrizLaj + offset(log(AreaTot)), 
                      data = Ipom_Completo, ziformula = ~ Geofaceis, family = nbinom2)
# Hipótese 2: A Vegetação Rupestre Aberta (MatrizVrup) no entorno é o driver
mod_vrup <- glmmTMB(AbundTot ~ MatrizVrup + offset(log(AreaTot)), 
                    data = Ipom_Completo, ziformula = ~ Geofaceis, family = nbinom2)
# Hipótese 3: A Vegetação Rupestre Arbustiva (MatrizVArb) no entorno é o driver
mod_varb <- glmmTMB(AbundTot ~ MatrizVArb + offset(log(AreaTot)), 
                    data = Ipom_Completo, ziformula = ~ Geofaceis, family = nbinom2)
# Hipótese 4: Modelo Nulo (A vizinhança não importa, apenas a Geoface interna)
mod_nulo <- glmmTMB(AbundTot ~ 1 + offset(log(AreaTot)), 
                    data = Ipom_Completo, ziformula = ~ 1, family = nbinom2)
mod_comparativo<-AICtab(mod_lajedo, mod_vrup, mod_varb, mod_nulo, weights = TRUE)

#predição fixada em Geofácies "VEGETAÇÃO RUPESTRE ABERTA"
# Isso simula o efeito do entorno APENAS para quem está dentro desse habitat
pred_perda <- ggpredict(mod_lajedo,terms = "MatrizLaj [0:100]", 
  condition = c(Geofaceis = "VEGETAÇÃO RUPESTRE ABERTA",MatrizVArb = 0,AreaTot = 10))

library(ggplot2)

# 1. Criar o dataframe com os valores EXATOS da sua tabela
df_simulacao_exato <- data.frame(
  MatrizLaj = c(0, 13, 25, 37, 50, 63, 75, 100),
  predicted = c(2.85, 3.15, 3.46, 3.79, 4.20, 4.64, 5.09, 6.18),
  conf.low  = c(2.63, 2.91, 3.17, 3.42, 3.69, 3.97, 4.25, 4.86),
  conf.high = c(3.09, 3.41, 3.77, 4.21, 4.77, 5.42, 6.10, 7.85)
) %>% 
  mutate(x_perdido = 100 - MatrizLaj) # Transforma presença em perda

# 2. Gerar o gráfico
ggplot(df_simulacao_exato, aes(x = x_perdido, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  labs(
    title = "Impacto da Perda de Matriz (Lajedo)",
    subtitle = "Abundância predita em Vegetação Rupestre Aberta",
    x = "% de Perda da Matriz de Lajedo (Entorno 10m)",
    y = "Abundância Predita (Contagem)"
  ) +
  theme_light() +
  scale_x_continuous(breaks = seq(0, 100, 10))




#====================================================================================================================
#3)O tipo de geoface determina efeito denso-dependência?
#Sim, Canga aberta é a área fonte
#====================================================================================================================
# Garantir que os pontos usem o mesmo CRS do Shapefile
Ipomsf <- st_as_sf(IpomoeaAbund, coords = c("X", "Y"), crs = st_crs(SnShape))
# Criar lista de possíveis fontes
fontes <- list(Lajedo = SnShape %>% filter(toupper(trimws(GeoF_2021)) == "LAJEDO"),
  Vrup_Aberta = SnShape %>% filter(toupper(trimws(GeoF_2021)) == "VEGETAÇÃO RUPESTRE ABERTA"),
  Vrup_Arbustiva = SnShape %>% filter(toupper(trimws(GeoF_2021)) == "VEGETAÇÃO RUPESTRE ARBUSTIVA"),
  Mata_Baixa = SnShape %>% filter(toupper(trimws(GeoF_2021)) == "MATA BAIXA"))
# Loop para calcular a distância mínima de cada plot até a borda de cada fisionomia
for(nome in names(fontes)) {
# Calcula a distância
  dist_temp <- st_distance(Ipomsf, fontes[[nome]])
  # Define o nome da coluna
  col_name <- paste0("Dist_", nome)
  # CORREÇÃO: Salva no dataframe que EXISTE (IpomoeaAbund)
 IpomoeaAbund[[col_name]] <- as.numeric(apply(dist_temp, 1, min))
  # Cria a versão padronizada (Z-score)
  IpomoeaAbund[[paste0(col_name, "_z")]] <- as.numeric(scale(IpomoeaAbund[[col_name]]))
}
# Modelando a abundância em função da distância de cada fonte
Lajedo<- glmmTMB(AbundTot ~ Dist_Lajedo_z + offset(log(AreaTot)), data =  IpomoeaAbund, ziformula = ~ Geofaceis, family = nbinom2)
Aberta<- glmmTMB(AbundTot ~ Dist_Vrup_Aberta_z + offset(log(AreaTot)), data = IpomoeaAbund, ziformula = ~ Geofaceis, family = nbinom2)
Arbustiva<- glmmTMB(AbundTot ~ Dist_Vrup_Arbustiva_z + offset(log(AreaTot)), data =  IpomoeaAbund, ziformula = ~ Geofaceis, family = nbinom2)
Mata<- glmmTMB(AbundTot ~ Dist_Mata_Baixa_z + offset(log(AreaTot)), data = IpomoeaAbund, ziformula = ~ Geofaceis, family = nbinom2)
ranking_modelos <- AICtab(Lajedo, Aberta, Arbustiva, Mata, weights = TRUE, base = TRUE)
print(ranking_modelos)
# Verificando se o coeficiente é negativo (indica decaimento com a distância)
# Mata Baixa funciona como um núcleo ou "área fonte". À medida que se afasta, a densidade de plantas cai.

summary(Mata)

mu_dist <- mean( IpomoeaAbund$Dist_Vrup_Aberta, na.rm = TRUE)
sd_dist <- sd( IpomoeaAbund$Dist_Vrup_Aberta, na.rm = TRUE)

preditAberta <- ggpredict(Aberta, terms = "Dist_Vrup_Aberta_z [all]")

# Converter escala Z de volta para Metros reais
preditAberta$x_met <- (preditAberta$x * sd_dist) + mu_dist

# Plotagem
ggplot(preditAberta, aes(x = x_met, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "darkred", alpha = 0.2) +
  geom_line(color = "darkred", size = 1.2) +
  geom_rug(data = IpomoeaAbund, aes(x = Dist_Vrup_Aberta), inherit.aes = FALSE, alpha = 0.3) +
  labs(
    title = "Efeito da Distância da Vegetação Rupestre Aberta",
    subtitle = "Hipótese de Área Fonte",
    x = "Distância até a mancha mais próxima (m)", 
    y = "Abundância Predita de Ipomoea"
  ) +
  theme_classic(base_size = 14)












