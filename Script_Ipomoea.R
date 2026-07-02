#'#################################################################
#' Script para análise da abundância de Ipomoea 
#' em diferentes geofácies e a influência da matriz circundante. 
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

#rm(list=ls()) #Clean the workspace



SnShape<-st_read("Data/PRJ112_DadosIpomoea_Densidades_v00.shp")
geof<-st_read("Data/Serra_Norte_2024_wLoc.shp")

SnShape$Geofacie%>%unique()

  
  
#Cria dados de abundância
IpomoeaAbund<- SnShape%>%
  st_drop_geometry()%>%  #Necessário para performance
# Arruma diferentes grafias para Vegetação rupestre aberta
# Junta Vegetação rupestre aberta e arbustiva em um único grupo de 
  mutate(Geofacie=
 ifelse(Geofacie %in% 
        c("vegetação rupestre aberta","Vegetação Rupestre Aberta","Vegetação Rupestre Arbustiva"), #Se a geofacie for uma dessas, então...
          "Vegetação Rupestre", Geofacie))%>%                                                       # agrupa em "Vegetação Rupestre"
      #filtra campo brejoso
  filter(Geofacie != "Campo Brejoso")
 

IpomoeaAbund$Geofacie%>%unique()

# Referências para análises
IpomoeaAbund$Geofacie<- relevel(as.factor(IpomoeaAbund$Geofacie), ref = "Vegetação Rupestre") #Define vegetação rupestre como referência
IpomoeaAbund$PLATO<- relevel(as.factor(IpomoeaAbund$PLATO), ref = "N1") #Define N1 como referência


# Verifica distribuição das manchas de habitat entre os platôs
# Note grande desigualdade marcada pelo baixo nº de lajedos e campos graminosos em N2 e N3, o que pode afetar a análise de abundância
IpomoeaAbund%>%dplyr::select(Geofacie,PLATO)%>%table()

IpomoeaAbund_noN23<-IpomoeaAbund%>%filter(!(PLATO %in% c("N2","N3"))) #Remove N2 e N3 para conseguir analisar lajedo e campo graminoso de forma adequada.

#==========================================================================================================
#1)A geoface determina a abundância?
#==========================================================================================================
# Esse modelo está com uma convergência ruim
# Portanto, foi descartado
#modAbundZI<- glmmTMB(N_Ind_Tota ~ Geofacie*PLATO ,data = IpomoeaAbund,
#                            ziformula = ~0,family = nbinom2)
#performance::check_zeroinflation(modAbundZI) #Ratio: 0.89 muito alta, é necessário zero-inflated
#ATENÇÂO:Haverá esse aviso!
# Warning: dropping columns from rank-deficient conditional model: Geofacie Campo graminoso:PLATON2 indica
# Indica que Geofacie * PLATO é a melhor solução visto que Geofacie não são comparáveis entre todos os Platôs e vice-versa!


modAbundZI<- glmmTMB(N_Ind_Tota ~ Geofacie*PLATO ,data = IpomoeaAbund,
                     ziformula = ~Geofacie,family = nbinom2)
modAbundZI$fit$convergence # Deve ser 0, Perfect!

performance::check_zeroinflation(modAbundZI) #Ratio: 1.00 Perfect! 

summary(modAbundZI)

options(na.action = "na.fail") # Necessário para o dredge

dredge_model<-dredge(modAbundZI) # Função maravilhosa, já faz as combinações de variáveis, testando assim diferentes modelos


best_model<-get.models(dredge_model, subset = 1)[[1]] # pega o melhor modelo dentre todos aqueles gerados em dredge_model (no caso é o modelo completo)

best_model%>%summary()
#modAbundZI%>%summary() Não vejo necessidade dessa linha, uma vez que já teremos e sabemos qual o melhor modelo

#Os nomes de legendas etc... deveram ser definidos para o painel do simpósio
plot(ggeffects::predict_response(best_model, terms = c("Geofacie", "PLATO")))+
  theme_classic(base_size=26)+theme(legend.position="bottom")

# Eu prefiro esse padrão de gráfico

# 1. Extrai as médias e erros
library(emmeans)
library(ggplot2)
library(dplyr)
library(multcomp)

# 1. Extrai as médias estimadas do modelo
em_regua <- emmeans(modAbundZI, ~ Geofacie | PLATO, type = "response")

# 2. Calcula as letras de Tukey (Compact Letter Display) DENTRO de cada Platô
letras_tukey <- cld(em_regua, Letters = letters, adjust = "tukey") %>% 
  as.data.frame()

# 3. Identificação dinâmica das colunas para evitar conflitos de versão
col_resposta <- intersect(c("response", "rate"), colnames(letras_tukey))
col_erro     <- intersect(c("SE", "std.error"), colnames(letras_tukey))

# Limpa espaços em branco criados pelo cld() na coluna de letras
letras_tukey$.group <- trimws(letras_tukey$.group)

# 4. Filtragem dos dados (remove o erro gigante do Lajedo no Platô 4 e NAs)
dados_grafico <- letras_tukey %>% 
  filter(!is.na(.data[[col_resposta]])) %>% 
  filter(.data[[col_erro]] < 100)

# 5. Criando o gráfico unificado com facetas e letras
plot_letras <- ggplot(dados_grafico, aes(x = Geofacie, y = .data[[col_resposta]], fill = Geofacie)) +
  # Barras de médias
  geom_bar(stat = "identity", color = "black", width = 0.6, alpha = 0.9) +
  # Barras de erro (garantindo mínimo zero)
  geom_errorbar(aes(
    ymin = pmax(0, .data[[col_resposta]] - .data[[col_erro]]), 
    ymax = .data[[col_resposta]] + .data[[col_erro]]
  ), width = 0.2, color = "black", linewidth = 0.6) +
  # ADICIONA AS LETRAS ACIMA DA BARRA DE ERRO
  geom_text(aes(
    y = .data[[col_resposta]] + .data[[col_erro]], 
    label = .group
  ), vjust = -0.5, size = 4, fontface = "bold") +
  # Divisão em 4 janelas (2x2)
  facet_wrap(~ PLATO, scales = "free_y", ncol = 2) +
  # Rótulos e estética
  labs(
    title = "Abundância de Ipomoea com Comparações de Tukey",
    subtitle = "Letras diferentes indicam diferença significativa (p < 0.05) dentro de cada platô",
    x = "Geofacie",
    y = "Número Médio de Indivíduos"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.8)

# Exibe o painel completo no RStudio
print(plot_letras)

# Para salvar o gráfico finalizado na sua pasta de trabalho:
# ggsave("Painel_Abundancia_Letras_Tukey.png", plot = plot_letras, width = 9, height = 8, dpi = 300)

# ggsave("Painel_Abundancia_Platos.png", plot = plot_final, width = 8, height = 7, dpi = 300)

#################################################################################################
# Agora retirando N2 e N3 pois a falta de lajedos e campos graminosos tornam esses platôs pouco comparáveis
############################################################################################################
IpomoeaAbund_noN23<-IpomoeaAbund %>% 
  filter(PLATO %in% c("N1", "N4/N5"))
modAbundZInoN23<- glmmTMB(N_Ind_Tota ~ Geofacie*PLATO ,data = IpomoeaAbund_noN23,
                      ziformula = ~Geofacie,family = nbinom2)

performance::check_zeroinflation(modAbundZInoN23)
modAbundZI$fit$convergence
summary(modAbundZInoN23)

# Test post-hoc pair to pair
medias <- emmeans(modAbundZInoN23, ~ Geofacie | PLATO, type = "response")

tukeyResults <- cld(medias, Letters = letters, adjust = "tukey")

print(tukeyResults)

plot(ggeffects::predict_response(modAbundZInoN23, terms = c("Geofacie", "PLATO")))+
  theme_classic(base_size=26)+theme(legend.position="bottom")

# Preparando os dados para o plot abaixo
pred_data <- predict_response(modAbundZInoN23, terms = c("Geofacie", "PLATO"))

medias <- emmeans(modAbundZInoN23, ~ Geofacie | PLATO, type = "response")
letras_data <- cld(medias, Letters = letters, adjust = "sidak") %>% 
  as.data.frame() %>% 
  mutate(.group = trimws(.group)) # Remove espaços extras das letras

plot_df <- as.data.frame(pred_data) %>% 
  rename(Geofacie = x, PLATO = group) %>% 
  left_join(letras_data, by = c("Geofacie", "PLATO"))

# Gostei desse tipo de gráfico também
# Preciso confirmar com os demais se será útil

# 1. Ajustando a tabela para remover o Infinito do Lajedo N4/N5
plot_df_limpo <- plot_df %>%
  mutate(
    conf.high = ifelse(is.infinite(conf.high), predicted, conf.high),
    conf.low = ifelse(conf.low < 0.001, 0, conf.low)
  )
teto_grafico <- max(plot_df_limpo$conf.high, na.rm = TRUE) + 3
# 2. Gerando o gráfico definitivo
ggplot(plot_df_limpo, aes(x = Geofacie, y = predicted, color = PLATO, group = PLATO)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.2, position = position_dodge(0.4), linewidth = 1) +
  geom_point(position = position_dodge(0.4), size = 5) +
  
  # Todas as letras posicionadas exatamente no topo (conf.high)
  geom_text(aes(y = conf.high, label = .group), 
            position = position_dodge(0.4), 
            vjust = -0.6,     # Empurra as letras levemente para cima das barras/pontos
            size = 6,         
            show.legend = FALSE) +
  
  # DEFINE O LIMITE MÁXIMO DO EIXO Y EM 25
  # O 'expand' garante que o título do eixo Y e as letras tenham espaço para respirar 
  scale_y_continuous(breaks = seq(0, teto_grafico + 5, by = 5))+
  theme_classic(base_size = 26) +
  theme(legend.position = "bottom",
        # Margens generosas nas laterais para nunca cortar os textos da imagem
        plot.margin = margin(t = 20, r = 20, b = 10, l = 30), 
        axis.text.x = element_text()) + 
  labs(x = "Geofacie", y = "Predicted abundance", color = "Platô")

IpomoeaAbund%>%glimpse()

#ATENÇÃO: A barra de erro padrão infinita em Lajedo N4/N5 é portque os 21 plots são todos 0, ou seja, em nehum á Ipomoea
IpomoeaAbund%>%dplyr::select(Geofacie,PLATO)%>%table()
IpomoeaAbund %>%
  filter(Geofacie == "Lajedo", PLATO == "N4/N5") %>%
  dplyr::select(N_Ind_Tota) %>%
  table()

 #N4/N5 têm uma alta porcentagem de zeros
tabela_zeros_N4N5 <- IpomoeaAbund %>%
  # Filtra apenas o platô de interesse
  filter(PLATO == "N4/N5") %>%
  # Agrupa por Geofacie
  group_by(Geofacie) %>%
  # Calcula as estatísticas de ausência
  summarise(
    Quantidade_Zeros = sum(N_Ind_Tota == 0),
    Total_Parcelas   = n(),
    Percentual_Zeros = (Quantidade_Zeros / Total_Parcelas) * 100
  )

# Visualiza o dataframe gerado
print(tabela_zeros_N4N5)

