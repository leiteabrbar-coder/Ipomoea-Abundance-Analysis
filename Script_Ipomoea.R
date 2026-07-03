#'#################################################################
#' Script para análise da abundância de Ipomoea 
#' em diferentes geofácies
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
  st_drop_geometry()%>%  #Necessário para performance, transforma o shp em dataframe
# Arruma diferentes grafias para Vegetação rupestre aberta
# Junta Vegetação rupestre aberta e arbustiva em um único grupo de 
  mutate(Geofacie=
 ifelse(Geofacie %in% 
        c("vegetação rupestre aberta","Vegetação Rupestre Aberta","Vegetação Rupestre Arbustiva"), #Se a geofacie for uma dessas, então...
          "Vegetação Rupestre", Geofacie))%>%                                                       # agrupa em "Vegetação Rupestre"
      #filtra campo brejoso, pois não é representativo em N2 e N3
  filter(Geofacie != "Campo Brejoso")
 

IpomoeaAbund$Geofacie%>%unique()

# Referências para análises
IpomoeaAbund$Geofacie<- relevel(as.factor(IpomoeaAbund$Geofacie), ref = "Vegetação Rupestre") #Define vegetação rupestre como referência
IpomoeaAbund$PLATO<- relevel(as.factor(IpomoeaAbund$PLATO), ref = "N1") #Define N1 como referência


# Verifica distribuição das manchas de habitat entre os platôs
# Note grande desigualdade marcada pelo baixo nº de lajedos e campos graminosos em N2 e N3, o que pode afetar a análise de abundância
IpomoeaAbund%>%dplyr::select(Geofacie,PLATO)%>%table()

IpomoeaAbund_noN23<-IpomoeaAbund%>%filter(!(PLATO %in% c("N2","N3"))) #Remove N2 e N3 para conseguir analisar lajedo e campo graminoso de forma adequada.

#########################################################################################
#1)A geoface determina a abundância?
#########################################################################################
#À partir daqui a modelagem está perfeita!
modAbundZI<- glmmTMB(N_Ind_Tota ~ Geofacie*PLATO ,data = IpomoeaAbund,
                     ziformula = ~Geofacie,family = nbinom2)
modAbundZI$fit$convergence # Deve ser 0, Perfect!

performance::check_zeroinflation(modAbundZI) #Ratio: 1.00 Perfect! 

summary(modAbundZI)

options(na.action = "na.fail") # Necessário para o dredge

dredge_model<-dredge(modAbundZI) # Função maravilhosa, já faz as combinações de variáveis, testando assim diferentes modelos


best_model<-get.models(dredge_model, subset = 1)[[1]] # pega o melhor modelo dentre todos aqueles gerados em dredge_model 
#(no caso é o modelo completo)

best_model%>%summary()
#modAbundZI%>%summary() Não vejo necessidade dessa linha, uma vez que já teremos e sabemos qual o melhor modelo

#Os nomes de legendas etc... deveram ser definidos para o painel do simpósio
plot(ggeffects::predict_response(best_model, terms = c("Geofacie", "PLATO")))+
  theme_classic(base_size=26)+theme(legend.position="bottom")

# Eu prefiro esse padrão de gráfico


em_values<- emmeans(modAbundZI, ~ Geofacie | PLATO, type = "response")
letras_tukey <- cld(em_values, Letters = letters, adjust = "tukey") %>% 
  as.data.frame()
col_resposta <- intersect(c("response", "rate"), colnames(letras_tukey))
col_erro     <- intersect(c("SE", "std.error"), colnames(letras_tukey))

letras_tukey$.group <- trimws(letras_tukey$.group)

# 4. Filtragem dos dados (remove o erro gigante do Lajedo no Platô 4 e NAs)
dados_grafico <- letras_tukey %>% 
  filter(!is.na(.data[[col_resposta]])) %>% 
  filter(.data[[col_erro]] < 100)


plot_letras <- ggplot(dados_grafico, aes(x = Geofacie, y = .data[[col_resposta]], fill = Geofacie)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, alpha = 0.9) +
  geom_errorbar(aes(
    ymin = pmax(0, .data[[col_resposta]] - .data[[col_erro]]), 
    ymax = .data[[col_resposta]] + .data[[col_erro]]
  ), width = 0.2, color = "black", linewidth = 0.6) +
  geom_text(aes(
    y = .data[[col_resposta]] + .data[[col_erro]], 
    label = .group
  ), vjust = -0.5, size = 4, fontface = "bold") +
  facet_wrap(~ PLATO, scales = "free_y", ncol = 2) +
 
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

print(plot_letras)

# ggsave("Painel_Abundancia_Platos.png", plot = plot_final, width = 8, height = 7, dpi = 300)

#################################################################################################
# Agora retirando N2 e N3 pois a falta de lajedos e campos graminosos tornam esses platôs pouco comparáveis
#1)A geoface determina a abundância?
#2)Há diferença na probabilidade ocorrência entre as geofácies?
############################################################################################################
IpomoeaAbund_noN23<-IpomoeaAbund %>% 
  filter(PLATO %in% c("N1", "N4/N5"))
modAbundZInoN23<- glmmTMB(N_Ind_Tota ~ Geofacie*PLATO ,data = IpomoeaAbund_noN23,
                      ziformula = ~Geofacie,family = nbinom2)

performance::check_zeroinflation(modAbundZInoN23)
modAbundZInoN23$fit$convergence
summary(modAbundZInoN23)
emm_abund <- emmeans(modAbundZInoN23, specs = ~ Geofacie * PLATO, type = "response")
summary(pairs(emm_abund, by = "PLATO"), adjust = "tukey")
cld_abund <- cld(emm_abund, by = "PLATO", adjust = "tukey", Letters = letters)

# Abaixo serão construídos os gráficos
#ATENÇÃO: A barra de erro padrão infinita em Lajedo N4/N5 é porque os 21 plots são todos 0, ou seja, em nenhum há Ipomoea
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

# 1. Ajuste fino e isolamento completo dos dados para o gráfico
plotAbund<- as.data.frame(resultado_pos_hoc) %>% 
  mutate(
    letra = trimws(.group),
    # Garante limites válidos para os intervalos de confiança
    LCL_adjusted = ifelse(is.na(asymp.LCL) | asymp.LCL < 0, 0, asymp.LCL),
    UCL_adjusted = ifelse(is.na(asymp.UCL) | asymp.UCL > 10000 | asymp.UCL == Inf, response, asymp.UCL),
    
    # Recalcula o intervalo de forma limpa para o Campo graminoso em N4/N5 usando o SE
    UCL_adjusted= ifelse(Geofacie == "Campo graminoso" & PLATO == "N4/N5", response + (1.96 * SE), UCL_adjusted),
    LCL_adjusted= ifelse(Geofacie == "Campo graminoso" & PLATO == "N4/N5", max(0, response - (1.96 * SE)), LCL_adjusted),
    
    # Define a posição vertical segura de cada letra para evitar sobreposição
    posicao_letra = ifelse(Geofacie == "Lajedo" & PLATO == "N4/N5", response + 1.2, UCL_adjusted + 1.2)
  )

# 2. Construção do Gráfico com sintaxe atualizada (linewidth)
ggplot(plotAbund, aes(x = Geofacie, y = response, color = factor(PLATO), group = factor(PLATO))) +
  geom_errorbar(aes(ymin = LCL_adjusted, ymax = UCL_adjusted), 
                width = 0.12, position = position_dodge(0.4), linewidth = 0.8) +
  geom_point(position = position_dodge(0.4), size = 4.5) +
  geom_text(aes(y = posicao_letra, label = letra), 
            vjust = 0, size = 5.5, fontface = "bold", 
            position = position_dodge(0.4), show.legend = FALSE) + 
  scale_color_manual(values = c("N1" = "#f3716d", "N4/N5" = "#1ebbc7"), name = "Platô") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 20, by = 10)) +
  labs(x = "Geofacie",y = "Predicted abund") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(color = "black", size = 12),
     axis.text.y = element_text(color = "black", size = 12),
     axis.line = element_line(color = "black", linewidth = 0.6),
     legend.position = "bottom")
ggsave("PlotAbun.png", width = 8, height = 7, dpi = 300)

#Gráfico de Probabilidade de ocorrência
# Campo graminoso é um filtro para a ocorrência de Ipomoea
# Gerando as médias para a inflação de zeros


# 1. Extrai as médias da inflação de zeros (component = "zi") na escala de probabilidade
emm_zi <- emmeans(modAbundZInoN23, specs = ~ Geofacie, component = "zi", type = "response")

# 2. Gera as letras do teste de Tukey para a ausência
cld_zi <- cld(emm_zi, adjust = "tukey", Letters = letters) %>% 
  as.data.frame() %>% 
  mutate(letra = trimws(.group))

# 3. Transforma chance de ausência em PROBABILIDADE DE OCORRÊNCIA (presença)
# Correção: Usando os nomes diretos das colunas sem a função any_of()
tabela_ocorrencia <- cld_zi %>% 
  mutate(
    Prob_Ocorrencia = 1 - response,
    LCL_limpo = 1 - asymp.UCL,  # Limite inferior da presença usa o teto da ausência
    UCL_limpo = 1 - asymp.LCL   # Limite superior da presença usa o piso da ausência
  ) %>% 
  dplyr::select(Geofacie, Prob_Ocorrencia, LCL_limpo, UCL_limpo, letra)

# 4. Printa o resultado final estruturado na tela
print(tabela_ocorrencia)

# 5. Printa os p-valores ajustados de Tukey para cada par
summary(pairs(emm_zi), adjust = "tukey")
cld_zi_preparado <- cld_zi %>% 
  mutate(
    Prob_Ocorrencia = 1 - response,
    IC_inferior     = 1 - asymp.UCL,  # Limite inferior da presença
    IC_superior     = 1 - asymp.LCL   # Limite superior da presença
  ) %>% 
  mutate(
    # Corrige matematicamente caso os erros passem dos limites lógicos de 0% e 100%
    IC_inferior = ifelse(IC_inferior < 0, 0, IC_inferior),
    IC_superior = ifelse(IC_superior > 1, 1, IC_superior),
    .group = trimws(.group) # Garante que as letras não tenham espaços
  )

# Gráfico de Probabilidade de Ocorrência
ggplot(cld_zi_preparado, aes(x = Geofacie, y = Prob_Ocorrencia, fill = Geofacie)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = IC_inferior, ymax = IC_superior), width = 0.2) +
  geom_text(aes(y = IC_superior, label = .group), vjust = -0.5, size = 5) + # Letras no topo
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) + # Transforma o eixo Y em % (0 a 100%)
  labs(
    x = "Geofácies",
    y = "Probabilidade de Ocorrência (%)",
    title = "Probabilidade de Ocorrência de Ipomoea por Geofácies"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")



# 1. Encontra o valor máximo para definir o último corte
max_val <- max(IpomoeaAbund$N_Ind_Tota, na.rm = TRUE)

# 2. Cria sequências de 10 em 10 a partir do 0 até ultrapassar o máximo
cortes_sequencia <- seq(0, ceiling(max_val / 10) * 10, by = 10)

# 3. Une o -1 (para isolar o zero) com a sequência de 10 em 10
todos_cortes <- c(-1, cortes_sequencia)

# 4. Cria os rótulos de forma automatizada (ex: "0", "1-10", "11-20"...)
rotulos <- c("0", paste0(cortes_sequencia[-length(cortes_sequencia)] + 1, "-", cortes_sequencia[-1]))

# 5. Aplica os cortes na base de dados
IpomoeaAbund$Faixa_10 <- cut(IpomoeaAbund$N_Ind_Tota, 
                             breaks = todos_cortes, 
                             labels = rotulos, 
                             include.lowest = TRUE)

# 6. Gera o gráfico finalizado com TODOS os números aparecendo
ggplot(IpomoeaAbund, aes(x = Faixa_10, fill = Faixa_10 == "0")) +
  geom_bar(width = 1.0, color = "white", show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#e63946", "FALSE" = "#4a90e2")) +
  scale_y_continuous(
    breaks = seq(0, 1000, by = 50), 
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(x= "Total de indivíduos",y= "Total de parcelas") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 9))+
 stat_count(geom = "text", aes(label = after_stat(count)), 
             vjust = -0.5, fontface = "bold", size = 3.5)






















