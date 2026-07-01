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
library(multcompView)
library(openxlsx)
library(pscl)
library(sandwich)
library(sf)
library(sjPlot)
library(terra)
library(tidyverse)

rm(list=ls()) #Clean the workspace



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
#Campo graminoso é um filtro para espécie, quando ocorre há uma alta densidade
#Pois ocorre nos gaps das gramineas
#Lajedo e cangas é abundante porém menos densa, seria devido ao tamanho da área?
#==========================================================================================================

modAbundZI<- glmmTMB(N_Ind_Tota ~ Geofacie*PLATO ,data = IpomoeaAbund,
                      ziformula = ~0,family = nbinom2)

# Warning: dropping columns from rank-deficient conditional model: GeofacieCampo graminoso:PLATON2 indica
# Indica que Geofacie * PLATO é a melhor solução visto que Geofacie não são comparáveis entre todos os Platôs e vice-versa!


performance::check_zeroinflation(modAbundZI) #Não parece necessitar de zero-inflado


library(MuMIn)
options(na.action = "na.fail") # Necessário para o dredge

dredge_model<-dredge(modAbundZI)


best_model<-get.models(dredge_model, subset = 1)[[1]] # pega o melhor modelo (no caso é o modelo completo)

best_model%>%summary()
modAbundZI%>%summary()

plot(ggeffects::predict_response(best_model, terms = c("Geofacie", "PLATO")))+
  theme_classic(base_size=26)+theme(legend.position="bottom")


# Agora retirando N2 e N3 pois a falta de lajedos e campos graminosos tornam esses platôs pouco comparáveis
modAbundZInoN23<- glmmTMB(N_Ind_Tota ~ Geofacie*PLATO ,data = IpomoeaAbund_noN23,
                      ziformula = ~0,family = nbinom2)

performance::check_zeroinflation(modAbundZInoN23)

summary(modAbundZInoN23)

plot(ggeffects::predict_response(modAbundZInoN23, terms = c("Geofacie", "PLATO")))+
  theme_classic(base_size=26)+theme(legend.position="bottom")



IpomoeaAbund%>%glimpse()




# FIM DA AJUDA ==========================================================================================================


