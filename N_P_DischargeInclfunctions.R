#a program to read N and P values of fish and feed anc compuate discharge to the environment
rm(list = ls(all = TRUE))    #delete variables

library("rfishbase")  #package to use fishbase in r
library("openxlsx")
library("rlang")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("data.table")
library("stringr")
library("openxlsx")



#-----------------------read dataset

#calculate C,N,P content of fish (based on energy and fat content ) and compute discharge based on Czamanski et al 2011, Marine Biology,
#Carbon, nitrogen and phosphorus elemental stoichiometry in aquacultured and wild-caught fish and consequences
#for pelagic nutrient dynamics
#And also a simplified NPZ method approach

#N and P content of fish in % of DM
#read Zach's data and make calculations using built functions
fishNutrition <-
  data.table(read.csv("FishGenus.csv", stringsAsFactors = FALSE))
fishNutrition <-fishNutrition %>%filter(Preparation=='raw')   #filter out all non raw observations (e.g. dried, cooked)

fishNutrition1<-fishNutrition %>% mutate(N=fishN(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                         ,P=fishP(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                        ,N_fat=fishN_viaFat(Water,Fat.total)
                                        ,P_fat=fishP_viaFat(Water,Fat.total)
                                        ,N_built_in=Nitrogen.total*100/(100-Water),P_built_in=Phosphorus/1000*100/(100-Water)) #conversion to units of percentages in DM. Nitrogen.total is in units of g/100 edible gram; Phosphorous in the dataset is in mg/100 edible gram


#N and P content of feed in % of DM (read from feed composition tables, USDA 1982)

feedNutrition <-read.xlsx("United-States-Canadian-Tables-of-Feed-1982-pages-68-921-with_CrudeProtein.xlsx",startRow = 2)
feedNutrition$`Phosphorus.(%)`<-as.numeric(feedNutrition$`Phosphorus.(%)`)
feedNutrition<-feedNutrition %>% mutate( Nitrogen = `Crude.protein.(%)`/ 6.25)
feedNutrition$Nitrogen=as.numeric(feedNutrition$Nitrogen)   

#define categorial feeds: soy, other products, animal and FM&O (take only 100% DM)
#soy: Entry numbers: 601-620
s=c(602,604,608,610,612,614,616,618,620)

#animal by products, Entry number: meat byproducts 385-392, poultry 479-481 (take only 100% DM)
a=c(386,388,390,392,480,482,484)

#other crops: cassava 141-142, peanut extr. 461-464, linseed extr. 345-348 , 
#corn gluten meal 223-228, pea (not protein concentrate) 455-456, rape 493-496, sunflower oil and meal 633-636, 
#wheat byproduct 683-688, Sorghum grain 573-580, triticale 667-668, Navy beans (instead of faba) 57-58,
cr=c(142,462,464,346,348,224,226,228,456,494,496,634,636,684,686,688,578,580,668,58)

#fishmeal and oil
# fish 301-304, fish alewife meal 305-306,anchovy meal 307-308, catfish meal 313-314, herring meal 317-318, 
# mackerel meal 322, manhaden meal 323-324, redfish meal 325-326, salmon meal 329-330, sardine meal 333-334
#tuna meal 337-338, white fish 341-342
f=c(302,304,306,308,314,318,322,324,326,330,334,338,342)

# calculate weighted average N and P content of feed
soy_N=mean(feedNutrition$Nitrogen[which(feedNutrition$Entry.Number %in% s)],na.rm=TRUE)   
soy_P=mean(feedNutrition$`Phosphorus.(%)`[which(feedNutrition$Entry.Number %in% s)],na.rm=TRUE) 
animalbyproduct_N=mean(feedNutrition$Nitrogen[which(feedNutrition$Entry.Number %in% a)],na.rm=TRUE)
animalbyproduct_P=mean(feedNutrition$`Phosphorus.(%)`[which(feedNutrition$Entry.Number %in% a)],na.rm=TRUE)
othercrop_N=mean(feedNutrition$Nitrogen[which(feedNutrition$Entry.Number %in% cr)],na.rm=TRUE)
othercrop_P=mean(feedNutrition$`Phosphorus.(%)`[which(feedNutrition$Entry.Number %in% cr)],na.rm=TRUE)
FMFO_N=mean(feedNutrition$Nitrogen[which(feedNutrition$Entry.Number %in% f)],na.rm=TRUE)
FMFO_P=mean(feedNutrition$`Phosphorus.(%)`[which(feedNutrition$Entry.Number %in% f)],na.rm=TRUE)

  
  
#Compute emission of N and P to the environment based on feed and fish N and P rations (dummy example)

R<-fishNutrition1[which(fishNutrition1$Scientific.Name=="Salmo salar"),]

FCR=1.3   #should be in DM
fish_N=mean(R$N,na.rm=TRUE)
fish_P=mean(R$P,na.rm=TRUE)

feed_N=0.3*soy_N+0.3*FMFO_N+0.2*othercrop_N+animalbyproduct_N*0.2# in % DM
feed_P=0.3*soy_P+0.3*FMFO_P+0.2*othercrop_P+animalbyproduct_P*0.2# in % DM

#N, P discharge based on Czamanski et al 2011, Marine Biology
RT<-N_P_discharge_model(fish_N,fish_P,feed_N,feed_P,FCR)  #fish_N,fish_P,feed_N,feed_P,FCR = N and P in % of dry matter, FCR in kg feed/1 kg
#RT=output (E_N,F_N,E_P,F_P) results in kg compared to FCR
#E are discharge to the enviroment and F are feces based on the above reference that implements the Sterner 1990 model
lost_to_environment_percentage_N=(RT[1]+RT[2])/(feed_N/100*FCR)*100
lost_to_environment_percentage_P=(RT[3]+RT[4])/(feed_P/100*FCR)*100

#simplified NPZ discharge model
RT_simp<-N_P_discharge_simple_model(fish_N,fish_P,feed_N,feed_P,FCR)  #fish_N,fish_P,feed_N,feed_P,FCR = N and P in % of dry matter, FCR in kg feed/1 kg
#RT=output (N_discharge,P discharge) results in kg (in accordance to FCR kg to kg) 

lost_to_environment_percentage_N_simp=(RT_simp[1])/(feed_N/100*FCR)*100   #percentage from feed input 
lost_to_environment_percentage_P_simp=(RT_simp[2])/(feed_P/100*FCR)*100   # "

#------------------------------------------------------------------------------------
#------------------------------------------functions---------------------------------
#------------------------------------------------------------------------------------

#first method: estimate from C content, Fig 1 of the above reference, linear regression values. Assuming metabolizable energy density is equal throughout whole fish

fishN <- function(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal){
  DM = (100 - Water) / 100     #convert to dry matter fraction

#In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
fish_C_percent <- 1/DM*Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal/11.4
N_slope <- -0.175   #from fig 1
N_intercept <- 18.506  #intercept (personal communication with authors)
fish_N <- fish_C_percent * N_slope + N_intercept #in % of DM
return(fish_N)}

fishP <- function(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal){
  
DM = (100 - Water) / 100     #convert to dry matter fraction

  #In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
  fish_C_percent <- 1/DM*Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal/11.4
  P_slope <- -0.083  #from fig 1
  P_intercept <- 6.311  #intercept (personal communication with authors)
  fish_P <- fish_C_percent * P_slope + P_intercept #in % of DM
  return(fish_P)}

#second method: estimate from total lipids, figure 2, assuming lipid density is equal throughout whole body

fishN_viaFat <- function(Water,Fat.total){
  DM = (100 - Water) / 100     #Dry Matter percentage
#Fat.total in units of g to 100 g
fat = Fat.total * 1  / DM       #conversion to DM
  
#fat percentage in whole body in % of DM. 
fish_C_via_fat = fat * 0.31 + 38  #linear regression from fig 2
fish_N_via_fat = fat * -0.08 + 12 #linear regression from fig 2
return(fish_N_via_fat)}

fishP_viaFat <- function(Water,Fat.total){
  DM = (100 - Water) / 100     #Dry Matter percentage

  #Fat.total in units of g to 100 g
  fat = Fat.total * 1  / DM      #conversion to DM
  
  #fat percentage in whole body in % of DM. Assuming Fat is contained only (mostly) in the edible portion.
  fish_C_via_fat = fat * 0.31 + 38  #linear regression from fig 2
  fish_P_via_fat = fat * -0.04 + 3.2 #linear regression from fig 2
  return(fish_P_via_fat)}

############################################################


N_P_discharge_model <- function(fish_N_percent,fish_P_percent,feed_N_percent,feed_P_percent,FCR){
  #G=I-[E+F]   from the above reference "Methods", where G is growth, I ingestion, E emission to enviroment (resupply) and F is feces
  #FCR = I/G  FCR assumes I=1 unit weight -->E+F=FCR-1
  #in terms of nutrient conservation: E_N+F_N=feed_N*(FCR-1)-fish_N*1; E_P+F_P=feed_P(FCR-1)-fish_P*1
  #convert N and P to molar weight 
  Lm=0.75;
  #browser()
  fish_N=fish_N_percent/14;feed_N=feed_N_percent/14; #stoichometric molar ratio
  fish_P=fish_P_percent/31;feed_P=feed_P_percent/31;
  beta_N<-0.91
  R_fish<-fish_N/fish_P;
  R_feed<-feed_N/feed_P;
  beta_P<- -0.19*R_fish/R_feed+0.8   #from fig 4 of reference
  if (R_feed*beta_N/beta_P<R_fish){
  alpha_N<-Lm;
  alpha_P<-Lm*R_feed*beta_N/beta_P/R_fish}
  else if (R_feed*beta_N/beta_P>R_fish){
    alpha_P<-Lm;
  alpha_N<-Lm*(R_fish)/(R_feed*beta_N/beta_P)}

  R_resupply<-((1-alpha_N)*beta_N)/((1-alpha_P)*beta_P)*R_feed #equation 7 from reference
  R_feces<-R_feed*(1-beta_N)/(1-beta_P)
  constant_N<-(feed_N_percent/100*(FCR)-fish_N_percent/100*1)/0.014; #conservation of nutrient: feed-fish_content=feces+resupply: units of molar weight: 0.014 kg/mole
  constant_P<-(feed_P_percent/100*(FCR)-fish_P_percent/100*1)/0.031; #conservation of nutrient: feed-fish_content=feces+resupply: units of molar weight: 0.031 kg/mole
  E_P<-(constant_N-R_feces*constant_P)/(R_resupply-R_feces);
  E_N<-R_resupply*E_P
  F_N<-constant_N-E_N;
  F_P<-constant_P-E_P;
  output<-c(E_N*0.014,F_N*0.014,E_P*0.031,F_P*0.031)
  return(output)
}
  
  
N_P_discharge_simple_model <- function(fish_N_percent,fish_P_percent,feed_N_percent,feed_P_percent,FCR){
  #simplified model of discharge, similar to many NPZ models 
  discharge_N<-(feed_N_percent/100*(FCR)-fish_N_percent/100*1); #conservation of nutrient: feed-fish_content=discharge to enviroment. units: weights as FCR
  discharge_P<-(feed_P_percent/100*(FCR)-fish_P_percent/100*1); #conservation of nutrient: feed-fish_content=discharge to enviromen.  weights as FCR 
  if(discharge_N<=0){discharge_N=0}else{discharge_N=discharge_N} #to avoid negative emissions
if(discharge_P<=0){discharge_P=0}else{discharge_P=discharge_P}
RT=c(discharge_N,discharge_P)  #in units of FCR e.g. kg N and P for 1 kg of live weight fish
 return(RT) }
