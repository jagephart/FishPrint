#a program to read N and P values of fish and feed and compute discharge to the environment
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
  data.table(read.csv("AFCD_live.csv", stringsAsFactors = FALSE))
fishNutrition <-fishNutrition %>%filter(Processing=='r') %>%  #filter out all non raw observations (e.g. dried, cooked)
  slice(c(1:1928))    #cut rows without names of species 
fishNutrition1<-fishNutrition %>% mutate(N_C=fishN_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                         ,P_C=fishP_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                         ,N_fat=fishN_viaFat(Water,Fat.total)
                                         ,P_fat=fishP_viaFat(Water,Fat.total)
                                         ,N_built_in=Nitrogen.total*100/(100-Water),P_built_in=Phosphorus/1000*100/(100-Water)) #conversion to units of percentages in DM. Nitrogen.total is in units of g/100 edible gram; Phosphorous in the dataset is in mg/100 edible gram

fishNutrition2<-fishNutrition1 %>% group_by(species) %>% summarise(P_byC=mean(P_C,na.rm = TRUE),N_byC=mean(N_C,na.rm = TRUE),
                                                      P_fat_1=mean(P_fat,na.rm = TRUE),N_fat_1=mean(N_fat,na.rm = TRUE),
                                                      P_built_in1=mean(P_built_in,na.rm = TRUE),N_built_in1=mean(N_built_in,na.rm = TRUE)) %>%                                 
                                                      rowwise() %>%                                                    
                                                     mutate(N_avg=mean(c(N_byC,N_fat_1,N_built_in1),na.rm = TRUE))%>%
                                                      mutate(P_avg=mean(c(P_byC,P_fat_1,P_built_in1),na.rm = TRUE))%>%
                                                      select(c("species","N_avg","P_avg"))

#read LCA db  
db_full<-data.table(read.csv("LCA_compiled_20201109.csv", stringsAsFactors = FALSE))


#merge with with full LCA
db_model_full<-merge(x=db_full, y=fishNutrition2, by.x = "Species.scientific.name",by.y = "species" ,all.x=TRUE)
write.xlsx(db_model_full,'LCA_full_data_compiled.xlsx')



#------------------------------------------------------------------------------------
#------------------------------------------functions---------------------------------
#------------------------------------------------------------------------------------

#first method: estimate from C content, Fig 1 of the above reference, linear regression values. Assuming metabolizable energy density is equal throughout whole fish

fishN_via_C <- function(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal){
  DM = (100 - Water) / 100     #convert to dry matter fraction
  
  #In aquacultured fish, whole-body C content was estimated from energy content using a converting factor of 1 g C = 11.4 kcal
  fish_C_percent <- 1/DM*Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal/11.4
  N_slope <- -0.175   #from fig 1
  N_intercept <- 18.506  #intercept (personal communication with authors)
  fish_N <- fish_C_percent * N_slope + N_intercept #in % of DM
  return(fish_N)}

fishP_via_C <- function(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal){
  
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

