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
fishNutrition <- data.table(read.csv("AFCD_live.csv", stringsAsFactors = FALSE))

#build taxa data 
fish_Nutrition_genus<-fishNutrition[,c("species","order","family","genus")] %>%distinct()


fishNutrition <-fishNutrition %>%filter(Processing=='r') %>%  #filter out all non raw observations (e.g. dried, cooked)
  slice(c(1:1928))    #cut rows without names of species 
fishNutrition1<-fishNutrition %>% mutate(N_C=fishN_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                         ,P_C=fishP_via_C(Water,Energy.total.metabolizable.calculated.from.the.energy.producing.food.components.original.as.from.source.kcal)
                                         ,N_fat=fishN_viaFat(Water,Fat.total)
                                         ,P_fat=fishP_viaFat(Water,Fat.total)
                                         ,N_built_in=Nitrogen.total*100/(100-Water),P_built_in=Phosphorus/1000*100/(100-Water)) #conversion to units of percentages in DM. Nitrogen.total is in units of g/100 edible gram; Phosphorous in the dataset is in mg/100 edible gram
fishNutrition1$N_built_in[which(fishNutrition1$N_built_in>5)]=NaN   #remove large unrealstic numbers

fishNutrition2<-fishNutrition1 %>% group_by(species) %>% summarise(P_byC=mean(P_C,na.rm = TRUE),N_byC=mean(N_C,na.rm = TRUE),
                                                      P_fat_1=mean(P_fat,na.rm = TRUE),N_fat_1=mean(N_fat,na.rm = TRUE),
                                                      P_built_in1=mean(P_built_in,na.rm = TRUE),N_built_in1=mean(N_built_in,na.rm = TRUE)) %>%                                 
                                                      rowwise() %>%                                                    
                                                     mutate(N_avg=mean(c(N_byC,N_fat_1,N_built_in1),na.rm = TRUE))%>%
                                                      mutate(P_avg=mean(c(P_byC,P_fat_1,P_built_in1),na.rm = TRUE))%>%
                                                      select(c("species","N_avg","P_avg"))



#read LCA db  
db_full<-data.table(read.csv("LCA_compiled_20201109.csv", stringsAsFactors = FALSE))
#filter those entries without scentific name
db_full<-db_full %>% filter(Species.scientific.name != "")
#add taxa data to our db
db_full<-merge(x=db_full, y=fish_Nutrition_genus, by.x="Species.scientific.name",by.y="species",all.x=TRUE) 
#merge with with full LCA
db_model_full<-merge(x=db_full, y=fishNutrition2, by.x="Species.scientific.name",by.y="species",all.x=TRUE)

# for N
db_model_full_try<-db_model_full

missing_N<-which(is.na(db_model_full$N_avg))
db_model_full$Species.common.name[missing_N]
db_model_full$Species.scientific.name[missing_N]
for (lm in 1:length(missing_N)){
 
   #match in descending order of importance until value is obtained
  #first match - by common name and alternative scientific names
mm<-agrep(db_model_full_try$Species.common.name[missing_N[lm]], fishNutrition1$alt.scinames, max = 5, ignore.case = TRUE)
N_built_2<-mean(fishNutrition1$N_built_in[mm],na.rm = TRUE);N_fat_2<-mean(fishNutrition1$N_fat[mm],na.rm = TRUE);N_C_2<-mean(fishNutrition1$N_C[mm],na.rm = TRUE);
dd<-mean(c(N_built_2,N_fat_2,N_C_2),na.rm = TRUE)

#second match - by common name and used name
if (is.na(dd)){mm<-agrep(db_model_full_try$Species.common.name[missing_N[lm]], fishNutrition1$Food.name.in.English, max=5, ignore.case = TRUE)
N_built_2<-mean(fishNutrition1$N_built_in[mm],na.rm = TRUE);N_fat_2<-mean(fishNutrition1$N_fat[mm],na.rm = TRUE);N_C_2<-mean(fishNutrition1$N_C[mm],na.rm = TRUE);
dd<-mean(c(N_built_2,N_fat_2,N_C_2),na.rm = TRUE)}

#third match - by genus
if (is.na(dd) & !is.na(db_model_full_try$genus[missing_N[lm]])){mm<-agrep(db_model_full_try$genus[missing_N[lm]], fishNutrition1$genus, ignore.case = TRUE,value = FALSE)
N_built_2<-mean(fishNutrition1$N_built_in[mm],na.rm = TRUE);N_fat_2<-mean(fishNutrition1$N_fat[mm],na.rm = TRUE);N_C_2<-mean(fishNutrition1$N_C[mm],na.rm = TRUE);
dd<-mean(c(N_built_2,N_fat_2,N_C_2),na.rm = TRUE)}

#fourth match - by family
if (is.na(dd) & !is.na(db_model_full_try$family[missing_N[lm]])){mm<-agrep(db_model_full_try$family[missing_N[lm]], fishNutrition1$family, ignore.case = TRUE,value = FALSE)
N_built_2<-mean(fishNutrition1$N_built_in[mm],na.rm = TRUE);N_fat_2<-mean(fishNutrition1$N_fat[mm],na.rm = TRUE);N_C_2<-mean(fishNutrition1$N_C[mm],na.rm = TRUE);
dd<-mean(c(N_built_2,N_fat_2,N_C_2),na.rm = TRUE)}

#fifth match - by order
if (is.na(dd) & !is.na(db_model_full_try$order[missing_N[lm]])){mm<-agrep(db_model_full_try$order[missing_N[lm]], fishNutrition1$order, ignore.case = TRUE,value=FALSE)
N_built_2<-mean(fishNutrition1$N_built_in[mm],na.rm = TRUE);N_fat_2<-mean(fishNutrition1$N_fat[mm],na.rm = TRUE);N_C_2<-mean(fishNutrition1$N_C[mm],na.rm = TRUE);
dd<-mean(c(N_built_2,N_fat_2,N_C_2),na.rm = TRUE)}

#sixth match - by loosely resembling words
if (is.na(dd) & !is.na(db_model_full_try$Species.common.name[missing_N[lm]])){mm<-agrep(db_model_full_try$Species.common.name[missing_N[lm]], fishNutrition1$alt.scinames, max.distance =0.5, ignore.case = TRUE,value=FALSE)
N_built_2<-mean(fishNutrition1$N_built_in[mm],na.rm = TRUE);N_fat_2<-mean(fishNutrition1$N_fat[mm],na.rm = TRUE);N_C_2<-mean(fishNutrition1$N_C[mm],na.rm = TRUE);
dd<-mean(c(N_built_2,N_fat_2,N_C_2),na.rm = TRUE)}

#seventh match - by loosely resembling words
if (is.na(dd) & !is.na(db_model_full_try$Species.common.name[missing_N[lm]])){mm<-agrep(db_model_full_try$Species.common.name[missing_N[lm]], fishNutrition1$Food.name.in.English, max.distance =0.5, ignore.case = TRUE,value=FALSE)
N_built_2<-mean(fishNutrition1$N_built_in[mm],na.rm = TRUE);N_fat_2<-mean(fishNutrition1$N_fat[mm],na.rm = TRUE);N_C_2<-mean(fishNutrition1$N_C[mm],na.rm = TRUE);
dd<-mean(c(N_built_2,N_fat_2,N_C_2),na.rm = TRUE)}

#if (is.na(dd)){print(lm)}
db_model_full_try$N_avg[missing_N[lm]]<-dd

rm(N_built_2,N_fat_2,N_C_2,mm,dd)
}


# for P
missing_P<-which(is.na(db_model_full_try$P_avg))
db_model_full_try$Species.common.name[missing_P]
db_model_full_try$Species.scientific.name[missing_P]
for (lm in 1:length(missing_P)){
  
  #match in descending order of importance until value is obtained
  #first match - by common name and alternative scientific names
  mm<-agrep(db_model_full_try$Species.common.name[missing_P[lm]], fishNutrition1$alt.scinames, max = 5, ignore.case = TRUE)
  P_built_2<-mean(fishNutrition1$P_built_in[mm],na.rm = TRUE);P_fat_2<-mean(fishNutrition1$P_fat[mm],na.rm = TRUE);P_C_2<-mean(fishNutrition1$P_C[mm],na.rm = TRUE);
  dd<-mean(c(P_built_2,P_fat_2,P_C_2),na.rm = TRUE)
  
  #second match - by common name and used name
  if (is.na(dd)){mm<-agrep(db_model_full_try$Species.common.name[missing_P[lm]], fishNutrition1$Food.name.in.English, max=5, ignore.case = TRUE)
  P_built_2<-mean(fishNutrition1$P_built_in[mm],na.rm = TRUE);P_fat_2<-mean(fishNutrition1$P_fat[mm],na.rm = TRUE);P_C_2<-mean(fishNutrition1$P_C[mm],na.rm = TRUE);
  dd<-mean(c(P_built_2,P_fat_2,P_C_2),na.rm = TRUE)}
  
  #third match - by genus
  if (is.na(dd) & !is.na(db_model_full_try$genus[missing_N[lm]])){mm<-agrep(db_model_full_try$genus[missing_P[lm]], fishNutrition1$genus, ignore.case = TRUE,value = FALSE)
  P_built_2<-mean(fishNutrition1$P_built_in[mm],na.rm = TRUE);P_fat_2<-mean(fishNutrition1$P_fat[mm],na.rm = TRUE);P_C_2<-mean(fishNutrition1$P_C[mm],na.rm = TRUE);
  dd<-mean(c(P_built_2,P_fat_2,P_C_2),na.rm = TRUE)}
  
  #fourth match - by family
  if (is.na(dd) & !is.na(db_model_full_try$family[missing_P[lm]])){mm<-agrep(db_model_full_try$family[missing_P[lm]], fishNutrition1$family, ignore.case = TRUE,value = FALSE)
  P_built_2<-mean(fishNutrition1$P_built_in[mm],na.rm = TRUE);P_fat_2<-mean(fishNutrition1$P_fat[mm],na.rm = TRUE);P_C_2<-mean(fishNutrition1$P_C[mm],na.rm = TRUE);
  dd<-mean(c(P_built_2,P_fat_2,P_C_2),na.rm = TRUE)}
  
  #fifth match - by order
  if (is.na(dd) & !is.na(db_model_full_try$order[missing_P[lm]])){mm<-agrep(db_model_full_try$order[missing_P[lm]], fishNutrition1$order, ignore.case = TRUE,value=FALSE)
  P_built_2<-mean(fishNutrition1$P_built_in[mm],na.rm = TRUE);P_fat_2<-mean(fishNutrition1$P_fat[mm],na.rm = TRUE);P_C_2<-mean(fishNutrition1$P_C[mm],na.rm = TRUE);
  dd<-mean(c(P_built_2,P_fat_2,P_C_2),na.rm = TRUE)}
  
  #sixth match - by loosely resembling words
  if (is.na(dd) & !is.na(db_model_full_try$Species.common.name[missing_P[lm]])){mm<-agrep(db_model_full_try$Species.common.name[missing_P[lm]], fishNutrition1$alt.scinames, max.distance =0.5, ignore.case = TRUE,value=FALSE)
  P_built_2<-mean(fishNutrition1$P_built_in[mm],na.rm = TRUE);P_fat_2<-mean(fishNutrition1$P_fat[mm],na.rm = TRUE);P_C_2<-mean(fishNutrition1$P_C[mm],na.rm = TRUE);
  dd<-mean(c(P_built_2,P_fat_2,P_C_2),na.rm = TRUE)}
  
  #seventh match - by loosely resembling words
  if (is.na(dd) & !is.na(db_model_full_try$Species.common.name[missing_P[lm]])){mm<-agrep(db_model_full_try$Species.common.name[missing_P[lm]], fishNutrition1$Food.name.in.English, max.distance =0.5, ignore.case = TRUE,value=FALSE)
  P_built_2<-mean(fishNutrition1$P_built_in[mm],na.rm = TRUE);P_fat_2<-mean(fishNutrition1$P_fat[mm],na.rm = TRUE);P_C_2<-mean(fishNutrition1$P_C[mm],na.rm = TRUE);
  dd<-mean(c(P_built_2,P_fat_2,P_C_2),na.rm = TRUE)}
  
  #if (is.na(dd)){print(lm)}
  db_model_full_try$P_avg[missing_P[lm]]<-dd
  
  rm(P_built_2,P_fat_2,P_C_2,mm,dd)
}

write.xlsx(db_model_full_try,'LCA_full_data_compiled.xlsx')


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

