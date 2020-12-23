#a program to read N and P values of aquatic plants
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

#calculate C,N,P content of fish (based on energy and fat content )
#N and P content of fish in % of DM
#read Zach's data, filter only aquatic plants and make calculations using built functions
aquatic_plants=c("Gracilaria chilensis", "Laminaria digitata", "Macrocystis pyrifera", "Saccharina latissima")
fishNutrition <- data.table(read.csv("AFCD_live.csv", stringsAsFactors = FALSE)) %>%filter(species %in% aquatic_plants)


fishNutrition$Nitrogen.total[which(fishNutrition$Nitrogen.total>100)]=NaN   #remove large unrealistic numbers
fishNutrition<- fishNutrition %>%select(c("species","Nitrogen.total","Phosphorus"))

#-----------------------compute N and P for aquatic plants in % of wet weight
fishNutrition1<-fishNutrition %>% mutate(N_wet=Nitrogen.total,P_wet=Phosphorus/1000) %>% #conversion to units of percentages in WET WEIGHT. Nitrogen.total is in units of g/100 edible gram; Phosphorous in the dataset is in mg/100 edible gram
group_by(species) %>% summarise(P_avg=mean(P_wet,na.rm = TRUE),N_avg=mean(N_wet,na.rm = TRUE))                          








