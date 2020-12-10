

#a program to add N, P and edible portions to excel using three different methods.


library("rfishbase")  #package to use fishbase in r
library("openxlsx")
library("rlang")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("data.table")
library("stringr")
library("openxlsx")
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

#load data from excel files
db <-data.table(read.csv("missing_NP_values.csv", stringsAsFactors = FALSE))
db_full<-data.table(read.csv("LCA_compiled_20201109.csv", stringsAsFactors = FALSE))


db$Scientific.Name<-as.factor(db$Scientific.Name)
peter_db <-read.xlsx("Peters edible muscle & protein portion of seafood.xlsx")
#leave only max edible fraction
peter_db<-peter_db %>%select(TaxonName,`Fraction.from.live.weight.%`)%>%mutate(Edible_portion_b=`Fraction.from.live.weight.%`/100) %>% select(TaxonName,Edible_portion_b)

#make sure fishNutition1 of fish content was computed using N_P_DischargeIncFunction
db_n<-fishNutrition1[,c("species","N","P","N_fat","P_fat","N_built_in","P_built_in","Edible.portion.coefficient","Parts")]

#aggregate based in species
r<-db_n %>% group_by(species) %>% summarise(P_byC=mean(P,na.rm = TRUE),N_byC=mean(N,na.rm = TRUE),
                                                    P_fat_1=mean(P_fat,na.rm = TRUE),N_fat_1=mean(N_fat,na.rm = TRUE),
                                                    P_built_in1=mean(P_built_in,na.rm = TRUE),N_built_in1=mean(N_built_in,na.rm = TRUE)) %>%
                                                    
  r_edible<-db_n %>% group_by(species)%>%filter(Parts=="f")%>%
                                                  summarize(Edible_portion=mean(Edible.portion.coefficient,na.rm = TRUE)) %>%
                                                 select(species,Edible_portion)

#merge edible and N/P
r_op<-merge(x=r, y=r_edible, by.x = "species",by.y = "species" ,all.x=TRUE)


#merge with Peter's excel
db1<-merge(x=r_op, y=peter_db, by.x = "species",by.y = "TaxonName" ,all.x=TRUE)

#merge with our missing data excel
db_model<-merge(x=db, y=db1, by.x = "Scientific.Name",by.y = "species" ,all.x=TRUE)
write.xlsx(db_model,'missing NP.xlsx')

#merge with with full LCA
db_model_full<-merge(x=db_full, y=db1, by.x = "Species.scientific.name",by.y = "species" ,all.x=TRUE)
write.xlsx(db_model_full,'LCA_full_data_compiled.xlsx')

#subset of the above for presentation
db_model_short<-db_model_full%>%select(Species.scientific.name,Species.common.name,Edible_portion,Edible_portion_b) %>%
  group_by(Species.scientific.name) %>% summarize(Edible_portion=mean(Edible_portion,na.rm = TRUE),
                                                  Edible_portion_b=mean(Edible_portion_b,na.rm = TRUE),
                                                  Species.common.name=unique(Species.common.name))
write.xlsx(db_model_short,'LCA_full_data__short_compiled.xlsx')
