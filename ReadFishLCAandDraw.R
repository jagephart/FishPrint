#a program to read and draw basic data from LCA inventory dataset
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



#-----------------------our dataset
Fish <-
  read.xlsx(
    "Calculator_20200428.xlsx",
    sheet = "Model_input"
  )
FishGroups <-
  read.xlsx(
    "Calculator_20200428.xlsx",
    sheet = "Diet and Taxonomy"
  )
Fishmerge <-
  merge(
    x = Fish,
    y = FishGroups,
    by.x = c("Species.common.name"),
    by.y = c("ASFIS.species.(ASFIS.species)"),
    all.x = TRUE,
    all.y = FALSE
  )

# remove lines with no scientific names which is the basis for merging
Fishmerge<-Fishmerge %>%
  drop_na(Species.scientific.name)

# -----------------------Load aquaculture data fishstatj 

fishstatj <- data.table(read.csv("globalAquaculutreFishStatJ.csv",stringsAsFactors = FALSE))
fishstatj$Species..ASFIS.species..2<-str_replace_all(fishstatj$Species..ASFIS.species..2,"X","0")
fishstatj$Species..ASFIS.species..2=sapply(fishstatj$Species..ASFIS.species..2,as.numeric)
fishstatj<-within(fishstatj,{Species = Species..ASFIS.species..1;rm(Species..ASFIS.species..1)})
fishstatj <- within(fishstatj, {ASFIS =Species..ASFIS.species..2;rm(Species..ASFIS.species..2)})
#reorder by column index and remove
fishstatj <- fishstatj[,c(7,1,2,3,4,6,8)]

## -----------------------add ecological traits through RFishBase and merge into a new dataset


ec<-data.table(ecology(species_list =, fields =, server =))
oxy<-data.table(oxygen(species_list =, fields =, server =))
mat<-data.table(maturity(species_list =, fields =, server =))

l <- merge(x=fishstatj, y=ec[,c("Species","EcologyRefNo","HabitatsRef","DietTroph","DietSeTroph","DietTLu","DietseTLu","FoodTroph","FoodSeTroph")], by.x=c("Species"),by.y =c("Species"), all.x=TRUE, all.y=FALSE) #merge with ecology
l1 <- merge(x=l, y=oxy[,c("Species","Temp","Weight","Salinity","Oxygen","OxygenCons","Oxygenmgl")], by.x=c("Species"),by.y =c("Species"), all.x=TRUE, all.y=FALSE) #merge with oxygen database
l2 <- merge(x=l1, y=mat[,c("Species","AgeMatMin","AgeMatMin2",
                           "LengthMatMin","LengthMatMin2","tm","Lm","SE_Lm","SD_Lm","LCL_Lm","UCL_Lm")], by.x=c("Species"),by.y =c("Species"), all.x=TRUE, all.y=FALSE,allow.cartesian=TRUE)    #merge with maturity life history

#more on merges: http://www.datasciencemadesimple.com/join-in-r-merge-in-r/

#remove unecessary columns
l2<-within(l2,{rm("DateModified",'DateEntered','Expert.y',"DateChecked","TS.y","Modified.y","Comment","Locality","C_Code","E_CODE","Entered.y","DietRemark","Entered.x","Modified.x","Datemodified","Dateentered","Expert.x","Datechecked","TS.x","Number","FoodRef","IHRemarks","SubstrateRef","OIRemarks","LengthMatRef","Sex","Aquaculture.area..FAO.major.fishing.area.")})
rm(l,l1)  #clear uncessary variables and files



#---------------merge only with ecology
Fishmerge <-
  merge(
    x = Fishmerge,
    y = ec[, c(
      "Species",
      "EcologyRefNo",
      "HabitatsRef",
      "DietTroph",
      "DietSeTroph",
      "DietTLu",
      "DietseTLu",
      "FoodTroph",
      "FoodSeTroph"
    )],
    by.x = c("Species.scientific.name"),
    by.y = c("Species"),
    all.x = TRUE,
    all.y = FALSE)

#---------------merge with FishStatJ
Fishmerge <-
  merge(
    x = Fishmerge,
    y = fishstatj[, c(
      "Species",
      "Species..Family.",
      "ASFIS",
      "Species..Main.grouping.",
      "Species..Order.",
      "Environment..Environment."
    )],
    by.x = c("Species.scientific.name"),
    by.y = c("Species"),
    all.x = TRUE,
    all.y = FALSE)


#--------------------------------

Fishmerge$Taxonomy <-
  ordered(Fishmerge$Taxonomy ,
          levels = c("Teleost", "Mollusca", "Crustacean", "Reptilia", ""))

Fishmerge$HighLevelISSCAAP = cut(Fishmerge$ISSCAAP_Group, c(10, 20, 30, 40, 50, 60))
levels(Fishmerge$HighLevelISSCAAP) = c("freshwater", "diadromous", "marine", "Crustacean", "Mollusca")



#----------store dataset
fwrite(Fishmerge, "LCA_data.csv")


#-----------linear regression between FCR and other variables

#remove duplicate observations (rows) that were added through the merge (join) and focus on the attributes

df<-Fishmerge %>% distinct(Species.scientific.name,FCR,FoodTroph,Feed_FMFO_percent,USD.per.tonne,ASFIS,.keep_all=TRUE)

linearModFCR <- lm(FCR~FoodTroph+Feed_FMFO_percent+USD.per.tonne+ASFIS, data=df)
linearModFCR <- lm(FCR~USD.per.tonne, data=df)
summary(linearModFCR) 


#-----------------
# draw 1
#-----------------

fig1 <- ggplot(df, aes(USD.per.tonne,FCR,color = HighLevelISSCAAP,label = Species.scientific.name)) +
geom_text(
  check_overlap = TRUE,
  color = "black",
  nudge_y = 0.1,
  nudge_x = 700,
  size = 2.5,
  alpha = 0.5
) +
geom_point(size = 2.0, aes(shape = Diet)) +
  ylim(0, 3)
fig1
ggsave("FCRbyPrice.jpg")

#-----------------
# draw 2
#-----------------

fig2 <-
  ggplot(
    Fishmerge,
    aes(FCR, ISSCAAP_Group, color = HighLevelISSCAAP, label = Species.scientific.name)
  ) +
  geom_boxplot(outlier.shape = 3) +
  geom_text(
    check_overlap = TRUE,
    color = "black",
    nudge_y = 2,
    size = 2.5,
    alpha = 0.5
  ) +
  geom_point(size = 2.0, aes(shape = Diet)) +
  xlim(0, 3) +
  ylab("ISSCAAP")
fig2
ggsave("FCRbyISSCAAP.jpg")

#-------------
#draw 3
#------------
#linear regression between Trophic level and FMFO
plot(density(Fishmerge$Feed_FMFO_percent,na.rm = TRUE))
linearMod <- lm(FoodTroph~Feed_FMFO_percent, data=df)
summary(linearMod) 
# Extract coefficients
beta = coef(linearMod)

fig3<-ggplot(Fishmerge,aes(Feed_FMFO_percent,FoodTroph,color=HighLevelISSCAAP,label=
                             Species.scientific.name))+
geom_point(aes(shape = Diet))+
  geom_text(
    check_overlap = TRUE,
    color = "black",
    nudge_y = 0.1,
    size = 2.5,
    alpha = 0.5
  )+
  geom_abline(intercept=beta[1],slope=beta[2])
fig3
ggsave("TrophicLevelbyFMFO.jpg")



  
  #intake$Nutrient <- ordered(intake$Nutrient, levels = c("protein", "iron", "zinc", 
  #                                    "vitamin A", "calcium"))
  #nutrition$key <- factor(nutrition$key, levels = c("protein", "iron", "zinc","vitA"), 
  #                       labels = c("Protein (g/100g)", "Iron (mg/100g)", "Zinc (mg/100g)",
  #                                   "Vit. A (RAE \u03bcg/100g)"))
  
  #nutrition<-data.frame(taxa=llg$Species..Main.grouping.,protein=llg$protein,zinc=llg$zinc,iron=llg$iron,vitA=llg$vita)
  #colnames(nutrition) <- c("specie", "protein","zinc","iron","vitA")
  #Fishmerge<-pivot_longer(Fishmerge,-Taxonomy, names_to = "key", values_to = "values")
  #nutrition<-nutrition %>% gather(key, value,protein:vitA
  