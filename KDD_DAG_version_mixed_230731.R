# Revision for KDD Paper
# Final Version Code
# 2023/07/09
# Written by Inyoung Jun

#0. Set up
library(tidyverse)
`%!in%` = Negate(`%in%`)
daydiff <- function(a,b,c,d){
  result<-(a-b)*365+c-d
  return(result)
} #a:Later Year, b:Former Year, c:Later Date, d:Former Date
options(scipen = 999)
setwd("/blue/prosperi/inyoungjun/KDD2023")

##### Import Dataset
d <-readRDS("/home/inyoungjun/DATA/data/dnew.RDS") %>%  #diagnosis
  rename(PID = Deidentified_Patient_ID, D_Year = Year, D_Date = Date_of_Service) %>%
  filter(D_Year>=2011 & D_Year<=2019) %>% select(-ICD_Code3)
a <-readRDS("/home/inyoungjun/DATA/data/a.rds") %>% #antibiogram
  rename(PID = Deidentified_Patient_ID, A_Year = Year, A_Date = Date_of_Service)
demo <-readRDS("/home/inyoungjun/DATA/data/demo.RDS") #demographics
ad <-readRDS("/home/inyoungjun/DATA/data/ad.RDS") %>%  #admission
  rename(PID = Deidentified_Patient_ID, 
         AD_Year = Year_of_Admit, AD_Date = admit_date, 
         DC_Year = Year_of_Discharge, DC_Date = discharge_date,
         ICU = icu_y_n)
m <-readRDS("/home/inyoungjun/DATA/data/m.RDS") #medication
l <-readRDS("/home/inyoungjun/DATA/data/l.RDS") #laboratory
p <-readRDS("/home/inyoungjun/DATA/data/p.RDS") #procedure

##### Study Population - Invasive MRSA
#Step1: Culture Test on body site
a_Invasive_MRSA <- a %>%
  filter(str_detect(tolower(organism),"staph") & str_detect(tolower(organism),"aureus")) %>%
  filter(str_detect(tolower(suscept),"resist")) %>% 
  filter(str_detect(tolower(specimen_type),"blood|fluid|bone|brain|heart|liver|kidney|pancreas|ovary|spleen|lung")|str_detect(tolower(specimen_source),"venipuncture|blood|fluid|bone|brain|heart|liver|kidney|pancreas|ovary|spleen|lung")) %>% 
  filter(str_detect(tolower(antibiotic),"oxacillin|methicillin")) %>% 
  arrange(PID,A_Year,A_Date) %>% 
  distinct(PID,.keep_all = TRUE) %>% 
  mutate(BC_Year = A_Year, BC_Date = A_Date) %>% 
  select(PID,BC_Year,BC_Date,specimen_type,specimen_source)#1433 unique patients
-sort(-table(droplevels(as.factor(a_Invasive_MRSA$BC_Year))))#2011 - 2019
-sort(-table(droplevels(as.factor(a_Invasive_MRSA$specimen_type))))#Blood/Fluid
-sort(-table(droplevels(as.factor(a_Invasive_MRSA$specimen_type))))#2011 - 2019
n_distinct(a_Invasive_MRSA)#1433
#Step2: One year prior medical history
d_Invasive_MRSA_oneyear <- d %>%
  filter(D_Year>=2011) %>% 
  filter(PID %in% a_Invasive_MRSA$PID) %>% 
  left_join(select(a_Invasive_MRSA,c("PID","BC_Year","BC_Date","specimen_type","specimen_source"))) %>% 
  group_by(PID) %>% 
  mutate(Diff_oneyear = daydiff(BC_Year,D_Year,BC_Date,D_Date)) %>% 
  filter(Diff_oneyear>=365)
d_Invasive_oneyear_ID<-unique(d_Invasive_MRSA_oneyear$PID)#914

#Step3: Attach Admission and Discharge
ad_Invasive_MRSA <- ad %>% 
  filter(PID %in% d_Invasive_oneyear_ID) %>% 
  left_join(select(d_Invasive_MRSA_oneyear,c("PID","BC_Year","BC_Date","specimen_type","specimen_source")), by="PID") %>%
  filter(daydiff(BC_Year,AD_Year,BC_Date,AD_Date)>=0 & daydiff(DC_Year,BC_Year,DC_Date,BC_Date)>=0) %>%
  group_by(PID) %>% 
  mutate(R_AD_Year = ifelse(daydiff(lag(DC_Year),AD_Year,lag(DC_Date),AD_Date)>=0,lag(AD_Year),AD_Year),
         R_AD_Date = ifelse(daydiff(lag(DC_Year),AD_Year,lag(DC_Date),AD_Date)>=0,lag(AD_Date),AD_Date)) %>% 
  mutate(AD_Year2 = ifelse(is.na(R_AD_Year)==TRUE,AD_Year,R_AD_Year),
         AD_Date2 = ifelse(is.na(R_AD_Date)==TRUE,AD_Date,R_AD_Date)) %>% 
  distinct(PID,AD_Year,AD_Date,.keep_all = TRUE) %>% 
  select(-AD_Year,-AD_Date,-R_AD_Year,-R_AD_Date) %>% 
  rename(AD_Year = AD_Year2, AD_Date = AD_Date2) %>% 
  arrange(PID,-DC_Year,-DC_Date) %>%
  distinct(PID,.keep_all = TRUE) %>% 
  mutate(HCI = ifelse(daydiff(BC_Year,AD_Year,BC_Date,AD_Date)>3,1,0)) %>% 
  select(PID,AD_Year,AD_Date,BC_Year,BC_Date,DC_Year,DC_Date,ICU,HCI,specimen_type,specimen_source) %>% 
  mutate(length_addc = daydiff(DC_Year,AD_Year,DC_Date,AD_Date))
head(ad_Invasive_MRSA)

n_distinct(ad_Invasive_MRSA$PID)#875 -> 878 (three patients showed duplicate AD-DCrecord) -> 872 (after adjusted)
table(ad_Invasive_MRSA$HCI,useNA = "ifany") #181 HCI, 691 CA; 160 HCI, 712 CA
round(160/872*100,2)#18.35%


##### Attach Demographic Records
demo_MRSA <- demo %>% filter(PID %in% ad_Invasive_MRSA$PID)

##### Attach Medication Records
m_MRSA <- m %>% filter(PID %in% demo_MRSA$PID) %>% 
  rename(MT_Date = Med_Taken_Date_of_Service) %>% 
  left_join(select(ad_Invasive_MRSA,c("PID","AD_Year","AD_Date","BC_Year","BC_Date","DC_Year","DC_Date",'specimen_source','specimen_type'))) %>% 
  select(PID,AD_Year,AD_Date,BC_Year,BC_Date,M_Year,M_Date,MT_Date,DC_Year,DC_Date,Simple_Generic_Med_Name,Med_Strength,Total_Dose,MAR_Action,Dispensed) %>%
  filter(daydiff(M_Year,AD_Year,MT_Date,AD_Date)>=0 & daydiff(DC_Year,M_Year,DC_Date,MT_Date)>=0) %>% 
  arrange(PID,M_Year,MT_Date) %>% 
  filter(MAR_Action %!in% c("MISSED","STOPPED","HELD")) %>% 
  distinct(PID,M_Year,M_Date,MT_Date,Simple_Generic_Med_Name,.keep_all = TRUE) %>% 
  mutate(M_Time2 = ifelse(daydiff(M_Year,BC_Year,MT_Date,BC_Date)==10, "Time3",NA),
         M_Time1 = ifelse(daydiff(M_Year,BC_Year,MT_Date,BC_Date)==3, "Time2",NA),
         M_Time0 = ifelse(daydiff(M_Year,BC_Year,MT_Date,BC_Date)<3,"Time1",NA)) %>%
  mutate(M_Group = coalesce(M_Time0,M_Time1,M_Time2)) %>% 
  mutate(Vanco = ifelse(str_detect(tolower(Simple_Generic_Med_Name),"vancomycin"),1,0)) %>% 
  filter(M_Group %in% c("Time1","Time2","Time3")) %>% 
  group_by(PID,M_Group) %>%
  filter(Vanco==max(Vanco)) %>% 
  distinct(Vanco,.keep_all = TRUE)
table(m_MRSA$M_Group)#Time0, Time1, Time2

n_distinct(m_MRSA$PID[m_MRSA$M_Group=="Time1"])#312 -> 275(MAR)-> 227 =>812 => 812
n_distinct(m_MRSA$PID[m_MRSA$M_Group=="Time2"])#973 -> 914(MAR) ->672 =>707 => 707
n_distinct(m_MRSA$PID[m_MRSA$M_Group=="Time3"])#592 -> 550 (MAR) ->734 =>663 => 427 (when we changed the criteria to day 10)
table(m_MRSA$M_Group,m_MRSA$Vanco)

##### make dataset
data_Invasive_MRSA <- m_MRSA %>%  
  left_join(select(ad_Invasive_MRSA,c("PID","ICU","HCI","specimen_type","specimen_source")),by="PID") %>% 
  left_join(demo_MRSA,by="PID") %>% 
  mutate(mortality30 = ifelse(daydiff(Death_Year,BC_Year,Death_Date,BC_Date)<=30,1,0)) %>% 
  arrange(PID) %>% 
  mutate(mortality30 = replace_na(mortality30, 0)) %>% 
  select(-M_Time0,-M_Time1,-M_Time2)
colnames(data_Invasive_MRSA)
n_distinct(data_Invasive_MRSA$PID)#872 -> 817

time_counts <- data_Invasive_MRSA %>% group_by(PID) %>% 
  summarise(num_time_values = n_distinct(M_Group))
patients_with_all_times <- time_counts %>% filter(num_time_values == 3) %>% 
  select(PID)#663 => 427 (patients who had all time1,2,3)
patients_with_two_times <- time_counts %>% filter(num_time_values == 2) %>% 
  select(PID)#663 => 427 (patients who had all time1,2,3)
patients_with_one_times <- time_counts %>% filter(num_time_values == 1) %>% 
  select(PID)#663 => 427 (patients who had all time1,2,3)

table(data_Invasive_MRSA[data_Invasive_MRSA$M_Group=="Time1",]$Vanco)

#Join the filtered patients back to the original data
data_Invasive_MRSA <- data_Invasive_MRSA %>% 
   filter(PID %in% patients_with_all_times$PID)
n_distinct(data_Invasive_MRSA$PID)#663->427

#make a time-varying variable1 (Ex. C-reactive) => Too little information
# CR <- l %>% 
#       filter(PID %in% data_Invasive_MRSA$PID) %>% 
#       filter(str_detect(tolower(Lab_Name),"c reactive|c-reactive")) %>% 
#       full_join(select(data_Invasive_MRSA_first, c("PID","AD_Year","AD_Date","BC_Year","BC_Date","DC_Year","DC_Date",)), by="PID") %>% 
#       filter(daydiff(L_Year,AD_Year,L_Date,AD_Date)>=0 & daydiff(DC_Year,L_Year,DC_Date,L_Date)>=0) %>% 
#       mutate(L_Time3 = ifelse(daydiff(L_Year,BC_Year,L_Date,BC_Date)>=10,"Time3",NA), #interval
#              L_Time2 = ifelse(daydiff(L_Year,BC_Year,L_Date,BC_Date)>=3 & daydiff(L_Year,BC_Year,L_Date,BC_Date)<10,"Time2",NA), #interval
#              L_Time1 = ifelse(daydiff(L_Year,BC_Year,L_Date,BC_Date)<3,"Time1",NA)) %>% #interval
#       mutate(L_Group = coalesce(L_Time1, L_Time2, L_Time3)) %>% 
#       filter(is.na(L_Group)==FALSE) %>% 
#       filter(Result_Unit == "mg/L") %>% 
#       mutate(C_Reactive = ifelse(is.na(as.numeric(Lab_Result)), NA, as.numeric(Lab_Result))) %>% 
#       arrange(PID,L_Year,L_Date) %>%
#       group_by(PID,L_Group) %>% 
#       mutate(C_Reactive_AVG = mean(C_Reactive)) %>% 
#       distinct(PID, L_Group,.keep_all = TRUE) %>% 
#       select(PID, L_Group, C_Reactive_AVG) %>% 
#       mutate(M_Group = L_Group)
#Select the first case of each patient
data_Invasive_MRSA_first <- data_Invasive_MRSA %>% filter(M_Group=="Time1")

#make a time-varying variable2 (EX. Nephrotoxicity using Creatine)
CREA <- l %>% 
  filter(PID %in% data_Invasive_MRSA$PID) %>% 
  filter(str_detect(tolower(Lab_Name),"creatinine")) %>% 
  full_join(select(data_Invasive_MRSA_first, c("PID","AD_Year","AD_Date","BC_Year","BC_Date","DC_Year","DC_Date",)), by="PID") %>% 
  filter(daydiff(L_Year,AD_Year,L_Date,AD_Date)>=0 & daydiff(DC_Year,L_Year,DC_Date,L_Date)>=0) %>% 
  mutate(L_Time3 = ifelse(daydiff(L_Year,BC_Year,L_Date,BC_Date)>=10,"Time3",NA), #interval
         L_Time2 = ifelse(daydiff(L_Year,BC_Year,L_Date,BC_Date)>=3 & daydiff(L_Year,BC_Year,L_Date,BC_Date)<10,"Time2",NA), #interval
         L_Time1 = ifelse(daydiff(L_Year,BC_Year,L_Date,BC_Date)<3,"Time1",NA)) %>% #interval
  mutate(L_Group = coalesce(L_Time1, L_Time2, L_Time3)) %>% 
  filter(is.na(L_Group)==FALSE) %>% 
  filter(Result_Unit == "mg/dL") %>% 
  mutate(Creatinine = ifelse(is.na(as.numeric(Lab_Result)), NA, as.numeric(Lab_Result))) %>% 
  arrange(PID,L_Year,L_Date) %>%
  group_by(PID,L_Group) %>% 
  mutate(Creatinine_AVG = mean(Creatinine)) %>% 
  distinct(PID, L_Group,.keep_all = TRUE) %>% 
  select(PID, L_Group, Creatinine_AVG) %>% 
  mutate(M_Group = L_Group)
library(data.table)
CREA_Wide <- CREA %>% mutate(Var = Creatinine_AVG) %>% dcast(PID~L_Group) %>% 
  mutate(Nephrotoxicity1 = ifelse((Time1-Time2>=Time1*0.5),1,0),
         Nephrotoxicity2 = ifelse((Time1-Time3>=Time1*0.5),1,0)) %>% 
  mutate(Nephrotoxicity_any = ifelse(Nephrotoxicity1==1|Nephrotoxicity2==1,1,0))
head(CREA_Wide)
table(CREA_Wide$Nephrotoxicity1,CREA_Wide$Nephrotoxicity2)
table(CREA_Wide$Nephrotoxicity1)#28
table(CREA_Wide$Nephrotoxicity2)#39
table(CREA_Wide$Nephrotoxicity_any)#44 or 352

#CCI Score
d_comorbidity<- d %>% 
  select(PID, D_Year, D_Date, ICD_Code, ICD_Code2, ICD_Description)%>% 
  filter(PID %in% data_Invasive_MRSA$PID) %>% 
  left_join(select(data_Invasive_MRSA,c(PID,BC_Year,BC_Date)),by="PID") %>% 
  filter(daydiff(BC_Year,D_Year,BC_Date,D_Date)>=0) 
n_distinct(d_comorbidity$PID)#427

CI_ID<-d_comorbidity[,c("PID", "ICD_Code2")]
CI<-CI_ID %>% 
  rename(
    id  = PID,
    code  = ICD_Code2
  )
str(CI)
#install.packages("comorbidity")
library(comorbidity)
#install.packages("checkmate")
library(checkmate)
CI_pop<-comorbidity(
  x=CI,
  id="id",
  code="code",
  map="charlson_icd9_quan",
  assign0=TRUE,
  labelled = TRUE,
  tidy.codes = TRUE
) %>% rename(PID = id) 

CI_score_pop <-as.vector(score(CI_pop, weights="charlson", assign0=FALSE))
#Based on the CCI score, the severity of comorbidity was categorized into three grades:
#mild, with CCI scores of 1-2; moderate, with CCI scores of 3-4; and severe, with CCI scores 5+.
summary(CI_score_pop)#min 0, max 23, median 3, mean 3.923
table(CI_score_pop)
CI_pop$CI_score = CI_score_pop
CI_pop #comorbiditiy dataset

#Antibiogram Resistance Result
top_res_rank_class <- read.csv("top_res_rank_class_220905.csv")#add gentamicin (aminoglycosides)
aminoglycosides<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="aminoglycosides"]
betalactams<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="betalactams"]
carbapenems<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="carbapenems"]
fluoroquinolones<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="fluoroquinolones"]
glycopeptides<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="glycopeptides"]
sulfonamides<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="sulfonamides"]
tetracyclines<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="tetracyclines"]
polypeptides<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="polypeptides"]
others<-top_res_rank_class$Antibiotics[top_res_rank_class$Antibiotics_class=="others"]

a_pop <- a %>%
  filter(PID %in% data_Invasive_MRSA$PID) %>% 
  left_join(select(data_Invasive_MRSA, c("PID","AD_Year","AD_Date","DC_Year","DC_Date","BC_Year","BC_Date"))) %>% 
  filter(daydiff(A_Year,AD_Year,A_Date,AD_Date)<0) %>% 
  filter(suscept == "Resistant") %>% 
  mutate(Class = ifelse(str_detect(tolower(antibiotic),paste(tolower(aminoglycosides),collapse="|")),"aminoglycosides",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(betalactams),collapse = "|")),"betalactams",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(carbapenems),collapse = "|")),"carbapenems",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(fluoroquinolones),collapse = "|")),"fluoroquinolones",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(glycopeptides),collapse="|")),"glycopeptides",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(sulfonamides),collapse="|")),"sulfonamides",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(tetracyclines),collapse = "|")),"tetracyclines",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(polypeptides),collapse="|")),"polypeptides",
                 ifelse(str_detect(tolower(antibiotic),paste(tolower(others),collapse="|")),"others",NA)))))))))) %>% 
  distinct(PID,antibiotic,Class) %>% 
  filter(is.na(Class)==FALSE)
table(a_pop$Class, useNA="ifany")
n_distinct(a_pop$PID)#237 out of 427 had resistance before starting this infection
a_pop_Wide <- a_pop %>% dcast(PID ~ Class) %>% 
              mutate(across(-1, ~ifelse(.>=1,1,.))) %>% 
              mutate(MDR = ifelse(rowSums(select(., -1))>=3,1,0))

##########################################################
#1. Dataset
data <- data_Invasive_MRSA %>% 
        mutate(Diff_BC = daydiff(M_Year,BC_Year,MT_Date,BC_Date),
               Diff_Death = daydiff(Death_Year,BC_Year,Death_Date,BC_Date)) %>% 
        mutate(Outcome = as.integer(Diff_Death <=30)) %>% 
        mutate(Outcome = if_else(is.na(Outcome),0,Outcome)) %>% 
        mutate(MRSA_Age = Current_Age - (2019 - BC_Year)) %>%
        left_join(CREA_Wide) %>% 
        left_join(CI_pop) %>% 
        left_join(a_pop_Wide) %>% 
        #left_join(ICD_variable) %>% 
        filter(Race == "BLACK"|Race == "WHITE") %>% #to make obvious comparison
        mutate(Sex = ifelse(Sex=="FEMALE",1,0),Race = ifelse(Race=="BLACK",1,0),ICU=ifelse(ICU=="Y",1,0))
        
#saveRDS(data, "data.RDS") #230725
###########################################
data_T1 <- data %>% filter(M_Group=="Time1")
data_T2 <- data %>% filter(M_Group=="Time2")
data_T3 <- data %>% filter(M_Group=="Time3")

data_combined <- data_T1 %>% 
                 select(PID, MRSA_Age, Sex, Race,
                        CI_score,ICU,MDR,
                        Nephrotoxicity_any,Vanco,Outcome) %>%
                 rename(Vanco_T1 = Vanco) %>% 
                 left_join(select(data_T2,c("PID","Vanco")),by="PID") %>% 
                 rename(Vanco_T2 = Vanco) %>% 
                 left_join(select(data_T3,c("PID","Vanco")),by="PID") %>% 
                 rename(Vanco_T3 = Vanco) %>% 
                 select(-M_Group.x, -M_Group.y, -M_Group) %>%  #427 obs with MDR NA 190 and Nephrotoxicity NA 31
                 mutate(across(c("MDR","Nephrotoxicity_any"), ~replace_na(., 0)))
#write.csv(data_combined,"data_combined.csv") #KDD version variables

########################################## Add newly implemented variables
# Attach variables of interest
Sources_Blood <- data_Invasive_MRSA %>%
  mutate(Source_blood = ifelse(str_detect(tolower(specimen_source),"venipuncture|blood"),1,0)) %>% 
  filter(M_Group=="Time1") %>% 
  select(PID,Source_blood)                        
table(Sources_Blood$Source_blood)#338 vs 89

Severity_of_Infection <- data_combined %>% 
  mutate(Severity = ifelse(ICU==1,1,0)) %>% 
  select(PID, Severity)
table(Severity_of_Infection$Severity)#330 vs 487

Kidney_Impairment <- d %>%
  filter(PID %in% data_Invasive_MRSA$PID) %>% 
  left_join(select(data_Invasive_MRSA, PID, BC_Year, BC_Date), by = "PID") %>% 
  group_by(PID) %>% 
  mutate(Diff_BC_D = daydiff(BC_Year, D_Year, BC_Date, D_Date)) %>% 
  filter(Diff_BC_D >= 0) %>% 
  mutate(Kidney_Impairment_before = ifelse(str_detect(ICD_Code2, "^5851|^5852|^5853|^5854|^5855|^5856|^5845|^5846|^586|^5859"),1,0)) %>%
  group_by(PID) %>% 
  arrange(PID,desc(Kidney_Impairment_before)) %>% 
  distinct(PID,.keep_all = TRUE)
table(Kidney_Impairment$Kidney_Impairment_before,useNA = "ifany") #407 vs 410

Previous_vanco <- m_MRSA <- m %>% filter(PID %in% demo_MRSA$PID) %>% 
  rename(MT_Date = Med_Taken_Date_of_Service) %>% 
  left_join(select(ad_Invasive_MRSA,c("PID","AD_Year","AD_Date","BC_Year","BC_Date","DC_Year","DC_Date"))) %>% 
  select(PID,AD_Year,AD_Date,BC_Year,BC_Date,M_Year,M_Date,MT_Date,DC_Year,DC_Date,Simple_Generic_Med_Name,Med_Strength,Total_Dose,MAR_Action,Dispensed) %>%
  filter(daydiff(AD_Year,M_Year,AD_Date,MT_Date)>0) %>% 
  arrange(PID,M_Year,MT_Date) %>% 
  filter(MAR_Action %!in% c("MISSED","STOPPED","HELD")) %>% 
  distinct(PID,M_Year,M_Date,MT_Date,Simple_Generic_Med_Name,.keep_all = TRUE) %>% 
  mutate(M_Time2 = ifelse(daydiff(M_Year,BC_Year,MT_Date,BC_Date)==10, "Time3",NA),
         M_Time1 = ifelse(daydiff(M_Year,BC_Year,MT_Date,BC_Date)==3, "Time2",NA),
         M_Time0 = ifelse(daydiff(M_Year,BC_Year,MT_Date,BC_Date)<3,"Time1",NA)) %>%
  mutate(M_Group = coalesce(M_Time0,M_Time1,M_Time2)) %>% 
  mutate(Prev_Vanco = ifelse(str_detect(tolower(Simple_Generic_Med_Name),"vancomycin"),1,0)) %>% 
  filter(M_Group %in% c("Time1","Time2","Time3")) %>% 
  group_by(PID,M_Group) %>% 
  arrange(PID,desc(Prev_Vanco)) %>% 
  distinct(PID,.keep_all = TRUE) %>% 
  select(PID, Prev_Vanco) %>% 
  mutate(Prev_Vanco=ifelse(Prev_Vanco==1,1,0))

data_combined_mix <- data_combined %>% 
                     left_join(select(Sources_Blood,c("PID","Source_blood")),by="PID") %>% 
                     left_join(select(Kidney_Impairment,c("PID","Kidney_Impairment_before")),by="PID") %>% 
                     left_join(select(Previous_vanco,c("PID","Prev_Vanco")),by="PID")
data_combined_mix_c <- na.omit(data_combined_mix)#400

write.csv(data_combined_mix,"data_combined_mix.csv")
write.csv(data_combined_mix_c,"data_combined_mix_c.csv")


#====================================================================================
###########################################################################################
# Attach Socio-Economic Variables from other sources
# Income, Crime, Area, Medical Assesibility, Insurance
###########################################################################################
# Income
setwd("/blue/prosperi/inyoungjun/PSB2023")
fips_county_FL<-read.csv("fips_county_FL.csv")
zip_county_FL<-read.csv("ZIP-COUNTY-FIPS_2017-06.csv") %>% 
  filter(STATE == "FL") %>% 
  mutate(ZIP_three = substr(ZIP,1,3)) %>% 
  group_by(COUNTYNAME) %>% 
  distinct(ZIP_three,.keep_all = TRUE) %>% 
  rename(FIPS_name = STCOUNTYFP)

library(tidycensus)
#https://cran.r-project.org/web/packages/tidycensus/tidycensus.pdf
#Ref: https://rpubs.com/kitadasmalley/getACSVars
census_api_key("35f8b75a23bef45fa36ac7e7bc67db771a4e1efd",install=TRUE, overwrite=TRUE) #deleted the key for protecting privacy
readRenviron("~/.Renviron")

Income<-NULL
year <- seq(2012,2019,by=1)
for (i in year){
  FL_Dem <- get_acs(geography = "county", year=i,
                    state = "FL", geometry = TRUE,
                    variables = c(population = "B02001_001", #Total
                                  median.household.income = "B19013_001")) #median household income in past 12 months
  FL_Dem_Pct<-as.data.frame(FL_Dem)[,c(1,3:4)]%>%
    spread(variable, estimate)%>%
    mutate(year = i)
  Income <- rbind(Income, FL_Dem_Pct)
}
table(Income$year)
head(Income)#

# Crime & Urban
crime_county_FL<-read.csv("combined_crime_1319.csv") %>% select(year,area,county,violent_crime)
head(crime_county_FL) #need to normalized upon population

# Medical Accessibility
medical_access_county_FL<-read.csv("combined_medical_accessibility_1219.csv") 
#%>% select(year,area,county,violent_crime)
head(medical_access_county_FL)

# Insurance
insurance_county_FL<-read.csv("combined_uninsured_1219.csv") %>% rename(year = Year,county=County) %>% select(-FIPS)
head(insurance_county_FL)

########## Merge
FL_Dem_match <- Income %>%
  mutate(FIPS_name = as.numeric(substr(GEOID,1,5))) %>% 
  left_join(fips_county_FL) %>% 
  select(year,FIPS_name, county_name, 
         median.household.income,population) %>%
  mutate(population = as.numeric(population)) %>% 
  rename(county=county_name) %>% 
  left_join(zip_county_FL) %>% 
  arrange(year,ZIP_three) %>% 
  select(year,county,ZIP, ZIP_three,FIPS_name,median.household.income,population) %>% 
  left_join(crime_county_FL,by=c("year","county")) %>% 
  mutate(violent_crime = as.numeric(violent_crime)) %>% 
  mutate(crime_rate = round(violent_crime/population*100,2)) %>% 
  left_join(medical_access_county_FL,by=c("year","county")) %>% 
  left_join(insurance_county_FL,by=c("year","county")) %>% 
  arrange(county) %>% 
  mutate(area = ifelse(area == "Urban",1,0))
head(FL_Dem_match,30)

FL_Dem_average <- FL_Dem_match %>% 
  group_by(year,ZIP_three) %>%
  summarise(Area = mean(area,na.rm = TRUE),
            Median_household_income = mean(median.household.income,na.rm=TRUE),
            Crime_rate = mean(crime_rate,na.rm=TRUE),
            Medical_rate = mean(medical_rate,na.rm=TRUE),
            Uninsured_perc = mean(Uninsured_perc,na.rm=TRUE)) %>% 
  mutate(across(c(Median_household_income:Uninsured_perc), ~round(.,2))) %>% 
  mutate(Area = ifelse(Area>=0.5,1,0)) %>% arrange(ZIP_three) %>% 
  rename(BC_Year = year,Grouped_Zip_Code=ZIP_three)
FL_Dem_average
#saveRDS(FL_Dem_average,"FL_Dem_average.RDS")

##### Attach socio variables 
head(data_combined_mix_c)
colnames(data_Invasive_MRSA)
data_combined_mix_c_socio <- data_combined_mix_c %>% 
                             left_join(select(data_Invasive_MRSA,c("PID","BC_Year","Grouped_Zip_Code")),by="PID") %>%
                             select(-M_Group.x,-M_Group.y) %>%
                             group_by(PID) %>% 
                             distinct(PID,.keep_all = TRUE) %>% 
                             left_join(FL_Dem_average,by=c("BC_Year","Grouped_Zip_Code"))#400
data_combined_mix_c_socio_c <- na.omit(data_combined_mix_c_socio) %>% 
                               mutate(MRSA_Age_bin = ifelse(MRSA_Age>=65,1,0)) %>% 
                               mutate(Income_bin = ifelse(Median_household_income>=48300,1,0)) %>% 
                               mutate(Facility_bin = ifelse(Medical_rate>=180,1,0)) %>% 
                               mutate(Insurance_bin = ifelse(Uninsured_perc >=21.8,1,0))#395

# Define a function for min-max scaling
min_max_scaling <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
#MRSA_Age, CI_score, Median_household_income, Medical_rate, Uninsured_perc
data_combined_mix_c_socio_c$MRSA_Age_s = min_max_scaling(data_combined_mix_c_socio_c$Median_household_income)
data_combined_mix_c_socio_c$Comorbidity = min_max_scaling(data_combined_mix_c_socio_c$CI_score)
data_combined_mix_c_socio_c$Income = min_max_scaling(data_combined_mix_c_socio_c$Median_household_income)
data_combined_mix_c_socio_c$Facility = min_max_scaling(data_combined_mix_c_socio_c$Medical_rate)
data_combined_mix_c_socio_c$Insurance = min_max_scaling(data_combined_mix_c_socio_c$Uninsured_perc)
summary(data_combined_mix_c_socio_c)
#write.csv(data_combined_mix_c_socio_c,"data_combined_mix_c_socio_c.csv")
colnames(data_combined_mix_c_socio_c)
table(data_combined_mix_c_socio_c$Outcome,useNA = "ifany")
365/395
colnames(data_combined_mix_c_socio_c)
summary(data_combined_mix_c_socio_c)
#Table 1
library(table1)
colnames(data_combined_mix_c_socio_c)
table1(~as.factor(ICU)+as.factor(MDR)+ as.factor(Nephrotoxicity_any)+as.factor(Vanco_T1)+as.factor(Vanco_T2)+as.factor(Vanco_T3)+
         as.factor(Source_blood)+as.factor(Kidney_Impairment_before)+as.factor(Prev_Vanco)+as.factor(Area),data=data_combined_mix_c_socio_c)

#Figure2 - ggplot by year

# Create the ggplot figure with theme_blank()#395

ggplot(data_combined_mix_c_socio_c, aes(x = factor(BC_Year))) +
  geom_bar(fill = "navy", color = "black") +
  geom_text(stat = "count", aes(label = stat(count)), vjust = -0.5, size = 4) +
  labs(x = "Year",
       y = "Count") +
  theme_void() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14,angle=90),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0),
    axis.text.y = element_text(size = 12, angle = 0, hjust = 1),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )
#Let's draw in Rate
# Calculate the total count first
total_count <- nrow(data_combined_mix_c_socio_c)

ggplot(data_combined_mix_c_socio_c, aes(x = factor(BC_Year))) +
  geom_bar(fill = "navy", color = "black") +
  geom_text(aes(label = scales::percent(after_stat(count/total_count))), stat = "count", vjust = -0.5, size = 4) +
  labs(x = "Year",
       y = "Count") +
  theme_void() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14, angle=90),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0),
    axis.text.y = element_text(size = 12, angle = 0, hjust = 1),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )

#revise Figure 2 with one bar graph and two line graphs
# Replace "my_data" with the actual name of your data frame
d_year_count <- d %>% group_by(D_Year) %>% count(PID) %>% mutate(value = 1) %>% group_by(D_Year) %>% summarize(sum=sum(value))

Year <- c(2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)
Invasive_MRSA_Case <- c(29, 45, 51, 53, 67, 55, 69, 26)#case
Patient_in_EHR <- c(29516, 31210, 32356, 32621, 32301, 31493, 30096, 23632)#EHR patient
alachua <-c(251596,252585,255606,259215,264127,266501,269190,269489)
duval <- c(880349,886277,897094,911106,927038,939167,950469,959853)
Population <- alachua + duval #alachua + duval county population

my_data <- data.frame(Year, 
                      Invasive_MRSA_Case, 
                      Population, 
                      Patient_in_EHR)

my_data$Year <- as.factor(my_data$Year)

ggplot(my_data, aes(x=Year)) +
  geom_bar(aes(y=Invasive_MRSA_Case, fill = "Invasive MRSA Case"), stat="identity") +
  geom_line(aes(y=Patient_in_EHR/350, color = "Patient in EHR", group=1), size=1, linetype="twodash") +
  geom_line(aes(y=Population/10000, color = "Population", group=1), size=1) +
  scale_fill_manual(values=c("Invasive MRSA Case"="#1B9E77")) +
  scale_color_manual(values=c("Patient in EHR"="#D95F02", "Population"="#7570B3")) +
  scale_y_continuous(limits = c(0, 70),
                     sec.axis = sec_axis(~.*350, name = "Number of Patients in EHR"), 
                     name = "Invasive MRSA Case") +
  scale_y_continuous(sec.axis = sec_axis(~.*10000, name = "Population in UF Health Regions")) +
  labs(x = "Year", 
       y = "Number of Patients in EHR",
       fill = NULL, 
       color = NULL) +
  theme_minimal() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom") 

########
colnames(data_combined_mix_c_socio_c)
#Table 2 - logistic regression and FACTS
model <- glm(Outcome ~ MRSA_Age_bin, data = data_combined_mix_c_socio_c, family=binomial())
model <- glm(Outcome ~ MRSA_Age_bin + Comorbidity + Prev_Vanco + MDR + Source_blood + ICU, data = data_combined_mix_c_socio_c, family=binomial())
model <- glm(Outcome ~ Race, data = data_combined_mix_c_socio_c, family=binomial())
model <- glm(Outcome ~ Race + Source_blood + ICU + Prev_Vanco + MDR + Comorbidity, data = data_combined_mix_c_socio_c, family=binomial())
model <- glm(Outcome ~ Income_bin , data = data_combined_mix_c_socio_c, family=binomial())
model <- glm(Outcome ~ Income_bin + Comorbidity + Race + Area + Sex, data = data_combined_mix_c_socio_c, family=binomial())
model <- glm(Outcome ~ Sex , data = data_combined_mix_c_socio_c, family=binomial())
model <- glm(Outcome ~ Sex + Comorbidity + Prev_Vanco + MDR + ICU + Source_blood, data = data_combined_mix_c_socio_c, family=binomial())

# Extract coefficient estimates
coef_estimates <- coef(model)
# Calculate odds ratio
odds_ratio <- exp(coef_estimates)
# Calculate confidence intervals for odds ratio
conf_intervals <- confint(model)
# Create a data frame with odds ratio and confidence intervals
result_df <- data.frame(
  Odds_Ratio = odds_ratio,
  CI_Lower = exp(conf_intervals[, "2.5 %"]),
  CI_Upper = exp(conf_intervals[, "97.5 %"])
)
round(result_df,2)
