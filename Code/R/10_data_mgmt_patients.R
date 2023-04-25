# selection of eligible individuals for analysis
# identifying first diagnosis of a mental illness per individual
# all left-truncation/right-censoring done here

library(data.table)
library(tictoc)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
#filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean_old"
#filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed_old"

load(file=file.path(filepath_read,"tblBAS.RDAta"))
load(file=file.path(filepath_read,"ICD10_F.RDAta"))
load(file=file.path(filepath_read,"tblFUPwide.RDAta"))
load(file=file.path(filepath_read,"tblLTFU.RDAta"))
load(file=file.path(filepath_read,"tblREGIMEN.RData"))

tic()

min_age <- 15

DTbas <- tblBAS[,.(patient,birth_d,enrol_d,sex)]

### Patient selection ###

DTbas <- merge(DTbas,tblFUPwide,by="patient",all.x=TRUE)
#DTbas <- DTbas[!last_scheme_code%in%c("GMS","GMT")]

N <- uniqueN(DTbas,"patient")
print(paste0("Starting number of individuals: ",N))

DTbas <- DTbas[!is.na(birth_d) & year(birth_d)<=2022]
N_prev <- N
N <- uniqueN(DTbas,"patient")
print(paste0("*after excluding invididuals with missing or invalid date of birth: ",N, " (",N-N_prev,")"))

DTbas <- DTbas[!is.na(start)]
N_prev <- N
N <- uniqueN(DTbas,"patient")
print(paste0("*after excluding individuals with no follow-up during study period: ",N, " (",N-N_prev,")"))

DTbas <- DTbas[,start:=pmax(birth_d+min_age*365.25,start)]  # left-truncation at age 15
DTbas <- DTbas[start<end]
N_prev <- N
N <- uniqueN(DTbas,"patient")
print(paste0("*after excluding individuals aged below ",min_age," at their end of follow-up: ",N, " (",N-N_prev,")"))

setorder(tblREGIMEN,"patient","moddate")
DTbas <- tblREGIMEN[art==1 & art_type!="Other",.(patient,first_art_d=moddate)][DTbas,on="patient",mult="first"]
DTbas <- DTbas[!is.na(first_art_d)]
DTbas[,start:=pmax(first_art_d,start)]  # left-truncation at ART initiation
N_prev <- N
N <- uniqueN(DTbas,"patient")
print(paste0("*after excluding individuals not on ART: ",N, " (",N-N_prev,")"))

### Identifying mental illnesses ###

setorder(tblICD10_F,"patient","icd10_date")
DTbas <- tblICD10_F[icd_23!=17,.(patient,mhd_d=icd10_date)][DTbas,on="patient",mult="first"]                  # any mental illness (excl. smoking)

### Saving patient-level dataset ###

DTbas <- DTbas[,.(patient,sex,birth_d,mhd_d,start,end)]

save(DTbas,file=file.path(filepath_write,"AfA_base.RData"))

toc()
