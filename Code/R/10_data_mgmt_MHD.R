# selection of eligible individuals for analysis
# identifying first diagnosis of a mental illness per individual

library(data.table)
library(tictoc)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"

load(file=file.path(filepath_read,"tblBAS.RDAta"))
load(file=file.path(filepath_read,"ICD10_F.RDAta"))
load(file=file.path(filepath_read,"tblFUPwide.RDAta"))
load(file=file.path(filepath_read,"tblLTFU.RDAta"))

tic()

min_age <- 15

DTmhd <- tblBAS[,.(patient,birth_d,enrol_d,sex)]

### Patient selection ###

DTmhd <- merge(DTmhd,tblFUPwide,by="patient",all.x=TRUE)
DTmhd <- DTmhd[!last_scheme_code%in%c("GMS","GMT")]

N <- uniqueN(DTmhd,"patient")
print(paste0("Starting number of individuals: ",N))

DTmhd <- DTmhd[!is.na(birth_d) & year(birth_d)<=2022]
N_prev <- N
N <- uniqueN(DTmhd,"patient")
print(paste0("*after excluding invididuals with missing or invalid date of birth: ",N, " (",N-N_prev,")"))

DTmhd <- DTmhd[!is.na(start)]
N_prev <- N
N <- uniqueN(DTmhd,"patient")
print(paste0("*after excluding individuals with no follow-up during study period: ",N, " (",N-N_prev,")"))

DTmhd <- DTmhd[as.numeric(end-birth_d)/365.25>=min_age]
N_prev <- N
N <- uniqueN(DTmhd,"patient")
print(paste0("*after excluding individuals aged below ",min_age," at their end of follow-up: ",N, " (",N-N_prev,")"))

### Identifying mental illnesses ###

setorder(tblICD10_F,"patient","icd10_date")
DTmhd <- tblICD10_F[icd_23!=17,.(patient,mhd_d=icd10_date)][DTmhd,on="patient",mult="first"]                  # any mental illness (excl. smoking)
DTmhd <- tblICD10_F[icd_23%in%0:9,.(patient,org_d=icd10_date)][DTmhd,on="patient",mult="first"]               # organic (F00-F09)
DTmhd <- tblICD10_F[icd_23%in%c(10:16,18,19),.(patient,sub_d=icd10_date)][DTmhd,on="patient",mult="first"]    # substance use excl. smoking (F10-F16, F18, F19)
DTmhd <- tblICD10_F[icd_23%in%20:29,.(patient,psy_d=icd10_date)][DTmhd,on="patient",mult="first"]             # psychotic (F20-F29)
DTmhd <- tblICD10_F[icd_23==31,.(patient,bip_d=icd10_date)][DTmhd,on="patient",mult="first"]                  # bipolar (F31)
DTmhd <- tblICD10_F[icd_23%in%c(32,33)|substr(icd10_code,1,5)=="F34.1",.(patient,dep_d=icd10_date)][DTmhd,on="patient",mult="first"]    # depression (F32, F33, F34.1)
DTmhd <- tblICD10_F[icd_23%in%40:49,.(patient,anx_d=icd10_date)][DTmhd,on="patient",mult="first"]             # anxiety (F40-F49)
DTmhd <- tblICD10_F[substr(icd10_code,1,5)=="F43.2",.(patient,adj_d=icd10_date)][DTmhd,on="patient",mult="first"]   # adjustment disorder (F43.2)

DTmhd <- DTmhd[,.(patient,sex,birth_d,mhd_d,org_d,sub_d,psy_d,bip_d,dep_d,anx_d,adj_d,start,end)]

save(DTmhd,file=file.path(filepath_write,"AfA_MHD.RData"))

toc()
