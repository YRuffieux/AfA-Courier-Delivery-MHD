# descriptive table by individual, overall and stratified by ART delivery type, and by baseline medical scheme

library(tictoc)
library(data.table)
library(writexl)
library(tableone)

filepath_load <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

tic()

load(file=file.path(filepath_load,"AfA_VL.RData"))

DTu <- DTrna[!is.na(rna_v)]               # removing 'fake' tests
setorder(DTu,"patient","rna_d")

# "ever" indicators
DTu[,`:=`(courier_ever=as.numeric(any(courier==1)),mhd_ever=as.numeric(any(mhd_ind==1)),
          last_rna_d=max(rna_d),N_rna=.N,delta_rna=as.numeric(rna_d-data.table::shift(rna_d,type="lag"))),
    by="patient"]

print("Median time between VL tests + IQR:") 
print(DTu[,quantile(delta_rna,p=c(0.25,0.5,0.75),na.rm=T)])

# deduplication, keeping first RNA record for each individual
DTu <- unique(DTu,by="patient")
setnames(DTu,"age_current","age_base")

print("Median follow-up time + IQR:")
print(DTu[,quantile(as.numeric(last_rna_d-rna_d),p=c(0.25,0.5,0.75))])

print("Median number of VL tests + IQR:")
print(DTu[,quantile(N_rna,p=c(0.25,0.5,0.75))])

# formatting
DTu[,sex:=as.character(sex)]
DTu[sex=="1",sex:="Male"]
DTu[sex=="2",sex:="Female"]
DTu[,`:=`(age_base_cat=cut(age_base,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
          calyear_base_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
          sex=factor(sex,levels=c("Male","Female")),
          courier_ever=as.character(courier_ever),
          mhd_ever=as.character(mhd_ever),
          art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")),
          scheme_code=NULL)]
DTu[courier_ever==0,courier_ever:="No"]
DTu[courier_ever==1,courier_ever:="Yes"]
DTu[,courier_ever:=factor(courier_ever,levels=c("No","Yes"))]
DTu[mhd_ever==0,mhd_ever:="No"]
DTu[mhd_ever==1,mhd_ever:="Yes"]
DTu[,mhd_ever:=factor(mhd_ever,levels=c("No","Yes"))]
DTu[!scheme_code_base%in%c("BON","PLM"),scheme_code_base:="Other"]
DTu[,scheme_code_base:=factor(scheme_code_base,levels=c("BON","PLM","Other"))]

remove_space <- function(x) gsub("\\( ","\\(",x)

# no stratification
overall_df <- CreateTableOne(vars = c("courier_ever","mhd_ever","sex","age_base_cat","age_base","calyear_base_cat","art_type_cf"),data=DTu)
overall_df <- print(overall_df,nonnormal="age_base",showAllLevels=TRUE,printToggle=FALSE)
overall_df <- data.table(data.frame(cbind(row.names(overall_df),overall_df)))
overall_df[,Overall:=remove_space(Overall)]
write_xlsx(overall_df,path=file.path(filepath_tables,"descriptive_individuals.xlsx"))

# by "ever" courier delivery yes/no
courier_df <- CreateTableOne(vars = c("mhd_ever","sex","age_base_cat","age_base","calyear_base_cat","art_type_cf","scheme_code_base"),
                             strata="courier_ever",
                             test=FALSE,addOverall=TRUE,includeNA=TRUE,data=DTu)
courier_df <- print(courier_df,nonnormal="age_base",showAllLevels=TRUE,printToggle=FALSE)
courier_df <- data.table(data.frame(cbind(row.names(courier_df),courier_df)))
cols <- c("Overall","No","Yes")
courier_df[,(cols):=lapply(.SD,remove_space),.SDcols=cols]
courier_df <- courier_df[,.(V1,level,No,Yes,Overall)]
write_xlsx(courier_df,path=file.path(filepath_tables,"descriptive_individuals_by_courier_status.xlsx"))

# by medical scheme
scheme_df <- CreateTableOne(vars = c("courier_ever","mhd_ever","sex","age_base_cat","age_base","calyear_base_cat","art_type_cf"),
                            strata="scheme_code_base",
                            test=FALSE,addOverall=TRUE,includeNA=TRUE,data=DTu)
scheme_df <- print(scheme_df,nonnormal="age_base",showAllLevels=TRUE,printToggle=FALSE)
scheme_df <- data.table(data.frame(cbind(row.names(scheme_df),scheme_df)))
cols <- c("Overall","BON","PLM","Other")
scheme_df[,(cols):=lapply(.SD,remove_space),.SDcols=cols]
scheme_df <- scheme_df[,.(V1,level,BON,PLM,Other,Overall)]
write_xlsx(scheme_df,path=file.path(filepath_tables,"descriptive_individuals_by_medical_scheme.xlsx"))


toc()