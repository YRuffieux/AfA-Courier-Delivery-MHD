# descriptive table by VL measurement, stratified by viral suppression status

library(tictoc)
library(data.table)
library(writexl)
library(tableone)

filepath_load <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

tic()

load(file=file.path(filepath_load,"AfA_VL.RData"))

setorder(DTrna,"patient","rna_d")

# formatting
DTrna[,sex:=as.character(sex)]
DTrna[sex=="1",sex:="Male"]
DTrna[sex=="2",sex:="Female"]
DTrna[,`:=`(age_current_cat=cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
            calyear_current_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
            sex=factor(sex,levels=c("Male","Female")),
            art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")),
            VLS_50=factor(ifelse(rna_v<50,"Suppressed","Unsuppressed"),levels=c("Suppressed","Unsuppressed")),
            VLS_400=factor(ifelse(rna_v<400,"Suppressed","Unsuppressed"),levels=c("Suppressed","Unsuppressed")),
            VLS_1000=factor(ifelse(rna_v<1000,"Suppressed","Unsuppressed"),levels=c("Suppressed","Unsuppressed")),
            courier=as.character(courier),
            mhd_ind=as.character(mhd_ind))]
DTrna[courier==0,courier:="No"]
DTrna[courier==1,courier:="Yes"]
DTrna[,courier:=factor(courier,levels=c("Yes","No"))]
DTrna[mhd_ind==0,mhd_ind:="No"]
DTrna[mhd_ind==1,mhd_ind:="Yes"]
DTrna[,mhd_ind:=factor(mhd_ind,levels=c("No","Yes"))]

# number of tests per patient
DTrna[,`:=`(N_tests=.N,N_tests_courier=sum(courier=="Yes"),N_tests_noncourier=sum(courier=="No")),by="patient"]
DTu <- unique(DTrna[,.(patient,N_tests,N_tests_courier,N_tests_noncourier)])
print("Median number of VL measurements per patient while on courier delivery: ")
print(DTu[,quantile(N_tests_courier,p=c(0.25,0.5,0.75))])
print("Median number of VL measurements per patient while not on courier delivery: ")
print(DTu[,quantile(N_tests_noncourier,p=c(0.25,0.5,0.75))])
print("Median number of VL measurements per patient overall: ")
print(DTu[,quantile(N_tests,p=c(0.25,0.5,0.75))])
rm(DTu)

# stratified by courier status (yes/no)
df_out <- CreateTableOne(vars = c("mhd_ind","sex","age_current_cat","age_current","calyear_current_cat","art_type_cf","VLS_400"),
                                     strata="courier",data=DTrna,test=FALSE,addOverall=TRUE,includeNA=TRUE)
df_out <- print(df_out,nonnormal="age_current",showAllLevels=TRUE,printToggle=FALSE)
df_out <- data.table(cbind(row.names(df_out),df_out))

write_xlsx(df_out,path=file.path(filepath_tables,"descriptive_by_courier_status.xlsx"))

# stratified by calendar period and courier status (yes/no)
df_out <- CreateTableOne(vars = c("mhd_ind","sex","age_current_cat","age_current","art_type_cf","VLS_400"),
                         strata=c("courier","calyear_current_cat"),data=DTrna,test=FALSE,addOverall=TRUE,includeNA=TRUE)
df_out <- print(df_out,nonnormal="age_current",showAllLevels=TRUE,printToggle=FALSE)
df_out <- data.table(cbind(row.names(df_out),df_out))

write_xlsx(df_out,path=file.path(filepath_tables,"descriptive_by_period_and_courier_status.xlsx"))

# stratified by courier pharmacy
df_out <- CreateTableOne(vars = c("mhd_ind","sex","age_current_cat","age_current","calyear_current_cat","art_type_cf","VLS_400"),
                         strata="courier_cat",data=DTrna,test=FALSE,addOverall=TRUE,includeNA=TRUE)
df_out <- print(df_out,nonnormal="age_current",showAllLevels=TRUE,printToggle=FALSE)
df_out <- data.table(cbind(row.names(df_out),df_out))

write_xlsx(df_out,path=file.path(filepath_tables,"descriptive_by_courier_pharmacy.xlsx"))

rm(df_out)

toc()