library(tictoc)
library(data.table)
library(writexl)
library(tableone)

filepath_load <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

tic()

load(file=file.path(filepath_load,"AfA_VL_courier_MHD.RData"))

setorder(DTrna,"patient","rna_d")

# formatting
DTrna[,sex:=as.character(sex)]
DTrna[sex=="1",sex:="Male"]
DTrna[sex=="2",sex:="Female"]
DTrna[,`:=`(age_current_cat=cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
            calyear_current_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
            sex=factor(sex,levels=c("Male","Female")),
            VL_400=as.numeric(rna_v<=400),
            courier=as.character(courier),
            mhd_ind=as.character(mhd_ind))]
DTrna[courier==0,courier:="No"]
DTrna[courier==1,courier:="Yes"]
DTrna[,courier:=factor(courier,levels=c("No","Yes"))]
DTrna[mhd_ind==0,mhd_ind:="No"]
DTrna[mhd_ind==1,mhd_ind:="Yes"]
DTrna[,mhd_ind:=factor(mhd_ind,levels=c("No","Yes"))]

df_out <- CreateTableOne(vars = c("courier","mhd_ind","sex","age_current_cat","age_current","calyear_current_cat","art_type"),
                                     strata="VL_400",data=DTrna,test=FALSE,addOverall=TRUE)
df_out <- print(df_out,nonnormal="age_current",showAllLevels=TRUE,printToggle=FALSE)
df_out <- data.table(cbind(row.names(df_out),df_out))

write_xlsx(df_out,path=file.path(filepath_tables,"descriptive_by_VL.xlsx"))

toc()