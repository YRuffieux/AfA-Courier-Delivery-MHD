# descriptive table by individual, overall and stratified by ART delivery type at first VL test
# default: removing people who are in Bonitas at first RNA test

library(tictoc)
library(data.table)
library(writexl)
library(tableone)

filepath_load <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

with_BON <- "no"                # "no", "yes", "only"

tic()

load(file=file.path(filepath_load,"AfA_VL.RData"))

DTu <- DTrna[!is.na(rna_v)]               # removing 'fake' tests
setorder(DTu,"patient","rna_d")

# "ever" indicators
DTu[,`:=`(courier_ever=as.numeric(any(courier==1)),mhd_ever=as.numeric(any(mhd_ind==1)),
          last_rna_d=max(rna_d),N_rna=.N,delta_rna=as.numeric(rna_d-data.table::shift(rna_d,type="lag"))),
    by="patient"]

# courier status at 50% or more tests
DTu[,N_courier:=sum(courier),by="patient"]
DTu[,courier_mode:=as.numeric(N_courier/N_rna>=0.5)]

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
          courier_base=as.character(courier),
          courier_mode=as.character(courier_mode),
          mhd_ever=as.character(mhd_ever),
          art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")),
          scheme_code=NULL)]
DTu[,courier:=NULL]
DTu[courier_ever==0,courier_ever:="No"]
DTu[courier_ever==1,courier_ever:="Yes"]
DTu[,courier_ever:=factor(courier_ever,levels=c("No","Yes"))]
DTu[courier_base==0,courier_base:="No"]
DTu[courier_base==1,courier_base:="Yes"]
DTu[,courier_base:=factor(courier_base,levels=c("No","Yes"))]
DTu[courier_mode==0,courier_mode:="No"]
DTu[courier_mode==1,courier_mode:="Yes"]
DTu[,courier_mode:=factor(courier_mode,levels=c("No","Yes"))]
DTu[mhd_ever==0,mhd_ever:="No"]
DTu[mhd_ever==1,mhd_ever:="Yes"]
DTu[,mhd_ever:=factor(mhd_ever,levels=c("No","Yes"))]
DTu[!scheme_code_base%in%c("BON","PLM"),scheme_code_base:="Other"]

if(with_BON=="no")
{
  # removing BON
  DTu <- DTu[scheme_code_base!="BON"]
  DTu[,scheme_code_base:=factor(scheme_code_base,levels=c("PLM","Other"))]
} else if(with_BON=="only")
{
  DTu <- DTu[scheme_code_base=="BON"]
} else if(with_BON=="yes")
{
  DTu[,scheme_code_base:=factor(scheme_code_base,levels=c("PLM","Other","BON"))]
}

remove_space <- function(x) gsub("\\( ","\\(",x)

# by base courier delivery yes/no
courier_df <- CreateTableOne(vars = c("mhd_ever","sex","age_base_cat","age_base","calyear_base_cat","art_type_cf","scheme_code_base"),
                             strata="courier_base",
                             test=FALSE,addOverall=TRUE,includeNA=TRUE,data=DTu)
courier_df <- print(courier_df,nonnormal="age_base",showAllLevels=TRUE,printToggle=FALSE)
courier_df <- data.table(data.frame(cbind(row.names(courier_df),courier_df)))
cols <- c("Overall","No","Yes")
courier_df[,(cols):=lapply(.SD,remove_space),.SDcols=cols]
courier_df <- courier_df[,.(V1,level,No,Yes,Overall)]
rm(cols)
if(with_BON=="no")
{
  write_xlsx(courier_df,path=file.path(filepath_tables,"descriptive_individuals_by_courier_status.xlsx"))
} else if(with_BON=="only")
{
  write_xlsx(courier_df,path=file.path(filepath_tables,"Scheme-specific","descriptive_individuals_by_courier_status_BON.xlsx"))
} else if(with_BON=="yes")
{
  write_xlsx(courier_df,path=file.path(filepath_tables,"Sensitivity","descriptive_individuals_by_courier_status_allschemes.xlsx"))
}

toc()