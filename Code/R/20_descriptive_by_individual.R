# descriptive table, by individual, stratified by MHD status
# breakdown of MHDs by type and sex

library(tictoc)
library(data.table)
library(writexl)
library(tableone)

filepath_load <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

tic()

load(file=file.path(filepath_load,"AfA_VL_courier_MHD.RData"))

which_mhd <- c("mhd","org","sub","psy","bip","dep","anx","adj")

DTu <- copy(DTrna)
setorder(DTu,"patient","rna_d")

# "ever" indicators
DTu[,`:=`(courier_ever=as.numeric(any(courier==1)),mhd_ever=as.numeric(any(mhd_ind==1)),org_ever=as.numeric(any(org_ind==1)),
          sub_ever=as.numeric(any(sub_ind==1)),psy_ever=as.numeric(any(psy_ind==1)),bip_ever=as.numeric(any(bip_ind==1)),
          dep_ever=as.numeric(any(dep_ind==1)),anx_ever=as.numeric(any(anx_ind==1)),adj_ever=as.numeric(any(adj_ind==1)),
          last_rna_d=max(rna_d),N_rna=.N,delta_rna=as.numeric(rna_d-data.table::shift(rna_d,type="lag"))),
    by="patient"]

time_between_rna <- DTu[,quantile(delta_rna,p=c(0.25,0.5,0.75),na.rm=T)]

# deduplication, keeping first RNA record for each individual
DTu <- unique(DTu,by="patient")
setnames(DTu,"age_current","age_base")

time_followup <- DTu[,quantile(as.numeric(last_rna_d-rna_d),p=c(0.25,0.5,0.75))]
N_measurements <- DTu[,quantile(N_rna,p=c(0.25,0.5,0.75))]

# formatting
DTu[,sex:=as.character(sex)]
DTu[sex=="1",sex:="Male"]
DTu[sex=="2",sex:="Female"]
DTu[,`:=`(age_base_cat=cut(age_base,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
          calyear_base_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
          sex=factor(sex,levels=c("Male","Female")),
          courier_ever=as.character(courier_ever),
          art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+NRTI","PI+2NRTI")))]
DTu[courier_ever==0,courier_ever:="No"]
DTu[courier_ever==1,courier_ever:="Yes"]
DTu[,courier_ever:=factor(courier_ever,levels=c("No","Yes"))]

# looping over types of mental disorders
MHD_count_df <- data.frame(NULL)

N_total <- DTu[,.N]
N_men <- DTu[sex=="Male",.N]
N_women <- DTu[sex=="Female",.N]

for(v in paste0(which_mhd,"_ever"))
{
  DTtemp <- copy(DTu)
  
  if(v=="mhd_ever")   # creating single table with/without mental illness
  {
    overall_df <- CreateTableOne(vars = c("sex","age_base_cat","age_base","calyear_base_cat","art_type_cf","courier_ever"),
                                 strata="mhd_ever",data=DTtemp,test=FALSE,addOverall=TRUE)
    overall_df <- print(overall_df,nonnormal="age_base",showAllLevels=TRUE,printToggle=FALSE)
    overall_df <- data.table(data.frame(cbind(row.names(overall_df),overall_df)))
    colnames(overall_df)[4] <- "No MHD"
    colnames(overall_df)[5] <- "MHD"
    write_xlsx(overall_df,path=file.path(filepath_tables,"descriptive_by_individual.xlsx"))
    rm(overall_df)
  }
  N_exp_total <- DTtemp[get(v)==1,.N]
  N_exp_men <- DTtemp[get(v)==1 & sex=="Male",.N]
  N_exp_women <- DTtemp[get(v)==1 & sex=="Female",.N]
  MHD_count_df <- rbind(MHD_count_df,data.frame(exp=v,
                                                men=paste0(N_exp_men," (",format(round(100*N_exp_men/N_men,digits=1),nsmall=1),"%)"),
                                                women=paste0(N_exp_women," (",format(round(100*N_exp_women/N_women,digits=1),nsmall=1),"%)"),
                                                total=paste0(N_exp_total," (",format(round(100*N_exp_total/N_total,digits=1),nsmall=1),"%)")))
  rm(DTtemp,N_exp_total,N_exp_men,N_exp_women)
}     

MHD_count_df <- data.table(MHD_count_df)
write_xlsx(MHD_count_df,path=file.path(filepath_tables,"descriptive_by_MHD.xlsx"))

toc()