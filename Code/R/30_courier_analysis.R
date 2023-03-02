# logistic regression on viral load suppression (no/yes)
# robust standard errors: https://www.r-bloggers.com/2021/05/clustered-standard-errors-with-r/

library(data.table)
library(lmtest)
library(sandwich)
library(tictoc)
library(writexl)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

rf_vect <- c("courier","mhd_ind","sex","age_current_cat","calyear_current_cat","art_type_cf")

tic()

load(file=file.path(filepath_read,"AfA_VL_courier_MHD.RData"))

VLS_threshold <- 50

# formatting
DTrna[,`:=`(vls_ind=as.numeric(rna_v<=VLS_threshold),
            age_current_cat=cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
            calyear_current_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
            art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+NRTI","PI+2NRTI")))]
DTrna[,age_current_cat:=relevel(age_current_cat,"[40,50)")]

stopifnot(DTrna[,all(!is.na(art_type_cf))])

df_out <- data.frame(NULL)

# unadjusted regressions
for(v in rf_vect)
{
  reg_formula <- as.formula(paste0("vls_ind ~",v))
  lreg <- glm(reg_formula,data=DTrna,family="binomial")
  lreg_coeffs <- coeftest(lreg,vcov=vcovCL,cluster=~patient)  # from sandwich package
  lreg_CIs <- coefci(lreg,vcov=vcovCL,cluster=~patient)       # from lmtest package
  rows <- rownames(lreg_coeffs)[grep(v,rownames(lreg_coeffs))]
  df_out <- rbind(df_out,data.frame(risk_factor=rows,OR_uni=paste0(round(exp(lreg_coeffs[rows,1]),digits=2)," (",
                                                            round(exp(lreg_CIs[rows,1]),digits=2),"-",
                                                            round(exp(lreg_CIs[rows,2]),digits=2),")")))
  rm(reg_formula,lreg,lreg_coeffs,lreg_CIs,rows)
}

# adjusted regression, no interaction
reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+")))
lreg <- glm(reg_formula, data=DTrna,family="binomial")
lreg_coeffs <- coeftest(lreg,vcov=vcovCL,cluster=~patient)
lreg_CIs <- coefci(lreg,vcov=vcovCL,cluster=~patient)

df_out$OR_multi <- paste0(round(exp(lreg_coeffs[-1,1]),2)," (",
                          round(exp(lreg_CIs[-1,1]),2),"-",
                          round(exp(lreg_CIs[-1,2]),2),")")
df_out <- rbind(df_out,c("courier:mhd",NA,NA))

# adjusted regression, with courier-MHD interaction
reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+"),"+ courier:mhd_ind"))
lreg_int <- glm(reg_formula, data=DTrna,family="binomial")
lreg_int_coeffs <- coeftest(lreg_int,vcov=vcovCL,cluster=~patient)
lreg_int_CIs <- coefci(lreg_int,vcov=vcovCL,cluster=~patient)

df_out <- cbind(df_out,paste0(round(exp(lreg_int_coeffs[-1,1]),2)," (",
                       round(exp(lreg_int_CIs[-1,1]),2),"-",
                       round(exp(lreg_int_CIs[-1,2]),2),")"))
colnames(df_out)[4] <- "OR_with_interaction"

write_xlsx(df_out,path=file.path(filepath_write,paste0("ORs_",VLS_threshold,".xlsx")))

# additional analyses by calendar period

reg_formula <- as.formula(paste0("vls_ind ~",paste0(setdiff(rf_vect,"calyear_current_cat"),collapse="+")))

lreg <- glm(reg_formula, data=DTrna,family="binomial",subset=(calyear_current_cat=="[2011,2014)"))
lreg_coeffs <- coeftest(lreg,vcov=vcovCL,cluster=~patient)
lreg_CIs <- coefci(lreg,vcov=vcovCL,cluster=~patient)
print(paste0("2011-2013: ",
             round(exp(lreg_coeffs[2,1]),2)," (",
             round(exp(lreg_CIs[2,1]),2),"-",
             round(exp(lreg_CIs[2,2]),2),")"))

lreg <- glm(reg_formula, data=DTrna,family="binomial",subset=(calyear_current_cat=="[2014,2017)"))
lreg_coeffs <- coeftest(lreg,vcov=vcovCL,cluster=~patient)
lreg_CIs <- coefci(lreg,vcov=vcovCL,cluster=~patient)
print(paste0("2014-2016: ",
             round(exp(lreg_coeffs[2,1]),2)," (",
             round(exp(lreg_CIs[2,1]),2),"-",
             round(exp(lreg_CIs[2,2]),2),")"))

lreg <- glm(reg_formula, data=DTrna,family="binomial",subset=(calyear_current_cat=="[2017,2020)"))
lreg_coeffs <- coeftest(lreg,vcov=vcovCL,cluster=~patient)
lreg_CIs <- coefci(lreg,vcov=vcovCL,cluster=~patient)
print(paste0("2017-2019: ",
             round(exp(lreg_coeffs[2,1]),2)," (",
             round(exp(lreg_CIs[2,1]),2),"-",
             round(exp(lreg_CIs[2,2]),2),")"))

lreg <- glm(reg_formula, data=DTrna,family="binomial",subset=(calyear_current_cat=="[2020,Inf)"))
lreg_coeffs <- coeftest(lreg,vcov=vcovCL,cluster=~patient)
lreg_CIs <- coefci(lreg,vcov=vcovCL,cluster=~patient)
print(paste0("2020-2022: ",
             round(exp(lreg_coeffs[2,1]),2)," (",
             round(exp(lreg_CIs[2,1]),2),"-",
             round(exp(lreg_CIs[2,2]),2),")"))

toc()
