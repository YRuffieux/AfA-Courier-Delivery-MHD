# logistic regression on viral load suppression (no/yes)
# using generalized estimating equations to produce odds ratios and 95% CIs: https://www.jstatsoft.org/article/view/v015i02
# exposure=courier pharmacy (none/A/B/other)

library(data.table)
library(geepack)
library(tictoc)
library(writexl)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

rf_vect <- c("courier_cat","mhd_ind","sex","age_current_cat","calyear_current_cat","art_type_cf")
correlation_structure <- "independence"

tic("Overall")

load(file=file.path(filepath_read,"AfA_VL_courier_MHD.RData"))

VLS_threshold <- 400

# formatting
DTrna[,`:=`(vls_ind=as.numeric(rna_v<VLS_threshold),
            age_current_cat=cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
            calyear_current_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
            art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")),
            patient=factor(patient))]
DTrna[,age_current_cat:=relevel(age_current_cat,"[40,50)")]

stopifnot(DTrna[,all(!is.na(art_type_cf))])

df_out <- data.frame(NULL)

tic("Unadjusted regressions")
for(v in rf_vect)
{
  reg_formula <- as.formula(paste0("vls_ind ~",v))
  lreg <- geeglm(reg_formula,data=DTrna,family="binomial",corstr=correlation_structure,id=patient)
  cc <- data.table(coef(summary(lreg)))
  cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
  out <- cc[rf!="(Intercept)",.(rf,OR_uni=paste0(round(estimate,digits=2)," (",round(lower,digits=2),"-",round(upper,digits=2),")"))]
  df_out <- rbind(df_out,out)
  rm(lreg,reg_formula,cc,out)
}
toc()

tic("Adjusted regression")
reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+")))
lreg <- geeglm(reg_formula,data=DTrna,family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
out <- cc[rf!="(Intercept)",.(OR_mult=paste0(round(estimate,digits=2)," (",round(lower,digits=2),"-",round(upper,digits=2),")"))]
toc()

df_out <- cbind(df_out,out)
rm(out,cc,lreg,reg_formula)

write_xlsx(df_out,path=file.path(filepath_write,paste0("ORs_pharmacy_",VLS_threshold,".xlsx")))

toc()
