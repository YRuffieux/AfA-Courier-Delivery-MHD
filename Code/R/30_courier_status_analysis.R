# logistic regression on viral load suppression (no/yes)
# using generalized estimating equations to produce odds ratios and 95% CIs: https://www.jstatsoft.org/article/view/v015i02
# exposure=courier status no/yes
# 4-5 minutes total runtime with exchangeable correlation structure

library(data.table)
library(geepack)
library(tictoc)
library(writexl)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

rf_vect <- c("courier","mhd_ind","sex","age_current_cat","calyear_current_cat","art_type_cf")
correlation_structure <- "exchangeable"
courier_lag <- 0                 # in months: 0, 6, or 12
first_line_art <- FALSE          # whether to censor tests occuring during 2nd+ line ART
VLS_threshold <- 400

tic("Overall")

format_CI <- function(estimate,lower,upper,digits=2)
  paste0(format(round(estimate,digits),nsmall=digits)," (",
         format(round(lower,digits),nsmall=digits),"-",
         format(round(upper,digits),nsmall=digits),")")

if(courier_lag==0)
{
  load(file=file.path(filepath_read,"AfA_VL_courier_MHD.RData"))
} else
{
  load(file=file.path(filepath_read,paste0("AfA_VL_courier_MHD_lag",courier_lag,".RData")))
}

if(first_line_art)
  DTrna <- DTrna[art_first_line==1]

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
  out <- cc[rf!="(Intercept)",.(rf,OR_uni=format_CI(estimate,lower,upper))]
  df_out <- rbind(df_out,out)
  rm(lreg,reg_formula,cc,out)
}
toc()

tic("Adjusted regression")
reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+")))
lreg <- geeglm(reg_formula,data=DTrna,family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
out <- cc[rf!="(Intercept)",.(OR_mult=format_CI(estimate,lower,upper))]
toc()

df_out <- cbind(df_out,out)
rm(out,cc,lreg,reg_formula)

savename <- "ORs_courier"
if(courier_lag!=0)
  savename <- paste0(savename,"_lag",courier_lag)
if(first_line_art==TRUE)
  savename <- paste0(savename,"_firstline")
savename <- paste0(savename,"_vls",VLS_threshold,".xlsx")

print(paste0("Save name: ",savename))
write_xlsx(df_out,path=file.path(filepath_write,savename))

tic("By calendar period")
reg_formula <- as.formula(paste0("vls_ind ~",paste0(setdiff(rf_vect,c("calyear_current_cat","art_type_cf")),collapse="+")))

lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat=="[2011,2014)"],family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
print(cc[rf=="courier",paste0("2011-2013: ",round(estimate,digits=2)," (",round(lower,digits=2),"-",round(upper,digits=2),")")])

lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat=="[2014,2017)"],family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
print(cc[rf=="courier",paste0("2014-2016: ",round(estimate,digits=2)," (",round(lower,digits=2),"-",round(upper,digits=2),")")])

lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat=="[2017,2020)"],family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
print(cc[rf=="courier",paste0("2017-2019: ",round(estimate,digits=2)," (",round(lower,digits=2),"-",round(upper,digits=2),")")])

lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat=="[2020,Inf)"],family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
print(cc[rf=="courier",paste0("2020-2022: ",round(estimate,digits=2)," (",round(lower,digits=2),"-",round(upper,digits=2),")")])

lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat!="[2011,2014)"],family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
print(cc[rf=="courier",paste0("2014-2022: ",round(estimate,digits=2)," (",round(lower,digits=2),"-",round(upper,digits=2),")")])

toc()




rm(lreg,cc)
gc(verbose=FALSE)

toc()
