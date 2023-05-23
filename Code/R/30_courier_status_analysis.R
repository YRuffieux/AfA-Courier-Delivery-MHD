# logistic regression on viral load suppression (no/yes)
# using generalized estimating equations to produce odds ratios and 95% CIs: https://www.jstatsoft.org/article/view/v015i02
# overall, and by calendar period
# exposure=courier status no/yes
# ~10 minute runtime

library(data.table)
library(geepack)
library(tictoc)
library(writexl)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

rf_vect <- c("courier","mhd_ind","sex","age_current_cat","calyear_current_cat","art_type_cf")
correlation_structure <- "exchangeable"
courier_lag <- 0                 # in months: 0, 6, or 12
VLS_threshold <- 400
include_untested <- FALSE         # whether to include untested follow-up - will be set to unsuppressed VL every six months

tic("Overall")

format_CI <- function(estimate,lower,upper,digits=2)
  paste0(format(round(estimate,digits),nsmall=digits)," (",
         format(round(lower,digits),nsmall=digits),"-",
         format(round(upper,digits),nsmall=digits),")")

if(courier_lag==0)
{
    load(file=file.path(filepath_read,"AfA_VL.RData"))
} else
{
  load(file=file.path(filepath_read,paste0("AfA_VL_lag",courier_lag,".RData")))
}

if(include_untested)
{
  DTrna[is.na(rna_v),rna_v:=VLS_threshold]
} else
{
  DTrna <- DTrna[!is.na(rna_v)]
}

# formatting
DTrna[,`:=`(vls_ind=as.numeric(rna_v<VLS_threshold),
            age_current_cat=cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
            calyear_current_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
            art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")),
            patient=factor(patient))]
DTrna[,`:=`(age_current_cat=relevel(age_current_cat,"[40,50)"),
            art_type_cf_2=as.character(art_type_cf))]
DTrna[art_type_cf%in%c("II+2NRTI","PI+2NRTI"),art_type_cf_2:="II/PI+2NRTI"]
DTrna[,art_type_cf_2:=factor(art_type_cf_2,levels=c("NNRTI+2NRTI","II/PI+2NRTI"))]    # need to group II and PI for pre-2018 period

stopifnot(DTrna[,all(!is.na(art_type_cf))])

df_out <- data.frame(NULL)

tic("Unadjusted regressions (overall)")
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

tic("Adjusted regression (overall)")
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
if(include_untested)
  savename <- paste0(savename,"_with_untested")


write_xlsx(df_out,path=file.path(filepath_write,paste0(savename,"_vls",VLS_threshold,"_overall.xlsx")))

rm(df_out)

df_list <- list()

for(cp in DTrna[,levels(calyear_current_cat)])
{
  tic(paste0("By calendar period: ",cp))
  
  df_out <- data.frame(NULL)
  rf_vect_cp <- setdiff(rf_vect,"calyear_current_cat")
  if(cp%in%c("[2011,2014)","[2014,2017)"))   # need to aggregate PI and II for those calendar periods
   rf_vect_cp[rf_vect_cp=="art_type_cf"] <- "art_type_cf_2"
  
  for(v in rf_vect_cp)   # univariate analyses
  {
    reg_formula <- as.formula(paste0("vls_ind ~",v))
    lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat==cp],family="binomial",corstr=correlation_structure,id=patient)
    cc <- data.table(coef(summary(lreg)))
    cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
    out <- cc[rf!="(Intercept)",.(rf,OR_uni=format_CI(estimate,lower,upper))]
    df_out <- rbind(df_out,out)
    rm(lreg,reg_formula,cc,out)
  }
  
  reg_formula_cp <- as.formula(paste0("vls_ind ~",paste0(rf_vect_cp,collapse="+")))
  lreg <- geeglm(reg_formula_cp,data=DTrna[calyear_current_cat==cp],family="binomial",corstr=correlation_structure,id=patient)
  cc <- data.table(coef(summary(lreg)))
  cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
  out <- cc[rf!="(Intercept)",.(OR_mult=format_CI(estimate,lower,upper))]
  
  df_out <- cbind(df_out,out)
  
  df_list[[gsub("[[)]","",cp)]] <- df_out
  rm(out,cc,lreg,reg_formula_cp,df_out)
  
  toc()
}

write_xlsx(df_list,path=file.path(filepath_write,paste0(savename,"_vls",VLS_threshold,"_by_period.xlsx")))

gc(verbose=FALSE)

toc()
