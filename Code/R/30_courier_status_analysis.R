# logistic regression on viral load suppression (no/yes)
# using generalized estimating equations to produce odds ratios and 95% CIs: https://www.jstatsoft.org/article/view/v015i02
# overall, and by calendar period
# exposure=courier status no/yes
# ~5-10 minute runtime
# primary analysis: without Bonitas

library(data.table)
library(geepack)
library(tictoc)
library(writexl)
library(stringr)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

which_scheme <- "noBON"         # current options: noBON (primary analysis), AllSchemes, BON, PLM (analysis left-truncated at start of 2016), Other

stopifnot(which_scheme%in%c("noBON","AllSchemes","BON","PLM","Other"))

rf_vect <- c("courier","mhd_ind","sex","age_current_cat","art_type","scheme_code","calyear_current_cat")
correlation_structure <- "exchangeable"
courier_lag <- 0              # in months: 0, 6, or 12
VLS_threshold <- 400

if(!which_scheme%in%c("AllSchemes","noBON"))
  filepath_write <- file.path(filepath_write,"Scheme-specific")
if(courier_lag>0 || VLS_threshold!=400 || which_scheme=="AllSchemes")
  filepath_write <- file.path(filepath_write,"Sensitivity")

if(which_scheme!="noBON" && (courier_lag!=0 || VLS_threshold!=400))
  stop("No need to do sensitivity analysis on scheme-specific data or data with Bonitas")

tic("Overall")

if(courier_lag==0)
{
    load(file=file.path(filepath_read,"AfA_VL.RData"))
} else
{
  load(file=file.path(filepath_read,paste0("AfA_VL_lag",courier_lag,".RData")))
}

DTrna <- DTrna[!is.na(rna_v)] # removing 'fake' tests
DTrna[!scheme_code%in%c("BON","PLM"),scheme_code:="Other"]
if(which_scheme=="noBON")
{
  DTrna <- DTrna[scheme_code_base!="BON"]
  DTrna <- DTrna[scheme_code!="BON"]
  DTrna[,`:=`(scheme_code=factor(scheme_code,levels=c("PLM","Other"),labels=c(1,2)))]
} else if(which_scheme=="AllSchemes")
{
  DTrna[,`:=`(scheme_code=factor(scheme_code,levels=c("BON","PLM","Other"),labels=1:3))]
}
setnames(DTrna,"art_type_cf","art_type")
DTrna[,scheme_code_base:=NULL]

if(which_scheme=="PLM")
  DTrna <- DTrna[scheme_code=="PLM" & year(rna_d)>=2016]
if(which_scheme=="Other")
  DTrna <- DTrna[scheme_code=="Other"]
if(which_scheme=="BON")
  DTrna <- DTrna[scheme_code=="BON"]
if(!which_scheme%in%c("AllSchemes","noBON"))
  rf_vect <- setdiff(rf_vect,"scheme_code")

savename <- "ORs_courier"
if(courier_lag!=0)
  savename <- paste0(savename,"_lag",courier_lag)
savename <- paste0(savename,"_",which_scheme)
savename <- paste0(savename,"_vls",VLS_threshold)
print(savename)

# formatting
DTrna[,`:=`(vls_ind=as.numeric(rna_v<VLS_threshold),
            courier=factor(courier+1),
            mhd_ind=factor(mhd_ind+1),
            sex=factor(sex),
            age_current_cat=factor(cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE,labels=FALSE)),
            art_type=factor(art_type,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI"),labels=1:3),
            patient=factor(patient))]
DTrna[,`:=`(age_current_cat=relevel(age_current_cat,3),
            art_type_agg=as.character(art_type))]
DTrna[art_type%in%c(2,3),art_type_agg:=2]
DTrna[,art_type_agg:=factor(art_type_agg,levels=c(1,2))]    # need to group II and PI for pre-2018 period

stopifnot(DTrna[,all(!is.na(art_type))])

# different calyear categorization when restricting to PLM
if(which_scheme=="PLM")
{
  DTrna[,calyear_current_cat:=factor(cut(year(rna_d),breaks=c(2016,2018,2020,Inf),right=FALSE,labels=FALSE))]
} else
{
  DTrna[,calyear_current_cat:=factor(cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE,labels=FALSE))]
}

df_out <- data.table(NULL)
format_CI <- function(a,b,c,digits=2) 
  paste0(trimws(format(round(a,digits),nsmall=digits))," (",trimws(format(round(b,digits),nsmall=digits)),"-",trimws(format(round(c,digits),nsmall=digits)),")")

tic("Unadjusted regressions (overall)")
for(v in rf_vect)
{
  reg_formula <- as.formula(paste0("vls_ind ~",v))
  lreg <- geeglm(reg_formula,data=DTrna,family="binomial",corstr=correlation_structure,id=patient)
  cc <- data.table(coef(summary(lreg)))
  cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
  out <- cc[rf!="(Intercept)",.(rf,uOR=format_CI(estimate,lower,upper))]
  df_out <- rbind(df_out,out)
  rm(lreg,reg_formula,cc,out)
}
toc()

tic("Adjusted regression (overall)")
reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+")))
lreg <- geeglm(reg_formula,data=DTrna,family="binomial",corstr=correlation_structure,id=patient)
cc <- data.table(coef(summary(lreg)))
cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
out <- cc[rf!="(Intercept)",.(aOR=format_CI(estimate,lower,upper))]
toc()

df_out <- cbind(df_out,out)
rm(out,cc,lreg,reg_formula)

df_small_overall <- c(1,df_out[rf=="courier2",uOR],"",1,df_out[rf=="courier2",aOR])

# setting up Excel tables of results
Z <- rbind(data.table(rf=paste0(rf_vect,"0"),uOR="",aOR=""),
           data.table(rf=paste0(rf_vect,"1"),uOR="1",aOR="1"))
df_out <- rbind(df_out,Z)
df_out[,`:=`(rf_char=factor(str_sub(rf,start=1,end=-2),levels=rf_vect),rf_num=str_sub(rf,start=-1))]
setorder(df_out,"rf_char","rf_num")
df_out[,`:=`(rf_char=NULL,rf_num=NULL)]
rm(Z)

write_xlsx(df_out,path=file.path(filepath_write,paste0(savename,".xlsx")))

rm(df_out)

# by calendar period

df_list <- list()
df_small <- data.table(NULL)

for(cp in DTrna[,levels(calyear_current_cat)])
{
  tic(paste0("By calendar period - ",cp))
  
  df_out <- data.table(NULL)
  rf_vect_cp <- setdiff(rf_vect,"calyear_current_cat")
  if(which_scheme!="PLM" & cp%in%c(1,2))   # need to aggregate PI and II for 2011-2016
   rf_vect_cp[rf_vect_cp=="art_type"] <- "art_type_agg"
  
  for(v in rf_vect_cp)   # univariate analyses
  {
    reg_formula <- as.formula(paste0("vls_ind ~",v))
    lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat==cp],family="binomial",corstr=correlation_structure,id=patient)
    cc <- data.table(coef(summary(lreg)))
    cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
    out <- cc[rf!="(Intercept)",.(rf,uOR=format_CI(estimate,lower,upper))]
    df_out <- rbind(df_out,out)
    rm(lreg,reg_formula,cc,out)
  }
  
  reg_formula_cp <- as.formula(paste0("vls_ind ~",paste0(rf_vect_cp,collapse="+")))
  lreg <- geeglm(reg_formula_cp,data=DTrna[calyear_current_cat==cp],family="binomial",corstr=correlation_structure,id=patient)
  cc <- data.table(coef(summary(lreg)))
  cc <- cc[,.(rf=names(coef(lreg)),estimate=exp(Estimate),lower=exp(Estimate-qnorm(0.975)*Std.err),upper=exp(Estimate+qnorm(0.975)*Std.err))]
  out <- cc[rf!="(Intercept)",.(rf,aOR=format_CI(estimate,lower,upper))]
  
  stopifnot(all(df_out[,rf]==out[,rf]))
  
  df_out <- merge(df_out,out,by="rf")
  
  df_small[,eval(cp):=c(1,df_out[rf=="courier2",uOR],"",1,df_out[rf=="courier2",aOR])]
  
  # setting up Excel tables of results
  Z <- rbind(data.table(rf=paste0(rf_vect_cp,"0"),uOR="",aOR=""),
             data.table(rf=paste0(rf_vect_cp,"1"),uOR="1",aOR="1"))
  df_out <- rbind(df_out,Z)
  df_out[,`:=`(rf_char=factor(str_sub(rf,start=1,end=-2),levels=rf_vect_cp),rf_num=str_sub(rf,start=-1))]
  setorder(df_out,"rf_char","rf_num")
  df_out[,`:=`(rf_char=NULL,rf_num=NULL)]
  rm(Z)
  
  df_list[[cp]] <- df_out
  rm(out,cc,lreg,reg_formula_cp,df_out,rf_vect_cp)
  
  toc()
}

df_small[,Overall:=df_small_overall]
rm(df_small_overall)

write_xlsx(df_list,path=file.path(filepath_write,paste0(savename,"_by_period.xlsx")))
write_xlsx(df_small,path=file.path(filepath_write,paste0(savename,"_by_period_small.xlsx")))

gc(verbose=FALSE)

toc()
