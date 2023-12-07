# unadjusted and adjusted hazard ratios for predictors associated with switching ART delivery type
# from two-state model with states "Courier" and "Retail"
# only done for main analysis (no Bonitas), restriction already done upstream (12_data_mgmt_mstate.R and 13_data_mgmt_mstate_split.R)

library(tictoc)
library(data.table)
library(survival)
library(mstate)
library(writexl)
library(stringr)

filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"


tic("Overall")

load(file=file.path(filepath_processed,"AfA_mstate_split.RData"))

covariates <- c("vls_ind","mhd_ind","sex","age_cat","calyear_cat","art_type")

DTms_split[!scheme_code_base%in%c("BON","PLM"),scheme_code_base:="Other"]
DTms_split[,scheme_code_base:=factor(scheme_code_base,levels=c("PLM","Other"))]

DTms_split[,`:=`(vls_ind=factor(vls_ind+1),mhd_ind=factor(mhd_ind+1))]

tic("Expanding risk factors by transition")

DTms_split[,art_type:=factor(as.numeric(art_type))]

DTms_split[,trans:=1]
DTms_split[from=="Courier",trans:=2]
tmat <- matrix(c(NA,1,2,NA),nrow=2,byrow=TRUE)
dimnames(tmat) <- list(from=c("Retail","Courier"),to=c("Retail","Courier"))
attr(DTms_split,"trans") <- tmat
class(DTms_split) <- c("msdata","data.frame")

DTms_split <- data.table(expand.covs(DTms_split,covs=covariates))

toc()

covnames <- names(DTms_split)
covariates_exp <- covnames[substr(covnames,nchar(covnames)-1,nchar(covnames))%in%c(".1",".2")]
main_formula <- paste0("Surv(start,end,status)~",paste0(covariates_exp,collapse="+"),"+strata(trans) + strata(scheme_code_base)")
main_formula <- as.formula(main_formula)
print(main_formula)

df_unadjusted <- data.frame(NULL)

tic("Unadjusted regressions")
for(v in covariates)
{
  V <- covnames[grepl(v,covnames) & substr(covnames,nchar(covnames)-1,nchar(covnames))%in%c(".1",".2")]
  uni_formula <- paste0("Surv(start,end,status)~",paste0(V,collapse="+"),"+strata(trans) + strata(scheme_code_base)")
  uni_formula <- as.formula(uni_formula)
  res_uni <- coxph(uni_formula,data=DTms_split)
  X <- summary(res_uni)$coefficients
  Y <- summary(res_uni)$conf.int
  Z <- cbind(X[,2],Y[,c(3,4)])
  df_unadjusted <- rbind(df_unadjusted,Z)
  rm(X,Y,Z,res_uni,uni_formula)
}
df_unadjusted <- data.table(cbind(covariates_exp,df_unadjusted))
df_unadjusted[,covariates_exp:=factor(covariates_exp,levels=covariates_exp)]
colnames(df_unadjusted) <- c("rf","uHR","uLCL","uUCL")
toc()

tic("Adjusted regression")
res <- coxph(main_formula,data=DTms_split)
X <- summary(res)$coefficients
Y <- summary(res)$conf.int
df_adjusted <- data.frame(cbind(X[,2],Y[,c(3,4)]))
df_adjusted <- data.table(cbind(covariates_exp,df_adjusted))
df_adjusted[,covariates_exp:=factor(covariates_exp,levels=covariates_exp)]
colnames(df_adjusted) <- c("rf","aHR","aLCL","aUCL")
rm(X,Y)
toc()

# setting up Excel table of results

format_CI <- function(a,b,c,digits=2) 
  paste0(trimws(format(round(a,digits),nsmall=digits))," (",trimws(format(round(b,digits),nsmall=digits)),"-",trimws(format(round(c,digits),nsmall=digits)),")")

stopifnot(all(df_unadjusted[,rf]==df_adjusted[,rf]))

df_all <- merge(df_unadjusted,df_adjusted,by="rf")
df_all <- df_all[,.(rf,uHR=format_CI(uHR,uLCL,uUCL),aHR=format_CI(aHR,aLCL,aUCL))]
df_t1 <- df_all[grepl(".1",rf,fixed=TRUE)]
df_t2 <- df_all[grepl(".2",rf,fixed=TRUE)]

Z <- rbind(data.table(rf=paste0(covariates,"0.x"),uHR="",aHR=""),
          data.table(rf=paste0(covariates,"1.x"),uHR="1",aHR="1"))
df_t1 <- rbind(df_t1,Z)
df_t1[,`:=`(rf_short=factor(str_sub(rf,1,-4),levels=covariates),rf_num=str_sub(rf,-3,-3))]
setorder(df_t1,"rf_short","rf_num")
df_t1[,`:=`(rf_short=NULL,rf_num=NULL)]
df_t2 <- rbind(df_t2,Z)
df_t2[,`:=`(rf_short=factor(str_sub(rf,1,-4),levels=covariates),rf_num=str_sub(rf,-3,-3))]
setorder(df_t2,"rf_short","rf_num")
df_t2[,`:=`(rf_short=NULL,rf_num=NULL)]
rm(Z)

df_out <- list("retail_to_courier"=df_t1,"courier_to_retail"=df_t2)

write_xlsx(df_out,path=file.path(filepath_tables,"mstate_HRs_courier_noBON.xlsx"))

X <- DTms_split[from=="Retail",.(events=sum(status),py=sum(end-start)/365.25),by=c("scheme_code_base","vls_ind")]
setorder(X,"scheme_code_base","vls_ind")
X[,rate:=events/py]

X[scheme_code_base=="PLM",rate[2]/rate[1]]
X[scheme_code_base=="Other",rate[2]/rate[1]]

Y <- DTms_split[from=="Retail",.(events=sum(status),py=sum(end-start)/365.25),by=c("vls_ind")]
Y[,rate:=events/py]
Y[,rate[2]/rate[1]]

toc()

