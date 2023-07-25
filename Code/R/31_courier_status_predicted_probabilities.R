# adjusted/unadjusted predicted probabilities of VL suppression by ART delivery method, by calendar period and overall
# using age on continuous scale, ART regimen and scheme code as a 0/1 continuous variable (NNRTI/non-NNRTI, PLM/non-PLM)
# 2-3 minutes

library(data.table)
library(geepack)
library(tictoc)
library(writexl)
library(ggplot2)
library(ggeffects)
library(miceadds)
library(stringr)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Plots"

which_scheme <- "All"         # current options: All, BON, PLM (analysis left-truncated at start of 2016), Other

rf_vect <- c("courier","mhd_ind","sex","age_current","art_type_num")
correlation_structure <- "exchangeable"
courier_lag <- 0                 # in months: 0, 6, or 12
VLS_threshold <- 400
include_untested <- FALSE     # whether to include untested follow-up - will be set to unsuppressed VL every six months
create_plots <- TRUE

tic("Overall")

if(courier_lag==0)
{
  load(file=file.path(filepath_read,"AfA_VL.RData"))
} else
{
  load(file=file.path(filepath_read,paste0("AfA_VL_lag",courier_lag,".RData")))
}

DTrna[!scheme_code%in%c("BON","PLM"),scheme_code:="Other"]
DTrna[,`:=`(scheme_code=factor(scheme_code,levels=c("BON","PLM","Other")),scheme_code_base=NULL)]
setnames(DTrna,"art_type_cf","art_type")

if(which_scheme=="PLM")
  DTrna <- DTrna[scheme_code=="PLM" & year(rna_d)>=2016]
if(which_scheme=="BON")
  DTrna <- DTrna[scheme_code=="BON"]
if(which_scheme=="Other")
  DTrna <- DTrna[scheme_code=="Other"]
if(which_scheme=="All")          # adjusting for scheme
  rf_vect <- c(rf_vect,"scheme_code_num")

savename <- "probs_courier"
if(courier_lag!=0)
  savename <- paste0(savename,"_lag",courier_lag)
if(include_untested)
  savename <- paste0(savename,"_with_untested")
if(which_scheme!="All")
  savename <- paste0(savename,"_",which_scheme)
savename <- paste0(savename,"_vls",VLS_threshold)
print(savename)

if(include_untested)
{
  DTrna[is.na(rna_v),rna_v:=VLS_threshold]
} else
{
  DTrna <- DTrna[!is.na(rna_v)]
}

# formatting
DTrna[,`:=`(vls_ind=as.numeric(rna_v<VLS_threshold),
            patient=factor(patient),
            calyear_num=as.numeric(rna_d-as.Date("2011-01-01"))/365.25,
            art_type_num=as.numeric(art_type!="NNRTI+2NRTI"),
            scheme_code_num=as.numeric(scheme_code!="PLM"))]
DTrna[,art_type:=NULL]

# different calyear cutoffs when restricting to PLM
if(which_scheme=="PLM")
{
  DTrna[,calyear_current_cat:=cut(year(rna_d),breaks=c(2016,2018,2020,Inf),right=FALSE)]
} else
{
  DTrna[,calyear_current_cat:=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE)]
}

DTrna[,sex:=as.numeric(sex)-1]

df_out <- data.table(col=c("retail_unadj","courier_unadj","","courier_adj","courier_unadj"))

# unadjusted probabilities, by calendar period and overall
probs_gee_crude_df <- data.table(NULL)
reg_formula <- as.formula("vls_ind ~ courier")
for(period in DTrna[,levels(calyear_current_cat)])
{
  lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat==period],family="binomial",corstr="exchangeable",id=patient)
  gres <- ggpredict(lreg,"courier")
  probs_gee_crude_df <- rbind(probs_gee_crude_df,
                              data.table(calyear_cat=rep(period,2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
  rm(lreg,gres)
}

lreg <- geeglm(reg_formula,DTrna,family="binomial",corstr="exchangeable",id=patient)
gres <- ggpredict(lreg,"courier")
probs_gee_crude_df <- rbind(probs_gee_crude_df,
                            data.table(calyear_cat=rep("Overall",2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
rm(reg_formula,lreg,gres)

# adjusted probabilities, by calendar period and overall
probs_gee_adj_df <- data.table(NULL)
reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+")))
print(reg_formula)
for(period in DTrna[,levels(calyear_current_cat)])
{
  lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat==period],family="binomial",corstr="exchangeable",id=patient)
  gres <- ggpredict(lreg,"courier")
  probs_gee_adj_df <- rbind(probs_gee_adj_df,
                            data.table(calyear_cat=rep(period,2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
  rm(lreg,gres)
}

rm(reg_formula)
reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+"),"+ calyear_num"))
lreg <- geeglm(reg_formula,data=DTrna,family="binomial",corstr="exchangeable",id=patient)
gres <- ggpredict(lreg,"courier")
probs_gee_adj_df <- rbind(probs_gee_adj_df,
                          data.table(calyear_cat=rep("Overall",2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
rm(reg_formula,lreg,gres)

probs_gee_crude_df[,courier:=factor(courier)]
probs_gee_adj_df[,courier:=factor(courier)]

if(which_scheme!="PLM")
{
  x_labels <- c("2011-2013","2014-2016","2017-2019","2020-2022")
} else
{
  x_labels <- c("2016-2017","2018-2019","2020-2022")
}

if(which_scheme!="BON")
{
  y_lim <- c(70,100)
} else
{
  y_lim <- c(30,100)
}

if(create_plots)
{
  pp_gee_crude <- ggplot(probs_gee_crude_df[calyear_cat!="Overall"],aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
    geom_point(size=2,position=position_dodge(0.5)) +
    geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
    labs(x="Time period",y="Percentage virally suppressed") +
    theme_bw() +
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
    scale_shape_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_color_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_x_discrete(labels=x_labels) +
    ylim(y_lim)
  ggsave(pp_gee_crude,filename=file.path(filepath_plot,paste0("crude_",savename,"_by_period.png")),height=4,width=6,dpi=600)
  
  pp_gee_adj <- ggplot(probs_gee_adj_df[calyear_cat!="Overall"],aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
    geom_point(size=2,position=position_dodge(0.5)) +
    geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
    labs(x="Time period",y="Percentage virally suppressed") +
    theme_bw() +
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
    scale_shape_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_color_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_x_discrete(labels=x_labels) +
    ylim(y_lim)
  ggsave(pp_gee_adj,filename=file.path(filepath_plot,paste0("adj_",savename,"_by_period.png")),height=4,width=6,dpi=600)
}

# Excel output
format_CI <- function(estimate,lower,upper,digits=2)
  paste0(format(round(estimate,digits),nsmall=digits)," (",
         format(round(lower,digits),nsmall=digits),"-",
         format(round(upper,digits),nsmall=digits),")")

probs_gee_adj_out <- probs_gee_adj_df[,.(calyear_cat=factor(calyear_cat,levels=c(DTrna[,levels(calyear_current_cat)],"Overall")),
                                         courier,out=format_CI(probability,lower,upper,digits=1))]
probs_gee_adj_out <- dcast(probs_gee_adj_out,courier~calyear_cat,value.var="out")
probs_gee_crude_out <- probs_gee_crude_df[,.(calyear_cat=factor(calyear_cat,levels=c(DTrna[,levels(calyear_current_cat)],"Overall")),
                                             courier,out=format_CI(probability,lower,upper,digits=1))]
probs_gee_crude_out <- dcast(probs_gee_crude_out,courier~calyear_cat,value.var="out")

probs_gee_out <- rbind(probs_gee_crude_out,probs_gee_adj_out)

write_xlsx(probs_gee_out,path=file.path(filepath_write,paste0(savename,"_by_period.xlsx")))

toc()
