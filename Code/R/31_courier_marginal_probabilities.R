# adjusted + unadjusted retail/courier marginal probabilities of VL suppression, by ART delivery method, by calendar period and overall
# adjusted analysis: using estimated marginal means (EMMs) conditioned on courier status
# 2-3 minutes
# plots produced only for primary analysis (without Bonitas)

library(data.table)
library(geepack)
library(tictoc)
library(writexl)
library(ggplot2)
library(emmeans)
library(stringr)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Plots"

which_scheme <- "noBON"         # current options: noBON, AllSchemes, BON, PLM (analysis left-truncated at start of 2016), Other

stopifnot(which_scheme%in%c("noBON","AllSchemes","BON","PLM","Other"))

rf_vect <- c("courier","mhd_ind","sex","age_current_cat","art_type","calyear_current_cat","scheme_code")
correlation_structure <- "exchangeable"
courier_lag <- 0                 # in months: 0, 6, or 12
VLS_threshold <- 400

if(!which_scheme%in%c("AllSchemes","noBON"))
{
  filepath_write <- file.path(filepath_write,"Scheme-specific")
}
if(courier_lag>0 || VLS_threshold!=400 || which_scheme=="AllSchemes")
{
  filepath_write <- file.path(filepath_write,"Sensitivity")
}

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
  DTrna[,`:=`(scheme_code=factor(scheme_code,levels=c("PLM","Other")))]
} else if(which_scheme=="AllSchemes")
{
  DTrna[,`:=`(scheme_code=factor(scheme_code,levels=c("BON","PLM","Other")))]
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

savename <- "probs_courier"
if(courier_lag!=0)
  savename <- paste0(savename,"_lag",courier_lag)
savename <- paste0(savename,"_",which_scheme)
savename <- paste0(savename,"_vls",VLS_threshold)
print(savename)

# formatting
DTrna[,`:=`(sex=sex-1,
            vls_ind=as.numeric(rna_v<VLS_threshold),
            patient=factor(patient),
            age_current_cat=factor(cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE,labels=FALSE)),
            art_type=factor(art_type,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")))]
DTrna[,`:=`(age_current_cat=relevel(age_current_cat,3),
            art_type_agg=art_type)]
levels(DTrna$art_type_agg) <- c("NNRTI+2NRTI","PI+2NRTI","PI+2NRTI")

# calendar year: different cutoffs when restricting to PLM
if(which_scheme=="PLM")
{
  DTrna[,calyear_current_cat:=cut(year(rna_d),breaks=c(2016,2018,2020,Inf),right=FALSE)]
} else
{
  DTrna[,calyear_current_cat:=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE)]
}

# unadjusted probabilities, by calendar period and overall
probs_gee_crude_df <- data.table(NULL)
reg_formula <- as.formula("vls_ind ~ courier")
for(period in DTrna[,levels(calyear_current_cat)])
{
  lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat==period],family="binomial",corstr="exchangeable",id=patient)
  emme <- summary(emmeans(lreg,"courier",type="response"))
  probs_gee_crude_df <- rbind(probs_gee_crude_df,
                              data.table(calyear_cat=rep(period,2),courier=c(0,1),probability=100*emme$prob,lower=100*emme$asymp.LCL,upper=100*emme$asymp.UCL))
  rm(lreg,emme)
}

lreg <- geeglm(reg_formula,DTrna,family="binomial",corstr="exchangeable",id=patient)
emme <- summary(emmeans(lreg,"courier",type="response"))
probs_gee_crude_df <- rbind(probs_gee_crude_df,
                            data.table(calyear_cat=rep("Overall",2),courier=c(0,1),probability=100*emme$prob,lower=100*emme$asymp.LCL,upper=100*emme$asymp.UCL))
rm(reg_formula,lreg,emme)

# adjusted probabilities, by calendar period and overall
probs_gee_adj_df <- data.table(NULL)
for(period in DTrna[,levels(calyear_current_cat)])
{
  rf_vect_cp <- setdiff(rf_vect,"calyear_current_cat")
  if(period%in%c("[2011,2014)","[2014,2017)"))
    rf_vect_cp[rf_vect_cp=="art_type"] <- "art_type_agg"
  reg_formula_cp <- as.formula(paste0("vls_ind ~",paste0(rf_vect_cp,collapse="+")))
  lreg <- geeglm(reg_formula_cp,data=DTrna[calyear_current_cat==period],family="binomial",corstr="exchangeable",id=patient)
  emme <- summary(emmeans(lreg,"courier",type="response"))
  probs_gee_adj_df <- rbind(probs_gee_adj_df,
                            data.table(calyear_cat=rep(period,2),courier=c(0,1),probability=100*emme$prob,lower=100*emme$asymp.LCL,upper=100*emme$asymp.UCL))
  rm(lreg,emme,rf_vect_cp,reg_formula_cp)
}

reg_formula <- as.formula(paste0("vls_ind ~",paste0(rf_vect,collapse="+")))
lreg <- geeglm(reg_formula,data=DTrna,family="binomial",corstr="exchangeable",id=patient)
emme <- summary(emmeans(lreg,"courier",type="response"))
probs_gee_adj_df <- rbind(probs_gee_adj_df,
                          data.table(calyear_cat=rep("Overall",2),courier=c(0,1),probability=100*emme$prob,lower=100*emme$asymp.LCL,upper=100*emme$asymp.UCL))
rm(reg_formula,lreg,emme)

probs_gee_crude_df[,courier:=factor(courier)]
probs_gee_adj_df[,courier:=factor(courier)]

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

# plotting the percentages: only for primary analysis
if(which_scheme=="noBON")
{
  pp_gee_crude <- ggplot(probs_gee_crude_df[calyear_cat!="Overall"],aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
    geom_point(size=2,position=position_dodge(0.5)) +
    geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
    labs(x="Time period",y="Percentage virally suppressed") +
    theme_bw() +
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
    scale_shape_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_color_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_x_discrete(labels=c("2011-2013","2014-2016","2017-2019","2020-2022")) +
    coord_cartesian(ylim=c(75,100))
  ggsave(pp_gee_crude,filename=file.path(filepath_plot,paste0("unadj_",savename,"_by_period.png")),height=4,width=6,dpi=600)
  
  pp_gee_adj <- ggplot(probs_gee_adj_df[calyear_cat!="Overall"],aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
    geom_point(size=2,position=position_dodge(0.5)) +
    geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
    labs(x="Time period",y="Percentage virally suppressed") +
    theme_bw() +
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
    scale_shape_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_color_discrete(labels=c("Retail","Courier"),name=NULL) +
    scale_x_discrete(labels=c("2011-2013","2014-2016","2017-2019","2020-2022")) +
    coord_cartesian(ylim=c(75,100))
  ggsave(pp_gee_adj,filename=file.path(filepath_plot,paste0("adj_",savename,"_by_period.png")),height=4,width=6,dpi=600)
}

toc()
