# adjusted/unadjusted predicted probabilities from VL suppression logistic model, by calendar period
# 1-2 minutes

library(data.table)
library(geepack)
library(tictoc)
library(writexl)
library(ggplot2)
library(ggeffects)
library(miceadds)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Plots"

rf_vect <- c("courier","mhd_ind","sex","age_current_cat","calyear_current_cat","art_type_cf")
correlation_structure <- "exchangeable"
courier_lag <- 0                 # in months: 0, 6, or 12
VLS_threshold <- 400

tic("Overall")

if(courier_lag==0)
{
  load(file=file.path(filepath_read,"AfA_VL.RData"))
} else
{
  load(file=file.path(filepath_read,paste0("AfA_VL_lag",courier_lag,".RData")))
}


# formatting
DTrna[,`:=`(vls_ind=as.numeric(rna_v<VLS_threshold),
            age_current_cat=cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
            calyear_current_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
            art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")),
            patient=factor(patient))]
DTrna[,age_current_cat:=relevel(age_current_cat,"[40,50)")]

stopifnot(DTrna[,all(!is.na(art_type_cf))])


# probs <- DTrna[,.(prob=sum(vls_ind==1)/.N),by=.(courier,calyear_current_cat)]
# probs[,courier:=factor(courier)]
# setorder(probs,"calyear_current_cat","courier")

DTrna[,sex:=as.numeric(sex)-1]

# crude probabilities
probs_crude_df <- data.frame(NULL)
reg_formula <- as.formula("vls_ind ~ courier")
for(period in DTrna[,levels(calyear_current_cat)])
{
  lreg <- glm(reg_formula,data=DTrna[calyear_current_cat==period],family="binomial")
  gres <- ggpredict(lreg,"courier")
  probs_crude_df <- rbind(probs_crude_df,data.frame(calyear_cat=rep(period,2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
}
rm(reg_formula)

# adjusted probabilities, no clustering
probs_adj_df <- data.frame(NULL)
reg_formula <- as.formula("vls_ind ~ courier + mhd_ind + sex + age_current_cat")
for(period in DTrna[,levels(calyear_current_cat)])
{
  lreg <- glm(reg_formula,data=DTrna[calyear_current_cat==period],family="binomial")
  gres <- ggpredict(lreg,"courier")
  probs_adj_df <- rbind(probs_adj_df,data.frame(calyear_cat=rep(period,2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
}
rm(reg_formula)

# crude probabilities, with gen. estim. equations
probs_gee_crude_df <- data.frame(NULL)
reg_formula <- as.formula("vls_ind ~ courier")
for(period in DTrna[,levels(calyear_current_cat)])
{
  lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat==period],family="binomial",corstr="exchangeable",id=patient)
  gres <- ggpredict(lreg,"courier")
  probs_gee_crude_df <- rbind(probs_gee_crude_df,data.frame(calyear_cat=rep(period,2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
}
rm(reg_formula)

# adjusted probabilities, with gen. estim. equations
probs_gee_adj_df <- data.frame(NULL)
reg_formula <- as.formula("vls_ind ~ courier + mhd_ind + sex + age_current_cat")
for(period in DTrna[,levels(calyear_current_cat)])
{
  lreg <- geeglm(reg_formula,data=DTrna[calyear_current_cat==period],family="binomial",corstr="exchangeable",id=patient)
  gres <- ggpredict(lreg,"courier")
  probs_gee_adj_df <- rbind(probs_gee_adj_df,data.frame(calyear_cat=rep(period,2),courier=c(0,1),probability=100*gres$predicted,lower=100*gres$conf.low,upper=100*gres$conf.high))
}
rm(reg_formula)

probs_crude_df <- data.table(probs_crude_df)
probs_adj_df <- data.table(probs_adj_df)
probs_gee_crude_df <- data.table(probs_gee_crude_df)
probs_gee_adj_df <- data.table(probs_gee_adj_df)

probs_crude_df$courier <- factor(probs_crude_df$courier)
probs_adj_df$courier <- factor(probs_adj_df$courier)
probs_gee_crude_df$courier <- factor(probs_gee_crude_df$courier)
probs_gee_adj_df$courier <- factor(probs_gee_adj_df$courier)

pp_crude <- ggplot(probs_crude_df,aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
  geom_point(size=2,position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
  labs(x="Time period",y="Percentage virally suppressed") +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
  scale_shape_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_color_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_x_discrete(labels=c("2011-2013","2014-2016","2017-2019","2020-2022")) +
  ylim(c(80,100))
ggsave(pp_crude,filename=file.path(filepath_plot,paste0("crude_probs_by_period_",VLS_threshold,".png")),height=4,width=6,dpi=600)

pp_adj <- ggplot(probs_adj_df,aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
  geom_point(size=2,position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
  labs(x="Time period",y="Percentage virally suppressed") +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
  scale_shape_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_color_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_x_discrete(labels=c("2011-2013","2014-2016","2017-2019","2020-2022")) +
  ylim(c(80,100))
ggsave(pp_adj,filename=file.path(filepath_plot,paste0("adjusted_probs_by_period_",VLS_threshold,".png")),height=4,width=6,dpi=600)

pp_gee_crude <- ggplot(probs_gee_crude_df,aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
  geom_point(size=2,position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
  labs(x="Time period",y="Percentage virally suppressed") +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
  scale_shape_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_color_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_x_discrete(labels=c("2011-2013","2014-2016","2017-2019","2020-2022")) +
  ylim(c(70,100))
ggsave(pp_gee_crude,filename=file.path(filepath_plot,paste0("gee_crude_probs_by_period_",VLS_threshold,".png")),height=4,width=6,dpi=600)

pp_gee_adj <- ggplot(probs_gee_adj_df,aes(x=calyear_cat,y=probability,color=courier,shape=courier)) +
  geom_point(size=2,position=position_dodge(0.5)) +
  geom_errorbar(aes(ymax=upper,ymin=lower),width=0.5,position=position_dodge(0.5)) +
  labs(x="Time period",y="Percentage virally suppressed") +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank()) +
  scale_shape_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_color_discrete(labels=c("Non-courier","Courier"),name=NULL) +
  scale_x_discrete(labels=c("2011-2013","2014-2016","2017-2019","2020-2022")) +
  ylim(c(70,100))
ggsave(pp_gee_adj,filename=file.path(filepath_plot,paste0("gee_adjusted_probs_by_period_",VLS_threshold,".png")),height=4,width=6,dpi=600)

format_CI <- function(estimate,lower,upper,digits=2)
  paste0(format(round(estimate,digits),nsmall=digits)," (",
         format(round(lower,digits),nsmall=digits),"-",
         format(round(upper,digits),nsmall=digits),")")

probs_crude_out <- probs_crude_df[,.(calyear_cat,courier,out=format_CI(probability,lower,upper,digits=1))]
probs_crude_out <- dcast(probs_crude_out,courier~calyear_cat,value.var="out")
probs_adj_out <- probs_adj_df[,.(calyear_cat,courier,out=format_CI(probability,lower,upper,digits=1))]
probs_adj_out <- dcast(probs_adj_out,courier~calyear_cat,value.var="out")
probs_gee_adj_out <- probs_gee_adj_df[,.(calyear_cat,courier,out=format_CI(probability,lower,upper,digits=1))]
probs_gee_adj_out <- dcast(probs_gee_adj_out,courier~calyear_cat,value.var="out")
probs_gee_crude_out <- probs_gee_crude_df[,.(calyear_cat,courier,out=format_CI(probability,lower,upper,digits=1))]
probs_gee_crude_out <- dcast(probs_gee_crude_out,courier~calyear_cat,value.var="out")

write_xlsx(probs_crude_out,path=file.path(filepath_write,paste0("probs_crude_by_period_",VLS_threshold,".xlsx")))
write_xlsx(probs_adj_out,path=file.path(filepath_write,paste0("probs_adj_by_period_",VLS_threshold,".xlsx")))
write_xlsx(probs_gee_crude_out,path=file.path(filepath_write,paste0("probs_gee_crude_by_period_",VLS_threshold,".xlsx")))
write_xlsx(probs_gee_adj_out,path=file.path(filepath_write,paste0("probs_gee_adj_by_period_",VLS_threshold,".xlsx")))

toc()
