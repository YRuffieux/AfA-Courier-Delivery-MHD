# descriptive table by VL measurement, stratified by various combinations of calendar period, courier status, and medical scheme
# various descriptive plots

library(tictoc)
library(data.table)
library(writexl)
library(tableone)
library(ggplot2)
library(scales)

filepath_load <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Plots"

tic()

load(file=file.path(filepath_load,"AfA_VL.RData"))

setorder(DTrna,"patient","rna_d")
DTrna <- DTrna[!is.na(rna_v)]

# formatting
DTrna[,sex:=as.character(sex)]
DTrna[sex=="1",sex:="Male"]
DTrna[sex=="2",sex:="Female"]
DTrna[,`:=`(age_current_cat=cut(age_current,breaks=c(15,30,40,50,60,70,Inf),right=FALSE),
            calyear_current_cat=cut(year(rna_d),breaks=c(2011,2014,2017,2020,Inf),right=FALSE),
            sex=factor(sex,levels=c("Male","Female")),
            art_type_cf=factor(art_type_cf,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI")),
            VLS_50=factor(ifelse(rna_v<50,"Suppressed","Unsuppressed"),levels=c("Suppressed","Unsuppressed")),
            VLS_400=factor(ifelse(rna_v<400,"Suppressed","Unsuppressed"),levels=c("Suppressed","Unsuppressed")),
            VLS_1000=factor(ifelse(rna_v<1000,"Suppressed","Unsuppressed"),levels=c("Suppressed","Unsuppressed")),
            courier=as.character(courier),
            mhd_ind=as.character(mhd_ind),
            scheme_code_base=NULL)]
DTrna[courier==0,courier:="No"]
DTrna[courier==1,courier:="Yes"]
DTrna[,courier:=factor(courier,levels=c("No","Yes"))]
DTrna[mhd_ind==0,mhd_ind:="No"]
DTrna[mhd_ind==1,mhd_ind:="Yes"]
DTrna[,mhd_ind:=factor(mhd_ind,levels=c("No","Yes"))]
DTrna[!scheme_code%in%c("BON","PLM"),scheme_code:="Other"]
DTrna[,scheme_code:=factor(scheme_code,levels=c("BON","PLM","Other"))]

# number of tests per patient
DTrna[,`:=`(N_tests=.N,N_tests_courier=sum(courier=="Yes"),N_tests_noncourier=sum(courier=="No")),by="patient"]
DTu <- unique(DTrna[,.(patient,N_tests,N_tests_courier,N_tests_noncourier)])
print("Median number of VL measurements per patient while on courier delivery: ")
print(DTu[,quantile(N_tests_courier,p=c(0.25,0.5,0.75))])
print("Median number of VL measurements per patient while not on courier delivery: ")
print(DTu[,quantile(N_tests_noncourier,p=c(0.25,0.5,0.75))])
print("Median number of VL measurements per patient overall: ")
print(DTu[,quantile(N_tests,p=c(0.25,0.5,0.75))])
rm(DTu)

remove_space <- function(x) gsub("\\( ","\\(",x)

# stratified by courier status (yes/no)
courier_df <- CreateTableOne(vars = c("mhd_ind","sex","age_current_cat","age_current","calyear_current_cat","art_type_cf","VLS_400","scheme_code"),
                                     strata="courier",data=DTrna,test=FALSE,addOverall=TRUE,includeNA=TRUE)
courier_df <- print(courier_df,nonnormal="age_current",showAllLevels=TRUE,printToggle=FALSE)
courier_df <- data.table(cbind(row.names(courier_df),courier_df))
cols <- c("Overall","No","Yes")
courier_df[,(cols):=lapply(.SD,remove_space),.SDcols=cols]
courier_df <- courier_df[,.(V1,level,No,Yes,Overall)]
write_xlsx(courier_df,path=file.path(filepath_tables,"descriptive_VL_by_courier_status.xlsx"))

# stratified by medical scheme
scheme_df <- CreateTableOne(vars = c("mhd_ind","sex","age_current_cat","age_current","calyear_current_cat","art_type_cf","VLS_400","courier"),
                         strata=c("scheme_code"),data=DTrna,test=FALSE,addOverall=TRUE,includeNA=TRUE)
scheme_df <- print(scheme_df,nonnormal="age_current",showAllLevels=TRUE,printToggle=FALSE)
scheme_df <- data.table(cbind(row.names(scheme_df),scheme_df))
cols <- c("Overall","BON","PLM","Other")
scheme_df[,(cols):=lapply(.SD,remove_space),.SDcols=cols]
scheme_df <- scheme_df[,.(V1,level,BON,PLM,Other,Overall)]
write_xlsx(scheme_df,path=file.path(filepath_tables,"descriptive_VL_by_scheme.xlsx"))

# colorblind-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Number of VL tests by year and medical scheme (BON/PLM/Other)
DTrna[,`:=`(scheme_code_cat=scheme_code,calyear=factor(year(rna_d)))]
DTrna[!scheme_code%in%c("BON","PLM"),scheme_code_cat:="Other"]
DTrna[,scheme_code_cat:=factor(scheme_code_cat,levels=c("BON","PLM","Other"))]
pp_stack_scheme <- ggplot(data=DTrna,aes(x=calyear,fill=scheme_code_cat)) +
  geom_bar(position="stack",stat="count") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.title=element_text(size=10), axis.text=element_text(size=10),legend.position="bottom") +
  scale_fill_manual(name="Scheme",labels=c("A","B","Other"),values=cbPalette[c(4,6,8)])+
  labs(x="Year of viral load test",y="Count",fill="Medical scheme")
ggsave(pp_stack_scheme,filename=file.path(filepath_plot,"number_VL_test_by_year_and_scheme.png"),height=4,width=6,dpi=600)

# Number of VL tests by year and ART regimen (NNRTI/II/PI)
pp_stack_regimen <- ggplot(data=DTrna,aes(x=calyear,fill=art_type_cf)) +
  geom_bar(position="stack",stat="count") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.title=element_text(size=10), axis.text=element_text(size=10),legend.position="bottom") +
  scale_fill_manual(values=cbPalette[c(3,5,7)]) +
  labs(x="Year of viral load test",y="Count",fill="ART regimen")
ggsave(pp_stack_regimen,filename=file.path(filepath_plot,"number_VL_test_by_year_and_regimen.png"),height=4,width=6,dpi=600)

# Courier/retail status by scheme and month/year
DTrna_temp <- copy(DTrna)
DTrna_temp[,rna_m_y:=paste0(month(rna_d),"-",year(rna_d))]
df_plot <- DTrna_temp[,.(p=sum(courier=="Yes")/.N,N=.N),by=.(scheme_code,rna_m_y,year(rna_d))]
df_plot[,rna_m_y:=factor(rna_m_y,levels= c(outer(paste0(1:12,"-"),2011:2022,FUN=paste0)))]
df_plot[,x:=as.numeric(rna_m_y)]
df_plot[,`:=`(lcl=p-qnorm(0.975)*sqrt(p*(1-p)/N),ucl=p+qnorm(0.975)*sqrt(p*(1-p)/N))]
df_plot[,ucl:=pmin(ucl,1)]
setorder(df_plot,"scheme_code","rna_m_y")

lab <- c(outer(c("01-","07-"),2011:2022,FUN=paste0))
lab <- lab[-length(lab)]
pp_line_courier <- ggplot(data=df_plot[scheme_code!="PLM" | (scheme_code=="PLM" & year>=2016) ],aes(x=x,y=p,color=scheme_code)) +
  geom_line() +
  geom_ribbon(aes(ymin=lcl,ymax=ucl,fill=scheme_code),linetype=0,alpha=0.2) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="bottom",
        axis.text=element_text(size=6),
        axis.text.x=element_text(angle=90,vjust=0.5)) +
  scale_color_manual(name="Scheme",labels=c("A","B","Other"),values=cbPalette[c(4,6,8)])+
  scale_fill_manual(name="Scheme",labels=c("A","B","Other"),values=cbPalette[c(4,6,8)])+
  scale_x_continuous(breaks=seq(1,133,by=6),labels=lab) +
  scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
  labs(x="Month and year",y="Percentage on courier delivery")
ggsave(pp_line_courier,filename=file.path(filepath_plot,"proportion_on_courier_over_time_by_scheme.png"),height=4,width=6,dpi=600)

toc()