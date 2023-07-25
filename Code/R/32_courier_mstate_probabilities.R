# plotting probabilities of courier/retail over follow-up

library(mstate)
library(data.table)
library(tictoc)
library(ggplot2)
library(writexl)
library(patchwork)

filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Plots"
filepath_tables <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Tables"

tic()

which_schemes <- c("All","BON","PLM","Other")

load(file=file.path(filepath_processed,"AfA_mstate.RData"))

DTms[!scheme_code_base%in%c("BON","PLM"),scheme_code_base:="Other"]
DTms[,scheme_code_base:=factor(scheme_code_base,levels=c("BON","PLM","Other"))]

# plotting % in each scheme (PLM, BON, Other) by entry year
X <- unique(DTms[from!="Entry"],by="patient")
X[,start_year:=factor(year(start))]
pp_stack_scheme <- ggplot(data=X,aes(x=start_year,fill=scheme_code_base)) +
  geom_bar(position="stack",stat="count") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom") +
  labs(x = "Year of entry",y="Count",fill="Medical scheme")
ggsave(pp_stack_scheme,filename=file.path(filepath_plot,"nb_entries_by_year_and_scheme.png"),height=4,width=6,dpi=600)

pp_list <- list(NULL)

for(scheme in which_schemes)
{
  DTms_temp <- copy(DTms)
  
  if(scheme!="All")
    DTms_temp <- DTms_temp[scheme_code_base==scheme]
  
  # plotting % on courier at cohort entry by calendar year
  X <- DTms_temp[status==1 & from=="Entry"]
  X[,start_year:=year(end)]
  Y <- X[,.(Courier=sum(to=="Courier")/.N),by="start_year"]
  rm(X)
  setorder(Y,"start_year")
  Y[,`:=`(start_year=factor(start_year),Retail=1-Courier)]
  Y <- melt(Y,id.vars="start_year",measure_vars=c("Courier","Retail"),variable.name="outcome",value.name="pct")
  Y[,`:=`(outcome=factor(outcome,levels=c("Retail","Courier")),start_year=factor(start_year))]
  
  pp_stack_courier <- ggplot(data=Y,aes(x=start_year,y=pct,fill=outcome)) +
    geom_bar(position="stack",stat="identity") +
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom") +
    scale_y_continuous(labels=scales::percent,limits=c(-0.005,1.005),expand=c(0.01,0)) +
    labs(x = "Year of entry",y="Percentage",fill=NULL)
  if(scheme=="All")
  {
    ggsave(pp_stack_courier,filename=file.path(filepath_plot,"starting_probabilities_by_year_and_courier.png"),height=4,width=6,dpi=600)
  } else
  {
    ggsave(pp_stack_courier,filename=file.path(filepath_plot,paste0("starting_probabilities_by_year_and_courier_",scheme,".png")),height=4,width=6,dpi=600)
  }
  rm(Y)
  
  
  # defining transitions, everything needs to be numerical!
  # state: 1=Entry, 2=Retail, 3=Courier
  # transitions: 1=Entry->Retail, 2=Entry->Courier, 3=Retail->Courier, 4=Courier->Retail
  DTms_temp[,trans:=0]
  DTms_temp[from=="Entry" & to=="Retail",trans:=1]
  DTms_temp[from=="Entry" & to=="Courier",trans:=2]
  DTms_temp[from=="Retail" & to=="Courier",trans:=3]
  DTms_temp[from=="Courier" & to=="Retail",trans:=4]
  stopifnot(DTms_temp[,all(trans!=0)])
  DTms_temp[from=="Entry",from:=1]
  DTms_temp[from=="Retail",from:=2]
  DTms_temp[from=="Courier",from:=3]
  DTms_temp[to=="Retail",to:=2]
  DTms_temp[to=="Courier",to:=3]
  DTms_temp[,`:=`(from=as.numeric(from),to=as.numeric(to))]
  
  tmat <- matrix(c(NA,1,2,NA,NA,3,NA,4,NA),nrow=3,byrow=TRUE)
  dimnames(tmat) <- list(from=c("Entry","Retail","Courier"),to=c("Entry","Retail","Courier"))
  
  # switching to follow-up scale, t=0 corresponds to first VL load, unless the individual has been left-truncated
  setorder(DTms_temp,"patient","start")
  DTms_temp <- DTms_temp[,.(patient,start_min=start+1)][DTms_temp,on="patient",mult="first"]
  DTms_temp[,`:=`(start=as.numeric(start-start_min),end=as.numeric(end-start_min))]
  DTms_temp[,start_min:=NULL]
  if(scheme=="PLM")  # if PLM only: right-censoring at six year mark
  {
    DTms_temp <- data.table(survSplit(Surv(start,end,status)~.,cut=6*365.25,episode="cens",data=DTms_temp))
    DTms_temp <- DTms_temp[cens==1]
    DTms_temp[,cens:=NULL]
  }
  
  # making mstate object
  attr(DTms_temp,"trans") <- tmat
  class(DTms_temp) <- c("msdata","data.frame")
  
  # transition probabilities, starting from Entry state
  creg <- coxph(Surv(start,end,status)~strata(trans),data=DTms_temp,method="breslow")
  msf <- msfit(object=creg,vartype="aalen",trans=tmat)
  pt <- data.table(probtrans(msf,predt=-1,method="aalen")[[1]])
  pt[,time:=time/365.25]
  setnames(pt,c("pstate2","pstate3","se2","se3"),c("Retail","Courier","Retail_se","Courier_se"))
  pt_melt_point <- melt(pt,id.vars="time",measure.vars=c("Retail","Courier"),variable.name="Method",value.name="Probability")
  pt_melt_point[,Method:=factor(Method,levels=c("Retail","Courier"))]
  pt_melt_point <- pt_melt_point[time>=0 & time<=10]
  
  if(scheme!="PLM")
  {
    pp_tp <- ggplot(pt_melt_point, aes(x=time,y=Probability,group=Method,fill=Method)) +
      geom_area(position="fill") +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom") +
      scale_x_continuous(limits=c(0,10),breaks=seq(0,10,by=2), expand = c(0, 0.1)) +
      scale_y_continuous(labels=scales::percent, limits=c(-0.005,1.005),expand=c(0.01, 0)) +
      ggtitle(scheme) +
      labs(x = "Year of follow-up",y="Percentage",fill=NULL)
  } else
  {
    pp_tp <- ggplot(pt_melt_point, aes(x=time,y=Probability,group=Method,fill=Method)) +
      geom_area(position="fill") +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom") +
      scale_x_continuous(limits=c(0,6),breaks=seq(0,6,by=1), expand = c(0, 0.1)) +
      scale_y_continuous(labels=scales::percent, limits=c(-0.005,1.005),expand=c(0.01, 0)) +
      ggtitle(scheme) +
      labs(x = "Year of follow-up",y="Percentage",fill=NULL) 
  }
  pp_list[[scheme]] <- pp_tp
  
  # numerical values + 95% CIs
  format_CI <- function(a,b,c,digits=2) 
    paste0(trimws(format(round(a,2),nsmall=2))," (",trimws(format(round(b,2),nsmall=2)),"-",trimws(format(round(c,2),nsmall=2)),")")
  
  pt <- pt[pstate1!=1]
  if(scheme!="PLM")
  {
    xxx <- pt[,abs(outer(time,(0:20)/2,"-"))]
  } else
  {
    xxx <- pt[,abs(outer(time,(0:12)/2,"-"))] 
  }
  yyy <- apply(xxx,2,FUN=which.min)
  
  df_out <- pt[yyy]
  df_out[,`:=`(Retail=format_CI(Retail,Retail+qnorm(0.025)*Retail_se,Retail-qnorm(0.025)*Retail_se),
               Courier=format_CI(Courier,Courier+qnorm(0.025)*Courier_se,Courier-qnorm(0.025)*Courier_se))]
  df_out <- df_out[,.(time,Retail,Courier)]
  
  if(scheme=="All")
  {
    write_xlsx(df_out,path=file.path(filepath_tables,"mstate_probabilities.xlsx"))
  } else
  {
    write_xlsx(df_out,path=file.path(filepath_tables,paste0("mstate_probabilities_",scheme,".xlsx")))
  }  
  
  # cumulative hazards
  msf_ch <- rbind(data.table(msf$Haz[msf$Haz$trans==3,]),
                  data.table(msf$Haz[msf$Haz$trans==4,]))
  msf_ch[,`:=`(trans=factor(trans),time=time/365.25)]
  msf_ch <- msf_ch[time>=0 & time<=10]
  if(scheme!="PLM")
  {
    pp_haz <- ggplot(msf_ch,aes(x=time,y=Haz,color=trans)) +
      geom_line() +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom")+
      labs(x="Year of follow-up",y="Cumulative hazard") +
      scale_x_continuous(limits=c(0,10),breaks=seq(0,10,by=2), expand = c(0, 0.1)) +
      scale_color_discrete(name = "Transition", labels = c("Retail->Courier","Courier->Retail"))
    if(scheme=="All")
    {
      ggsave(pp_haz,filename=file.path(filepath_plot,"mstate_hazards.png"),height=4,width=6,dpi=600)
    } else
    {
      ggsave(pp_haz,filename=file.path(filepath_plot,paste0("mstate_hazards_",scheme,".png")),height=4,width=6,dpi=600)
    }
  } else
  {
    pp_haz <- ggplot(msf_ch,aes(x=time,y=Haz,color=trans)) +
      geom_line() +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom")+
      labs(x="Year of follow-up",y="Cumulative hazard") +
      scale_x_continuous(limits=c(0,6),breaks=seq(0,10,by=1), expand = c(0, 0.1)) +
      scale_color_discrete(name = "Transition", labels = c("Retail->Courier","Courier->Retail"))
    ggsave(pp_haz,filename=file.path(filepath_plot,"mstate_hazards_PLM.png"),height=4,width=6,dpi=600)
  }
}

pp_grid <- pp_list[["All"]] + pp_list[["BON"]] + pp_list[["PLM"]] + pp_list[["Other"]]
ggsave(pp_grid,filename=file.path(filepath_plot,"mstate_probabilities.png"),width=8,height=6)

toc()