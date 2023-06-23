library(mstate)
library(data.table)
library(tictoc)
library(ggplot2)

filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_plot <- "C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Output/Plots"

tic()

load(file=file.path(filepath_processed,"AfA_mstate.RData"))

# plotting % on courier at cohort entry by calendar year
X <- DTms[status==1 & from=="Entry"]
X[,start_year:=year(end)]
Y <- X[,.(Courier=sum(to=="Courier")/.N),by="start_year"]
rm(X)
setorder(Y,"start_year")
Y[,`:=`(start_year=factor(start_year),Retail=1-Courier)]
Y <- melt(Y,id.vars="start_year",measure_vars=c("Courier","Retail"),variable.name="outcome",value.name="pct")
Y[,`:=`(outcome=factor(outcome,levels=c("Retail","Courier")),start_year=factor(start_year))]

pp_stack <- ggplot(data=Y,aes(x=start_year,y=pct,fill=outcome)) +
  geom_bar(position="stack",stat="identity") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom") +
  scale_y_continuous(labels=scales::percent,limits=c(-0.005,1.005),expand=c(0.01,0)) +
  labs(x = "Year of entry",y="Percentage",fill=NULL)
ggsave(pp_stack,filename=file.path(filepath_plot,"starting_probabilities_by_calyear.png"),height=4,width=6,dpi=600)
  

# defining transitions, everything needs to be numerical!
# state: 1=Entry, 2=Retail, 3=Courier
# transitions: 1=Entry->Retail, 2=Entry->Courier, 3=Retail->Courier, 4=Courier->Retail
DTms[,trans:=0]
DTms[from=="Entry" & to=="Retail",trans:=1]
DTms[from=="Entry" & to=="Courier",trans:=2]
DTms[from=="Retail" & to=="Courier",trans:=3]
DTms[from=="Courier" & to=="Retail",trans:=4]
stopifnot(DTms[,all(trans!=0)])
DTms[from=="Entry",from:=1]
DTms[from=="Retail",from:=2]
DTms[from=="Courier",from:=3]
DTms[to=="Retail",to:=2]
DTms[to=="Courier",to:=3]
DTms[,`:=`(from=as.numeric(from),to=as.numeric(to))]

tmat <- matrix(c(NA,1,2,NA,NA,3,NA,4,NA),nrow=3,byrow=TRUE)
dimnames(tmat) <- list(from=c("Entry","Retail","Courier"),to=c("Entry","Retail","Courier"))

# switching to follow-up scale, t=0 corresponds to first VL load, unless the individual has been left-truncated
setorder(DTms,"patient","start")
DTms <- DTms[,.(patient,start_min=start+1)][DTms,on="patient",mult="first"]
DTms[,`:=`(start=as.numeric(start-start_min),end=as.numeric(end-start_min))]
DTms[,start_min:=NULL]

# making mstate object
attr(DTms,"trans") <- tmat
class(DTms) <- c("msdata","data.frame")

print("Transition frequencies (pretty meaningless here):")
print(events(DTms))

# transition probabilities, starting from Entry state
creg <- coxph(Surv(start,end,status)~strata(trans),data=DTms,method="breslow")
msf <- msfit(object=creg,vartype="aalen",trans=tmat)
pt <- data.table(probtrans(msf,predt=-1,method="aalen")[[1]])
pt[,time:=time/365.25]
setnames(pt,c("pstate2","pstate3"),c("Retail","Courier"))
pt_melt_point <- melt(pt,id.vars="time",measure.vars=c("Retail","Courier"),variable.name="Method",value.name="Probability")
pt_melt_point[,Method:=factor(Method,levels=c("Retail","Courier"))]
pt_melt_point <- pt_melt_point[time>=0 & time<=10]

pp_tp <- ggplot(pt_melt_point, aes(x=time,y=Probability,group=Method,fill=Method)) +
  geom_area(position="fill") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.title=element_text(size=10), axis.text=element_text(size=10), legend.position="bottom") +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,by=2), expand = c(0, 0.1)) +
  scale_y_continuous(labels=scales::percent, limits=c(-0.005,1.005),expand=c(0.01, 0)) +
  labs(x = "Year of follow-up",y="Percentage",fill=NULL)
ggsave(pp_tp,filename=file.path(filepath_plot,"transition_probabilities.png"),height=4,width=6,dpi=600)

toc()