library(data.table)
library(tictoc)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"
filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"

tic()

load(file=file.path(filepath_processed,"AfA_VL.RData"))
load(file=file.path(filepath_read,"tblARV.RData"))

DTrna <- DTrna[!is.na(rna_v)]               # removing 'fake' tests
setorder(DTrna,"patient","rna_d")

DTrna <- DTrna[,.(patient,rna_d,start,end)]

# censoring 6 months after the last VL measurement
DTrna <- DTrna[,.(patient,last_rna_d=rna_d)][DTrna,on="patient",mult="last"]
DTrna[,end:=pmin(end,last_rna_d+365.25/2)]
stopifnot(DTrna[,all(start<end)])
DTrna[,`:=`(rna_d=NULL,last_rna_d=NULL)]

close_date <- DTrna[,max(end)]

DTrna <- unique(DTrna)

# courier status information, AFA0800228, AFA0800248, mistake with AFA0800267, probably because of both courier and retail on samy day
DTms <- tblARV[,.(patient,med_sd,practice_number,courier)]
DTms <- DTms[courier!=9]
setorder(DTms,"patient","med_sd")
DTms[,practice_number:=NULL]

DTms[,delta:=abs(courier-data.table::shift(courier,type="lag")),by="patient"]
DTms[is.na(delta),delta:=0]
DTms[,n_switch:=cumsum(delta),by="patient"]
DTms <- unique(DTms,by=c("patient","n_switch"))
DTms[,`:=`(delta=NULL,n_switch=NULL)]
DTms[,med_ed:=data.table::shift(med_sd,type="lead"),by="patient"]
DTms[is.na(med_ed),med_ed:=close_date]
rm(close_date)

DTms <- merge(DTms,DTrna,by="patient")
DTms[,`:=`(start=pmax(start,med_sd),end=pmin(end,med_ed,na.rm=TRUE))]
DTms <- DTms[start<end]
DTms[,`:=`(med_sd=NULL,med_ed=NULL)]

DTms[courier==1,`:=`(from="Courier",to="Retail")]
DTms[courier==0,`:=`(from="Retail",to="Courier")]
DTms[,status:=1]
DTms[,`:=`(n=1:.N,N=.N),by="patient"]
DTms[n==N,status:=0]
DTms[,`:=`(n=NULL,N=NULL,courier=NULL)]

# appending 'dummy' starting states (necessary for transition probabilities)
X <- unique(DTms[,.(patient,start,from)],by="patient")   # starting state for each patient (courier/retail)
X[,`:=`(end=start,from="Entry",to=from,status=1)]
X[,start:=start-1]
Y <- unique(DTms[,.(patient,start,to)],by="patient")     # non-starting state (courier/retail)
Y[,`:=`(end=start,from="Entry",status=0)]
Y[,start:=start-1]

DTms <- rbind(DTms,X[,.(patient,start,end,from,to,status)])
DTms <- rbind(DTms,Y[,.(patient,start,end,from,to,status)])
setorder(DTms,"patient","start","to")
rm(X,Y)

toc()