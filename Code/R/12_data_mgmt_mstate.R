# data preparation for multi-state analysis of courier delivery

library(data.table)
library(tictoc)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"
filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"

set.seed(1512)

tic()

load(file=file.path(filepath_processed,"AfA_VL.RData"))
load(file=file.path(filepath_read,"tblARV.RData"))
load(file=file.path(filepath_read,"tblCOVERPERIODS.RData"))

DTrna <- DTrna[!is.na(rna_v)]               # removing 'fake' tests
setorder(DTrna,"patient","rna_d")

DTrna <- DTrna[,.(patient,rna_d,start,end)]

# fetching start/end dates from RNA table: censoring 6 months after the last VL measurement, left-truncating at first VL measurement
DTrna <- DTrna[,.(patient,first_rna_d=rna_d)][DTrna,on="patient",mult="first"]
DTrna <- DTrna[,.(patient,last_rna_d=rna_d)][DTrna,on="patient",mult="last"]
DTrna[,`:=`(start=pmax(start,first_rna_d),end=pmin(end,last_rna_d+365.25/2))]
DTrna <- DTrna[start<end]
DTrna[,`:=`(rna_d=NULL,first_rna_d=NULL,last_rna_d=NULL)]

close_date <- DTrna[,max(end)]

DTrna <- unique(DTrna)

# insurance scheme at baseline
setorder(tblCOVERPERIODS,"patient","coverfrom_date")
DTrna <- tblCOVERPERIODS[,.(patient,coverfrom_date,scheme_code_base=scheme_code)][DTrna,on=.(patient,coverfrom_date<=start),mult="last"]
DTrna <- DTrna[,.(patient,start=coverfrom_date,end,scheme_code_base)]
stopifnot(DTrna[,all(!is.na(scheme_code_base))])

# courier status information from ART table
DTms <- tblARV[,.(patient,med_sd,courier)]
DTms <- DTms[courier!=9]
setorder(DTms,"patient","med_sd")

# resolving cases where an individual has both courier and retail on same day, see for ex. AFA0800267 and AFA0800526
DTms[,`:=`(n_courier=sum(courier==1),n_retail=sum(courier==0)),by=.(patient,med_sd)]
DTms[n_courier>n_retail,courier:=1]   # more courier records on the day than retail records -> courier wins
DTms[n_courier<n_retail,courier:=0]   # more retail records on the day than courier records -> retail wins
DTms[,`:=`(n_courier=NULL,n_retail=NULL)]
DTms <- DTms[sample(1:.N)]            # otherwise choosing randomly by shuffling entire dataset and picking first records each day
DTms <- unique(DTms,by=c("patient","med_sd"))
setorder(DTms,"patient","med_sd")

# identifying switches in dispensing method and throwing away all intermediate records
DTms[,delta:=abs(courier-data.table::shift(courier,type="lag")),by="patient"]
DTms[is.na(delta),delta:=0]
DTms[,n_switch:=cumsum(delta),by="patient"]
DTms <- unique(DTms,by=c("patient","n_switch"))
DTms[,`:=`(delta=NULL,n_switch=NULL)]
DTms[,med_ed:=data.table::shift(med_sd,type="lead"),by="patient"]
DTms[is.na(med_ed),med_ed:=close_date]
rm(close_date)

# left-truncation and right-censoring based on follow-up time
DTms <- merge(DTms,DTrna,by="patient")
DTms[,`:=`(start=pmax(start,med_sd),end=pmin(end,med_ed,na.rm=TRUE))]
DTms <- DTms[start<end]
DTms[,`:=`(med_sd=NULL,med_ed=NULL)]

# multi-state format
DTms[courier==1,`:=`(from="Courier",to="Retail")]
DTms[courier==0,`:=`(from="Retail",to="Courier")]
DTms[,status:=1]
DTms[,`:=`(n=1:.N,N=.N),by="patient"]
DTms[n==N,status:=0]
DTms[,`:=`(n=NULL,N=NULL,courier=NULL)]

# appending 'dummy' starting states (necessary for transition probabilities)
X <- unique(DTms[,.(patient,start=start-1,end=start,scheme_code_base,from="Entry",to=from,status=1)],by="patient")
Y <- unique(DTms[,.(patient,start=start-1,end=start,scheme_code_base,from="Entry",to,status=0)],by="patient")
DTms <- rbind(DTms,X)
DTms <- rbind(DTms,Y)
setorder(DTms,"patient","start","to")
rm(X,Y)

# saving multi-state dataset
save(DTms,file=file.path(filepath_processed,"AfA_mstate.RData"))

toc()