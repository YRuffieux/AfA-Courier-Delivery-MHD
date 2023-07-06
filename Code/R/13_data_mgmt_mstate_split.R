# setting up data for time-to-event analysis in mstate, with time-updated variables

library(data.table)
library(tictoc)
library(zoo)
library(survival)
library(mstate)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"
filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"

tic("Overall")

tic("Loading data")
load(file=file.path(filepath_processed,"AfA_mstate.RData"))         # mstate data with no time-updating
load(file=file.path(filepath_processed,"AfA_VL.RData"))
load(file=file.path(filepath_processed,"AfA_base.RData"))
load(file=file.path(filepath_read,"tblREGIMEN.RData"))
source("C:/ISPM/HomeDir/HIV-mental disorders/AfA_Courier_Delivery/Repository/Code/R/utils/timeSplit_DT.R")
toc()

age_breaks <- c(15,30,40,50,60,70)
calyear_breaks <- c(2011,2014,2017,2020)

DTms <- DTms[from!="Entry",.(patient,start,end,from,to,status,scheme_code_base)]         # no need for entry state in this analysis

DTrna <- DTrna[!is.na(rna_v),.(patient,rna_d,rna_v,art_type_cf)]        # removing 'fake' tests
setorder(DTrna,"patient","rna_d")
DTrna[,`:=`(vls_ind=as.numeric(rna_v<400),rna_v=NULL)]

DTrna[,delta:=vls_ind-data.table::shift(vls_ind,type="lag"),by="patient"]
DTrna[is.na(delta),delta:=0]
DTrna[,N_vls_switches:=sum(abs(delta)),by="patient"]

# values below shouldn't change after incorporating the time-updated variables
N_indiv <- uniqueN(DTms,"patient")
N_status <- DTms[,sum(status)]
N_courier_to_retail <- DTms[status==1 & from=="Courier" & to=="Retail",.N]
N_retail_to_courier <- DTms[status==1 & from=="Retail" & to=="Courier",.N]
N_py <- DTms[,sum(as.numeric(end-start))]
N_py_courier <- DTms[from=="Courier",sum(as.numeric(end-start))]
N_py_retail <- DTms[from=="Retail",sum(as.numeric(end-start))]
N_status_BON <- DTms[scheme_code_base=="BON",sum(status)]
N_py_BON <- DTms[scheme_code_base=="BON",sum(as.numeric(end-start))]

# reminder: follow-up is already left-truncated at first VL measurement and right-censored six months after last VL measurement in 'AfA_mstate.RData'
DTms_split <- DTms[,.(patient,start,end,from,to,status,scheme_code_base)]
DTms_split <- DTms_split[,.(patient,enter_d=start)][DTms_split,on="patient",mult="first"]   # date patient enters the study, need this for later
DTms_split <- DTms_split[,.(patient,leave_d=end)][DTms_split,on="patient",mult="last"]      # date patient leaves the study, need this for later
rm(DTms)

tic("Splitting, VL suppression and ART regimen")

#### Splitting at change in viral suppression status ####

# for clarity, removing intermediate VL measurements between switches
DTrna[,delta:=abs(vls_ind-data.table::shift(vls_ind)),by="patient"]
DTrna[is.na(delta),delta:=0]
DTrna[,n_switch:=cumsum(delta),by="patient"]
X <- unique(DTrna,by=c("patient","n_switch"))

# merging mstate and VL datasets
DTms_split <- merge(DTms_split,X[,.(patient,start=rna_d,vls_ind)],by=c("patient","start"),all=TRUE)
setorder(DTms_split,"patient","start")

#### Splitting at change in ART regimen ####

DTreg <- tblREGIMEN[art_type%in%c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI"),.(patient,moddate,art_type)]
DTreg[,art_type:=factor(art_type,levels=c("NNRTI+2NRTI","II+2NRTI","PI+2NRTI"))]

setorder(DTreg,"patient","moddate")

# restricting to ARV records collected during follow-up
DTreg <- DTms_split[,.(patient,enter_d,leave_d)][DTreg,on="patient",nomatch=0,mult="first"]
DTreg <- DTreg[enter_d<=moddate & moddate<leave_d]
DTreg[,`:=`(enter_d=NULL,leave_d=NULL)]

# appending ART regimen at first RNA VL measurement, i.e. cohort entry
X <- unique(DTrna,by="patient")
DTreg <- rbind(DTreg,X[,.(patient,moddate=rna_d,art_type=art_type_cf)])
rm(X)
setorder(DTreg,"patient","moddate")

# for clarity, removing intermediate ART records between switches
DTreg[,delta:=abs(art_type!=data.table::shift(art_type)),by="patient"]
DTreg[is.na(delta),delta:=0]
DTreg[,n_switch:=cumsum(delta),by="patient"]
X <- unique(DTreg,by=c("patient","n_switch"))

# merging mstate and ART datasets
DTms_split <- merge(DTms_split,X[,.(patient,start=moddate,art_type)],by=c("patient","start"),all=TRUE)
setorder(DTms_split,"patient","start")
rm(X)

toc()

# filling in empty gaps
tic("Cleanup")

cols <- c("art_type","vls_ind","from","to","leave_d","status","scheme_code_base")
DTms_split[,(cols):=lapply(.SD,na.locf,na.rm=FALSE),.SDcols=cols,by="patient"]
rm(cols)
DTms_split[,end:=data.table::shift(start,type="lead"),by="patient"]
DTms_split[to==data.table::shift(to,type="lead"),status:=0]
DTms_split[is.na(end),end:=leave_d]
DTms_split[,`:=`(enter_d=NULL,leave_d=NULL)]

# DTms_split[,`:=`(art_type=na.locf(art_type,na.rm=FALSE),
#                  vls_ind=na.locf(vls_ind,na.rm=FALSE),
#                  from=na.locf(from,na.rm=FALSE),
#                  to=na.locf(to,na.rm=FALSE),
#                  leave_d=na.locf(leave_d,na.rm=FALSE),
#                  status=na.locf(status,na.rm=FALSE),
#                  scheme_code_base=na.locf(status,na.rm=FALSE),
#                  end=data.table::shift(start,type="lead")),
#            by="patient"]
# DTms_split[to==data.table::shift(to,type="lead"),status:=0]
# DTms_split[is.na(end),end:=leave_d]
# DTms_split[,`:=`(enter_d=NULL,leave_d=NULL)]

toc()

#### Appending baseline data + other time-updated variables #####

DTms_split <- merge(DTms_split,DTbas[,.(patient,sex,birth_d,mhd_d)],by="patient")
DTms_split[,sex:=factor(sex)]

# history of mental illness, time-updated
DTms_split <- timeSplit_DT(DTms_split,vars_date="mhd_d",vars_tu="mhd_ind",event="status",start_date="start",stop_date="end",id_var="patient")
DTms_split[,mhd_d:=NULL]

# age, time-updated
DTms_split[,`:=`(start=as.numeric(start-birth_d),end=as.numeric(end-birth_d))]
DTms_split <- data.table(survSplit(Surv(start,end,status)~.,data=DTms_split,cut=365.25*age_breaks,episode="age_cat"))
DTms_split[,age_cat:=factor(age_cat)]

# calendar year, time-updated
DTms_split <- copy(DTms_split)
DTms_split[,`:=`(start=as.numeric(birth_d+start),end=as.numeric(birth_d+end))]
calyear_cutoffs <- sort(as.numeric(as.Date(c(paste0(calyear_breaks,"-01-01"),"2022-07-01"))))
DTms_split <- data.table(survSplit(Surv(start,end,status)~.,data=DTms_split,cut=calyear_cutoffs,episode="calyear_cat"))
rm(calyear_cutoffs)
DTms_split[,`:=`(calyear_cat=factor(calyear_cat),birth_d=NULL)]

# switching to follow-up scale for analysis, in days since first RNA measurement
DTms_split <- DTms_split[,.(patient,enter_d=start)][DTms_split,on="patient",mult="first"]
DTms_split[,`:=`(start=start-as.numeric(enter_d),end=end-as.numeric(enter_d))]
DTms_split[,enter_d:=NULL]

##### checking that we didn't lose anything, that everyone starts at time 0, and saving dataset #######

stopifnot(N_indiv==uniqueN(DTms_split,"patient"))
stopifnot(N_status==DTms_split[,sum(status)])
stopifnot(N_courier_to_retail==DTms_split[status==1 & from=="Courier" & to=="Retail",.N])
stopifnot(N_retail_to_courier==DTms_split[status==1 & from=="Retail" & to=="Courier",.N])
stopifnot(N_py==DTms_split[,sum(as.numeric(end-start))])
stopifnot(N_py_courier==DTms_split[from=="Courier",sum(as.numeric(end-start))])
stopifnot(N_py_retail==DTms_split[from=="Retail",sum(as.numeric(end-start))])
stopifnot(N_status_BON==DTms_split[scheme_code_base=="BON",sum(status)])
stopifnot(N_py_BON==DTms_split[scheme_code_base=="BON",sum(as.numeric(end-start))])
stopifnot(DTms_split[,all(from!=to)])
stopifnot(all(unique(DTms_split,by="patient")[,start]==0))

save(DTms_split,file=file.path(filepath_processed,"AfA_mstate_split.RData"))

toc()
