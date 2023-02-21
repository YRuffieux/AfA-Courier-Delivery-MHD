# create table with one line per VL measurement, with all relevant exposures
# takes 1-2 minutes

library(data.table)
library(tictoc)
library(foreign)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"
filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_stata <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/Stata/clean"

courier_lag <- 0             # how many months to 'lag' the courier delivery and ART regimen
min_age <- 15
close_date <- as.Date("2022-07-01")
which_mhd <- c("mhd","org","sub","psy","bip","dep","anx","adj")

tic("Overall")

tic("Loading data")
load(file=file.path(filepath_processed,"AfA_MHD.RData"))
load(file=file.path(filepath_read,"tblLAB_RNA.RData"))
load(file=file.path(filepath_read,"tblREGIMEN.RData"))
load(file=file.path(filepath_read,"tblARV.RData"))
toc()

DTrna <- tblLAB_RNA[,.(patient,rna_d,rna_v)]
rm(tblLAB_RNA)
gc()

DTrna <- merge(DTrna,DTmhd[,.(patient,birth_d,start,end)],by="patient")

N <- uniqueN(DTmhd,"patient")

print(paste0("Starting number of patients: ",N))

# merging with RNA viral load

setorder(DTrna,"patient","rna_d")
DTrna <- DTrna[rna_d>start & rna_d<=end]
DTrna[,`:=`(start=NULL,end=NULL)]

DTrna <- DTrna[as.numeric(rna_d-birth_d)/365.25>=min_age]

#### identifying ART regimen at time of each viral load count, test patient: AFA0803914

setorder(tblREGIMEN,"patient","moddate")

DTreg <- copy(tblREGIMEN)
DTreg[,moddate:=moddate+courier_lag*365.25/12]
DTreg <- DTreg[moddate<=close_date]

setnames(DTreg,"moddate","art_sd")
DTreg[,art_ed:=data.table::shift(art_sd,fill=close_date,type="lead"),by="patient"]

# merging with VL table
DTrna <- DTreg[DTrna[,.(patient,rna_d,rna_v)],on=.(patient,art_sd<rna_d,art_ed>=rna_d)]
setnames(DTrna,"art_sd","rna_d")
DTrna[,art_ed:=NULL]
DTrna[is.na(art),`:=`(art=0,drug="",art_type="None")]

DTrna <- DTrna[art!=0]
N_prev <- N
N <- uniqueN(DTrna,"patient")
print(paste0("*after excluding invididuals with no RNA VL while on ART and when above age ",min_age-1,": ",N, " (",N-N_prev,")"))

#### identifying ARV delivery type (courier/non-courier) at time of each viral load count ####

# MED_ATC_J dataset has been restricted to ARVs upstream
DTarv <- tblARV[,.(patient,med_sd,courier)]
DTarv <- DTarv[courier!=9]
setorder(DTarv,"patient","med_sd","courier")
rm(tblARV)
gc()

DTarv[,med_sd:=med_sd+courier_lag*365.25/12]

# switches
DTarv[,delta:=courier-data.table::shift(courier,type="lag"),by="patient"]
DTarv[is.na(delta),delta:=0]
DTarv[,N_switches:=sum(abs(delta)),by="patient"]
DTarv[,med_ed:=data.table::shift(med_sd,fill=close_date,type="lead"),by="patient"]

# merging with VL table
DTrna <- DTarv[DTrna,on=.(patient,med_sd<rna_d,med_ed>=rna_d)]
DTrna <- DTrna[!is.na(courier)]

N_prev <- N
N <- uniqueN(DTrna,"patient")
print(paste0("*after excluding invididuals not in ARV table or with no viral load while receiving ARVs: " ,N, " (",N-N_prev,")"))

DTrna <- DTrna[,.(patient,rna_d=med_sd,rna_v,drug,art_type,courier,N_switches_arv=N_switches)]
DTrna[,delta:=courier-data.table::shift(courier,type="lag"),by="patient"]
DTrna[is.na(delta),delta:=0]
DTrna[,N_switches_vl:=sum(abs(delta)),by="patient"]
DTrna[,delta:=NULL]

tab_arv <- unique(DTrna,by="patient")[,prop.table(table(N_switches_arv))]
tab_vl <- unique(DTrna,by="patient")[,prop.table(table(N_switches_vl))]

##### appending MHD indicators #####

DTrna <- merge(DTrna,DTmhd,by="patient")

for(mhd in which_mhd)
{
  setnames(DTrna,paste0(mhd,"_d"),"exp_d")
  DTrna[,eval(paste0(mhd,"_ind")):=0]
  DTrna[!is.na(exp_d) & exp_d<rna_d,eval(paste0(mhd,"_ind")):=1]
  DTrna[,exp_d:=NULL]
}

DTrna[,`:=`(age_current=as.numeric(rna_d-birth_d)/365.25,start=NULL,end=NULL)]

#### excluding 'repeated' VL measurements (i.e. 4 months or less between them), good test patient: AFA0836808

print(paste0("Total number of VL tests: ",DTrna[,.N]))
setorder(DTrna,"patient","rna_d")

time_window <- 122         # max days for repeated VL measurement

tic("Excluding repeated VL measurements")
remove_repeated_vl <- function(d,tw)
{
  i <- which(diff(d)<tw)
  if(length(i)==0)
    return(d)
  else
    return(remove_repeated_vl(d[-(i[1]+1)],tw))
}

DTrna_reduced <- DTrna[,.(rna_d=remove_repeated_vl(rna_d,tw=time_window)),by="patient"]
DTrna <- merge(DTrna,DTrna_reduced,by=c("patient","rna_d"))
rm(DTrna_reduced)

DTrna[,delta:=as.numeric(rna_d)-data.table::shift(as.numeric(rna_d),type="lag",fill=-Inf),by="patient"]
stopifnot(all(DTrna[,delta>=time_window]))

print(paste0("*after removing repeated tests: ",DTrna[,.N]))

toc()

##### saving full VL table #############

if(courier_lag==0)
{
  save(DTrna,file=file.path(filepath_processed,"AfA_VL_courier_MHD.RData"))
  write.dta(DTrna,file=file.path(filepath_stata,"AfA_VL_courier_MHD.dta"))
} else
{
  save(DTrna,file=file.path(filepath_processed,paste0("AfA_VL_courier_MHD_lag",courier_lag,".RData")))
  write.dta(DTrna,file=file.path(filepath_stata,paste0("AfA_VL_courier_MHD_lag",courier_lag,".dta")))
}

toc()
