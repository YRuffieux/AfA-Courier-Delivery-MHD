# create table with one line per VL measurement, with all relevant exposures
# "filling in" stretches of untested follow-up with missing RNA values
# takes 1-2 minutes

library(data.table)
library(tictoc)
library(foreign)
library(zoo)

filepath_read <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"
filepath_processed <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/processed"

courier_lag <- 0
min_age <- 15
close_date <- as.Date("2022-07-01")

tic("Overall")

tic("Loading data")
load(file=file.path(filepath_processed,"AfA_base.RData"))
load(file=file.path(filepath_read,"tblLAB_RNA.RData"))
load(file=file.path(filepath_read,"tblREGIMEN.RData"))
load(file=file.path(filepath_read,"tblARV.RData"))
toc()

DTrna <- tblLAB_RNA[,.(patient,rna_d,rna_v)]

DTrna <- merge(DTrna,DTbas,by="patient")

N <- uniqueN(DTbas,"patient")

print(paste0("Starting number of patients: ",N))

# merging with RNA viral load

setorder(DTrna,"patient","rna_d")

DTrna <- DTrna[rna_d>start & rna_d<=end]

DTbas_untested <- DTbas[!DTrna,on="patient"]         # dataset with indiviudals having no VL test during follow-up

#### excluding 'repeated' VL measurements recursively (i.e. 4 months or less between them), good test patient: AFA0836808

print(paste0("Total number of VL tests: ",DTrna[,.N]))
setorder(DTrna,"patient","rna_d")

DTrna <- unique(DTrna,by=c("patient","rna_d"))         # removing tests on same day

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

print(paste0("*after removing repeated tests: ",DTrna[,.N]))

toc()

###### filling in stretches of follow-up with no testing, inserting dummy VLs every six months

tic("Filling in RNA tests")

# input dates assumed to be sorted already, with: start < all rna dates < end
fill_rna <- function(d=c(),start,stop,tw=365.25)
{
  x <- c(start,d,stop)
  i <- which(diff(x)>tw)
  if(length(i)==0)
    return(x[-c(1,length(x))])
  i <- i[1]
  x <- c(x[1:i],x[i]+tw/2,x[(i+1):length(x)])
  return(fill_rna(d=x[-c(1,length(x))],start=start,stop=stop,tw=tw))
}

DTrna_untested <- DTbas_untested[,.(rna_d=as.Date(fill_rna(start=as.numeric(start)[1],stop=as.numeric(end)[1])),
                                    rna_v=NA,sex,birth_d,mhd_d),
                                 by="patient"]
DTrna_untested <- DTrna_untested[!is.na(rna_d)]   # removing those with less than one year of follow-up

print(paste0("Number of untested persons added: ",uniqueN(DTrna_untested,"patient")))
print(paste0("*including ",DTrna_untested[,.N]," tests"))

# some good test cases:
# AFA0800342 has huge gap in the middle of follow-up
# AFA0800494 gap at the start
# AFA0800672 gap at the end

X <- DTrna[,.(rna_d=as.Date(fill_rna(d=as.numeric(rna_d),start=as.numeric(start)[1],stop=as.numeric(end)[1]))),
              by="patient"]
X <- merge(DTrna[,.(patient,rna_d,rna_v)],X,by=c("patient","rna_d"),all.y=TRUE)    # newly added RNA tests -> missing values
DTrna <- DTrna[,.(patient,sex,birth_d,mhd_d)][X,on="patient",mult="first"]
rm(X)
setorder(DTrna,"patient","rna_d")
gc(verbose=FALSE)

print(paste0("Number of additional tests in tested individuals: ", sum(DTrna[,sum(is.na(rna_v))])))

toc()

# appending untested individuals to dataset
DTrna_untested[,tested:=0]
DTrna <- rbind(DTrna[,.(patient,rna_d,rna_v,sex,birth_d,mhd_d,tested=1)],DTrna_untested)
rm(DTrna_untested)

#### identifying ART regimen at time of each viral load count, test patient: AFA0803914

setorder(tblREGIMEN,"patient","moddate")

DTreg <- copy(tblREGIMEN)
DTreg[,moddate:=moddate+courier_lag*365.25/12]
DTreg <- DTreg[moddate<=close_date]
DTreg[art_type=="Other",art:=0]

# filling ART interruptions with previous ART regimen (see for ex. AFA0800214, AFA1130893)
DTreg[,art_type_cf:=art_type]
DTreg[art==0,art_type_cf:=NA]
DTreg[,art_type_cf:=na.locf(art_type_cf,na.rm=FALSE),by="patient"]

# identifying 1st, 2nd, 3rd line ART
DTreg <- DTreg[!is.na(art_type_cf)]
DTreg[,switch_ind:=as.numeric(art_type_cf!=data.table::shift(art_type_cf,type="lag")),by="patient"]
DTreg[is.na(switch_ind),switch_ind:=0]
DTreg[,art_line:=cumsum(switch_ind)+1,by="patient"]

# merging with VL table
setnames(DTreg,"moddate","art_sd")
DTreg[,art_ed:=data.table::shift(art_sd,fill=close_date,type="lead"),by="patient"]
DTrna <- DTreg[DTrna[,.(patient,rna_d,rna_v,sex,birth_d,mhd_d)],on=.(patient,art_sd<rna_d,art_ed>=rna_d)]
setnames(DTrna,"art_sd","rna_d")
DTrna[,art_ed:=NULL]
DTrna[is.na(art),`:=`(art=0,drug="",art_type="None")]
if(courier_lag==0)
  stopifnot(DTrna[,all(!is.na(art_type_cf))])
DTrna <- DTrna[!is.na(art_type_cf)]

#### identifying ARV delivery type (courier/non-courier) at time of each viral load count ####

# MED_ATC_J dataset has been restricted to ARVs upstream
DTarv <- tblARV[,.(patient,med_sd,practice_number,courier)]

# manual correction, Optipharm is a courier delivery
DTarv[practice_number=="0197440",courier:=1]

# three categories for the courier pharmacy, two main ones and the rest
DTarv[,courier_cat:="None"]
DTarv[courier==1,courier_cat:="Other"]
DTarv[courier==1 & practice_number=="0126225",courier_cat:="A"]
DTarv[courier==1 & practice_number=="6065732",courier_cat:="B"]
DTarv[,courier_cat:=factor(courier_cat,levels=c("None","A","B","Other"))]

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
stopifnot(DTrna[,all(!is.na(courier))])

DTrna <- DTrna[,.(patient,sex,birth_d,rna_d=med_sd,mhd_d,rna_v,drug,art_type_cf,courier,courier_cat,
                  N_switches_arv=N_switches,art_first_line=as.numeric(art_line==1))]
DTrna[,delta:=courier-data.table::shift(courier,type="lag"),by="patient"]
DTrna[is.na(delta),delta:=0]
DTrna[,N_switches_vl:=sum(abs(delta)),by="patient"]
DTrna[,delta:=NULL]

tab_arv <- unique(DTrna,by="patient")[,prop.table(table(N_switches_arv))]
tab_vl <- unique(DTrna,by="patient")[,prop.table(table(N_switches_vl))]

##### appending MHD indicators #####

DTrna[,mhd_ind:=0]
DTrna[!is.na(mhd_d) & mhd_d<rna_d,mhd_ind:=1]

DTrna[,`:=`(age_current=as.numeric(rna_d-birth_d)/365.25)]

##### saving full VL table #############

DTrna[,`:=`(drug=NULL,birth_d=NULL,mhd_d=NULL)]

if(courier_lag==0)
{
  save(DTrna,file=file.path(filepath_processed,"AfA_VL.RData"))
} else
{
  save(DTrna,file=file.path(filepath_processed,paste0("AfA_VL_lag",courier_lag,".RData")))
}

toc()
