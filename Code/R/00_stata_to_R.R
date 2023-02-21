library(readstata13)
library(data.table)
library(tictoc)


filepath_source <- "C:/ISPM/Data/AfA-AfAc/Stata/AfA_20230109"
filepath_clean <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/Stata/clean"
filepath_write <- "C:/ISPM/Data/HIV-mental disorders/AfA_Courier_Delivery/R/clean"

tic()

tblICD10_F <- data.table(read.dta13(file.path(filepath_source,"ICD10_F.dta"),convert.factors=TRUE,nonint.factors=TRUE))
save(tblICD10_F,file=file.path(filepath_write,"ICD10_F.RData"))

tblMED_ATC_J <- data.table(read.dta13(file.path(filepath_source,"MED_ATC_J.dta"),convert.factors=TRUE,nonint.factors=TRUE))
save(tblMED_ATC_J,file=file.path(filepath_write,"MED_ATC_J.RData"))
# restricting dataset to ARVs
tblARV <-  tblMED_ATC_J[med_id%in% c("J05AE09","J05AE01","J05AR10","J05AE10","J05AE08","J05AE07","J05AF09","J05AF05","J05AF01","J05AF02","J05AF07", 
                                     "J05AF06","J05AF13","J05AG05","J05AG04","J05AG06","J05AG03","J05AG01","J05AJ01","J05AJ03","J05AX09","J05AX07","J05AX23") |
                        substr(med_id,1,5)=="J05AR"]
save(tblARV,file=file.path(filepath_write,"tblARV.RData"))
rm(tblMED_ATC_J,tblARV)
gc()

# tblART <- data.table(read.dta13(file.path(filepath_source,"tblART.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblART,file=file.path(filepath_write,"tblART.RData"))

# tblARVClaims <- data.table(read.dta13(file.path(filepath_source,"tblARVClaims.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblARVClaims,file=file.path(filepath_write,"tblARVClaims.RData"))

tblBAS <- data.table(read.dta13(file.path(filepath_source,"tblBAS.dta"),convert.factors=TRUE,nonint.factors=TRUE))
save(tblBAS,file=file.path(filepath_write,"tblBAS.RData"))

tblCOVERPERIODS <- data.table(read.dta13(file.path(filepath_source,"tblCOVERPERIODS.dta"),convert.factors=TRUE,nonint.factors=TRUE))
save(tblCOVERPERIODS,file=file.path(filepath_write,"tblCOVERPERIODS.RData"))

# tblDIS <- data.table(read.dta13(file.path(filepath_source,"tblDIS.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblDIS,file=file.path(filepath_write,"tblDIS.RData"))

# tblHOS <- data.table(read.dta13(file.path(filepath_source,"tblHOS.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblHOS,file=file.path(filepath_write,"tblHOS.RData"))

# tblLAB <- data.table(read.dta13(file.path(filepath_source,"tblLAB.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblLAB,file=file.path(filepath_write,"tblLAB.RData"))

# tblLAB_CD4 <- data.table(read.dta13(file.path(filepath_source,"tblLAB_CD4.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblLAB_CD4,file=file.path(filepath_write,"tblLAB_CD4.RData"))

tblLAB_RNA <- data.table(read.dta13(file.path(filepath_source,"tblLAB_RNA.dta"),convert.factors=TRUE,nonint.factors=TRUE))
save(tblLAB_RNA,file=file.path(filepath_write,"tblLAB_RNA.RData"))

tblLTFU <- data.table(read.dta13(file.path(filepath_source,"tblLTFU.dta"),convert.factors=TRUE,nonint.factors=TRUE))
save(tblLTFU,file=file.path(filepath_write,"tblLTFU.RData"))

# tblPROGRAM <- data.table(read.dta13(file.path(filepath_source,"tblPROGRAM.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblPROGRAM,file=file.path(filepath_write,"tblPROGRAM.RData"))

# tblREV_VITAL_STATUS <- data.table(read.dta13(file.path(filepath_source,"tblREV_VITAL_STATUS.dta"),convert.factors=TRUE,nonint.factors=TRUE))
# save(tblREV_VITAL_STATUS,file=file.path(filepath_write,"tblREV_VITAL_STATUS.RData"))

tblFUPwide <- data.table((read.dta13(file.path(filepath_clean,"FUPwide.dta"),convert.factors=TRUE,nonint.factors=TRUE)))
save(tblFUPwide,file=file.path(filepath_write,"tblFUPwide.RData"))

tblREGIMEN <- data.table((read.dta13(file.path(filepath_clean,"regimen.dta"),convert.factors=TRUE,nonint.factors=TRUE)))
save(tblREGIMEN,file=file.path(filepath_write,"tblREGIMEN.RData"))

toc()