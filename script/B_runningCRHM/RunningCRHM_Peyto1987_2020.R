# Running CRHM
library(CRHMr)
# CHRM set up
# Running CRHM
install.packages("usethis")
library(usethis)
install.packages("backports")
library(backports)
install.packages("devtools")
library(devtools)
install_github("CentreForHydrology/CRHMr")
library(CRHMr)

# Just formatting flow data to also use a validation
sr50<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/sr50_water/sr50_water.obs','etc/GMT+7')
result <- writeObsFile(sr50, 'D:/PeytoCRHM_1990_2020/data_raw/sr50_water/sr50_water.obs')

##########################################################################################
## ## Running Peyto 1987-2020 simulations

####
prjname <- 'D:/PeytoCRHM_1990_2020/chrm/projectfile/Peyto1987_2020.prj'

## 2 basinflow
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/2_Peyto1987_2020_basinflow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/2_basinflow.txt')

## 3 glacier h20
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/3_Peyto1987_2020_glacierh20.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/3_glacierh20.txt')


## 1 SWE
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/1_Peyto1987_2020_SWE.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/1_SWE.txt')

## 4 melt
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/4_Peyto1987_2020_melt.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/4_melt.txt')

## 5 T Rain Snow
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/5_Peyto1987_2020_TRainSnow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/5_TRainSnow.txt')

## 6 blowing snow
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/6_Peyto1987_2020_blowingsnow.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/6_BS.txt')

## 7 evap
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/7_Peyto1987_2020_evap.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/7_evap.txt')

## 8 SWEslope
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/8_Peyto1987_2020_Sweslope.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/8_Sweslope.txt')

## 9 Meltrunoff Infil
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/9_Peyto1987_2020_MeltrunoffInfil.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/9_MeltrunoffInfil.txt')

## 10 Runoff Snowinfil
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/10_Peyto1987_2020_RunoffSnowinfil.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/10_RunoffSnowinfil.txt')

## 11 ET
filename <- 'D:/PeytoCRHM_1990_2020/chrm/setoutputvar/11_Peyto1987_2020_ET.prj'
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)
setPrjOutputVariables(prjname, variables)
result <- automatePrj(prjname)
result<-runCRHM(CRHMfile = 'D:/PeytoCRHM_1990_2020/chrm/crhm_20200702/CRHM.exe',prjname, outFile='D:/PeytoCRHM_1990_2020/chrm/output/11_ET.txt')
