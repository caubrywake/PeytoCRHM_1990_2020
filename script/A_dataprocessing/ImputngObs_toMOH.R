## Peyto Observation
## inpute LOngwave from ERA to peyto os
library (CRHMr)
library(WISKIr)

## REvverything is regress to MOH!

# Create a new obs file
obs <- createObsDataframe(
  "1987-09-30",
  "2020-10-02",
  timestep = 1,
  variables = c("t", "ea", "u", "Qsi", "Qli"),
  reps = 1,
  timezone = "MST")

## load MOH
MOH<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_met/PeytoMainOld_tRHuQsiQli_hrly_6Sept1987_31July2018.obs', 'etc/GMT+7')
MOH <- changeRHtoEa(MOH)
MOH <- insertMissing(MOH)
## Ipmpute to obs
obs.MOH<- impute(obs, c(1,2,3,4,5), MOH, c(1,2,3,4,5))
plotObs(MOH)

##  Import ERA (mostly for LW) 
era<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_reanalysis/ERAI_tRHuQsiQlip_BiasCorrect2PeytoBow_Jan1979_Aug2019.obs', 'etc/GMT+7')
era<-changeRHtoEa(era)
## Regress ERA to obs
regera <- regress(obs.MOH, c(1,2,3,4,5), era,c(1,2,3,4,5), TRUE)
## Impute ERA to obs
obs.era<- impute(obs.MOH, c(1,2,3,4,5), era, c(1,2,3,4,5), regera$slope)
plotObs(obs.era)

## add i my obs already formatted
mnh<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_met/PeytoMNH_1987_2020.obs', 'etc/GMT+7')
mnh<-changeRHtoEa(mnh)
regmnh<- regress(obs.era, c(1,2,3,4,5), mnh,c(1,2,3,4,5), TRUE)
obs.mnh<- impute(obs.MOH, c(1,2,3,4,5), mnh, c(1,2,3,4,5), regmnh$slope)
plotObs(obs.mnh)

## Import precipitation
p<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_met/BowSummit_p_hrly_21Nov2008_01Oct2020.obs', 'etc/GMT+7')
obs.p <- assembleObs(obs.era, p)

# regress era P to bow summit P
regp <- regress(obs.p, c(6), era,c(6), TRUE)
## impute the data
pera<- impute(obs.p, c(6), era, c(6))
obs.pera <- assembleObs(obs.mnh, pera)

## Interpolate over larger gap
obs.filled<- interpolate(obs.pera, c(1,2,3,4,5, 6),methods = "linear",maxlength = 12)

## ea to rh
obs.filled <- changeEatoRH(obs.filled)


## Qa Qc -
# qsi cant be below 0
obs.filled <- minObs(obs.filled)
max.filled <- maxObs(obs.filled)
## Export obs
result <- writeObsFile(obs.filled, 'D:/PeytoCRHM_1990_2020/data_raw/peyto_met/PeytoMNH_1987_2020_MOH.obs')

