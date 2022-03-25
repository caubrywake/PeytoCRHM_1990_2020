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

# load MNH
MNH<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_met/PeytoMain_tRHuQsiQli_hrly_18July2013_11Aug2018.obs', 'etc/GMT+7')
MNH <- changeRHtoEa(MNH)
MNH <- insertMissing(MNH)
# Impute MNH in obs data frame
obs.MNH<- impute(obs, c(1,2,3,4,5), MNH, c(1,2,3,4,5))

## Obtain 2018-2020 data MNH data
Timeseries <- findWISKItimeseries('Peyto Main New') # find the time series numbers, and copy into next section
ta <- getWISKIvalues('57841042','2018-08-12','2020-10-05', timezone='MST')#
qli <- getWISKIvalues('57843042','2018-08-12','2020-10-05', timezone='MST')#
qsi <- getWISKIvalues('57845042','2018-08-12','2020-10-05', timezone='MST')#
rh <- getWISKIvalues('57853042','2018-08-12','2020-10-05', timezone='MST')#
u <-  getWISKIvalues('57871042','2018-08-12','2020-10-05', timezone='MST')#
# write to obs file
ta <- WISKItoObs(ta, timezone =  'Etc/GMT+7')
qli <- WISKItoObs(qli, timezone = 'Etc/GMT+7')
qsi <- WISKItoObs(qsi, timezone ='Etc/GMT+7')
rh <- WISKItoObs(rh, timezone = 'Etc/GMT+7')
u <- WISKItoObs(u, timezone = 'Etc/GMT+7')
# Assemble in 1 obs file
mnh <- assembleObs(ta, rh)
mnh<- assembleObs(mnh, u)
mnh <- assembleObs(mnh, qsi)
mnh <- assembleObs(mnh, qli)
mnh <- changeRHtoEa(mnh)
# impute in obs files
obs.mnh<- impute(obs.MNH, c(1,2,3,4,5), mnh, c(1,2,3,4,5))

## load data from server
## Big gap 2020-09-24 to 2020-10-03
mnhserver <- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_met/MNHfromserver.obs', 'etc/GMT+7')
mnhserver <- changeRHtoEa(mnhserver)

regserv<- regress(obs.mnh, c(1,2,3,4,5), mnhserver,c(1,2,5,3,4))
obs.serv<- impute(obs.mnh, c(1,2,3,4,5), mnhserver, c(1,2,5, 3,4), regserv$slope, regserv$intercept)
obs.serv<- interpolate(obs.serv, c(1,2,3,4,5),methods = "linear",maxlength = 6)

## load MOH
MOH<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_met/PeytoMainOld_tRHuQsiQli_hrly_6Sept1987_31July2018.obs', 'etc/GMT+6')
MOH <- changeRHtoEa(MOH)
MOH <- insertMissing(MOH)
## Regress MOH to MNH
regmoh <- regress(obs.serv, c(1,2,3,4,5), MOH,c(1,2,3,4,5), TRUE)
## Impute MOH into MNH
obs.MOH<- impute(obs.serv, c(1,2,3,4,5), MOH, c(1,2,3,4,5), regmoh$slope)
## Put Qli back in the obs
obs.MOH$qli.1 <- obs.serv$qli.1

##  Import ERA (mostly for LW) 
era<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_reanalysis/ERAI_tRHuQsiQlip_BiasCorrect2PeytoBow_Jan1979_Aug2019.obs', 'etc/GMT+6')
era<-changeRHtoEa(era)
## Regress ERA to obs
regera <- regress(obs.MOH, c(1,2,3,4,5), era,c(1,2,3,4,5), TRUE)
## Impute ERA to obs
obs.era<- impute(obs.MOH, c(1,2,3,4,5), era, c(1,2,3,4,5), regera$slope)

## Load athabasca 
atha<- readObsFile('D:/FireandIce/data_process/ice_aws/iceAWS_TRHUQsiQliPalb.obs', 'etc/GMT+7')
atha<-changeRHtoEa(atha)
regatha <- regress(obs.era, c(1,2,3,4,5),atha,c(1,2,3,4,5), TRUE)
obs.atha<- impute(obs.era, c(1,2,3,4,5), atha, c(1,2,3,4,5), regatha$slope)

## Import precipitation
p<- readObsFile('D:/PeytoCRHM_1990_2020/data_raw/peyto_met/BowSummit_p_hrly_21Nov2008_01Oct2020.obs', 'etc/GMT+7')
obs.p <- assembleObs(obs.era, p)

# regress era P to bow summit P
regp <- regress(obs.p, c(6), era,c(6), TRUE)
## impute the data
pera<- impute(obs.p, c(6), era, c(6))
obs.pera <- assembleObs(obs.atha, pera)

## Interpolate over larger gap
obs.filled<- interpolate(obs.pera, c(1,2,3,4,5, 6),methods = "linear",maxlength = 12)

## ea to rh
obs.filled <- changeEatoRH(obs.filled)


## Qa Qc -
# qsi cant be below 0
obs.filled <- minObs(obs.filled)
max.filled <- maxObs(obs.filled)
## Export obs
result <- writeObsFile(obs.filled, 'D:/PeytoCRHM_1990_2020/data_raw/peyto_met/PeytoMNH_1987_2020.obs')

