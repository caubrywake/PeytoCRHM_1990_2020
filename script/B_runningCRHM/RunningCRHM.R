# Running CRHM
library(CRHMr)
# CHRM set up
# Observation at Peyto

Peyto20002018 <- readObsFile('D:/InputDataComparison/data_processed/PeytoCompiled/MNH_filledwithAtha.obs','etc/GMT+7')
Peyto20002018 <- readObsFile('D:/InputDataComparison/data_processed/PeytoCompiled/MNH_filledwithAtha.obs','etc/GMT+7')
Peyto_obs <- plotObs(Peyto20002018)

filename = 'D:/Model_June2018/ProjFiles/Current/Peyto_Oct22.prj'

# Results

result <- setPrjBasinNamggsave("D:/Model_June2018/figs/R/Peyto_obs.png")

e(filename, 'Peyto2010')
result <- setPrjDates(filename,'2012 10 01', '2014 09 30')
a<- c('Firn1', 'FirnIce2', 'FirnIce3', 'Ice4', 'IceTongue5',
  'IceRock6' ,'IceRock7', 'Outlet8', 'Firn9', 'FirnIce10', 
  'Firnice11' , 'Ice12' ,'IceRock13', 'Debris14' , 'DebrisMoraine15',
  'Rock16' , 'FirnIce17' ,'Firnice18', 'IceRock19' , '20Moraine20',
  'RockT21' , 'RockM22' ,'RockL23', 'RockT24' , 'RockM25',
  'RockL26' , 'RockH27' ,'RockM28', 'Moraine29', 'IceMid30',
  'DebrisM31', 'RockH32', '33Nunatak33', 'RockH34')

#setPrjParameters()
variables <- readPrjOutputVariables(filename, asDataframe=FALSE)

#variables <- variables[-1]
setPrjOutputVariables(filename, variables)

result <- automatePrj(filename)

# Run CHRM
result<-runCRHM('D:/Model_June2018/CRHM.exe', filename)

# Read output file


