
#
# CVS - $Id: ERMCProbOfEffect_S3BL.R,v 1.8 2009/12/11 13:56:29 lmilazzo Exp $
#




#
#    --          function to evaluate the Sound Effects         --
#    --   taking into account the Sound Exposure Levels (SELs)  --
#
#
# the function calculates:
#
# a) the probability of no effect for each grid cell in the AOC
#    for each species and two effects TTS, PTS
#
# b) the expected numbers effected for each grid cell in the AOC
#    for each species and two effects TTS, PTS
#
# c) the uncertainties associated with a) and b)
#
getSoundEffects_SEL<- function(InitData, SystemData){ 

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  # (cell size: 0.5 degrees)
  GridLats<- InitData$InitialSELsObject$Locations[,2]
  GridLongs<- InitData$InitialSELsObject$Locations[,1]

  # no. of grid cells within the AOC
  NumberOfGridCells<- length(GridLats)

  # the grid cells within the AOC are identified by a `grid cell reference'
  # (an ID number); they are identified sequentially column by column,
  # starting from the bottom left corner;
  #
  # col 1: Lat;  col 2: Long;  col 3: ID
  GridKey<- cbind(GridLats, GridLongs, 1:NumberOfGridCells)

  # grid cell references;
  # grid location ID, per animal, per species
  GridCellReference<- InitData$GridCellReference 

  # species
  SpeciesList<- InitData$WorkingSpecies
  NumberOfSpecies<- length(SpeciesList)

  # SELs per animal, per species, accumulated during the `simulation time'
  # (duration of a simulation run);
  # the invalid simulated animals (animals out of sea) are marked with `NA'
  SELs<- SystemData$SELs
  # invalid simulated animals
  InvalidAnimals<- SystemData$AnimalOutOfSea

  # animal density and uncertainties
  GridDensity<- InitData$WorkingDensityArray
  GridDensityUncertainty<- InitData$WorkingDensityUncertainty


  # effects: TTS, PTS
  NumberOfEffects<- 2


  # define data structures ...
  ProportionEffectArray<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  NumberEffectArray<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  NumberEffectArrayUpper<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  NumberEffectArrayLower<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  TotalNumberEffected<- array(dim=c(NumberOfEffects, NumberOfSpecies))

  #--
  # ... species ...
  for(j in 1:NumberOfSpecies){

    cat(paste("data processing for species:", SpeciesList[j], "\n"))

    CurrentSpeciesInfo<- BAESpeciesInfo[which(BAESpeciesInfo$Species_ID==SpeciesList[j]),]
    EP_Params<- getEffectProbs_Params(CurrentSpeciesInfo)

    # trim the invalid simulated animals (i.e. the animals on land)
    # conventions: 1 = animals out of sea; 0 = animals in the sea
    ValidSELs<- SELs[which(InvalidAnimals[,j]!=1),j]
    ValidGridCellReference<- GridCellReference[which(InvalidAnimals[,j]!=1),j]

    # TTS and PTS probabilities associated with the SELs for the animals of a given species
    # matrix (Nx2) ... col 1: TTS;  col 2: PTS;  N = no. of animals 
    ProbabilityMatrix<- getEffectProbs_SEL(ValidSELs, EP_Params$SEL_DRParams)

    # average these for each grid cell and put into matrix form
    ExpectedProbs<- by(ProbabilityMatrix, ValidGridCellReference, function(x){apply(x,2,mean)})

    ProbMatrix<- matrix(unlist(ExpectedProbs), ncol=NumberOfEffects, byrow=T)

    GridOrder<- sort(unique(ValidGridCellReference))  
    # fill a results array allowing for the possibility not all cells were populated
    # e.g. low density animal may by chance have no animals randomly assigned to a particular cell
    # so GridOrder may not be a list of each number 1:k, some rows are therefore not filled (but will contain 0's from definition above)
    for(k in 1:length(GridOrder)){ProportionEffectArray[GridOrder[k], ,j]<- ProbMatrix[k,]} 

    # calculate expected value of affected MF in each cell for each effect type
    NumberEffectArray[,,j]<- ProportionEffectArray[,,j]*GridDensity[,j]
    # place-holder for uncertainties
    CoefficientOfVariation<- ifelse(GridDensityUncertainty[,j]==0, 1, GridDensityUncertainty[,j])
    # use coefficient of variation (uncertainty) and mean (abund ests) to get lognormal distribution parameters
    LMean<- log(GridDensity[,j])-0.5*log(CoefficientOfVariation^2+1)  # mean of distn on log scale
    LStDev<- sqrt(log(CoefficientOfVariation^2+1))  # sd of distn on log scale
    LowerMultiplier<- qlnorm(0.025, LMean, LStDev, lower.tail=F)  # lower 2.5% bound
    UpperMultiplier<- qlnorm(0.975, LMean, LStDev, lower.tail=F)  # upper 97.5% bound
    NumberEffectArrayUpper[,,j]<- ProportionEffectArray[,,j]*UpperMultiplier
    NumberEffectArrayLower[,,j]<- ProportionEffectArray[,,j]*LowerMultiplier
  }
    
  # calculate P(X=0) from Poisson dist with mean = NumberEffectArray elements/gridcells
  ProbNoEffectArray<- exp(-NumberEffectArray)
  ProbNoEffectArrayUpper<- exp(-NumberEffectArrayUpper)
  ProbNoEffectArrayLower<- exp(-NumberEffectArrayLower)

  # tally for the total calculation region for each effect type and species
  for(j in 1:NumberOfSpecies){
    for(k in 1:NumberOfEffects){
      TotalNumberEffected[k,j]<- sum(NumberEffectArray[,k,j])
    }
  }

  # calculate P(X=0) from Poisson dist with mean = tot.num.affect
  TotalProbNoConsequence<- exp(-TotalNumberEffected)

  SoundEffects_SEL<- list(ProbNoEffect=ProbNoEffectArray,
                          ProbNoEffectUpper=ProbNoEffectArrayUpper,
                          ProbNoEffectLower=ProbNoEffectArrayLower,
                          NumbersEffected=NumberEffectArray,
                          NumbersEffectedUpper=NumberEffectArrayUpper,
                          NumbersEffectedLower=NumberEffectArrayLower,
                          TotalNumberEffected=TotalNumberEffected,
                          TotalProbNoConsequence=TotalProbNoConsequence,
                          SpeciesList=SpeciesList,
                          GridKey=GridKey)

  return(SoundEffects_SEL)

}


#
#    --          function to evaluate the Sound Effects         --
#    --   taking into account the Sound Exposure Levels (SELs)  --
#
#
# the function calculates:
#
# a) the probability of no effect for each grid cell in the AOC
#    for each species and two effects TTS, PTS
#
# b) the expected numbers effected for each grid cell in the AOC
#    for each species and two effects TTS, PTS
#
# c) the uncertainties associated with a) and b)
#
getSoundEffects_SEL_m1<- function(InitData, SystemData){ 

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  # (cell size: 0.5 degrees)
  GridLats<- InitData$InitialSELsObject$Locations[,2]
  GridLongs<- InitData$InitialSELsObject$Locations[,1]

  # no. of grid cells within the AOC
  NumberOfGridCells<- length(GridLats)

  # the grid cells within the AOC are identified by a `grid cell reference'
  # (an ID number); they are identified sequentially column by column,
  # starting from the bottom left corner;
  #
  # col 1: Lat;  col 2: Long;  col 3: ID
  GridKey<- cbind(GridLats, GridLongs, 1:NumberOfGridCells)

  # grid cell references;
  # grid location ID, per animal, per species
  GridCellReference<- InitData$GridCellReference 

  # species
  SpeciesList<- InitData$WorkingSpecies
  NumberOfSpecies<- length(SpeciesList)

  # SELs per animal, per species, accumulated during the `simulation time'
  # (duration of a simulation run);
  # the invalid simulated animals (animals out of sea) are marked with `NA'
  SELs<- SystemData$SELs
  # invalid simulated animals
  InvalidAnimals<- SystemData$AnimalOutOfSea

  # animal density and uncertainties
  GridDensity<- InitData$WorkingDensityArray
  GridDensityUncertainty<- InitData$WorkingDensityUncertainty


  # effects: TTS, PTS
  NumberOfEffects<- 2


  # define data structures ...
  ProportionEffectArray<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  NumberEffectArray<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  NumberEffectArrayUpper<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  NumberEffectArrayLower<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  TotalNumberEffected<- array(dim=c(NumberOfEffects, NumberOfSpecies))

  #--
  # ... species ...
  for(j in 1:NumberOfSpecies){

    cat(paste("data processing for species:", SpeciesList[j], "\n"))

    CurrentSpeciesInfo<- BAESpeciesInfo[which(BAESpeciesInfo$Species_ID==SpeciesList[j]),]
    EP_Params<- getEffectProbs_Params_m1(CurrentSpeciesInfo)

    # trim the invalid simulated animals (i.e. the animals on land)
    # conventions: 1 = animals out of sea; 0 = animals in the sea
    ValidSELs<- SELs[which(InvalidAnimals[,j]!=1),j]
    ValidGridCellReference<- GridCellReference[which(InvalidAnimals[,j]!=1),j]

    # TTS and PTS probabilities associated with the SELs for the animals of a given species
    # matrix (Nx2) ... col 1: TTS;  col 2: PTS;  N = no. of animals 
    ProbabilityMatrix<- getEffectProbs_SEL(ValidSELs, EP_Params$SEL_DRParams)

    # average these for each grid cell and put into matrix form
    ExpectedProbs<- by(ProbabilityMatrix, ValidGridCellReference, function(x){apply(x,2,mean)})

    ProbMatrix<- matrix(unlist(ExpectedProbs), ncol=NumberOfEffects, byrow=T)

    GridOrder<- sort(unique(ValidGridCellReference))  
    # fill a results array allowing for the possibility not all cells were populated
    # e.g. low density animal may by chance have no animals randomly assigned to a particular cell
    # so GridOrder may not be a list of each number 1:k, some rows are therefore not filled (but will contain 0's from definition above)
    for(k in 1:length(GridOrder)){ProportionEffectArray[GridOrder[k], ,j]<- ProbMatrix[k,]} 

    # calculate expected value of affected MF in each cell for each effect type
    NumberEffectArray[,,j]<- ProportionEffectArray[,,j]*GridDensity[,j]
    # place-holder for uncertainties
    CoefficientOfVariation<- ifelse(GridDensityUncertainty[,j]==0, 1, GridDensityUncertainty[,j])
    # use coefficient of variation (uncertainty) and mean (abund ests) to get lognormal distribution parameters
    LMean<- log(GridDensity[,j])-0.5*log(CoefficientOfVariation^2+1)  # mean of distn on log scale
    LStDev<- sqrt(log(CoefficientOfVariation^2+1))  # sd of distn on log scale
    LowerMultiplier<- qlnorm(0.025, LMean, LStDev, lower.tail=F)  # lower 2.5% bound
    UpperMultiplier<- qlnorm(0.975, LMean, LStDev, lower.tail=F)  # upper 97.5% bound
    NumberEffectArrayUpper[,,j]<- ProportionEffectArray[,,j]*UpperMultiplier
    NumberEffectArrayLower[,,j]<- ProportionEffectArray[,,j]*LowerMultiplier
  }
    
  # calculate P(X=0) from Poisson dist with mean = NumberEffectArray elements/gridcells
  ProbNoEffectArray<- exp(-NumberEffectArray)
  ProbNoEffectArrayUpper<- exp(-NumberEffectArrayUpper)
  ProbNoEffectArrayLower<- exp(-NumberEffectArrayLower)

  # tally for the total calculation region for each effect type and species
  for(j in 1:NumberOfSpecies){
    for(k in 1:NumberOfEffects){
      TotalNumberEffected[k,j]<- sum(NumberEffectArray[,k,j])
    }
  }

  # calculate P(X=0) from Poisson dist with mean = tot.num.affect
  TotalProbNoConsequence<- exp(-TotalNumberEffected)

  SoundEffects_SEL<- list(ProbNoEffect=ProbNoEffectArray,
                          ProbNoEffectUpper=ProbNoEffectArrayUpper,
                          ProbNoEffectLower=ProbNoEffectArrayLower,
                          NumbersEffected=NumberEffectArray,
                          NumbersEffectedUpper=NumberEffectArrayUpper,
                          NumbersEffectedLower=NumberEffectArrayLower,
                          TotalNumberEffected=TotalNumberEffected,
                          TotalProbNoConsequence=TotalProbNoConsequence,
                          SpeciesList=SpeciesList,
                          GridKey=GridKey)

  return(SoundEffects_SEL)

}


#
#    --           function to evaluate the Sound Effects           --
#    --  taking into account the Peak Sound Pressure (PSP) levels  --
#
#
# the function calculates:
#
# a) the probability of no effect for each grid cell in the AOC
#    for each species and two effects TTS, PTS
#
# b) the expected numbers effected for each grid cell in the AOC
#    for each species and two effects TTS, PTS
#
getSoundEffects_PSP<- function(InitData, SystemData){ 

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  # (cell size: 0.5 degrees)
  GridLats<- InitData$InitialSELsObject$Locations[,2]
  GridLongs<- InitData$InitialSELsObject$Locations[,1]

  # no. of grid cells within the AOC
  NumberOfGridCells<- length(GridLats)

  # the grid cells within the AOC are identified by a `grid cell reference'
  # (an ID number); they are identified sequentially column by column,
  # starting from the bottom left corner;
  #
  # col 1: Lat;  col 2: Long;  col 3: ID
  GridKey<- cbind(GridLats, GridLongs, 1:NumberOfGridCells)

  # grid cell references;
  # grid location ID, per animal, per species
  GridCellReference<- InitData$GridCellReference 


  # species
  SpeciesList<- InitData$WorkingSpecies
  NumberOfSpecies<- length(SpeciesList)

  # SPLmaxs per animal, per species, recorded from the start of the simulation
  # the invalid simulated animals (animals out of sea) are marked with `NA'
  SPLmaxs<- SystemData$SPLmaxs
  # invalid simulated animals
  InvalidAnimals<- SystemData$AnimalOutOfSea

  # animal density and uncertainties
  GridDensity<- InitData$WorkingDensityArray
  GridDensityUncertainty<- InitData$WorkingDensityUncertainty


  # effects: TTS, PTS
  NumberOfEffects<- 2


  # define data structures ...
  ProportionEffectArray<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  NumberEffectArray<- array(0, dim=c(NumberOfGridCells, NumberOfEffects, NumberOfSpecies)) 
  TotalNumberEffected<- array(dim=c(NumberOfEffects, NumberOfSpecies))

  #--
  # ... species ...
  for(j in 1:NumberOfSpecies){

    cat(paste("data processing for species:", SpeciesList[j], "\n"))

    CurrentSpeciesInfo<- BAESpeciesInfo[which(BAESpeciesInfo$Species_ID==SpeciesList[j]),]
    EP_Params<- getEffectProbs_Params(CurrentSpeciesInfo)

    # trim the invalid simulated animals (i.e. the animals on land)
    # conventions: 1 = animals out of sea; 0 = animals in the sea
    ValidSPLmaxs<- SPLmaxs[which(InvalidAnimals[,j]!=1),j]
    ValidGridCellReference<- GridCellReference[which(InvalidAnimals[,j]!=1),j]

    # TTS and PTS probabilities associated with the PSPs for the animals of a given species
    # matrix (Nx2) ... col 1: TTS;  col 2: PTS;  N = no. of animals 
    ProbabilityMatrix<- getEffectProbs_PSP(ValidSPLmaxs, EP_Params$PSP_TCParams)

    # average these for each grid cell and put into matrix form
    ExpectedProbs<- by(ProbabilityMatrix, ValidGridCellReference, function(x){apply(x,2,mean)})

    ProbMatrix<- matrix(unlist(ExpectedProbs), ncol=NumberOfEffects, byrow=T)

    GridOrder<- sort(unique(ValidGridCellReference))  
    # fill a results array allowing for the possibility not all cells were populated
    # e.g. low density animal may by chance have no animals randomly assigned to a particular cell
    # so GridOrder may not be a list of each number 1:k, some rows are therefore not filled (but will contain 0's from definition above)
    for(k in 1:length(GridOrder)){ProportionEffectArray[GridOrder[k], ,j]<- ProbMatrix[k,]} 

    # calculate expected value of affected MF in each cell for each effect type
    NumberEffectArray[,,j]<- ProportionEffectArray[,,j]*GridDensity[,j]

  }
    
  # calculate P(X=0) from Poisson dist with mean = NumberEffectArray elements/gridcells
  ProbNoEffectArray<- exp(-NumberEffectArray)

  # tally for the total calculation region for each effect type and species
  for(j in 1:NumberOfSpecies){
    for(k in 1:NumberOfEffects){
      TotalNumberEffected[k,j]<- sum(NumberEffectArray[,k,j])
    }
  }

  # calculate P(X=0) from Poisson dist with mean = tot.num.affect
  TotalProbNoConsequence<- exp(-TotalNumberEffected)

  SoundEffects_PSP<- list(ProbNoEffect=ProbNoEffectArray,
                          NumbersEffected=NumberEffectArray,
                          TotalNumberEffected=TotalNumberEffected,
                          TotalProbNoConsequence=TotalProbNoConsequence,
                          SpeciesList=SpeciesList,
                          GridKey=GridKey)

  return(SoundEffects_PSP)

}
