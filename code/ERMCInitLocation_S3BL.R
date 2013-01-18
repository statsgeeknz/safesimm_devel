
#
# CVS - $Id: ERMCInitLocation_S3BL.R,v 1.8 2011/06/03 15:16:02 lmilazzo Exp $
#




#
#     -- function to set the Initial Conditions --
#
initializeSystem<- function(PrimaryShipStats, MultiShipStats, UserSettings, Seeds){
 
  #
  # this function set the Initial Conditions, in particular ...
  #   a) actual/working values for Species, Animal Densities
  #   b) spatial locations
  #   c) Sound Exposure Levels (SELs) and
  #      maximum values of the Sound Propagation (SP) levels
  #

  #--
    
  # RNG - assign a seed from the predetermined seed collection
  # by using RNG kind 1 = Marsaglia-Multicarry;
  rs<- c(1, Seeds[1,])
  rs_i<- as.integer(rs)
  .Random.seed<<- rs_i


  #--

  #
  # - initial conditions - actual/working values for Species, Animal Densities
  #


  # Area of Calculation (AOC);
  # calculate the geographic region of interest AOC;
  #
  # `AOC' is a matrix; col 1: lower then upper Longs [degrees]
  #                    col 2: lower then upper Lats  [degrees]
  #
  AOC<- getAreaOfCalc(PrimaryShipStats, MultiShipStats, UserSettings)

  # Animal Densities for the AOC;
  # calculate the density estimates for all animals over the AOC;
  # density estimates are indexed by centre of 0.5 degree cells;
  #
  # note that, in the current version of SAFESIMM, the spatial distribution of
  # animal densities is considered constant; it is evaluated only at the
  # beginning of the first simulation run;
  #
  # `DensityData' is a list containing: `Density', `Uncertainty', `Species',
  #                                     `GridLong', `GridLat'
  #
  if(exists("DensityData")==F){
    DensityData<- getAnimalDensity(AOC)
    DensityData<<- DensityData
  }

  # no. of species
  NumberOfSpecies<- ncol(DensityData$Density)
  # no. of simulated animals per species
  NumberOfAnimals<- UserSettings$NumberOfSimAnimals

  SumSpeciesDensity<- vector(length=NumberOfSpecies)
  for(j in 1:NumberOfSpecies){SumSpeciesDensity[j]<- sum(DensityData$Density[,j])}

  # minimum animal density;
  # animal densities below this value are considered negligible
  tolerance<- UserSettings$ZeroTolerance
  WorkingDensityArray<- DensityData$Density[, which(SumSpeciesDensity>tolerance)]
  if(NumberOfSpecies==1){
    WorkingDensityArray<- array(WorkingDensityArray, dim=c(length(WorkingDensityArray),1))
  }
  WorkingDensityUncertainty<- DensityData$Uncertainty[, which(SumSpeciesDensity>tolerance)]
  if(NumberOfSpecies==1){
    WorkingDensityUncertainty<- array(WorkingDensityUncertainty, dim=c(length(WorkingDensityUncertainty),1))
  }
  WorkingSpecies<- DensityData$Species[which(SumSpeciesDensity>tolerance)]
  WorkingNoOfSpecies<- ncol(WorkingDensityArray)


  #--

  #
  # - initial conditions - spatial locations
  #

 
  # observation data;
  # to take into account animals that have fixed (not random) start positions    
  ObserverObject<- UserSettings$ObserverObject

  if(is.null(ObserverObject)==F){
    ObserverProportion<- ObserverObject$ObservedNumber/SumSpeciesDensity
    ObserverPortion<- ifelse(ObserverProportion>1, NumberOfAnimals, round(NumberOfAnimals*ObserverProportion,0))
    ObserverPortion<- ObserverPortion[which(SumSpeciesDensity>tolerance)]
    ObserverObject$ObservationLocation<- ObserverObject$ObservationLocation[which(SumSpeciesDensity>tolerance)]
  }
  # if no observation data ...
  if(is.null(ObserverObject)){
    ObserverPortion<- rep(0, WorkingNoOfSpecies)
  }

  WorkingSumSpeciesD<- vector(length=WorkingNoOfSpecies)
  for(j in 1:WorkingNoOfSpecies){WorkingSumSpeciesD[j]<- sum(WorkingDensityArray[,j])}

  # probability of assigning a simulated animal to a given cell;
  # the allocation of simulated animals is in proportion to estimated
  # animal density per cell
  GridProbs<- array(dim=dim(WorkingDensityArray))
  for(j in 1:WorkingNoOfSpecies){
    GridProbs[,j]<- WorkingDensityArray[,j]/WorkingSumSpeciesD[j]
  }

  #
  # the grid cells within the AOC are identified by a `grid cell reference'
  # (an ID number); they are identified sequentially column by column,
  # starting from the bottom left corner;
  #
  # define a data structure for the grid cell references (the IDs for the
  # grid locations within the AOC)
  GridCellReference<- array(dim=c(NumberOfAnimals, WorkingNoOfSpecies))

  # no. of grid cells within the AOC
  NumberOfGridCells<- nrow(WorkingDensityArray)
    
  # define a data structure for the `initial locations within the cells';
  # (later this data structure will be used for the `initial starting positions')
  XYSimLocation<- array(dim=c(NumberOfAnimals, 2, WorkingNoOfSpecies))
  # define a data structure for the `initial locations of the cells'
  # (the initial coordinates Lat, Long of the centres of the cells within the AOC);
  InitialXYGrid<- array(dim=c(NumberOfAnimals, 2, WorkingNoOfSpecies))

  
  # ... species ...
  for(j in 1:WorkingNoOfSpecies){

    # random locations within 0.5 degree grid about centre;
    # random allocation to Lat/Long cell centres;
    # the sample function allocations probabilities of selection to outcomes
    # the required number of random draws are taken based on these probabilities

    # allocate random locations ...
    if(ObserverPortion[j]!=NumberOfAnimals){
      # simulated animals with unknown/random initial starting positions
      AnimalAllocation<- (NumberOfAnimals-ObserverPortion[j])
      XYSimLocation[1:AnimalAllocation,1,j]<- runif(AnimalAllocation, -0.25, 0.25)
      XYSimLocation[1:AnimalAllocation,2,j]<- runif(AnimalAllocation, -0.25, 0.25)
      Sampler<- sample(1:NumberOfGridCells, AnimalAllocation, replace=T, prob=GridProbs[,j])
      InitialXYGrid[1:AnimalAllocation,1,j]<- DensityData$GridLong[Sampler]
      InitialXYGrid[1:AnimalAllocation,2,j]<- DensityData$GridLat[Sampler]
      GridCellReference[1:AnimalAllocation,j]<- Sampler

      # initial starting positions;
      # random location = random location within the cell +
      #                   random location of the cell
      XYSimLocation[,,j]<- XYSimLocation[,,j]+InitialXYGrid[,,j]

    }
        
    # allocate remaining positions based on observations ...        
    if(is.null(ObserverObject)==F){
      empty_loc_index<- which(is.na(XYSimLocation[,1,j]))
      for(i in empty_loc_index){
        Sampler<- sample(1:nrow(ObserverObject$ObservationLocation[[j]]), 1, replace=T)
        # inflate number of `observed' animals to correct proportion of total sim animals
        XYSimLocation[i,,j]<- ObserverObject$ObservationLocation[[j]][Sampler,]
        GridCellReference[empty_loc_index,j]<- Sampler
      }
    }
    
  }
  

  #--
  
  #
  # - initial conditions - Sound Exposure Levels (SELs) and
  #                        maximum values of the Sound Propagation (SP) levels
  #


  NumberOfLoc<- length(DensityData$GridLat)

  # no. of bins (binning procedure)
  NumberOfBins<- UserSettings$SummaryBinNumber

  # the `History' object is associated with the data from the previous
  # simulation runs;

  # `HistoricalSELs' is a data structure containing:
  # a) summary of SELs, b) species, c) locations
  HistoricalSELs<- History$HistSELsObject

  # define a data structure for the summary (histograms) of SELs;
  # SELs values per bins, per species, per location
  SELsSummary<- array(0, dim=c(NumberOfBins, WorkingNoOfSpecies, NumberOfLoc))

  # for the first simulation run ...
  if(is.null(HistoricalSELs)==T){
    InitialSELsObject<- list(SummaryValues=SELsSummary,
                             Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                             Species=WorkingSpecies)
    History$HistSELsObject<<- InitialSELsObject
  }
  # for the following simulation runs ...
  if(is.null(HistoricalSELs)==F){
    for(k in 1:NumberOfLoc){
      loc_index<- which(HistoricalSELs$Locations[,1]==DensityData$GridLong[k] &
                            HistoricalSELs$Locations[,2]==DensityData$GridLat[k]) 
      for(m in 1:WorkingNoOfSpecies){
        species_index<- which(HistoricalSELs$Species==WorkingSpecies[m])
        if(length(loc_index)>0 & length(species_index)>0){
          SELsSummary[, m, k]<- HistoricalSELs$SummaryValues[, species_index, loc_index]
        }
      }
    }
    InitialSELsObject<- list(SummaryValues=SELsSummary,
                             Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                             Species=WorkingSpecies)
  }


  # `HistoricalSPLmaxs' is a data structure containing:
  # a) summary of SPLmaxs, b) species, c) locations
  HistoricalSPLmaxs<- History$HistSPLmaxsObject

  # define a data structure for the summary (histograms) of SPLmaxs;
  # SPLmaxs values per bins, per species, per location
  SPLmaxsSummary<- array(0, dim=c(NumberOfBins, WorkingNoOfSpecies, NumberOfLoc))

  # for the first simulation run ...
  if(is.null(HistoricalSPLmaxs)==T){
    InitialSPLmaxsObject<- list(SummaryValues=SPLmaxsSummary,
                                Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                                Species=WorkingSpecies)
    History$HistSPLmaxsObject<<- InitialSPLmaxsObject
  }
  # for the following simulation runs ...
  if(is.null(HistoricalSPLmaxs)==F){
    for(k in 1:NumberOfLoc){
      loc_index<- which(HistoricalSPLmaxs$Locations[,1]==DensityData$GridLong[k] &
                            HistoricalSPLmaxs$Locations[,2]==DensityData$GridLat[k]) 
      for(m in 1:WorkingNoOfSpecies){
        species_index<- which(HistoricalSPLmaxs$Species==WorkingSpecies[m])
        if(length(loc_index)>0 & length(species_index)>0){
          SPLmaxsSummary[, m, k]<- HistoricalSPLmaxs$SummaryValues[, species_index, loc_index]
        }
      }
    }
    InitialSPLmaxsObject<- list(SummaryValues=SPLmaxsSummary,
                                Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                                Species=WorkingSpecies)
  }


  #--

  InitData<- list(AOC=AOC,
                  GridCellReference=GridCellReference,
                  InitialXYGrid=InitialXYGrid,
                  XYSimLocation=XYSimLocation,
                  WorkingDensityArray=WorkingDensityArray,
                  WorkingDensityUncertainty=WorkingDensityUncertainty,
                  WorkingSpecies=WorkingSpecies,
                  InitialSELsObject=InitialSELsObject,
                  InitialSPLmaxsObject=InitialSPLmaxsObject)
  
  return(InitData)

}


#
#     -- function to set the Initial Conditions --
#
initializeSystem_m1<- function(PrimaryShipStats, MultiShipStats, UserSettings, Seeds){
 
  #
  # this function set the Initial Conditions, in particular ...
  #   a) actual/working values for Species, Animal Densities
  #   b) spatial locations
  #   c) Sound Exposure Levels (SELs) and
  #      maximum values of the Sound Propagation (SP) levels
  #

  #--
    
  # RNG - assign a seed from the predetermined seed collection
  # by using RNG kind 1 = Marsaglia-Multicarry;
  rs<- c(1, Seeds[1,])
  rs_i<- as.integer(rs)
  .Random.seed<<- rs_i


  #--

  #
  # - initial conditions - actual/working values for Species, Animal Densities
  #


  # Area of Calculation (AOC);
  # calculate the geographic region of interest AOC;
  #
  # `AOC' is a matrix; col 1: lower then upper Longs [degrees]
  #                    col 2: lower then upper Lats  [degrees]
  #
  AOC<- getAreaOfCalc(PrimaryShipStats, MultiShipStats, UserSettings)

  # Animal Densities for the AOC;
  # calculate the density estimates for all animals over the AOC;
  # density estimates are indexed by centre of grid  cells;
  #
  # note that, in the current version of SAFESIMM, the spatial distribution of
  # animal densities is considered constant; it is evaluated only at the
  # beginning of the first simulation run;
  #
  # `DensityData' is a list containing: `Density', `Uncertainty', `Species',
  #                                     `GridLong', `GridLat'
  #
  if(exists("DensityData")==F){
    DensityData<- getAnimalDensity_m1(AOC)
    DensityData<<- DensityData
  }

  # no. of species
  NumberOfSpecies<- ncol(DensityData$Density)
  # no. of simulated animals per species
  NumberOfAnimals<- UserSettings$NumberOfSimAnimals

  SumSpeciesDensity<- vector(length=NumberOfSpecies)
  for(j in 1:NumberOfSpecies){SumSpeciesDensity[j]<- sum(DensityData$Density[,j])}

  # minimum animal density;
  # animal densities below this value are considered negligible
  tolerance<- UserSettings$ZeroTolerance
  WorkingDensityArray<- DensityData$Density[, which(SumSpeciesDensity>tolerance)]
  if(NumberOfSpecies==1){
    WorkingDensityArray<- array(WorkingDensityArray, dim=c(length(WorkingDensityArray),1))
  }
  WorkingDensityUncertainty<- DensityData$Uncertainty[, which(SumSpeciesDensity>tolerance)]
  if(NumberOfSpecies==1){
    WorkingDensityUncertainty<- array(WorkingDensityUncertainty, dim=c(length(WorkingDensityUncertainty),1))
  }
  WorkingSpecies<- DensityData$Species[which(SumSpeciesDensity>tolerance)]
  WorkingNoOfSpecies<- ncol(WorkingDensityArray)


  #--

  #
  # - initial conditions - spatial locations
  #

 
  # observation data;
  # to take into account animals that have fixed (not random) start positions    
  ObserverObject<- UserSettings$ObserverObject

  if(is.null(ObserverObject)==F){
    ObserverProportion<- ObserverObject$ObservedNumber/SumSpeciesDensity
    ObserverPortion<- ifelse(ObserverProportion>1, NumberOfAnimals, round(NumberOfAnimals*ObserverProportion,0))
    ObserverPortion<- ObserverPortion[which(SumSpeciesDensity>tolerance)]
    ObserverObject$ObservationLocation<- ObserverObject$ObservationLocation[which(SumSpeciesDensity>tolerance)]
  }
  # if no observation data ...
  if(is.null(ObserverObject)){
    ObserverPortion<- rep(0, WorkingNoOfSpecies)
  }

  WorkingSumSpeciesD<- vector(length=WorkingNoOfSpecies)
  for(j in 1:WorkingNoOfSpecies){WorkingSumSpeciesD[j]<- sum(WorkingDensityArray[,j])}

  # probability of assigning a simulated animal to a given cell;
  # the allocation of simulated animals is in proportion to estimated
  # animal density per cell
  GridProbs<- array(dim=dim(WorkingDensityArray))
  for(j in 1:WorkingNoOfSpecies){
    GridProbs[,j]<- WorkingDensityArray[,j]/WorkingSumSpeciesD[j]
  }

  #
  # the grid cells within the AOC are identified by a `grid cell reference'
  # (an ID number); they are identified sequentially column by column,
  # starting from the bottom left corner;
  #
  # define a data structure for the grid cell references (the IDs for the
  # grid locations within the AOC)
  GridCellReference<- array(dim=c(NumberOfAnimals, WorkingNoOfSpecies))

  # no. of grid cells within the AOC
  NumberOfGridCells<- nrow(WorkingDensityArray)
    
  # define a data structure for the `initial locations within the cells';
  # (later this data structure will be used for the `initial starting positions')
  XYSimLocation<- array(dim=c(NumberOfAnimals, 2, WorkingNoOfSpecies))
  # define a data structure for the `initial locations of the cells'
  # (the initial coordinates Lat, Long of the centres of the cells within the AOC);
  InitialXYGrid<- array(dim=c(NumberOfAnimals, 2, WorkingNoOfSpecies))

  
  # ... species ...
  for(j in 1:WorkingNoOfSpecies){

    # random locations within each grid cell;
    # random allocation to Lat/Long cell centres;
    # the sample function allocations probabilities of selection to outcomes
    # the required number of random draws are taken based on these probabilities

    # allocate random locations ...
    if(ObserverPortion[j]!=NumberOfAnimals){
      # simulated animals with unknown/random initial starting positions
      AnimalAllocation<- (NumberOfAnimals-ObserverPortion[j])
      # note that the range for the lower and upper limits of the distr. 
      # depends on the grid size step
      XYSimLocation[1:AnimalAllocation,1,j]<- runif(AnimalAllocation, -0.042, 0.042)
      XYSimLocation[1:AnimalAllocation,2,j]<- runif(AnimalAllocation, -0.042, 0.042)
      Sampler<- sample(1:NumberOfGridCells, AnimalAllocation, replace=T, prob=GridProbs[,j])
      InitialXYGrid[1:AnimalAllocation,1,j]<- DensityData$GridLong[Sampler]
      InitialXYGrid[1:AnimalAllocation,2,j]<- DensityData$GridLat[Sampler]
      GridCellReference[1:AnimalAllocation,j]<- Sampler

      # initial starting positions;
      # random location = random location within the cell +
      #                   random location of the cell
      XYSimLocation[,,j]<- XYSimLocation[,,j]+InitialXYGrid[,,j]

    }
        
    # allocate remaining positions based on observations ...        
    if(is.null(ObserverObject)==F){
      empty_loc_index<- which(is.na(XYSimLocation[,1,j]))
      for(i in empty_loc_index){
        Sampler<- sample(1:nrow(ObserverObject$ObservationLocation[[j]]), 1, replace=T)
        # inflate number of `observed' animals to correct proportion of total sim animals
        XYSimLocation[i,,j]<- ObserverObject$ObservationLocation[[j]][Sampler,]
        GridCellReference[empty_loc_index,j]<- Sampler
      }
    }
    
  }
  

  #--
  
  #
  # - initial conditions - Sound Exposure Levels (SELs) and
  #                        maximum values of the Sound Propagation (SP) levels
  #


  NumberOfLoc<- length(DensityData$GridLat)

  # no. of bins (binning procedure)
  NumberOfBins<- UserSettings$SummaryBinNumber

  # the `History' object is associated with the data from the previous
  # simulation runs;

  # `HistoricalSELs' is a data structure containing:
  # a) summary of SELs, b) species, c) locations
  HistoricalSELs<- History$HistSELsObject

  # define a data structure for the summary (histograms) of SELs;
  # SELs values per bins, per species, per location
  SELsSummary<- array(0, dim=c(NumberOfBins, WorkingNoOfSpecies, NumberOfLoc))

  # for the first simulation run ...
  if(is.null(HistoricalSELs)==T){
    InitialSELsObject<- list(SummaryValues=SELsSummary,
                             Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                             Species=WorkingSpecies)
    History$HistSELsObject<<- InitialSELsObject
  }
  # for the following simulation runs ...
  if(is.null(HistoricalSELs)==F){
    for(k in 1:NumberOfLoc){
      loc_index<- which(HistoricalSELs$Locations[,1]==DensityData$GridLong[k] &
                            HistoricalSELs$Locations[,2]==DensityData$GridLat[k]) 
      for(m in 1:WorkingNoOfSpecies){
        species_index<- which(HistoricalSELs$Species==WorkingSpecies[m])
        if(length(loc_index)>0 & length(species_index)>0){
          SELsSummary[, m, k]<- HistoricalSELs$SummaryValues[, species_index, loc_index]
        }
      }
    }
    InitialSELsObject<- list(SummaryValues=SELsSummary,
                             Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                             Species=WorkingSpecies)
  }


  # `HistoricalSPLmaxs' is a data structure containing:
  # a) summary of SPLmaxs, b) species, c) locations
  HistoricalSPLmaxs<- History$HistSPLmaxsObject

  # define a data structure for the summary (histograms) of SPLmaxs;
  # SPLmaxs values per bins, per species, per location
  SPLmaxsSummary<- array(0, dim=c(NumberOfBins, WorkingNoOfSpecies, NumberOfLoc))

  # for the first simulation run ...
  if(is.null(HistoricalSPLmaxs)==T){
    InitialSPLmaxsObject<- list(SummaryValues=SPLmaxsSummary,
                                Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                                Species=WorkingSpecies)
    History$HistSPLmaxsObject<<- InitialSPLmaxsObject
  }
  # for the following simulation runs ...
  if(is.null(HistoricalSPLmaxs)==F){
    for(k in 1:NumberOfLoc){
      loc_index<- which(HistoricalSPLmaxs$Locations[,1]==DensityData$GridLong[k] &
                            HistoricalSPLmaxs$Locations[,2]==DensityData$GridLat[k]) 
      for(m in 1:WorkingNoOfSpecies){
        species_index<- which(HistoricalSPLmaxs$Species==WorkingSpecies[m])
        if(length(loc_index)>0 & length(species_index)>0){
          SPLmaxsSummary[, m, k]<- HistoricalSPLmaxs$SummaryValues[, species_index, loc_index]
        }
      }
    }
    InitialSPLmaxsObject<- list(SummaryValues=SPLmaxsSummary,
                                Locations=cbind(DensityData$GridLong, DensityData$GridLat), 
                                Species=WorkingSpecies)
  }


  #--

  InitData<- list(AOC=AOC,
                  GridCellReference=GridCellReference,
                  InitialXYGrid=InitialXYGrid,
                  XYSimLocation=XYSimLocation,
                  WorkingDensityArray=WorkingDensityArray,
                  WorkingDensityUncertainty=WorkingDensityUncertainty,
                  WorkingSpecies=WorkingSpecies,
                  InitialSELsObject=InitialSELsObject,
                  InitialSPLmaxsObject=InitialSPLmaxsObject)
  
  return(InitData)

}


