
#
# CVS - $Id: ERMCAnimalMover_S3BL.R,v 1.9 2011/11/21 12:30:32 lmilazzo Exp $
#





#
#     -- function to simulate the Movement of the Animals --
#
#
simulateAnimalMovement<- function(InitData, Seeds){

  #--
    
  # RNG - assign a seed from the predetermined seed collection
  # by using RNG kind 1 = Marsaglia-Multicarry;
  rs<- c(1, Seeds[3,])
  rs_i<- as.integer(rs)
  .Random.seed<<- rs_i


  #--

  # simulation time (duration of a simulation run)
  SimulRunTime<- UserSettings$AlgorithmRunTime    # [sec]
  
  # length of the time step
  TimeBlocks<- UserSettings$TimeStepEval          # [sec]
  # simulation time over sucessive runs
  GlobalSimulationTime<- UserSettings$ElapsedSimulationTime

  # 'centre' of the AOC, Long then Lat;
  CentreOfAOC<- c(mean(InitData$AOC[,1]), mean(InitData$AOC[,2]))
  WorkingBathyLat<- (Bathy$RowLatCoords-CentreOfAOC[2])*60*1.852235*1000    # [metres]

  NumberOfSimAnimals<- UserSettings$NumberOfSimAnimals
  SpeciesList<- InitData$WorkingSpecies
  NumberOfSpecies<- length(SpeciesList)

  # initial starting positions [degrees]
  # data format for `InitXY': (no. of animals) x XY x (no. of species);
  InitXY<- InitData$XYSimLocation
  # initial coordinates of the centres of the cells [degrees]
  InitialXYGrid<- InitData$InitialXYGrid

  NumberOfPlatforms<- 1+length(MultiShipStats)

  # ... for the primary platform ...
  DutyCycle<- PrimaryShipStats$DutyCycle
  SonarFrequencies<- PrimaryShipStats$SonarFrequency
  MultiplatSourceStrength<- PrimaryShipStats$SourceStrength
  # define a data structure for the locations of the platforms
  CombinedShipLocations<-array(dim=c(NumberOfPlatforms, 2))
  CombinedShipLocations[1,]<- PrimaryShipStats$ShipLocation
  # define a data structure for the sector depths;
  # `AllSectorDepths' is a data structure containing the depths of
  # the 16 sectors for each of the three columns, per platform;
  AllSectorDepths<- array(dim=c(NumberOfPlatforms, 16, 3))
  AllSectorDepths[1, , ]<- PrimaryShipStats$SectorDepths
  # define a data structure for the relative locations of the platforms
  PlatformsLoc<- array(dim=c(NumberOfPlatforms,2))
  # define a data structure for the Sound Propagation Loss (SP Loss)
  MultiSPLossData<- list()
  MultiSPLossData[[1]]<- PrimaryShipStats$SPLossData
  # define a data structure for the transmissions times of the platforms
  CombinedTransmissionTimes<-array(dim=c(NumberOfPlatforms, 2))
  CombinedTransmissionTimes[1,]<- c(PrimaryShipStats$TransmissionStart, PrimaryShipStats$TransmissionFinish)
  # ... for the others platforms ...
  if(is.null(MultiShipStats)==F){
    for(i in 2:NumberOfPlatforms){
      DutyCycle<- c(DutyCycle, MultiShipStats[[i-1]]$DutyCycle)
      SonarFrequencies<- c(SonarFrequencies, MultiShipStats[[i-1]]$SonarFrequency)
      CombinedShipLocations[i,]<- MultiShipStats[[i-1]]$ShipLocation
      AllSectorDepths[i, , ]<- MultiShipStats[[i-1]]$SectorDepths
      MultiSPLossData[[i]]<- MultiShipStats[[i-1]]$SPLossData
      MultiplatSourceStrength[i]<- MultiShipStats[[i-1]]$SourceStrength
      CombinedTransmissionTimes[i,]<- c(MultiShipStats[[i-1]]$TransmissionStart, MultiShipStats[[i-1]]$TransmissionFinish)
    }
  }

  # determine the locations of the simulated animals relative to the centre
  # of the AOC; 
  # transform absolute XY coordinates [degrees] in relative XY coordinates [metres] 
  for(j in 1:NumberOfSpecies){
    # units of `InitXY': [metres]
    InitXY[,,j]<- evaluateRelativePositions(InitXY[,,j], CentreOfAOC)
  }
  # determine the locations of the platforms relative to the centre of the AOC;
  # transform absolute XY coordinates [degrees] in relative XY coordinates [metres]
  PlatformsLoc<- evaluateRelativePositions(CombinedShipLocations, CentreOfAOC)

  # in the case of multiple simulation runs ...
  # take into account the SELs accumulated from the start of the simulation
  InitialSELsObject<- InitData$InitialSELsObject
  InitSELs<- getInitialValues(InitialSELsObject, SpeciesList, InitialXYGrid, Seeds)
  # take into account the SPLmaxs recorded from the start of the simulation
  InitialSPLmaxsObject<- InitData$InitialSPLmaxsObject
  InitSPLmaxs<- getInitialValues(InitialSPLmaxsObject, SpeciesList, InitialXYGrid, Seeds)

  # define a data structure for the simulated animals that are out of sea
  AnimalOutOfSea<- array(0, dim=c(NumberOfSimAnimals, NumberOfSpecies))
  # define a data structure for the final locations of the simulated animals
  FinalXY<- array(dim=c(NumberOfSimAnimals, 2, NumberOfSpecies))
  # define a data structure for the accumulated Sound Exposure Levels (SELs)
  AccumulatedSELs<- array(dim=c(NumberOfSimAnimals, NumberOfSpecies))
  # define a data structure for the maximum values of the Sound Propagation (SP) levels 
  SPLmaxs<- array(0, dim=c(NumberOfSimAnimals, NumberOfSpecies))

  
  #--
  # ... species ...
  for(j in 1:NumberOfSpecies){

    cat(paste("data processing for species:", SpeciesList[j], "\n"))

    CurrentSpeciesInfo<- BAESpeciesInfo[which(BAESpeciesInfo$Species_ID==SpeciesList[j]),]
    AudiogramAdjustment<- getAudiogramAdjustment(CurrentSpeciesInfo$GuildID, SonarFrequencies)
    HorizontalResponseToSound<- substring(CurrentSpeciesInfo$GuildID, 13, 13)
    VerticalResponseToSound<- substring(CurrentSpeciesInfo$GuildID, 14, 14)
    # GuildDiveCode<- substring(BAESpeciesInfo$GuildID, 4, 6)
    SpeciesJMaxDiveDepth<- CurrentSpeciesInfo$MaxDiveDepth        # [metres]
    SpeciesJMaxSpeed<- CurrentSpeciesInfo$MaxSwimSpeed            # [metres/sec]
    SpeciesJMaxDiveDuration<- CurrentSpeciesInfo$MaxDiveDuration  # [min]
    SpeciesJMinDiveDuration<- CurrentSpeciesInfo$MinDiveDuration  # [min]
    SpeciesJMinSurfaceTime<- CurrentSpeciesInfo$MinSurfaceTime    # [min]
    SpeciesJMaxSurfaceTime<- CurrentSpeciesInfo$MaxSurfaceTime    # [min]
    
    #--
    # ... individuals ...

    for(i in 1:NumberOfSimAnimals){
      # animals in the sea: depths>0; animals out of sea: depths<=0
      # conventions: 1 = animals out of sea; 0 = animals in the sea
      AnimalOutOfSea[i,j]<- ifelse(getDepth(InitXY[i,1,j], InitXY[i,2,j], WorkingBathyLat, CentreOfAOC)<=0, 1, 0)
    }

    # locations of the simulated animals relative to the centre of the AOC;
    SumXYPath<- InitXY[,,j]
	
    for(i in (1:NumberOfSimAnimals)[which(AnimalOutOfSea[,j]!=1)]){
      PreviousSimulRunTimes<- GlobalSimulationTime
      SimulTime<- 0

      # accumulated SEL; initial value
      SEL<- InitSELs[i,j]
      # maximum value of the SP levels; initial value
      SPLmax<- InitSPLmaxs[i,j]

      # Sound Propagation (SP) levels; initial values
      # `IntensityArray' is a matrix: a SP level per time step (dive depth), per platform
      IntensityArray<- getSoundPropLevels(SumXYPath[i,1], SumXYPath[i,2], 0,
                                          PlatformsLoc, MultiSPLossData,
                                          MultiplatSourceStrength, AllSectorDepths)


      #--
      # ... run the simulation over the `simulation time' ...

      # units of `SimulTime': [min];
      # units of `SimulRunTime': [sec];
      while(SimulTime<(SimulRunTime/60)){


        #--
        # ... dive phase (diving and re-emerging) ...
        
        #
        # v = dx/dt => dx = v*dt ;
        # dx[metres] = v[metres/sec.]*dt[min.] = v*dt*60
        #
        # dx = dx_d (diving) + dx_r (re-emerging)
        # dx_d = dx*1/2
        #
      
        LocalDepth<- getDepth(SumXYPath[i,1], SumXYPath[i,2], WorkingBathyLat, CentreOfAOC)  # [metres]
        CurrentDiveDuration<- runif(1, SpeciesJMinDiveDuration, SpeciesJMaxDiveDuration)        # [min]
        SpeedDepthLimit<- SpeciesJMaxSpeed*CurrentDiveDuration/2*60             # [metres]
        CurrentDepth<- min(LocalDepth, SpeedDepthLimit, SpeciesJMaxDiveDepth)   # [metres]

        # discretize the time interval `CurrentDiveDuration';
        # `EvaluationTime' is a series of discrete time steps (during the dive phase)
        EvaluationTime<- seq(0, CurrentDiveDuration, by=(TimeBlocks/60))
        NoOfTimeSteps<- length(EvaluationTime)
        FineScaleTime<- PreviousSimulRunTimes+SimulTime+EvaluationTime
        TotalNoTimeSteps<- length(FineScaleTime)

        # `z' is a vector containing the dive depths;
        # if the dive is deep and long enough, then calculate the depth at each time step
        if(CurrentDepth>UserSettings$ShallowDiveLimit&length(EvaluationTime)>1){
          z<- evaluateDiveDepths(EvaluationTime, CurrentDepth, SpeciesJMaxSpeed, VerticalResponseToSound)
        }
        # if the dive is too shallow, or it is less than the time step, then set to be surface i.e. depth=0
        # note this implies each `dive phase' is a time step at least, even if spent at surface
        if(CurrentDepth<=UserSettings$ShallowDiveLimit | length(EvaluationTime)==1){
          z<- rep(0, NoOfTimeSteps)
        }

        # Sound Propagation (SP) levels
        IntensityArray<- getSoundPropLevels(SumXYPath[i,1], SumXYPath[i,2], z,
                                            PlatformsLoc, MultiSPLossData,
                                            MultiplatSourceStrength, AllSectorDepths)

        # define a data structure for transmitting platforms
        Transmitting<- array(dim=c(TotalNoTimeSteps, NumberOfPlatforms))

        # determine which platform is transmitting over the current time
        for(k in 1:NumberOfPlatforms){
          Transmitting[,k]<- ifelse(FineScaleTime<CombinedTransmissionTimes[k,1] |
                                    FineScaleTime>CombinedTransmissionTimes[k,2],F,T)
        }
        # redefine the data structure `Transmitting';
        # now `Transmitting' stores the SP levels;
        # now `Transmitting' is a matrix: a SP level per time step,
        # per platform (only during diving);
        # the platforms that are not transmitting are marked with `0'
        Transmitting<- ifelse(Transmitting, IntensityArray, 0)


        #--
        # ... surface phase ...

        CurrentSurfaceDuration<- runif(1, SpeciesJMinSurfaceTime, SpeciesJMaxSurfaceTime)  # [min]
        # `CurrentSurfaceDuration' must be greater or at least equal to the time step
        CurrentSurfaceDuration<- ifelse((CurrentSurfaceDuration*60)<=TimeBlocks, TimeBlocks/60, CurrentSurfaceDuration)

        JumpLocation<- getNextMove(Transmitting, SumXYPath[i,], CurrentSurfaceDuration,
                                   HorizontalResponseToSound, CurrentSpeciesInfo,
                                   LocalDepth, PlatformsLoc,
                                   WorkingBathyLat, CentreOfAOC, Seeds)

        HorizontalPathX<- seq(SumXYPath[i,1], JumpLocation[1], length=(CurrentSurfaceDuration/(TimeBlocks/60)))
        HorizontalPathY<- seq(SumXYPath[i,2], JumpLocation[2], length=(CurrentSurfaceDuration/(TimeBlocks/60)))
	    if(CurrentSurfaceDuration == 1){
          HorizontalPathX<- seq(SumXYPath[i,1], JumpLocation[1], length=CurrentSurfaceDuration + 1)
          HorizontalPathY<- seq(SumXYPath[i,2], JumpLocation[2], length=CurrentSurfaceDuration + 1)
        }


        #--
        # ... dive phase & surface phase ...

        # `z' is a vector containing the dive depths;
        # note that `z' refers to both dive phase and surface phase - for the
        # surface movement, dive depth = 0;
        z<- c(z, rep(0, length(HorizontalPathX)))
        # `HorizontalPathX' and `HorizontalPathY' are vectors containing the
        # coordinates of the locations along a path during the surface movement,
        # between two follow-up dive phases;
        HorizontalPathX<- c(rep(SumXYPath[i,1], length(EvaluationTime)), HorizontalPathX)
        HorizontalPathY<- c(rep(SumXYPath[i,2], length(EvaluationTime)), HorizontalPathY)
        
        # Sound Propagation (SP) levels;
        # SP levels are root-mean-square (RMS) values
        IntensityArray<- getSoundPropLevels(HorizontalPathX, HorizontalPathY, z,
                                            PlatformsLoc, MultiSPLossData,
                                            MultiplatSourceStrength, AllSectorDepths)

        # `EvaluationTime' is a series of discrete time steps (dive phase & surface phase)
        EvaluationTime<- seq(0, ceiling(CurrentDiveDuration)+CurrentSurfaceDuration, by=(TimeBlocks/60))
        FineScaleTime<- PreviousSimulRunTimes+SimulTime+EvaluationTime
        TotalNoTimeSteps<- length(FineScaleTime)

        # define a data structure for transmitting platforms
        Transmitting<- array(dim=c(TotalNoTimeSteps, NumberOfPlatforms))

        # determine which platform is transmitting over the current time
        for(k in 1:NumberOfPlatforms){
          Transmitting[,k]<-
        ifelse(FineScaleTime<CombinedTransmissionTimes[k,1] |
               FineScaleTime>CombinedTransmissionTimes[k,2],F,T)
        }

        # redefine the data structure `Transmitting';
        # now `Transmitting' stores the SP levels;
        # now `Transmitting' is a matrix: a SP level per time step,
        # per platform (during diving and surface movement);
        # the platforms that are not transmitting are marked with `NA'
        Transmitting<- ifelse(Transmitting, IntensityArray, NA)

        # calculate/update the value of the accumulated Sound Exposure Level (SEL)
        SEL<- getSEL(SEL, Transmitting, AudiogramAdjustment, DutyCycle)
        # calculate/update the maximum value of the Sound Propagation (SP) levels
        SPLmax = getMaxSPLevel(SPLmax, Transmitting)

        # move horizontally and prepare for next dive
        SumXYPath[i,]<- JumpLocation

	    # truncate to the nearest second
        STime<- trunc((SimulTime+CurrentDiveDuration+CurrentSurfaceDuration)*60)
        SimulTime<- STime/60       # [min]

      }  # ... end of while loop

      # store the accumulated SEL for individual i, species j, after the `SimulRunTime'
      AccumulatedSELs[i,j]<- SEL
      # exclude the invalid simulated animals (animals out of sea);
      # the invalid simulated animals are marked with `NA'
      AccumulatedSELs[which(AnimalOutOfSea[,j]==1),j]<- NA
      # store the maximum value of the SP levels for individual i, species j, after the `SimulRunTime'
      SPLmaxs[i,j]<- SPLmax
      # exclude the invalid simulated animals (animals out of sea);
      # the invalid simulated animals are marked with `NA'
      SPLmaxs[which(AnimalOutOfSea[,j]==1),j]<- NA

      # store the final locations for individual i, species j,
      # after the `SimulRunTime'
      FinalXY[i,,j]<- SumXYPath[i,]
    
    }  # ... end of i loop

  }  # ... end of j loop

  UserSettings$ElapsedSimulationTime<<-GlobalSimulationTime+SimulTime

  # in the case of multiple simulation runs ...
  # store the SELs accumulated
  SELsToStore<- getSummary(AccumulatedSELs, InitialSELsObject, SpeciesList, InitialXYGrid) 
  History$HistSELsObject<<- updateHistory(History$HistSELsObject, SELsToStore)
  # store the maximum values of the SP levels
  SPLmaxsToStore<- getSummary(SPLmaxs, InitialSPLmaxsObject, SpeciesList, InitialXYGrid) 
  History$HistSPLmaxsObject<<- updateHistory(History$HistSPLmaxsObject, SPLmaxsToStore)

  SystemDataObject<- list(SELs=AccumulatedSELs, SPLmaxs=SPLmaxs,
                          FinalXY=FinalXY, AnimalOutOfSea=AnimalOutOfSea)

  return(SystemDataObject)


}


#
#     -- function to simulate the Movement of the Animals --
#
#
simulateAnimalMovement_m1<- function(InitData, Seeds){

  #--
    
  # RNG - assign a seed from the predetermined seed collection
  # by using RNG kind 1 = Marsaglia-Multicarry;
  rs<- c(1, Seeds[3,])
  rs_i<- as.integer(rs)
  .Random.seed<<- rs_i


  #--

  # simulation time (duration of a simulation run)
  SimulRunTime<- UserSettings$AlgorithmRunTime    # [sec]
  	
  # length of the time step
  TimeBlocks<- UserSettings$TimeStepEval          # [sec]
  # simulation time over sucessive runs
  GlobalSimulationTime<- UserSettings$ElapsedSimulationTime

  # 'centre' of the AOC, Long then Lat;
  CentreOfAOC<- c(mean(InitData$AOC[,1]), mean(InitData$AOC[,2]))
  WorkingBathyLat<- (Bathy$RowLatCoords-CentreOfAOC[2])*60*1.852235*1000    # [metres]

  NumberOfSimAnimals<- UserSettings$NumberOfSimAnimals
  SpeciesList<- InitData$WorkingSpecies
  NumberOfSpecies<- length(SpeciesList)

  # initial starting positions [degrees]
  # data format for `InitXY': (no. of animals) x XY x (no. of species);
  InitXY<- InitData$XYSimLocation
  # initial coordinates of the centres of the cells [degrees]
  InitialXYGrid<- InitData$InitialXYGrid

  NumberOfPlatforms<- 1+length(MultiShipStats)

  # ... for the primary platform ...
  DutyCycle<- PrimaryShipStats$DutyCycle
  SonarFrequencies<- PrimaryShipStats$SonarFrequency
  MultiplatSourceStrength<- PrimaryShipStats$SourceStrength
  # define a data structure for the locations of the platforms
  CombinedShipLocations<-array(dim=c(NumberOfPlatforms, 2))
  CombinedShipLocations[1,]<- PrimaryShipStats$ShipLocation
  # define a data structure for the relative locations of the platforms
  PlatformsLoc<- array(dim=c(NumberOfPlatforms,2))
  # define a data structure for the Sound Propagation Loss (SP Loss)
  MultiMWSELData<- list()
  MultiMWSELData[[1]]<- PrimaryShipStats$MWSELData
  # define a data structure for the transmissions times of the platforms
  CombinedTransmissionTimes<-array(dim=c(NumberOfPlatforms, 2))
  CombinedTransmissionTimes[1,]<- c(PrimaryShipStats$TransmissionStart, PrimaryShipStats$TransmissionFinish)
  # ... for the others platforms ...
  if(is.null(MultiShipStats)==F){
    for(i in 2:NumberOfPlatforms){
      DutyCycle<- c(DutyCycle, MultiShipStats[[i-1]]$DutyCycle)
      SonarFrequencies<- c(SonarFrequencies, MultiShipStats[[i-1]]$SonarFrequency)
      CombinedShipLocations[i,]<- MultiShipStats[[i-1]]$ShipLocation
      MultiMWSELData[[i]]<- MultiShipStats[[i-1]]$MWSELData
      MultiplatSourceStrength[i]<- MultiShipStats[[i-1]]$SourceStrength
      CombinedTransmissionTimes[i,]<- c(MultiShipStats[[i-1]]$TransmissionStart, MultiShipStats[[i-1]]$TransmissionFinish)
    }
  }

  # determine the locations of the simulated animals relative to the centre
  # of the AOC; 
  # transform absolute XY coordinates [degrees] in relative XY coordinates [metres] 
  for(j in 1:NumberOfSpecies){
    # units of `InitXY': [metres]
    InitXY[,,j]<- evaluateRelativePositions(InitXY[,,j], CentreOfAOC)
  }
  # determine the locations of the platforms relative to the centre of the AOC;
  # transform absolute XY coordinates [degrees] in relative XY coordinates [metres]
  PlatformsLoc<- evaluateRelativePositions(CombinedShipLocations, CentreOfAOC)

  
  # +++++ DEBUG +++++
  cat(paste("CentreOfAOC: ", CentreOfAOC[1], ", ", CentreOfAOC[2], "\n"))
  cat(paste("PlatformsLoc: ", PlatformsLoc[,1], ", ", PlatformsLoc[,2], "\n"))

		  
		  
  # in the case of multiple simulation runs ...
  # take into account the SELs accumulated from the start of the simulation
  InitialSELsObject<- InitData$InitialSELsObject
  InitSELs<- getInitialValues(InitialSELsObject, SpeciesList, InitialXYGrid, Seeds)

  # define a data structure for the simulated animals that are out of sea
  AnimalOutOfSea<- array(0, dim=c(NumberOfSimAnimals, NumberOfSpecies))
  # define a data structure for the final locations of the simulated animals
  FinalXY<- array(dim=c(NumberOfSimAnimals, 2, NumberOfSpecies))
  # define a data structure for the accumulated Sound Exposure Levels (SELs)
  AccumulatedSELs<- array(dim=c(NumberOfSimAnimals, NumberOfSpecies))
  # define a data structure for the number of sound responses for individual i, species j
  no_sresponses<- array(0, dim=c(NumberOfSimAnimals, NumberOfSpecies))


  
  
  # +++++ DEBUG +++++
  offset_debug<- (SimulRunTime/60)/10
  idxd<- 0
  stimes_debug<- array(0, dim=c(NumberOfSimAnimals*11, NumberOfSpecies)) 
  paths_debug<- array(0, dim=c(NumberOfSimAnimals*11, 2, NumberOfSpecies)) 
  sel_debug<- array(0, dim=c(NumberOfSimAnimals*11, NumberOfSpecies)) 
  
  
  
  
  
  #--
  # ... species ...
  for(j in 1:NumberOfSpecies){

    cat(paste("data processing for species:", SpeciesList[j], "\n"))

    CurrentSpeciesInfo<- BAESpeciesInfo[which(BAESpeciesInfo$Species_ID==SpeciesList[j]),]
    AudiogramAdjustment<- getAudiogramAdjustment(CurrentSpeciesInfo$GuildID, SonarFrequencies)
    HorizontalResponseToSound<- substring(CurrentSpeciesInfo$GuildID, 13, 13)
    VerticalResponseToSound<- substring(CurrentSpeciesInfo$GuildID, 14, 14)
    # GuildDiveCode<- substring(BAESpeciesInfo$GuildID, 4, 6)
    SpeciesJMaxDiveDepth<- CurrentSpeciesInfo$MaxDiveDepth        # [metres]
    SpeciesJMaxSpeed<- CurrentSpeciesInfo$MaxSwimSpeed            # [metres/sec]
    SpeciesJMaxDiveDuration<- CurrentSpeciesInfo$MaxDiveDuration  # [min]
    SpeciesJMinDiveDuration<- CurrentSpeciesInfo$MinDiveDuration  # [min]
    SpeciesJMinSurfaceTime<- CurrentSpeciesInfo$MinSurfaceTime    # [min]
    SpeciesJMaxSurfaceTime<- CurrentSpeciesInfo$MaxSurfaceTime    # [min]
    
    #--
    # ... individuals ...

    for(i in 1:NumberOfSimAnimals){
      # animals in the sea: depths>0; animals out of sea: depths<=0
      # conventions: 1 = animals out of sea; 0 = animals in the sea
      AnimalOutOfSea[i,j]<- ifelse(getDepth(InitXY[i,1,j], InitXY[i,2,j], WorkingBathyLat, CentreOfAOC)<=0, 1, 0)
    }

    # locations of the simulated animals relative to the centre of the AOC;
    SumXYPath<- InitXY[,,j]

    for(i in (1:NumberOfSimAnimals)[which(AnimalOutOfSea[,j]!=1)]){
      PreviousSimulRunTimes<- GlobalSimulationTime
      SimulTime<- 0

      # accumulated SEL; initial value
      SEL<- InitSELs[i,j]

      # Sound Exposure Levels (SELs), M-weighted; initial values
      # `mw_sels' is a 2D array: a SEL per time step (dive depth), per platform
      mw_sels<- getMWeightedSELs(SumXYPath[i,1], SumXYPath[i,2], 0,
                                 PlatformsLoc, MultiMWSELData,
                                 MultiplatSourceStrength)

      #--
      # ... run the simulation over the `simulation time' ...

	  
	  
	  # +++++ DEBUG +++++
      tstep_debug<- 0

  
  
      # units of `SimulTime': [min];
      # units of `SimulRunTime': [sec];
      while(SimulTime<(SimulRunTime/60)){


        #--
        # ... dive phase (diving and re-emerging) ...
        
        #
        # v = dx/dt => dx = v*dt ;
        # dx[metres] = v[metres/sec.]*dt[min.] = v*dt*60
        #
        # dx = dx_d (diving) + dx_r (re-emerging)
        # dx_d = dx*1/2
        #
      
        LocalDepth<- getDepth(SumXYPath[i,1], SumXYPath[i,2], WorkingBathyLat, CentreOfAOC)  # [metres]
        CurrentDiveDuration<- runif(1, SpeciesJMinDiveDuration, SpeciesJMaxDiveDuration)     # [min]
        SpeedDepthLimit<- SpeciesJMaxSpeed*CurrentDiveDuration/2*60             # [metres]
        CurrentDepth<- min(LocalDepth, SpeedDepthLimit, SpeciesJMaxDiveDepth)   # [metres]

        # discretize the time interval `CurrentDiveDuration';
        # `EvaluationTime' is a series of discrete time steps (during the dive phase)
        EvaluationTime<- seq(0, CurrentDiveDuration, by=(TimeBlocks/60))
        NoOfTimeSteps<- length(EvaluationTime)
        FineScaleTime<- PreviousSimulRunTimes+SimulTime+EvaluationTime
        TotalNoTimeSteps<- length(FineScaleTime)

        # `z' is a vector containing the dive depths;
        # if the dive is deep and long enough, then calculate the depth at each time step
        if(CurrentDepth>UserSettings$ShallowDiveLimit&length(EvaluationTime)>1){
          z<- evaluateDiveDepths(EvaluationTime, CurrentDepth, SpeciesJMaxSpeed, VerticalResponseToSound)
        }
        # if the dive is too shallow, or it is less than the time step, then set to be surface i.e. depth=0
        # note this implies each `dive phase' is a time step at least, even if spent at surface
        if(CurrentDepth<=UserSettings$ShallowDiveLimit | length(EvaluationTime)==1){
          z<- rep(0, NoOfTimeSteps)
        }

        # Sound Exposure Levels (SELs), M-weighted
        # `mw_sels' is a 2D array: a SEL per time step (dive depth), per platform
        mw_sels<- getMWeightedSELs(SumXYPath[i,1], SumXYPath[i,2], z,
                                   PlatformsLoc, MultiMWSELData,
                                   MultiplatSourceStrength)
	
        # define a data structure for transmitting platforms
        Transmitting<- array(dim=c(TotalNoTimeSteps, NumberOfPlatforms))

        # determine which platform is transmitting over the current time
        for(k in 1:NumberOfPlatforms){
          Transmitting[,k]<- ifelse(FineScaleTime<CombinedTransmissionTimes[k,1] |
                                    FineScaleTime>CombinedTransmissionTimes[k,2],F,T)
        }
        # redefine the data structure `Transmitting';
        # now `Transmitting' stores the SELs;
        # now `Transmitting' is a 2D array: a SEL per time step,
        # per platform (only during diving);
        # the platforms that are not transmitting are marked with `0'
        Transmitting<- ifelse(Transmitting, mw_sels, 0)


        #--
        # ... surface phase ...

        CurrentSurfaceDuration<- runif(1, SpeciesJMinSurfaceTime, SpeciesJMaxSurfaceTime)  # [min]
        # `CurrentSurfaceDuration' must be greater or at least equal to the time step
        CurrentSurfaceDuration<- ifelse((CurrentSurfaceDuration*60)<=TimeBlocks, TimeBlocks/60, CurrentSurfaceDuration)

        NextMove<- getNextMove_m1(Transmitting, SumXYPath[i,], CurrentSurfaceDuration,
                                  HorizontalResponseToSound, CurrentSpeciesInfo, no_sresponses[i,j],
                                  LocalDepth, PlatformsLoc,
                                  WorkingBathyLat, CentreOfAOC, Seeds)
        JumpLocation<- NextMove$Jump
        no_sresponses[i,j]<- NextMove$no_sresponses

        HorizontalPathX<- seq(SumXYPath[i,1], JumpLocation[1], length=(CurrentSurfaceDuration/(TimeBlocks/60)))
        HorizontalPathY<- seq(SumXYPath[i,2], JumpLocation[2], length=(CurrentSurfaceDuration/(TimeBlocks/60)))
	    if(CurrentSurfaceDuration == 1){
          HorizontalPathX<- seq(SumXYPath[i,1], JumpLocation[1], length=CurrentSurfaceDuration + 1)
          HorizontalPathY<- seq(SumXYPath[i,2], JumpLocation[2], length=CurrentSurfaceDuration + 1)
        }


        #--
        # ... dive phase & surface phase ...

        # `z' is a vector containing the dive depths;
        # note that `z' refers to both dive phase and surface phase - for the
        # surface movement, dive depth = 0;
        z<- c(z, rep(0, length(HorizontalPathX)))
        # `HorizontalPathX' and `HorizontalPathY' are vectors containing the
        # coordinates of the locations along a path during the surface movement,
        # between two follow-up dive phases;
        HorizontalPathX<- c(rep(SumXYPath[i,1], length(EvaluationTime)), HorizontalPathX)
        HorizontalPathY<- c(rep(SumXYPath[i,2], length(EvaluationTime)), HorizontalPathY)

        # Sound Exposure Levels (SELs), M-weighted
        mw_sels<- getMWeightedSELs(HorizontalPathX, HorizontalPathY, z,
                                   PlatformsLoc, MultiMWSELData,
                                   MultiplatSourceStrength)

        # `EvaluationTime' is a series of discrete time steps (dive phase & surface phase)
        EvaluationTime<- seq(0, ceiling(CurrentDiveDuration)+CurrentSurfaceDuration, by=(TimeBlocks/60))
        FineScaleTime<- PreviousSimulRunTimes+SimulTime+EvaluationTime
        TotalNoTimeSteps<- length(FineScaleTime)

        # define a data structure for transmitting platforms
        Transmitting<- array(dim=c(TotalNoTimeSteps, NumberOfPlatforms))

        # determine which platform is transmitting over the current time
        for(k in 1:NumberOfPlatforms){
          Transmitting[,k]<- ifelse(FineScaleTime<CombinedTransmissionTimes[k,1] |
                                    FineScaleTime>CombinedTransmissionTimes[k,2],F,T)
        }

        # redefine the data structure `Transmitting';
        # now `Transmitting' stores the SELs;
        # now `Transmitting' is a 2D array: a SEL per time step,
        # per platform (during diving and surface movement);
        # the platforms that are not transmitting are marked with `NA'
        Transmitting<- ifelse(Transmitting, mw_sels, NA)

        # calculate/update the value of the accumulated Sound Exposure Level (SEL)
        SEL<- getSEL_m1(SEL, Transmitting, DutyCycle)

        # move horizontally and prepare for next dive
        SumXYPath[i,]<- JumpLocation

	    # truncate to the nearest second
        STime<- trunc((SimulTime+CurrentDiveDuration+CurrentSurfaceDuration)*60)
        SimulTime<- STime/60       # [min]

		
		
			  
	    # +++++ DEBUG +++++
        if (SimulTime> tstep_debug) {
          tstep_debug<- tstep_debug + offset_debug
		  cat(paste("tstep_debug: ", tstep_debug, "\n"))
		  if (idxd == NumberOfSimAnimals*11) idxd<- 0
		  idxd<- idxd + 1
		  stimes_debug[idxd,j]<- SimulTime
          #cat(paste("SimulTime: ", SimulTime, "\n"))
          paths_debug[idxd,1,j]<- SumXYPath[i,1]
          paths_debug[idxd,2,j]<- SumXYPath[i,2]
		  #cat(paste(SumXYPath[i,1], ", ", SumXYPath[i,2], "\n"))
		  sel_debug[idxd,j]<- SEL
	      #cat(paste(SEL, "\n"))
		}
        #--

	  
	  
	  
	  
      }  # ... end of while loop

      # store the accumulated SEL for individual i, species j, after the `SimulRunTime'
      AccumulatedSELs[i,j]<- SEL
      # exclude the invalid simulated animals (animals out of sea);
      # the invalid simulated animals are marked with `NA'
      AccumulatedSELs[which(AnimalOutOfSea[,j]==1),j]<- NA

      # store the final locations for individual i, species j,
      # after the `SimulRunTime'
      FinalXY[i,,j]<- SumXYPath[i,]
	
    }  # ... end of i loop

  }  # ... end of j loop
 


  # +++++ DEBUG +++++
  dump('stimes_debug', file=paste('stimes_debug_', RunIdentifier, '_', DateStamp, '.dat', sep=''))
  dump('paths_debug', file=paste('paths_debug_', RunIdentifier, '_', DateStamp, '.dat', sep=''))
  dump('sel_debug', file=paste('sel_debug_', RunIdentifier, '_', DateStamp, '.dat', sep=''))




 
  UserSettings$ElapsedSimulationTime<<-GlobalSimulationTime+SimulTime

  # in the case of multiple simulation runs ...
  # store the SELs accumulated
  SELsToStore<- getSummary(AccumulatedSELs, InitialSELsObject, SpeciesList, InitialXYGrid) 
  History$HistSELsObject<<- updateHistory(History$HistSELsObject, SELsToStore)

  SystemDataObject<- list(SELs=AccumulatedSELs, SPLmaxs=NULL, SoundResponses= no_sresponses,
                          FinalXY=FinalXY, AnimalOutOfSea=AnimalOutOfSea)


  return(SystemDataObject)


}

