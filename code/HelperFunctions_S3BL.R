
#
# CVS - $Id: HelperFunctions_S3BL.R,v 1.16 2011/11/21 13:00:28 lmilazzo Exp $
#





#
#       -- function to evaluate the Area of Calculation (AOC) --
#
# this function draws a bounding box about the source that captures a given radius;
# only requires the location of the source, the desired radius and the radius of the globe;
# uses a modified Haversine formula;
# gives back a matrix of decimal coords (Long/Lat) for the upper-left and lower-right (in rows)
#
getAreaOfCalc<-function(PrimaryShipStats, MultiShipStats, UserSettings){

  Distance<- UserSettings$AOCRadius         # [km]
  EarthRadius<- UserSettings$EarthRadius    # [km]

  ExtremeLongLat<- matrix(nrow=2, ncol=2)   # col 1: lower then upper Longs
                                            # col 2: lower then upper Lats

  # ... multiple platforms ...
  if(is.null(MultiShipStats)==F){
    NumberOfPlatforms<- 1+length(MultiShipStats)
    # define a data structure for the locations of the platforms
    CombinedShipLocations<-array(dim=c(NumberOfPlatforms, 2))
    CombinedShipLocations[1,]<- PrimaryShipStats$ShipLocation
    for(i in 2:NumberOfPlatforms){CombinedShipLocations[i,]<- MultiShipStats[[i-1]]$ShipLocation}
    ExtremeLongLat[,1]<- range(CombinedShipLocations[,1]) # lower and upper Longs
    ExtremeLongLat[,2]<- range(CombinedShipLocations[,2]) # lower and upper Lats
  }
  # ... single platform ...
  if(is.null(MultiShipStats)){
    ExtremeLongLat[,1]<- rep(PrimaryShipStats$ShipLocation[1],2)  # lower and upper Longs are equal
    ExtremeLongLat[,2]<- rep(PrimaryShipStats$ShipLocation[2],2)  # lower and upper Lats are equal
  }

  LatBounds<- c(-Distance/1.85/60/360*2*pi, Distance/1.85/60/360*2*pi)   # [radians]
  RadianLatBounds<- ExtremeLongLat[,2]/360*2*pi+LatBounds                # [radians]
 
  # modified Haversine formula for the Long distance (at both extreme Lats of AOC)
  LongDist<- 2*asin(sin(Distance/(2*EarthRadius))/cos(RadianLatBounds))  # [radians]
  LongDist<- LongDist/(2*pi)*360                                         # [degrees]
  LongBounds<- c(-LongDist[1], LongDist[2])

  # define a data structure for the size of AOC
  RelativeBounds<- matrix(nrow=2, ncol=2)
  # most extreme Longs of AOC
  RelativeBounds[,1]<- ExtremeLongLat[,1]+LongBounds   # [degrees]
  # most extreme Lats of AOC
  RelativeBounds[,2]<- RadianLatBounds/(2*pi)*360      # [degrees]
  AOC<- RelativeBounds

  return(AOC)

}



#
#    -- function to evaluate the Animal Density   --
#
# this function determines the density of the animals for a given AOC;
#
getAnimalDensity<- function(AOC){

  # `AOC' is a matrix; col 1: lower then upper Longs [degrees]
  #                    col 2: lower then upper Lats  [degrees]
 
  # consider a grid that covers the entire globe
  # (grid spacing: 0.5 degrees)
  LongSeq<- seq(359.75, 0.25, by=-0.5)
  LatSeq<- seq(0.25, 179.75, by=0.5)

  UpperLeft<- c(floor(AOC[1,1] * 2) / 2, ceiling(AOC[2,2] * 2) / 2) 
  LowerRight<- c(ceiling(AOC[2,1] * 2) / 2, floor(AOC[1,2] * 2) / 2)
  # convert to 0:360 for the Long and to 0:180 for the Lats
  TempUpperLeft<- UpperLeft+c(180,90)   # increase to capture partial cells in AOC
  TempLowerRight<- LowerRight+c(180,90) # increase to capture partial cells in AOC

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  long_index<- which(LongSeq<=TempLowerRight[1] & LongSeq>=TempUpperLeft[1])
  lat_index<- which(LatSeq>=TempLowerRight[2] & LatSeq<=TempUpperLeft[2])
  # store the coord. (Lat, Long) scanning the grid column by column,
  # starting from the bottom left corner
  Centres<- expand.grid(LatSeq[lat_index], LongSeq[long_index])
  # col 1: Lat; col 2: Long
  names(Centres)<- c('Latitude', 'Longitude')
  Centres[,1]<- Centres[,1]-90
  Centres[,2]<- Centres[,2]-180

  # no. of grid cells within the AOC
  NumberOfGridCells<- nrow(Centres)

  DensityFiles<- list.files(AbundanceDataPath)
  SpeciesCodes<- sapply(DensityFiles, function(q){substring(q, 1, 5)}, USE.NAMES=F)
  # no. of species
  NumberOfSpecies<- length(DensityFiles)

  # define a data structure for animal density and uncertainty;
  # animal density and uncertainty, for each grid location and species;
  #
  # the animal densities are stored in the same sequence as the coordinates
  # of the grid locations - scanning the grid column by column, starting
  # from the bottom left corner
  # 
  # (note that, at this stage, uncertainty data are not available)
  #
  DensityArray<- array(NA, dim=c(NumberOfGridCells, NumberOfSpecies, 2))

  # evaluate AOC cell abundances ...
  for(i in 1:NumberOfSpecies){

    data_from_file<- read.csv(paste(AbundanceDataPath, DensityFiles[i], sep=''), header=T)
    cat(paste("species:", SpeciesCodes[i], "\n"))

    for(j in 1:NumberOfGridCells){
      species_index<- which(round(data_from_file$Longitude, digits=3)==round(Centres$Longitude[j], digits=3) &
                            round(data_from_file$Latitude, digits=3)==round(Centres$Latitude[j], digits=3))
      if(length(species_index)>=1){
        DensityArray[j,i,1]<- data_from_file$Abundance[species_index]
        # (not used until uncertainty data available)
        # DensityArray[j,i,2]<- data_from_file$Uncertainty[species_index]
      }
    }
    DensityArray[,i,]<- ifelse(is.na(DensityArray[,i,]),0,DensityArray[,i,])

  }

  # condense the data down to only the predicted to be present species
  if(NumberOfSpecies!=1){
    wspecies_index<- which(apply(DensityArray[,,1], 2, sum)!=0)
  }
  else{
    wspecies_index<- NumberOfSpecies
  }
  WorkingDensity<- as.data.frame(DensityArray[, wspecies_index,1])
  names(WorkingDensity)<- SpeciesCodes[wspecies_index]
  WorkingUncertainty<- as.data.frame(DensityArray[, wspecies_index,2])
  names(WorkingUncertainty)<- SpeciesCodes[wspecies_index]

  ADensityData<- list(Density=WorkingDensity, Uncertainty=WorkingUncertainty,
                      Species=names(WorkingDensity), GridLong=Centres[,2], GridLat=Centres[,1])

  return(ADensityData)

}



#
#    -- function to evaluate the Animal Density   --
#
# this function determines the density of the animals for a given AOC;
#
getAnimalDensity_m1<- function(AOC){

  # `AOC' is a matrix; col 1: lower then upper Longs [degrees]
  #                    col 2: lower then upper Lats  [degrees]
 
  # consider a grid that covers the entire globe ...

  # (grid spacing: 0.5 degrees)
  #LongSeq<- seq(359.75, 0.25, by=-0.5)
  #LatSeq<- seq(0.25, 179.75, by=0.5)
 
  # (grid spacing: 0.1 degrees)
  #LongSeq<- seq(359.95, 0.05, by=-0.1)
  #LatSeq<- seq(0.05, 179.95, by=0.1)

  # (grid spacing: 0.0833 degrees)
  LongSeq<- seq(360.001, 0.062, by=-0.0833)
  LatSeq<- seq(0.012, 179.94, by=0.0833)

  UpperLeft<- c(floor(AOC[1,1] * 2) / 2, ceiling(AOC[2,2] * 2) / 2) 
  LowerRight<- c(ceiling(AOC[2,1] * 2) / 2, floor(AOC[1,2] * 2) / 2)
  # convert to 0:360 for the Long and to 0:180 for the Lats
  TempUpperLeft<- UpperLeft+c(180,90)   # increase to capture partial cells in AOC
  TempLowerRight<- LowerRight+c(180,90) # increase to capture partial cells in AOC

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  long_index<- which(LongSeq<=TempLowerRight[1] & LongSeq>=TempUpperLeft[1])
  lat_index<- which(LatSeq>=TempLowerRight[2] & LatSeq<=TempUpperLeft[2])
  # store the coord. (Lat, Long) scanning the grid column by column,
  # starting from the bottom left corner
  Centres<- expand.grid(LatSeq[lat_index], LongSeq[long_index])
  # col 1: Lat; col 2: Long
  names(Centres)<- c('Latitude', 'Longitude')
  Centres[,1]<- Centres[,1]-90
  Centres[,2]<- Centres[,2]-180
  # no. of grid cells within the AOC
  NumberOfGridCells<- nrow(Centres)

  # memory management
  rm(LongSeq)
  rm(LatSeq)
  gc()

  DensityFiles<- list.files(AbundanceDataPath)
  SpeciesCodes<- sapply(DensityFiles, function(q){substring(q, 1, 5)}, USE.NAMES=F)
  # no. of species
  NumberOfSpecies<- length(DensityFiles)

  # define a data structure for animal density and uncertainty;
  # animal density and uncertainty, for each grid location and species;
  #
  # the animal densities are stored in the same sequence as the coordinates
  # of the grid locations - scanning the grid column by column, starting
  # from the bottom left corner
  # 
  # (note that, at this stage, uncertainty data are not available)
  #
  DensityArray<- array(NA, dim=c(NumberOfGridCells, NumberOfSpecies, 2))

  # evaluate AOC cell abundances ...
  for(i in 1:NumberOfSpecies){

    data_from_file<- read.csv(paste(AbundanceDataPath, DensityFiles[i], sep=''),
                              header=T)
    cat(paste("species:", SpeciesCodes[i], "\n"))

    for(j in 1:NumberOfGridCells){
      species_index<- which(round(data_from_file$Longitude, digits=2)==round(Centres$Longitude[j], digits=2) &
                            round(data_from_file$Latitude, digits=2)==round(Centres$Latitude[j], digits=2))
      if(length(species_index)>=1){
        DensityArray[j,i,1]<- data_from_file$Abundance[species_index]
        # (not used until uncertainty data available)
        # DensityArray[j,i,2]<- data_from_file$Uncertainty[species_index]
      }
    }
    DensityArray[,i,]<- ifelse(is.na(DensityArray[,i,]),0,DensityArray[,i,])

  }

  # condense the data down to only the predicted to be present species
  if(NumberOfSpecies!=1){
    wspecies_index<- which(apply(DensityArray[,,1], 2, sum)!=0)
  }
  else{
    wspecies_index<- NumberOfSpecies
  }
  WorkingDensity<- as.data.frame(DensityArray[, wspecies_index,1])
  names(WorkingDensity)<- SpeciesCodes[wspecies_index]
  WorkingUncertainty<- as.data.frame(DensityArray[, wspecies_index,2])
  names(WorkingUncertainty)<- SpeciesCodes[wspecies_index]

  ADensityData<- list(Density=WorkingDensity, Uncertainty=WorkingUncertainty,
                      Species=names(WorkingDensity), GridLong=Centres[,2], GridLat=Centres[,1])

  return(ADensityData)

}



#
#    -- function to evaluate the Relative Positions respect the centre of the AOC  --
#
# this function determines the locations (of simulated animals, platforms, ...)
# relative to the 'source'; the 'source' is the centre of the AOC (in the case
# of single platform, the source is the platform);
# transform absolute XY coordinates [degrees] in relative XY coordinates [metres];
#
evaluateRelativePositions<- function(XYLocations, ReferenceLocation){

  # units of `XYLocations' and `ReferenceLocation': # [degrees]

  Lat<- XYLocations[,2]

  Long<- (XYLocations[,1]-ReferenceLocation[1])   # degrees from source's Longitude line
  # apply the great circle distance conversion
  XYLocations[,1]<- sign(Long)*acos(sin(Lat/360*2*pi)*sin(Lat/360*2*pi)+cos(Lat/360*2*pi)^2*cos(abs(Long)/360*2*pi))
  XYLocations[,1]<- trunc(XYLocations[,1]*360/(2*pi)*60*1.852235*1000)     # [metres]

  Lat<- (Lat-ReferenceLocation[2])                # degrees from source's Latitude line
  XYLocations[,2]<- trunc(Lat*60*1.852235*1000)                            # [metres]

  # units of `XYLocations': # [metres]
  return(XYLocations)

}



#
#            -- function to evaluate the Depth --
#
getDepth<- function(RelativeX, RelativeY, BathyLat, AOCReference){

  # units of `RelativeX', `RelativeY', and `BathyLat': [metres]
  # units of `AOCReference': [degrees]

  Difference<- abs(BathyLat-RelativeY)
  nearestY_index<- which(Difference==min(Difference))
  Lat<- Bathy$RowLatCoords[nearestY_index]

  TempLong<- Bathy$ColLongCoords-AOCReference[1]
  TempLong<- sign(TempLong)*acos(sin(Lat/360*2*pi)*sin(Lat/360*2*pi)+cos(Lat/360*2*pi)^2*cos(abs(TempLong)/360*2*pi))
  TempLong<- TempLong*360/(2*pi)*60*1.852235*1000   # [metres]
  Difference<- abs(TempLong-RelativeX)
  nearestX_index<- which(Difference==min(Difference))

  return(Bathy$DepthMatrix[nearestY_index, nearestX_index])     # [metres]

}



#
#       -- function to evaluate the Dive Depths --
#
# this function determines the dive depths at each time step;
#
evaluateDiveDepths<- function(Time, DepthLimit, Speed, ResponseToSound){
  # note that `ResponseToSound' is not used;
  # at this stage, the vertical response is zero;

  # define a vector `z' for the dive depths;
  DepthVector<- rep(DepthLimit, length(Time))

  Temp<- Time*Speed*UserSettings$TimeStepEval
  Temp<- Temp[which(Temp<=DepthLimit)]
  DepthVector[1:length(Temp)]<- Temp
  DepthVector[length(Time):(length(Time)-length(Temp)+1)]<- Temp

  return(DepthVector)

}



#
#      -- function to evaluate the Animal Movement (next location) --
#
# this function determines the next location for the simulated animal to move to;
#
# movements are: 
# - random bivariate U(-1,1) in the XY plane - scaled for animal speed;
# - random yet directed e.g. towards or away source, movement based on a wrapped;
# normal distribution with parameters based on source location and recieved sound;
#
#
getNextMove<- function(SoundIntensities, StartXY, SurfTime,
                       ResponseToSound, SpeciesInfo,
                       PresentDepth, PlatformsLoc,
                       RelativeBathyLats, AOCCentre, Seeds){

  #--
    
  # RNG - assign a seed from the predetermined seed collection
  # by using RNG kind 1 = Marsaglia-Multicarry;
  rs<- c(1, Seeds[7,])
  rs_i<- as.integer(rs)
  .Random.seed<<- rs_i


  #--


  # data structure for all locations of the simulated
  # animals relative to the platforms
  PositionToAllShips<- t(StartXY-t(PlatformsLoc))

  # define a data structure for the directions
  MagnitudeVectors<- array(dim=dim(PositionToAllShips))
  
  SoundIntensities<- ifelse(is.na(SoundIntensities), 0, SoundIntensities)
  # get loudest sound for each platform
  SoundIntensities<- apply(SoundIntensities, 2, max)

  #
  # horizontal movement in response to sound (HM)
  #
  # no known response to sound source:        HM = 0
  # move horizontally away from sound source: HM = 1
  # move horizontally towards sound source:   HM = 2
  #

  # in the case of no sound, reset `ResponseToSound' to zero
  if(sum(SoundIntensities)==0){ResponseToSound<- '0'}

  # move randomly
  if(ResponseToSound=='0'){
    MoveDirection<- runif(1, 0, (2*pi))    # [radians]
  }
  # move horizontally away or towards
  if(ResponseToSound=='1' | ResponseToSound=='2'){
    DirectionToShips<- atan2(PositionToAllShips[,1], PositionToAllShips[,2])
    NormalisedPositionToShips<-cbind(sin(DirectionToShips), cos(DirectionToShips))    
    MagnitudeVectors[,1]<- SoundIntensities*NormalisedPositionToShips[,1]
    MagnitudeVectors[,2]<- SoundIntensities*NormalisedPositionToShips[,2]
    NetVectors<- apply(MagnitudeVectors, 2, sum)
    # angle from source to simulated animal
    SourceDirection<- atan2(NetVectors[1], NetVectors[2])
    ReactionIntensity<- UserSettings$ReactionIntensity
    # note that the angle is from source to simulated animal,
    # hence is a default away movement
    MoveDirection<- rwrpnorm(1, SourceDirection, sd=ReactionIntensity)    # [radians]
    # adjust for away or towards movement
    MoveDirection<- ifelse(ResponseToSound=='2', (MoveDirection - pi), MoveDirection)
  }
  
  # note that the signature of the function `atan2()' is: `atan2(y,x)';
  # here we call the function `atan2(x,y)' to follow the convention adopted in
  # the Sound Propagation Loss (SP Loss) model;


  # distance travelled
  MoveDistance<- runif(1, 0, SpeciesInfo$MaxSwimSpeed)*60*SurfTime

  # next location (potential, to be assessed)
  Temp<- StartXY+MoveDistance*c(sin(MoveDirection), cos(MoveDirection))
  # depth associated with the next location
  TempDepth<- getDepth(Temp[1], Temp[2], RelativeBathyLats, AOCCentre)

  # animals in the sea: depths>0; animals out of sea: depths<=0
  # 1 = in the sea; 0 = out of sea
  NotToLand<- TempDepth>0

  # current location
  # 1 = in favoured habitat; 0 = outside the favoured habitat
  InHabitat<- PresentDepth>=abs(SpeciesInfo$HabitatMinDepth)&PresentDepth<=abs(SpeciesInfo$HabitatMaxDepth)

  # next location (potential, to be assessed)
  # 1 = outside the favoured habitat; 0 = in favoured habitat
  OutOfHabitat<- TempDepth<abs(SpeciesInfo$HabitatMinDepth)|TempDepth>abs(SpeciesInfo$HabitatMaxDepth)

  # if the current location is in favoured habitat and
  # the animals are in the sea ...
  if(InHabitat&NotToLand){
    # if moving into favoured habitat, accept the next location
    if(OutOfHabitat==F){Jump<- Temp}
    else{
      # if moving into a not favoured habitat, accept the next location
      # with a given probability
      if(OutOfHabitat&(runif(1)>UserSettings$HabitatAffinity)){Jump<- Temp}
      else{Jump<- StartXY}
    }
  }
  # otherwise, when the current location is not in favoured habitat
  # and the animals are out of sea ...
  else{
    # if moving into the sea, accept the next location
    if(NotToLand){Jump<- Temp}
    else{Jump<- StartXY}
  }
  
  return(Jump)

}



#
#      -- function to evaluate the Animal Movement (next location) --
#
# this function determines the next location for the simulated animal to move to;
#
# movements are: 
# - random bivariate U(-1,1) in the XY plane - scaled for animal speed;
# - random yet directed e.g. towards or away source, movement based on a wrapped;
# normal distribution with parameters based on source location and recieved sound;
#
#
getNextMove_m1<- function(mw_sels, StartXY, SurfTime,
                          ResponseToSound, SpeciesInfo, no_sresponses,
                          PresentDepth, PlatformsLoc,
                          RelativeBathyLats, AOCCentre, Seeds){

  #--
    
  # RNG - assign a seed from the predetermined seed collection
  # by using RNG kind 1 = Marsaglia-Multicarry;
  #rs<- c(1, Seeds[7,])
  #rs_i<- as.integer(rs)
  #.Random.seed<<- rs_i


  #--


  # data structure for all locations of the simulated
  # animals relative to the platforms
  PositionToAllShips<- t(StartXY-t(PlatformsLoc))

  # define a data structure for the directions
  MagnitudeVectors<- array(dim=dim(PositionToAllShips))
  
  mw_sels<- ifelse(is.na(mw_sels), 0, mw_sels)
  # max. value of M-weighted SEL for each platform
  mw_sel_max<- apply(mw_sels, 2, max)

  #
  # horizontal movement in response to sound (HM)
  #
  # no known response to sound source:        HM = 0
  # move horizontally away from sound source: HM = 1
  # move horizontally towards sound source:   HM = 2
  #

  # in the case of no sound, reset `ResponseToSound' to zero
  if(sum(mw_sel_max)==0){ResponseToSound<- '0'}

  # if the animal doesn't respond to sound ...
  if(ResponseToSound=='0'){
    # move randomly ...
    MoveDirection<- runif(1, 0, (2*pi))    # [radians]
  }
  # if the animal responds to sound ...
  if(ResponseToSound=='1' | ResponseToSound=='2'){

    # probability of movement (horizontally away or towards)
    prob_mov<- getMovProb_SEL(mw_sel_max)

    if (runif(1, 0, 1)<max(prob_mov)){
      no_sresponses<- no_sresponses + 1

      # move horizontally away or towards ...
      DirectionToShips<- atan2(PositionToAllShips[,1], PositionToAllShips[,2])
      NormalisedPositionToShips<-cbind(sin(DirectionToShips), cos(DirectionToShips))    
      MagnitudeVectors[,1]<- mw_sel_max*NormalisedPositionToShips[,1]
      MagnitudeVectors[,2]<- mw_sel_max*NormalisedPositionToShips[,2]
      NetVectors<- apply(MagnitudeVectors, 2, sum)
      # angle from source to simulated animal
      SourceDirection<- atan2(NetVectors[1], NetVectors[2])
      ReactionIntensity<- UserSettings$ReactionIntensity
      # note that the angle is from source to simulated animal,
      # hence is a default away movement
      MoveDirection<- rwrpnorm(1, SourceDirection, sd=ReactionIntensity)    # [radians]
      # adjust for away or towards movement
      MoveDirection<- ifelse(ResponseToSound=='2', (MoveDirection - pi), MoveDirection)
    }
    else{
      # move randomly ...
      MoveDirection<- runif(1, 0, (2*pi))    # [radians]
    }
  
  }

  # distance travelled
  MoveDistance<- runif(1, 0, SpeciesInfo$MaxSwimSpeed)*60*SurfTime

  # next location (potential, to be assessed)
  Temp<- StartXY+MoveDistance*c(sin(MoveDirection), cos(MoveDirection))
  # depth associated with the next location
  TempDepth<- getDepth(Temp[1], Temp[2], RelativeBathyLats, AOCCentre)

  # animals in the sea: depths>0; animals out of sea: depths<=0
  # 1 = in the sea; 0 = out of sea
  NotToLand<- TempDepth>0

  # current location
  # 1 = in favoured habitat; 0 = outside the favoured habitat
  InHabitat<- PresentDepth>=abs(SpeciesInfo$HabitatMinDepth)&PresentDepth<=abs(SpeciesInfo$HabitatMaxDepth)

  # next location (potential, to be assessed)
  # 1 = outside the favoured habitat; 0 = in favoured habitat
  OutOfHabitat<- TempDepth<abs(SpeciesInfo$HabitatMinDepth)|TempDepth>abs(SpeciesInfo$HabitatMaxDepth)

  # if the current location is in favoured habitat and
  # the animals are in the sea ...
  if(InHabitat&NotToLand){
    # if moving into favoured habitat, accept the next location
    if(OutOfHabitat==F){Jump<- Temp}
    else{
      # if moving into a not favoured habitat, accept the next location
      # with a given probability
      if(OutOfHabitat&(runif(1)>UserSettings$HabitatAffinity)){Jump<- Temp}
      else{Jump<- StartXY}
    }
  }
  # otherwise, when the current location is not in favoured habitat
  # and the animals are out of sea ...
  else{
    # if moving into the sea, accept the next location
    if(NotToLand){Jump<- Temp}
    else{Jump<- StartXY}
  }

  NextMove<- list(Jump=Jump, no_sresponses=no_sresponses)
  
  return(NextMove)

}



#
#    --    function to evaluate the Probability of Movement       --
#    --       for an animal of a given species (from SEL)         --
#
# this function calculates the probability of movement, taking into account the SEL;
#
# the equation for the dose-response curve is the `logistic function';
#
getMovProb_SEL<- function(SEL){

  #
  # dose-response curve
  #   dose     = sound exposure level
  #   response = probability of movement
  #
  # logistic function
  FunctionEval<- exp(DRParamsMov[1]+SEL*DRParamsMov[2])
  MovP_SEL<- FunctionEval/(1+FunctionEval)
  
  return(MovP_SEL)

}



#
#    -- function to evaluate the Sound Propagation (SP) levels --
# 
getSoundPropLevels<- function(InputX, InputY, z,
                              PlatformsLoc, SPLossData,
                              MultiplatSourceStrength, SectorDepth){

  #
  # units of `InputX', `InputY', `z': [metres]
  # units of `PlatformsLoc': [metres]
  #
  # `InputX', `InputY', `z' refer to the locations of the simulated
  # animal relative to the centre of the AOC;
  # `InputX' and ``InputY' are vectors containing the coordinates
  # of the locations along a path during the surface movement,
  # between two follow-up dive phases;
  # `z' is a vector containing the dive depths;
  # note that `z' refers to both dive phase and surface phase - for
  # the surface movement, dive depth = 0;
  # 
  # `PlatformsLoc' refers to the locations of the platforms
  # relative to the centre of the AOC;
  #
  # the coordinate system is centred on the 'source'; the 'source' is the centre
  # of the AOC (in the case of single platform, the source is the platform);
  #

  #
  # within the Sound Propagation (SP) model, a given region of sea is
  # divided into sectors:
  #   3 cylinders are associated with a given sound source;
  #   50 distance ranges are defined for each cylinder;
  #   16 slides are selected for each cylinder;
  #   20 depth levels are selected for each slice, in each cylinder;
  #

  #
  # `SectorDepth' is a data structure containing the depths of
  # the 16 sectors, per column, per platform;
  #

  NoOfTimeSteps<- length(z)
  NoOfPlatforms<- nrow(PlatformsLoc)

  # define a data structure for Sound Propagation (SP) levels;
  # `SoundPropL' is a matrix;
  # the rows are associated with the time steps (dive depths)
  # the cols are associated with the platforms
  SoundPropL<- array(dim=c(NoOfTimeSteps, NoOfPlatforms))

  # ... platforms ...
  for(j in 1:NoOfPlatforms){

    # locations and distances of the simulated animal
    # respect a given platform
    x<- InputX-PlatformsLoc[j,1]; y<- InputY-PlatformsLoc[j,2]
    if(length(z) != length(x)){x<- rep(x, length(z)); y<- rep(y, length(z))}
    Temp<- data.frame(Distance=sqrt((x)^2+(y)^2), xdist=x, ydist=y)

    # ... dive phase & surface phase  ...
    for(i in 1:NoOfTimeSteps){

      # determine the slice
      # (note that there are 16 slices, defined in clockwise order
      # from the North-South direction)
      ObsSlice<- atan2(Temp$xdist[i], Temp$ydist[i])
  
      # note that the signature of the function `atan2()' is: `atan2(y,x)';
      # here we call the function `atan2(x,y)' to follow the convention adopted in
      # the Sound Propagation Loss (SP Loss) model;

      Slice<- round(ObsSlice/(pi/8))
      Slice<- ifelse(Slice>=0, Slice+1, 16-abs(Slice)+1)

      # first of the spatial resolutions 10km - 100km (cylinder #3)
      if(Temp$Distance[i]>=10000 & Temp$Distance[i]<100000){
        Range<- floor(Temp$Distance[i]/2000)+1 # fifty range steps = 2000m (round towards source)
        # SP sector max != bathy data and fudge for when dive hits seafloor which gives div by 0
        z[i]<- ifelse(z[i]>=SectorDepth[j,Slice,3], SectorDepth[j,Slice,]*0.99, z[i])
        Depth<- floor(z[i]/(SectorDepth[j,Slice,3]/20))+1 # depth is over 20 steps - assume +ve depth recordings
        PropLoss<- SPLossData[[j]][Range, Depth, Slice, 3]
      }
      # second of the spatial resolutions 1km - 10km (cylinder #2)
      if(Temp$Distance[i]>=1000 & Temp$Distance[i]<10000){
        Range<- floor(Temp$Distance[i]/200)+1 # fifty range steps = 200m (round towards source)
        # SP sector max != bathy data and fudge for when dive hits seafloor which gives div by 0
        z[i]<- ifelse(z[i]>=SectorDepth[j,Slice,2], SectorDepth[j,Slice,]*0.99, z[i])
        Depth<- floor(z[i]/(SectorDepth[j,Slice,2]/20))+1 # depth is over 20 steps - assume +ve depth recordings
        PropLoss<- SPLossData[[j]][Range, Depth, Slice, 2]
      }
      # third of the spatial resolutions 0m - 1000m (cylinder #1)
      if(Temp$Distance[i]<1000){
        Range<- floor(Temp$Distance[i]/20)+1 # fifty range steps = 20m (round towards source)
        # SP sector max != bathy data and fudge for when dive hits seafloor which gives div by 0
        z[i]<- ifelse(z[i]>=SectorDepth[j,Slice,1], SectorDepth[j,Slice,]*0.99, z[i])
        Depth<- floor(z[i]/(SectorDepth[j,Slice,1]/20))+1 # depth is over 20 steps - assume +ve depth recordings
        PropLoss<- SPLossData[[j]][Range, Depth, Slice, 1]
      }
      if(Temp$Distance[i]>=100000){
        PropLoss<- MultiplatSourceStrength[j]
      }
      PropLoss<- ifelse(PropLoss<0, MultiplatSourceStrength[j], PropLoss)
      SoundPropL[i,j]<- ifelse((MultiplatSourceStrength[j]-PropLoss)<=0,0,(MultiplatSourceStrength[j]-PropLoss))

    }  # end i loop

  }  # end j loop

  return(SoundPropL)

}


#
#    -- function to evaluate the Sound Exposure Levels (SELs), M-weighted --
# 
getMWeightedSELs<- function(InputX, InputY, z,
                            PlatformsLoc, mw_sel_data,
                            MultiplatSourceStrength){

  #
  # units of `InputX', `InputY', `z': [metres]
  # units of `PlatformsLoc': [metres]
  #
  # `InputX', `InputY', `z' refer to the location of a simulated
  # animal relative to the centre of the AOC;
  # `InputX' and ``InputY' are vectors containing the coordinates
  # of the locations along a path during the surface movement,
  # between two follow-up dive phases;
  # `z' is a vector containing the dive depths;
  # note that `z' refers to both dive phase and surface phase - for
  # the surface movement, dive depth = 0;
  # 
  # `PlatformsLoc' refers to the locations of the platforms
  # relative to the centre of the AOC;
  #
  # the coordinate system is centred on the 'source'; the 'source' is the centre
  # of the AOC (in the case of single platform, the source is the platform);
  #

  #
  # within the Sound Propagation (SP) model, a given region of sea is
  # divided into sectors:
  #   - 96 slides, defined in clockwise order from the North-South direction;
  #   - `N' distance ranges for each slide; each distance range is equal to 100 [m];
  #

  NoOfTimeSteps<- length(z)
  NoOfPlatforms<- nrow(PlatformsLoc)

  # define a data structure for the Sound Exposure Levels (SELs), M-weighted;
  # `mw_sel' is a 2D array;
  # the rows are associated with the time steps;
  # the cols are associated with the platforms;
  mw_sels<- array(dim=c(NoOfTimeSteps, NoOfPlatforms))

  # ... platforms ...
  for(j in 1:NoOfPlatforms){

    # normalizing the column #1 (distances)
    mw_sel_data[[j]][,1]<- mw_sel_data[[j]][,1]/100
	
    ndist_max<- max(mw_sel_data[[j]][,1])

    # locations and distances of the simulated animal
    # respect a given platform
    loc_x<- InputX-PlatformsLoc[j,1]
    loc_y<- InputY-PlatformsLoc[j,2]
    if(length(z) != length(loc_x)){
      loc_x<- rep(loc_x, length(z))
      loc_y<- rep(loc_y, length(z))
    }
    saloc<- data.frame(saloc_dist=sqrt((loc_x)^2+(loc_y)^2), saloc_x=loc_x, saloc_y=loc_y)

    # ... dive phase & surface phase  ...
    for(i in 1:NoOfTimeSteps){

      # determining the slice ID
      # (note that there are 96 slices, defined in clockwise order
      # from the North-South direction)
      saloc_angle<- atan2(saloc$saloc_x[i], saloc$saloc_y[i])
      slice_id<- round(saloc_angle/(pi/48))
      slice_id<- ifelse(slice_id>=0, slice_id+1, 96-abs(slice_id)+1)
	  # determining the normalized distance ID
      ndist<- round(saloc$saloc_dist[i]/100)

      if(ndist <= ndist_max){
        ndist_id<- which(mw_sel_data[[j]][,1]==ndist)
        mw_sels[i,j]<- mw_sel_data[[j]][ndist_id, (slice_id+1)]
      }
      else{
        mw_sels[i,j]<- 0
      }

    }  # end i loop

  }  # end j loop
  
  
  return(mw_sels)

}


#
#    -- function to evaluate/update the Sound Exposure Level --
#
getSEL<- function(SELcurrent, IntensityArray, AudiogramAdjustment, DutyCycle){

  # `IntensityArray' contains the SP levels;
  # `IntensityArray' is a matrix:
  #    - rows refer to the time steps (during diving and surface movement)
  #    - cols refer to the platforms (the platforms that are not
  #      transmitting are marked with `NA')

  Adjustment<- IntensityArray-matrix(rep(AudiogramAdjustment, nrow(IntensityArray)), ncol=ncol(IntensityArray), byrow=T)
  # mark with `NA' any zero (or negative) values for exclusion
  Adjustment<- ifelse(Adjustment<=0, NA, Adjustment)

  # apply the equation, accumulating current SEL with latest contribution (if any)
  # note the term within the log10 function is at least 1 as currentSEL>=0
  TempSum<- apply(Adjustment, 2, function(q){sum(10^((q)/10), na.rm=T)})

  # Sound Exposure Level
  SELvalue<- 10*log10(sum(DutyCycle*UserSettings$TimeStepEval*TempSum) + 10^(SELcurrent/10)) 

  return(SELvalue)

}



#
#    -- function to evaluate the accumulated Sound Exposure Level --
#
getSEL_m1<- function(mw_sel_old, mw_sels, DutyCycle){

  # `mw_sels' contains the SELs;
  # `mw_sels' is a matrix:
  #    - rows refer to the time steps (during diving and surface movement)
  #    - cols refer to the platforms (the platforms that are not
  #      transmitting are marked with `NA')

  # mark with `NA' any zero (or negative) values for exclusion
  mw_sels<- ifelse(mw_sels<=0, NA, mw_sels)

  # apply the equation, accumulating current SEL with latest contribution (if any)
  # note the term within the log10 function is at least 1 as mw_sel_old>=0
  sum_tmp<- apply(mw_sels, 2, function(q){sum(10^((q)/10), na.rm=T)})

  # updated value of the accumulated Sound Exposure Level
  mw_sel_new<- 10*log10(sum(DutyCycle*UserSettings$TimeStepEval*sum_tmp) + 10^(mw_sel_old/10)) 

  return(mw_sel_new)

}


#
#    -- function to evaluate the Adjustment to the Sound Intensity   --
#
# this function determines the appropriate adjustment for each species
# sensitivity to the frequency used;
#
getAudiogramAdjustment<- function(SpeciesInfo, SonarFreq){

  AudiogramGuildCode<- substring(SpeciesInfo, 2, 3)
  AudiogramGuildCode<- paste("Audiogram", AudiogramGuildCode, sep='')
  audiogram_index<- which(names(BAEAudiograms)==AudiogramGuildCode)
  cat(paste("audiogram:", AudiogramGuildCode, "\n"))

  WorkingAudiogram<- as.data.frame(BAEAudiograms[[audiogram_index]])
  WorkingAudiogram[,1]<- log10(WorkingAudiogram[,1])

  # evaluate the desired frequency by linear interpolation
  Adjustment<- approx(WorkingAudiogram, xout=log10(SonarFreq))$y
  
  # note in the case of multiple platform, `SonarFreq' is a vector;
  # as a result, `Adjustment' is a vector; the order of the elements respects
  # the order of the platforms;

  return(Adjustment)

}


#
#    --    function to evaluate the parameters for:    --
#    --          a) Dose-Response curves (SEL)         --
#    --          b) Threshold Criteria (PSP)           --
#
getEffectProbs_Params<- function(SpeciesInfo){

  #
  # > SEL
  #

  #
  # DRXXXXYY where ...
  #                XXXX = Guild ID for TTS & PTS dose-response
  #                 YY  = sonar operation mode
  #
  # DRXXXXYY is a matrix [2x2]
  # col 1: Temporary Threshold Shift (TTS)
  # col 2: Permanent Threshold Shift (PTS)
  # row 1: `n' (intercept), parameter for the logistic function
  # row 2: `m' (slope), parameter for the logistic function
  #
  # DR010100 = Cetaceans, non-pulsed mode
  # DR010101 = Cetaceans, pulsed mode
  # DR020200 = Pinnipeds, non-pulsed mode
  # DR020201 = Pinnipeds, pulsed mode
  #

  # Guild ID for TTS & PTS dose-response
  GuildID<- substring(SpeciesInfo$GuildID, 7, 10)
  # sonar opration mode
  # (all platforms operate with the same operation mode)
  OperationMode<- PrimaryShipStats$OperationMode

  DRCode<- paste("DR", GuildID, "0", OperationMode, sep='')
  cat(paste("dose-response curve:", DRCode, "\n"))
  params_index<- which(names(DRParamsSEffect)==DRCode)
  SEL_DRParams<- DRParamsSEffect[[params_index]]

  #--

  #
  # > PSP
  #

  #
  # TCXXXX where ...
  #            XXXX = Guild ID for TTS & PTS threshold criteria
  #
  # TCXXXX is a matrix [1x2]
  # col 1: Temporary Threshold Shift (TTS)
  # col 2: Permanent Threshold Shift (PTS)
  # row 1: threshold values
  #
  # TC0101 = Cetaceans
  # TC0202 = Pinnipeds
  #

  #
  # note that the Guild ID for TTS & PTS threshold criteria coincides with the
  # Guild ID for TTS & PTS dose-response 
  #

  TCCode<- paste("TC", GuildID, sep='')
  cat(paste("threshold criteria:", TCCode, "\n"))
  params_index<- which(names(ThresholdCriteriaParams)==TCCode)
  PSP_TCParams<- ThresholdCriteriaParams[[params_index]]


  #--

  # `EP_Params' is a data structure containing:
  #   a) parameters for the Dose-Response curves (SEL)
  #   b) parameters for the Threshold Criteria (PSP)
  EP_Params<- list(SEL_DRParams=SEL_DRParams, PSP_TCParams=PSP_TCParams)

  return(EP_Params)

}


#
#    --    function to evaluate the parameters for:    --
#    --          a) Dose-Response curves (SEL)         --
#
getEffectProbs_Params_m1<- function(SpeciesInfo){

  #
  # > SEL
  #

  #
  # DRXXXXYY where ...
  #                XXXX = Guild ID for TTS & PTS dose-response
  #                 YY  = sonar operation mode
  #
  # DRXXXXYY is a matrix [2x2]
  # col 1: Temporary Threshold Shift (TTS)
  # col 2: Permanent Threshold Shift (PTS)
  # row 1: `n' (intercept), parameter for the logistic function
  # row 2: `m' (slope), parameter for the logistic function
  #
  # DR010100 = Cetaceans, non-pulsed mode
  # DR010101 = Cetaceans, pulsed mode
  # DR020200 = Pinnipeds, non-pulsed mode
  # DR020201 = Pinnipeds, pulsed mode
  #

  # Guild ID for TTS & PTS dose-response
  GuildID<- substring(SpeciesInfo$GuildID, 7, 10)
  # sonar opration mode
  # (all platforms operate with the same operation mode)
  OperationMode<- PrimaryShipStats$OperationMode

  DRCode<- paste("DR", GuildID, "0", OperationMode, sep='')
  cat(paste("dose-response curve:", DRCode, "\n"))
  params_index<- which(names(DRParamsSEffect)==DRCode)
  SEL_DRParams<- DRParamsSEffect[[params_index]]


  #--

  # `EP_Params' is a data structure containing:
  #   a) parameters for the Dose-Response curves (SEL)
  #   b) parameters for the Threshold Criteria (PSP)
  EP_Params<- list(SEL_DRParams=SEL_DRParams, PSP_TCParams=NULL)

  return(EP_Params)

}


#
#    --    function to evaluate the Probabilities of the Effect     --
#    --       for the animals of a given species (from SELs)        --
#
# this function calculates the probabilities of the effect for a given species,
# taking into account the accumulated SELs;
#
# the equation for the dose-response curve is the `logistic function';
#
getEffectProbs_SEL<- function(SELs, DRParams){
  
  #
  # SELs = sound exposure levels for the animals of a given species
  #

  # no. of simulated animals per species
  NumberOfAnimals<- length(SELs)

  # effects: TTS, PTS
  NumberOfEffects<- ncol(DRParams)

  # define a data structure for the probability of the effect
  # col 1: TTS;  col 2: PTS
  EffectP_SEL<- array(dim=c(NumberOfAnimals, NumberOfEffects))

  #
  # dose-response curve
  #   dose     = sound exposure levels
  #   response = probability of the effect
  #
  for(i in 1:NumberOfEffects){
    # logistic function
    FunctionEval<- exp(DRParams[1,i]+SELs*DRParams[2,i])
    EffectP_SEL[,i]<- FunctionEval/(1+FunctionEval)
  }
    
  EffectP_SEL<- as.matrix(EffectP_SEL)
  
  return(EffectP_SEL)

}


#
#     --  function to evaluate/update the maximum value of  --
#     --         the Sound Propagation (SP) levels          --
#
getMaxSPLevel<- function(spl_max, SPLs){

  # `spl_max' is the current value of the maximum SP level
  #
  # `SPLs' contains the SP levels (SP levels are RMS values);
  # `SPLs' is a matrix:
  #    - rows refer to the time steps (during diving and surface movement)
  #    - cols refer to the platforms (the platforms that are not
  #      transmitting are marked with `NA')

  spl_max_tmp = max(SPLs, na.rm=TRUE)
  
  if(spl_max_tmp>spl_max){
    spl_max = spl_max_tmp
  }  

  return(spl_max)

}


#
#    --    function to evaluate the Probabilities of the Effect     --
#    --       for the animals of a given species (from PSPs)        --
#
# this function calculates the probabilities of the effect for a given species,
# taking into account the PSPs;
#
getEffectProbs_PSP<- function(SPLmaxs, TCParams){

  #
  # SPLmaxs = maximum values of the Sound Propagation (SP) levels for the
  #           animals of a given species
  #

  # no. of simulated animals per species
  NumberOfAnimals<- length(SPLmaxs)

  # define a data structure for the Peak Sound Pressure (PSP) levels
  PSPs<- array(0, dim=(NumberOfAnimals))

  # to calculate the PSP levels, add 3dB to the the maximum values of Sound
  # Propagation (SP) levels (note that the SP levels are RMS values);
  for(i in 1:NumberOfAnimals){
    if(SPLmaxs[i]>0){PSPs[i] = SPLmaxs[i] + 3}
  }

  # effects: TTS, PTS
  NumberOfEffects<- ncol(TCParams)

  # define a data structure for the probability of the effect
  # col 1: TTS;  col 2: PTS
  EffectP_PSP<- array(0, dim=c(NumberOfAnimals, NumberOfEffects))

  # compare the PSP levels with the threshold values
  for(i in 1:NumberOfAnimals){
      if(PSPs[i]>=TCParams[1,1]){EffectP_PSP[i,1] = 1}
      if(PSPs[i]>=TCParams[1,2]){EffectP_PSP[i,2] = 1}
  }
    
  EffectP_PSP<- as.matrix(EffectP_PSP)
  
  return(EffectP_PSP)

}



#
#    --  function to evaluate the Initial Values of X (SELs, SPLmaxs)  --
#
# in the case of multiple simulation runs, a sampling procedure is carried out
# to evaluate the initial values for the quantities X;
# the function samples from the approximate PDF (histogram) of the quantities X
# for each cell;
#
# (X are the Sound Exposure Levels or the maximum values of the Sound
#  Propagation levels) 
# 
getInitialValues<- function(HistXObject, Species, InitCoordCellC, Seeds){

  #--
    
  # RNG - assign a seed from the predetermined seed collection
  # by using RNG kind 1 = Marsaglia-Multicarry;
  rs<- c(1, Seeds[5,])
  rs_i<- as.integer(rs)
  .Random.seed<<- rs_i


  #--

  NumberOfSimAnimals<- UserSettings$NumberOfSimAnimals
  NumberOfSpecies<- length(Species)

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  CellCentresCoord<- HistXObject$Locations
  NumberOfLoc<- nrow(HistXObject$Locations)

  # no. of bins (binning procedure)
  NumberOfBins<- UserSettings$SummaryBinNumber

  # define a data structure for X
  XArray<- array(0, dim=c(NumberOfSimAnimals, NumberOfSpecies))

  # in the case of multiple simulation runs, msr !=0
  msr=0
  for(i in 1:length(HistXObject$SummaryValues)){
    if(HistXObject$SummaryValues[i]!=0){
      msr=1
    }
  }
    
  if(msr!=0){
    # sampling procedure
    for(j in 1:NumberOfSpecies){
      for(k in 1:NumberOfLoc){
        cell_index<- which(InitCoordCellC[,1,j]==CellCentresCoord[k,1] &  
                           InitCoordCellC[,2,j]==CellCentresCoord[k,2])
        if(length(cell_index)!=0){
          random_bin_index<- sample(1:(NumberOfBins-1), length(cell_index), replace=T)
          BinValues<- runif(length(random_bin_index),
                            round(HistXObject$SummaryValues[random_bin_index,j,k],3),
                            round(HistXObject$SummaryValues[random_bin_index+1,j,k],3))
          XArray[cell_index,j]<- BinValues
        }
      }
    }
  }

  return(XArray)

}



#
#    -- function to evaluate the summary (histograms) of X (SELs, SPLmaxs)  --
#
# let us consider X for all simulated animals starting their movement
# paths from a given cell within the AOC;
#
# this function estimates the approximate probability density function (PDF)
# for these X; it calculates the histogram for X for each cell within the AOC;
#
# (X are the Sound Exposure Levels or the maximum values of the Sound
#  Propagation levels)
#
getSummary<- function(XData, HistXObject, Species, InitCoordCellC){

  NumberOfSpecies<- length(Species)

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  CellCentresCoord<- HistXObject$Locations
  NumberOfLoc<- nrow(HistXObject$Locations)
  
  # if the X value is not available, initialize to 0
  for(i in 1:length(XData)){
    if(is.na(XData[i])==T){
      XData[i]=0
    }
  }

  # no. of bins (binning procedure)
  NumberOfBins<- UserSettings$SummaryBinNumber

  # define a data structure for the summary (histograms) of X;
  # X values per bins, per species, per location
  XSummary<- array(0, dim=c(NumberOfBins, NumberOfSpecies, NumberOfLoc))

  # each simulated animal has a starting location (AOC cell);
  # a given cell within the AOC is the starting location for one or more
  # simulated animals;
  #
  # summaries are made for each cell (i) and species (j)
  for(j in 1:NumberOfSpecies){
    UniqueCells<- unique(InitCoordCellC[,,j])
    for(i in 1:nrow(UniqueCells)){
      XData_tmp<- XData[which(InitCoordCellC[,1,j]==UniqueCells[i,1] &
                              InitCoordCellC[,2,j]==UniqueCells[i,2]),j]       
      Breaks<- quantile(XData_tmp, probs=seq(0, 1, length=NumberOfBins), na.rm=T)
      cell_index<- which(CellCentresCoord[,1]==UniqueCells[i,1] &
                         CellCentresCoord[,2]==UniqueCells[i,2])      
      XSummary[,j,cell_index]<- Breaks
    }

  }

  # `XSummaryData' is a data structure containing:
  # a) summary of X, b) species, c) locations
  XSummaryData<- list(SummaryValues=XSummary, Species=Species, Locations=HistXObject$Locations)

  return(XSummaryData)

}



#
#    -- function to update the summary (histograms) of X (SELs, SPLmaxs) --
#
# the inputs (`HistXObject' and `NewXObject') and the output are data
# structures containing summaries (histograms) of X;
#
# (X are the Sound Exposure Levels or the maximum values of the Sound
#  Propagation levels)
#
updateHistory<- function(HistXObject, NewXObject){

  HistSpeciesList<- unique(c(HistXObject$Species, NewXObject$Species))
  NumberOfSpecies<- length(HistSpeciesList)

  # coordinates (Lat, Long) of the centres of the cells within the AOC
  HistCellCoord<- unique(rbind(HistXObject$Locations, NewXObject$Locations))
  NumberOfLocations<- nrow(HistCellCoord)

  # no. of bins (binning procedure)
  NumberOfBins<- UserSettings$SummaryBinNumber

  # define a data structure for the summary (histograms) of X
  # X values per bins, per species, per location
  HistXSummary<- array(0, dim=c(NumberOfBins, NumberOfSpecies, NumberOfLocations))

  for(j in 1:NumberOfSpecies){
    write_species_idx<- which(HistXObject$Species==HistSpeciesList[j])
    overwrite_species_idx<- which(NewXObject$Species==HistSpeciesList[j])
    for(k in 1:NumberOfLocations){
      write_loc_idx<- which(HistXObject$Locations[,1]==HistCellCoord[k,1] &
                            HistXObject$Locations[,2]==HistCellCoord[k,2])
      overwrite_loc_idx<- which(NewXObject$Locations[,1]==HistCellCoord[k,1] &  
                                NewXObject$Locations[,2]==HistCellCoord[k,2])
      # if the latest simul. run has data for a species/locations similar to previous ...
      if(length(write_species_idx)>0 & length(write_loc_idx)>0){
        InsertVector<- HistXObject$SummaryValues[, write_species_idx, write_loc_idx]
        HistXSummary[,j,k]<- InsertVector
      }
      if(length(overwrite_species_idx)>0 & length(overwrite_loc_idx)>0){
        InsertVector<- NewXObject$SummaryValues[, overwrite_species_idx, overwrite_loc_idx]
        HistXSummary[,j,k]<- InsertVector
      }
    }
  }

  # `UpdatedHistXObject' is a data structure containing:
  # a) summary of X, b) species, c) locations
  UpdatedHistXObject<- list(SummaryValues=HistXSummary, Species=HistSpeciesList, Locations=HistCellCoord)

  return(UpdatedHistXObject)

}
