
#
# CVS - $Id: ERMCMain_S3BL.R,v 1.8 2009/12/11 13:53:01 lmilazzo Exp $
#





ERMCmain<- function(PrimaryShipStats, MultiShipStats, UserSettings){

  # Random Number Generator - settings
  RNGkind(kind="Marsaglia-Multicarry", normal.kind='Box-Muller')


  #--
  ##                  -- Initialization --              ##
  #
  
  # the `History' object is associated with the data from the previous
  # simulation runs;
  if(exists("History")==F){
    History<<- NULL
  }

  # initialize the system
  #
  # `InitData' is a list containing: `AOC', `GridCellReference',
  #                                  `InitialXYGrid', `XYSimLocation',
  #                                  `WorkingDensityArray', `WorkingDensityUncertainty',
  #                                  `WorkingSpecies',
  #                                  `InitialSELsObject', InitialSPLmaxsObject
  #
  InitData<- initializeSystem(PrimaryShipStats, MultiShipStats, UserSettings, SeedSet)
  cat("initialization ... done\n")


  #--
  ##               -- Surface and Dive Movement  --       ##
  #
  # simulate the movement of the animals
  #
  # `SystemData' is a list containing: `SELs', `SPLmaxs', 
  #                                    `FinalXY', `AnimalOutOfSea'
  #
  SystemData<- simulateAnimalMovement(InitData, SeedSet)
  cat("calculation of SELs and SPL (max. values) ... done\n")

  #--
  ##                 -- Effects of Sound --               ##
  #
  #
  # `SE_SEL' is a list containing: `ProbNoEffect', `ProbNoEffectUpper', `ProbNoEffectLower',
  #                                `NumbersEffected', `NumbersEffectedUpper`, `NumbersEffectedLower',
  #                                `TotalNumberEffected', `TotalProbNoConsequence',
  #                                `SpeciesList', `GridKey'
  #
  # `SE_PSP' is a list containing: `ProbNoEffect',
  #                                `NumbersEffected',
  #                                `TotalNumberEffected', `TotalProbNoConsequence',
  #                                `SpeciesList', `GridKey'
  #
  SE_SEL<- getSoundEffects_SEL(InitData, SystemData)
  SE_PSP<- getSoundEffects_PSP(InitData, SystemData)

  SoundEffects<- list(SE_SEL=SE_SEL, SE_PSP=SE_PSP)  
  cat("calculation of the Sound Effects ... done\n")
 

  return(SoundEffects)

}



ERMCmain_m1<- function(PrimaryShipStats, MultiShipStats, UserSettings){

  # Random Number Generator - settings
  RNGkind(kind="Marsaglia-Multicarry", normal.kind='Box-Muller')


  #--
  ##                  -- Initialization --              ##
  #
  
  # the `History' object is associated with the data from the previous
  # simulation runs;
  if(exists("History")==F){
    History<<- NULL
  }

  # initialize the system
  #
  # `InitData' is a list containing: `AOC', `GridCellReference',
  #                                  `InitialXYGrid', `XYSimLocation',
  #                                  `WorkingDensityArray', `WorkingDensityUncertainty',
  #                                  `WorkingSpecies',
  #                                  `InitialSELsObject', InitialSPLmaxsObject
  #
  InitData<- initializeSystem_m1(PrimaryShipStats, MultiShipStats, UserSettings, SeedSet)
  cat("initialization ... done\n")


  #--
  ##               -- Surface and Dive Movement  --       ##
  #
  # simulate the movement of the animals
  #
  # `SystemData' is a list containing: `SELs', `SPLmaxs', `SoundResponses',
  #                                    `FinalXY', `AnimalOutOfSea'
  #
  SystemData<- simulateAnimalMovement_m1(InitData, SeedSet)
  cat("calculation of SELs and SPL (max. values) ... done\n")

  #--
  ##                 -- Effects of Sound --               ##
  #
  #
  # `SE_SEL' is a list containing: `ProbNoEffect', `ProbNoEffectUpper', `ProbNoEffectLower',
  #                                `NumbersEffected', `NumbersEffectedUpper`, `NumbersEffectedLower',
  #                                `TotalNumberEffected', `TotalProbNoConsequence',
  #                                `SpeciesList', `GridKey'
  #
  SE_SEL<- getSoundEffects_SEL_m1(InitData, SystemData)
  numberValidAnimals<- apply(SystemData$AnimalOutOfSea, 2, function(q){length(which(q==0))})
  
  SoundEffects<- list(SE_SEL=SE_SEL, SE_PSP=NULL, SRES=SystemData$SoundResponses, Num_valid_animals=numberValidAnimals)  
  cat("calculation of the Sound Effects ... done\n")
 

  return(SoundEffects)

}
