# EMC Growth Function with dependent functions
#
# Functions in this file

#  1. Growth_EMC             : main function          
#  2. Cwt                    : Carbon from dry weight (mg) to mg C 
                               # (Hofmann & Lascara 2000) (Growth_Fn_EMC_Cwt)
#  3. DL_EMC                 : Day length (DL) 
                               # website: http://www.gandraxa.com/length_of_day.xml 
                               # (Growth_Fn_EMC_DayLength)
#  4. GonadDW                : Max Gonad DW (mg) of gravid females at length 
                               # (Morris etal 1988) (Growth_Fn_EMC_GonadDW)
#  5. IGRmax                 : maximum instantaneous growth rate (mm) 
                               # (Candy & Kawaguchi 2006) (Growth_Fn_EMC_IGRmax)
#  6. IMP_EMC                : Intermoult Period 
                               # (Kawaguchi etal 2006)(Growth_Fn_EMC_IMP)
#  7. Fcf_EMC                : Feeding rate per year (m3 per year) 
                               # (Hofmann & Lascara 2000)(Growth_Fn_EMC_Fcf)
#  8. MetC_EMC               : Metabolic requirements as food (mgC.yr-1) from dry weight
                               # (Ikeda & Dixon 1982)(Growth_Fn_EMC_MetC)
#  9. MetC_EMC_from_CarbonWT : Metabolic requirements as food (mgC.yr-1) from 
                               # Carbon Wt rather than dry weight 
                               #(Hofmann & Lascara 2000)(Growth_Fn_EMC_MetC_from_CarbonWT)
# 10. DW_from_Length         : DW from length (mm)
                               #(Growth_Fn_EMC_DW_from_Length)
# 11. CW_from_Length         : Carbon weight from length (mm)
                               # (Growth_Fn_EMC_CW_from_Length)
# 12. Growth_EMC_allocate    : routine to allocate carbon between reproduction 
                               # and growth based on proportional allocation 
                               # (Growth_Fn_EMC_allocate)
# 13. IGR_Atkinson_MaxFood   : Atkinson etal 2006 model(accompanied by parameters - Atkinson_B)
# 14. DWfromC                : Dry Weight (mg) from Carbon mgC 
                               #(Hofmann & Lascara 2000)(Growth_Fn_EMC_DWfromC)
# 15. Length_from_DW         : Length (mm) from DW (mg)
                               #(Growth_Fn_EMC_Length_from_DW)
# 16. Growth.VB              : von Bertalanffy growth 
                               # (Growth_Fn_VB)
# 17. Growth.PVB             : von Bertalanffy growth function - part year
                               # (Growth_Fn_PVB)
# 18. Growth.PVB_F           : von Bertalanffy growth function - part year 
                               #(Candy & Kawaguchi 2006) (Growth_Fn_PVB_F)
# 19. Growth.SVB             : von Bertalanffy growth function - seasonal
                               # (Siegel 1986) (eqn 4 in Candy & Kawaguchi 2006)
                               # (Growth_Fn_SVB)
# 20. Reprod_EggsMax         : Minimum reproductive output (size of gonads) 
                               # (Morris etal 1988; Nicol etal 1995)(Reprod_Fn_Eggs_max)


###############################################################################
###############################################################################

# 1. Growth_EMC

###############################################################################


#  from Constable & Kawaguchi 2012  

# Implementation Notes
#
#  1. Spatial units of Env$Food should be the same as those in Comp_Filtration_Rate
#      e.g. if food is density per m2 then the population density shuold be by area.
#           if food density is by volumen then so should the population density etc.


# Returned Data
  
#   A list of data is returned:  
#              Kstate      - active state for growth, reproduction, carbon, logical for whether Female
#            , Total       - matrix of the following data at the last moult or totals over the projection period
#                              Age
#                              Length
#                              DryWeight
#                              MoultsN       - number of full moults during projection
#                              Eggs          - total reproductive output (eggs)
#                              Ci            - Ingested carbon (filtered carbon)
#                              Ce            - egested carbon (faeces)(mg) (adjusted as assimilated carbon plus carbon not used in met, moult, G, R)
#                              C_Lost_moults - carbon lost in moulted exoskeletons (mg)
#                              Cr            - remaining carbon not consumed from filtered water (mg) 
#            , Res         - results per moult

# Input Data
  
# 1. Age         - (yrs)  at last moult
# 2. Length      - (mm)
# 3. DryWeight     - (mg dry weight)
# 4. StartProj   - proportion of calendar year passed to start this projection period (yrs)
# 5. EndProj     - length of time for projection of growth and reproduction (yrs)
# 6. Kstate      - list of State variables
#      Gro         - growth state vector  
#         * iG                   - investment in growth (mgC) within moult (magnitude of growth state)  
#         * IMP                  - prop of Inter Moult Period completed  
#         * IMPtime               - time of IMP so far (years)  
#         * DayLength_LastMoult  - the length of the day at the last moult  
#      Rep         reproductive state vector  
#         * iR         - investment in reproduction (mgC) (magnitude of reproductive state)  
#         * Moults     - how many moults have passed in stage of cycle
#         * Cycle      - which stage of the cycle is currently active
#      Carbon
#         * Cingest 
#         * Cassim  
#         * Cegest  
#         * Cexcess 
#      Female     - logical as to whether krill is female  
# old 7. Pop_Filtration_Rate  - total volume the population has filtered over one year
                            # e.g. total population filtration rate

# new 7. Filt_Eff  - Filtration efficiency (arising from competition - Proportion of maximum filtration rate possible)
# 8. Local       list of variables for local environment
#    * AreaLat     latitude of location to give day length (degrees)  
#    * Env         matrix of records of environmental conditions, where columns are  
#        1 "dt"              prop of year for which environmental conditions apply (yrs)  
#        2 "Temp"            temperature in oC  
#        3 "Food"            density of carbon in mgC.m-3 (for filtering) 
#        4 "CalPropYrStart"  prop of calendar year to start of period (yrs)  
# 9. EMC_p       list of krill EMC growth parameters


# old Growth_EMC<-function(Age,Length,DryWeight,Kstate,Pop_Filtration_Rate,StartProj,EndProj,Local,EMC_p){
Growth_EMC<-function(Age,Length,DryWeight,Kstate,Filt_Eff,StartProj,EndProj,Local,EMC_p){

# 1. Save initial conditions at beginning of projection period ####  

Age0<-Age
L0<-Length
DW0<-DryWeight


# 2. Truncate environment matrix for this projection period ####
  
  EnvMat<-Local$Env
  EMcols<-ncol(EnvMat)
  FirstRow<-max(which(EnvMat[,"CalPropYrStart"]<=StartProj))
  LastRow<-max(which(EnvMat[,"CalPropYrStart"]<EndProj))
  EnvMat<-matrix(EnvMat[c(FirstRow:LastRow),],ncol=EMcols,dimnames=dimnames(EnvMat))
  Records<-nrow(EnvMat)

  # when the start and/or end of the projection period occurs during a record then
  # the proportion of the year that a record accounts for needs to be adjusted so
  # that the sum of proportions of the year across all records equals the duration of the projection period

  # Also, need to adjust the amount of food available assuming same consumption/time

  # thus, adjust proportion of year in first and last record
  
  if (EnvMat[1,"CalPropYrStart"]<StartProj){
    dtOld<-EnvMat[1,"dt"]
    EnvMat[1,"dt"]<-EnvMat[1,"dt"]-(StartProj-EnvMat[1,"CalPropYrStart"])
    EnvMat[1,"CalPropYrStart"]<-StartProj
    EnvMat[1,"Food"]<-EnvMat[1,"Food"]/dtOld*EnvMat[1,"dt"]
  }
  if (EnvMat[Records,"CalPropYrEnd"]>EndProj){
    dtOld<-EnvMat[Records,"dt"]
    EnvMat[Records,"dt"]<-EnvMat[Records,"dt"]-(EnvMat[Records,"CalPropYrEnd"]-EndProj)
    EnvMat[Records,"CalPropYrEnd"]<-EndProj
    EnvMat[Records,"Food"]<-EnvMat[Records,"Food"]/dtOld*EnvMat[Records,"dt"]

  }
#    print(EnvMat)
#    print(dimnames(EnvMat)[[2]])
#    for (wi in 1:nrow(EnvMat)) cat(round(EnvMat[wi,],3),"\n",sep = "\t")  
  
# 3. Loop through environment matrix to do moults in each record ####

Res <- NULL   # results for each moult

LastMoult_Age                 <- NA  # at last full moult
LastMoult_Length              <- NA  # at last full moult
LastMoult_DryWeight           <- NA  # DryWeight   # at last full moult
LastMoult_GonadDW             <- 0  # Gonad DryWeight   # at last full moult
LastMoult_Total_MoultsN       <- 0   # Number of full moults during projection
LastMoult_Total_Eggs          <- 0   # total reproductive output (eggs)
LastMoult_Total_Ci            <- 0   # Ingested carbon (filtered carbon)
LastMoult_Total_Ca            <- 0   # assimilated carbon
LastMoult_Total_Ce            <- 0   # egested carbon (faeces)(mg) (adjusted as assimilated carbon plus carbon in excess)
LastMoult_Total_C_Lost_moults <- 0   # carbon lost in moulted exoskeletons (mg)


for (em in 1:nrow(EnvMat)){
  
# 3.1 calculate the number of moults in the environment record ####
      # according to the length of IMP arising from the temperature and the proportion of year

  MoultsPerRec<-EnvMat[em,"dt"]/IMP_EMC(EnvMat[em,"Temp"])
  
  # determine number of moults in record (note that more than one record may be be required for a moult)
  
  MoultsByRecordEnd<-MoultsPerRec+Kstate$Gro$IMP

  Nmoults<-floor(MoultsByRecordEnd)
  
# 4. determine a Moult-Matrix for calculating ration and metabolic demand per moult ####
  
#   The matrix has one record per moult.
#   For each moult, the total ration will be the sum of ration across all records minus the proportion of the first record not in the moult (i.e. before the moult) and the proportion of the last record not in the moult (i.e. after the moult)
#   
#   Columns for each moult:
#   1. the proportion of IMP covered by the record (not including the proportion of IMP at start of this routine)  
#   2. proportion of the record prior to the moult  
#   3. proportion of the record after the moult
  
  MoultMatrix<-NULL
  
  # If no moult is completed in the record then accumulate all data without moult
  if (Nmoults<1) {  # moult is not completed in this interval
    MoultMatrix<-matrix(c(MoultsPerRec,0,0),nrow=1)
    Nmoults<-1
    DoMoult<-c(FALSE) # logical vector on whether to moult for each record of the moult matrix
    LastMoultIncomplete<-TRUE
    
  } else {
    # If at least one moult is completed then have one record in MoultMatrix for each moult plus an additional record for
    # any time in the projection period remaing as an incomplete moult
    
    # what proportion of the IMP is at the end of the previous environment record
    IMPendPrevInterval<-Kstate$Gro$IMP

    LastMoultIncomplete<-FALSE
    DoMoult<-rep(TRUE,Nmoults)
    
    if (MoultsByRecordEnd>Nmoults) {
      Nmoults<-Nmoults+1  # last "moult" is only a partial IMP
      DoMoult<-c(DoMoult,FALSE)
      LastMoultIncomplete<-TRUE
    }
    for (m in 1:Nmoults){
      if (MoultsByRecordEnd>=m) {
        IMPprop<-1    #stage of IMP at end of record
        PropRecGTm<-(MoultsByRecordEnd-m)/MoultsPerRec # proportion of record left over
      } else { # last moult is incomplete
        IMPprop<-MoultsByRecordEnd-(Nmoults-1)  # proportion of the IMP completed by end of record for incomplete moult
        PropRecGTm<-0    # proportion of period in next IMP
      }
      
      if(m==1) IMPprop<-IMPprop-Kstate$Gro$IMP
      PropRecLTm<-(m-1-IMPendPrevInterval)/MoultsPerRec
      if(PropRecLTm<0) PropRecLTm<-0
    
      
      MoultMatrix<-rbind(MoultMatrix,c(IMPprop,PropRecLTm,PropRecGTm))
    } # end for Nmoults
  } # end if Nmoults>=1


# 5. loop through each moult and then do any leftover assignment ####

   for (m in 1:Nmoults){

    
# 5.1 determine time of the moult plus time of each record in the moult (yr) ####

  Trec     <- EnvMat[em,"dt"]*(1-MoultMatrix[m,2]-MoultMatrix[m,3])
 
## 4.3 update attributes of the IMP in readiness for the moult (noting that the moult may not have been completed in last record)
IMPtime<-Kstate$Gro$IMPtime+Trec

if(DryWeight>0){ # only carry this out if the krill have not starved
  
# 5.2 Volume of water filtered (m3 for time of record) ####
    
    Fvol <- Fcf_EMC(DryWeight)*Trec

# 5.3 determine food assimilated, metabolic demand and cost of moult ####

    # Carbon ingested, assimilated and egested (mgC)
    # routine to allow for maximum ingestion rate given food density
    # NB - reduction in assimilation efficiency with increased ingestion is not included here.  This is back calculated
    #      later by egesting mC assimilated here but not used in growth, reproduction and metabolism


 # old   PopFiltration<-Pop_Filtration_Rate*EnvMat[em,"dt"]  # convert per annum population filtration to the amount filtered within the time period of the food available
# old    Cingest <-  {IngestAvail <- EnvMat[em,"Food"]/PopFiltration*Fvol

# EnvMat[em,"Food"] should be mgC.m-3

   Cingest <- EnvMat[em,"Food"]*Fvol*Filt_Eff

 #   if(EMC_p$Ingest$Max){
        # in order to take account of competition, a simple function is used here
 #       CingestMax <- EMC_p$Ingest$MaxIngest*Fvol
 #       if (Cingest>CingestMax) Cingest <- CingestMax
 #       }
 
    Cassim  <- EMC_p$AssimilationEfficiency*Cingest
    Cegest  <- Cingest-Cassim
    Cexcess <- Kstate$Carbon$Cexcess   # gathers up excess that can be used for moult - need to recall from previous incomplete moult 

# 5.4 initialise variables for moult

    Tsteps  <- ceiling(Trec/EMC_p$Met$dt)  # number of time steps
    dst     <- Trec/Tsteps                  # length of time step (years)

    #         determine maximum increment for growth and reproduction given length and maximum investment in carbon to achieve it
    #         (noting that carbon is calculated based on ideal weight)  

    dGmax   <- IGRmax(Length)             # maximum growth increment in a moult (mm)
    iGmax   <- CW_from_Length((Length+dGmax),EMC_p$LW,2)-CW_from_Length(Length,EMC_p$LW,2) # (mgC)

    Rmax    <- GonadDW(Length,EMC_p$LW)   # maximum dry weight investment in gonad over 2 IMPs (mg)
    iRmax   <- Rmax*EMC_p$GonadConversion # max resource requirements for gonad (mgC) (note that this could be achieved in 1 IMP)

    CwIdeal <- CW_from_Length((Length),EMC_p$LW,2) # ideal weight given the current length
#old    Cw      <- Cwt(DryWeight)

    BodyCarbon<-Kstate$Carbon$BodyCarbon #new   not differentiating between body and gonad at this point
    BodyCarbonMax<-CwIdeal+iGmax+iRmax  #new

    Ca_dst  <- Cassim/Tsteps

#old    iG      <- Kstate$Gro$iG     # investment in growth as the difference between current weight and ideal weight (CwIdeal)
#old    iR      <- Kstate$Rep$iR

# 5.5 integrate metabolism and body carbon (body + gonad) weight through time steps

for (ts in 1:Tsteps){
  
   # Metabolic demand given Temperature (mgC for moult time in record)
   Cw_prop_ideal<-(BodyCarbon/CwIdeal)
  Cmet<-MetC_EMC_from_CarbonWT(EnvMat[em,"Temp"],BodyCarbon,EMC_p$Met,Cw_prop_ideal)*dst
     if(EMC_p$Met$Male$HiRate){ # if males have higher metabolic rate than females, this only applies when males are mature (length) and the body is not starved or hibernating
       if(!Kstate$Female & Length<=EMC_p$Met$Male$Mature_Length & Cw_prop_ideal>EMC_p$Met$Starving_Prop_Ideal_Cw) Cmet<-Cmet*EMC_p$Met$Male$Multiplier
       }
  BodyCarbon<-BodyCarbon+Ca_dst-Cmet

  if(BodyCarbon>BodyCarbonMax){
    Cexcess<-Cexcess+BodyCarbon-BodyCarbonMax
    BodyCarbon<-BodyCarbonMax
  }
  
  } # end timestep


# 4.6.1  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#       If do Moult

if(DoMoult[m]){
  
  
  # calculate day length at this moult as this can affect different parts of the reproduction cycle
  
  TimeOfYear <- EnvMat[em,"CalPropYrStart"]+EnvMat[em,"dt"]*(1-MoultMatrix[m,3])
  DayLength<-DL_EMC(Local$AreaLat,TimeOfYear)
  
 
  # Food ration available for growth and reproduction = C assimilated less metabolic demand
  
  Ration  <- BodyCarbon-CwIdeal
  
  
  # Carbon lost with the moulted exoskeleton and needs replacing - this is removed from carbon budget in a series of step
  
  Carbon_ThisMoult<-DW_from_Length(Length,EMC_p$LW,2)*EMC_p$Moult$BodyDwProp*EMC_p$Moult$PropC
  
  
  # Step 1 : if Cexcess greater than zero then use excess
  
  Cexcess <- Cexcess - Carbon_ThisMoult
  
  # Step 2 : if insufficient in Excess then take from Ration
  
  if (Cexcess<0) {Ration <- Ration+Cexcess; Cexcess<-0}
  

    ### 4.4.2  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     Determine actual investment (mgC) in growth and reproduction based on ration, taking account of maxima for investment.
#     Options are the following:
#       1. If Mature, Female and reproductive cycle is in stage of egg production i.e. growth is influenced by investment in egg production
#     2. if immature or male i.e. growth is not influenced by egg production
#     3. if ration =< 0, i.e. reproductive reserves have been depleted to 0 and shrinkage will occur


if(Ration > 0) {
      
      # 4.4.2.2 ************************************************
      #         If Mature, Female and reproductive cycle is egg production 
      #         i.e. growth is influenced by investment in egg production
  
      if(((EMC_p$MaturityByAge & Age>=EMC_p$AgeMature) | 
            (!EMC_p$MaturityByAge & Length>=EMC_p$LengthMature)) & 
           Kstate$Female & 
           !(Kstate$Rep$Cycle<3)){
        
        # determine investment in each of growth and reproduction (mgC)
        res<-Growth_EMC_allocate(Ration,iRmax,iGmax,EMC_p$PropReprodBeforeGrowth)
        iG<-res[1]
        iR<-res[2]
        
        
        # 4.4.2.3 ************************************************
        #         If immature or male i.e. growth not influenced by egg production
      } else { 
        if (Ration > iGmax) {iG<-iGmax} else {iG<-Ration}
        iR<-0
      }

      # 4.4.2.4  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
    } else {  #  If Ration <= 0
      
      iG<-Ration
      iR<-0   # this assumes that any resources in reproduction are resorbed for use in metabolism
    }



# 5.6 Update growth and reproduction for moult

     Gonad_DW_this_moult<-0  # set to zero unless changed in moult
     Rout <- 0  # set to zero unless changed in moult
    AgeReprod<-NA  # only have !NA when reproducing

      
      # 4.6.1.1 ************************************************
      #         growth increment
      Kstate$Carbon$BodyCarbon<-  CwIdeal+iG
      if(Kstate$Carbon$BodyCarbon<=0)  Kstate$Carbon$BodyCarbon<-0
        
      Weight_new<-DWfromC(Kstate$Carbon$BodyCarbon)
      Length_new<-Length_from_DW(Weight_new,EMC_p$LW,2)

    
      # for animals that shrink, this will be limited by daily shrink rate

      if(iG < 0 & Weight_new>0) {
        if(!EMC_p$Shrink) {
          Length_new<-Length 
        } else {
          # calculate new length based on maximum shrinkage in carbon 
          IMPdays<-(IMPtime*365)
          Length_new_MaxShrink<-Length-(EMC_p$ShrinkParams$slope*Length + EMC_p$ShrinkParams$int)*IMPdays
         if(Length_new_MaxShrink>Length_new) Length_new<-Length_new_MaxShrink
         }
      } # end shrinkage

      # 4.6.1.2 
      #  reproduction - only for mature females
      #  given a new moult, update the stage in the annual reproductive cycle and, if in correct stage, the investment in the gonad
      
            
      if(((EMC_p$MaturityByAge & Age>EMC_p$AgeMature) | 
            (!EMC_p$MaturityByAge & Length>EMC_p$LengthMature)) & 
           Kstate$Female){
#        
        # add moult
        Kstate$Rep$Moults<-Kstate$Rep$Moults+1
        
        # action depends on the stage during the year in the reproductive cycle
        #    Stage 1: Regressing
        #    Stage 2: Maturing
        #    Stage 3: Egg production
        
        if(Kstate$Rep$Cycle<3) {
          
          # 4.6.1.2.1 if regressing or body maturing then no investment in reproduction
          
          if(Kstate$Rep$Cycle==1) {
            # if regressing then check day length and see if regression should
            #         be accelerated (not more than 3 moults after critical day length
            #        if moults <2 then moults equals 2
            if(DayLength<EMC_p$CritDayLength & Kstate$Rep$Moults<2) Kstate$Rep$Moults<-2
          }
          if(Kstate$Rep$Cycle==1 & floor(Kstate$Rep$Moults)==EMC_p$nMoultsRegress) {
            # if moults = nMoultsRegress and regressing then stage set to body maturing and moults set to zero
            Kstate$Rep$Cycle  <- 2
            Kstate$Rep$Moults <- 0
          }
          
          if(Kstate$Rep$Cycle==2 & Kstate$Rep$Moults==EMC_p$nMoultsMature) {
            # if moults = nMoultsMature and body maturing then stage set to egg production and moults set to zero
            Kstate$Rep$Cycle  <- 3
            Kstate$Rep$Moults <- 0
          }
          
        } else {
          # 4.6.1.2.2 if producing eggs then update investment in reproduction and state
          Gonad_DW_this_moult<-iR/EMC_p$GonadConversion # dry weight for reporting
          
          if(Kstate$Rep$Moults>=2 & iR>0){ # begin egg production
            iRminForSpawn<-EMC_p$MinReprodProp*iRmax
            if (iR >=iRminForSpawn) {
              Rout<-Gonad_DW_this_moult
              AgeReprod<-Age+IMPtime
              iR<-0
              Kstate$Rep$Moults<-0
            } else {
              if (Kstate$Rep$Moults>2 & 
                    DayLength < EMC_p$DLreprod &
                    Kstate$Gro$DayLength_LastMoult > DayLength) {
                # failed to reproduce in 2 moult cycles, regression is entered if DeclineDayLength & Daylength< DLreprod
                Kstate$Rep$Moults<-0
                Kstate$Rep$Cycle<-1
              }
            }
          } # end egg production
          
          # if time of year has passed critical day length and krill is still producing eggs
          # then test to see if regression stage (Cycle = 1) should begin
          
          if(Kstate$Rep$Cycle==3
             & DayLength < EMC_p$CritDayLength
             & (Kstate$Gro$DayLength_LastMoult > DayLength)){
            Kstate$Rep$Moults<-0
            Kstate$Rep$Cycle<-1
            # note that iR is not set to zero here and the energy reserves can be used during next IMP for growth
          }
        } # end reproduction
      } # end mature female


      # 4.6.1.3 ************************************************
      #         update state of krill in projection period
 
if(Weight_new>0){
  dG<-Length_new-Length
  Length<-Length_new
  DryWeight<-Weight_new
  GonadDW<-Gonad_DW_this_moult
  Age<-Age+IMPtime
  Kstate$Gro$IMP <- 0  # update IMP for Kstate$Gro here
  Kstate$Gro$IMPtime<-0
  Kstate$Gro$DayLength_LastMoult<-DayLength
} else { # cohort died from starvation
  dG<-0
  Length<-0
  DryWeight<-0
  GonadDW<-0
  Age<-Age
  Kstate$Gro$IMP <- 0 
  Kstate$Gro$IMPtime<-IMPtime
  Kstate$Gro$DayLength_LastMoult<-Kstate$Gro$DayLength_LastMoult
}

# note the condition relative to ideal weight

StarveCondition<-(DryWeight+GonadDW-DW_from_Length(Length,EMC_p$LW,(if(Kstate$Female) 2 else 1)))/Length

# 4.7 -------------------------------------------------------------------------------------------------------------------------------------------------
#     accumulate results for moult (as a vector and converted to matrix at end of projection period
Res<-c(Res      # results so far
       ,m        # number of the moult
       ,Age
       ,if(Kstate$Female) 2 else 1
       ,IMPtime
       ,Ration
       ,Cingest
       ,Cmet
       ,Kstate$Carbon$Cegest+Cegest+Cexcess
       ,Carbon_ThisMoult
       ,Length
       ,DryWeight
       ,Gonad_DW_this_moult
       ,StarveCondition
       ,AgeReprod #RepAge
       ,Rout*EMC_p$GonadDW_to_eggs # if Gonad_DW above critical size
       )

LastMoult_Age                 <- Age  
LastMoult_Length              <- Length  
LastMoult_DryWeight           <- DryWeight
LastMoult_StarveCondition     <- StarveCondition
LastMoult_GonadDW             <- Gonad_DW_this_moult
LastMoult_Total_MoultsN       <- if(Weight_new>0) (LastMoult_Total_MoultsN+1) else LastMoult_Total_MoultsN
LastMoult_Total_Eggs          <- LastMoult_Total_Eggs+Rout*EMC_p$GonadDW_to_eggs
LastMoult_Total_Ci            <- LastMoult_Total_Ci+Kstate$Carbon$Cingest+Cingest
LastMoult_Total_Ca            <- LastMoult_Total_Ca+Kstate$Carbon$Cassim+Cassim
LastMoult_Total_Ce            <- LastMoult_Total_Ce+Kstate$Carbon$Cegest+Cegest+Cexcess
LastMoult_Total_C_Lost_moults <- LastMoult_Total_C_Lost_moults+Carbon_ThisMoult
  
Kstate$Carbon$Cingest    <- 0
Kstate$Carbon$Cassim     <- 0
Kstate$Carbon$Cegest     <- 0
Kstate$Carbon$Cexcess    <- 0
Kstate$Carbon$BodyCarbon <- Kstate$Carbon$BodyCarbon+iR


      # 4.6.2  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #       If NOT do Moult
    } else { # end DoMoult  ELSE
      
      # need to put portion in iG and iR along with prop IMP and IMPtime so far
           
      Kstate$Gro$IMP    <- (Kstate$Gro$IMP+MoultMatrix[m,1])
      Kstate$Gro$IMPtime <- IMPtime
      Kstate$Carbon$Cingest    <- Cingest
      Kstate$Carbon$Cassim     <- Cassim
      Kstate$Carbon$Cegest     <- Cegest
      Kstate$Carbon$Cexcess    <- Cexcess
      Kstate$Carbon$BodyCarbon <- BodyCarbon
      
    } # end else not do moult

} else { # end if krill had not already starved 
  # else do accounting of the time
  
  Kstate$Gro$IMPtime <- IMPtime
  
} # end else if starved   
    
  } # end loop m


} # end records of Environment Matrix

  
  #####################################################################################################################################################
    
  # 6. Close function

Res_Total <-list(Vars = 12 # vector of results accumulated across projection period
                ,ByRow = TRUE
                ,Data  = c(LastMoult_Age
                           ,if(Kstate$Female) 2 else 1
                           ,LastMoult_Length
                           ,LastMoult_DryWeight
                           ,LastMoult_GonadDW
                           ,LastMoult_StarveCondition
                           ,LastMoult_Total_MoultsN
                           ,LastMoult_Total_Eggs
                           ,LastMoult_Total_Ci
                           ,LastMoult_Total_Ca
                           ,LastMoult_Total_Ce
                           ,LastMoult_Total_C_Lost_moults
                          ) # end Data
                ,Names = c("Age"   # age at last moult
                           ,"Sex"
                           ,"Length"
                           ,"DryWeight"
                           ,"GonadDW"
                           ,"StCond"        # starvation condition
                           ,"MoultsN"       # number of full moults during projection
                           ,"Eggs"          # total reproductive output (eggs)
                           ,"Ci"            # Ingested carbon (filtered carbon)
                           ,"Ca"            # assimilated carbon
                           ,"Ce"            # egested carbon (faeces)(mg) (adjusted as assimilated carbon plus carbon not used in met, moult, G, R)
                           ,"C_Lost_moults" # carbon lost in moulted exoskeletons (mg)
                )# end names
       ) # end Res_Total

Res_Moults<- list(Vars = 15 # vector of results accumulated across projection period
                 ,ByRow = TRUE
                 ,Data  = Res # results per moult
                 ,Names = c("Moult"
                            ,"Age"
                            ,"Sex"
                            ,"IMPtime"
                            ,"Ration"
                            ,"Cingest"
                            ,"Cmet"
                            ,"Cegest"
                            ,"Carbon_ThisMoult"
                            ,"Length"
                            ,"DryWeight"
                            ,"GonadDW"    # >0 even if no eggs when only first moult, or insufficient to reproduce
                            ,"StCond" # starvation condition of krill
                            ,"AgeReprod"
                            ,"Eggs"
                 )# end names
          ) # end results per moult


  return(list(Age = Age
            , Length    = Length
            , DryWeight = DryWeight
            , GonadDW   = LastMoult_GonadDW
            , Kstate = Kstate    # active state for growth, reproduction, carbon & logical for whether Female
            , Total  = Res_Total # vector of results accumulated across projection period
            , Moults  = Res_Moults # vector of results accumulated across projection period
            ))
  
} # end EMC growth function


###############################################################################
###############################################################################

# 2. Growth_Fn_EMC_Cwt


###############################################################################

# Carbon from dry weight (mg) to mg C (Hofmann & Lascara 2000)

Cwt<-function(DW){0.366*DW^1.037}



###############################################################################
###############################################################################

# 3. source(paste(dir_fns,"Growth_Fn_EMC_DayLength.r",sep=""))       # checked


###############################################################################


# 7. Day length (DL)
# website: http://www.gandraxa.com/length_of_day.xml

#  when DL < 0  then 24 hr darkness, DL>2 then 24 hr light

DL_EMC<-function(Lat,PropCalYr){
#  Axis<-23.439/180*pi  # degrees - does change over 1000s of years so assume constant
#  m<-1-tan(Lat/180*pi)*tan(Axis*cos(0.0172*JulianDay))
#  m<-1-tan(Lat/180*pi)*tan(Axis*cos(0.0172*365*PropCalYr))


  m<-1-tan(0.01745329*Lat)*tan(0.4090877*cos(6.278*PropCalYr))

  m[(m<0)] <-0
  m[(m>2)] <-2

#  return(24*acos(1-m)/pi)
  return(7.639437*acos(1-m))
  }

# tdays<-c(1:365)/365
# tLat<--55
# tDL<-DL_EMC(tLat,tdays)
# plot(tdays,tDL)
#


###############################################################################
###############################################################################

# 4.source(paste(dir_fns,"Growth_Fn_EMC_GonadDW.r",sep=""))         # checked


###############################################################################

# Max Gonad DW (mg) of gravid females at length
# derived from Morris etal 1988 equations

GonadDW<-function(L,LW){  # DW
      return(DW_from_Length(L,LW,3)-DW_from_Length(L,LW,2))
          }


###############################################################################
###############################################################################

# 5. source(paste(dir_fns,"Growth_Fn_EMC_IGRmax.r",sep=""))          # checked


###############################################################################

# 5. IGR (mm)
# derived from Candy & Kawaguchi 2006, Table 2 December predicted IGR at size
# method to estimate parameters - max likelihood
# IGR = m*Length   Length<= Length at maturity
# IGR = m*Length_at_maturity*(1/(Length-Length_at_maturity+1))^z    Length>Length_at_maturity

IGRmax<- function(Length){
#       IGR <- 0.22498*Length
#       IGR[Length>36.173]<-(8.138202)*(1/(Length[Length>36.173]-35.173))^0.6364


  IGR <- 0.09326*Length
  IGR[Length>39.28832]<-(3.664029)*(1/(Length[Length>39.28832]-38.28832))^0.37142

     #  very old  IGR<-8.63-0.14*Length
    return(IGR)
    }



###############################################################################
###############################################################################

# 6. source(paste(dir_fns,"Growth_Fn_EMC_IMP.r",sep=""))             # checked


###############################################################################

# 4. IMP
# Kawaguchi etal 2006 multiplied by (1/365) [inserted into log domain]
# IMP=exp(3.5371-0.5358*log(Temp+2))/365
# returns proportion of the year
IMP_EMC<-function(Temp) exp(-2.362797-0.5358*log(Temp+2))



###############################################################################
###############################################################################

# 7. source(paste(dir_fns,"Growth_Fn_EMC_Fcf.r",sep=""))             # checked


###############################################################################

# 3. Feeding rate per year (m3 per year)
# assumes no energetic cost in this model - excess consumption is not used and added back to egested carbon.

# Hofmann & Lascara 2000 regression: Feeding rate (in their model it is m3 per day)
# note that the Figure 3 cannot be replicated with the  parameters of Table 3 in that paper.
# the coefficients here are  altered  x 10-1 in order to replicate the figure)

# converted to m3 per year

Fcf_EMC<-function(dw){  # dw = dry weight (mg)
     FR<-dw*0
     FR[dw<26]<- (0.000085*dw[dw<26]^0.825) # DW<26 mg
     FR[dw>84]<- (0.000343*dw[dw>84]^0.514) # DW>84 mg
     Wt26<-(84-dw)/58
     Wt84<-(dw-26)/58
     FR[dw>=26 & dw<=84]<-(Wt26[dw>=26 & dw<=84]*(0.000085*dw[dw>=26 & dw<=84]^0.825)+
                           Wt84[dw>=26 & dw<=84]*(0.000343*dw[dw>=26 & dw<=84]^0.514))
     FR*365  # converts to m3 per year  ### note that 5*  would be used to make outcome comparable with Quetin & Ross review in Nicol etal 1995
    }



###############################################################################
###############################################################################

# 8. source(paste(dir_fns,"Growth_Fn_EMC_MetC.r",sep=""))            # checked (takes vector of temperature but single value of Length


###############################################################################

# 2. Metabolic requirements as food (mgC.yr-1)

# metabolism microlitre per hour
# metabolic requirements in Carbon per day
# Hofmann & Lascara 2000 based on Ikeda & Dixon 1982  - conversion ratio of 0.5357 ml 02 consumed per mg C.
# multiply hour by 24 to give per day

# returns metabolic requirements (mgC) of krill at length over a year for each temperature
#    16.35244=365*(24/535.7)  using carbon conversion of Hofmann & Lascara
MetC_EMC<-function(Temp   # vector of temperature
                  ,DW     # dry weight (mg)
                  ) return(16.35244*10^(0.02438*Temp-0.1838)*DW^(0.8918-0.0109*Temp))

###############################################################################
###############################################################################

# 9. source(paste(dir_fns,"Growth_Fn_EMC_MetC_from_CarbonWT.r",sep=""))


###############################################################################

# 2. Metabolic requirements as food (mgC.yr-1)

# metabolism microlitre per hour
# metabolic requirements in Carbon per day
# Hofmann & Lascara 2000 based on Ikeda & Dixon 1982  - conversion ratio of 0.5357 ml 02 consumed per mg C.
# multiply hour by 24 to give per day

#formula was for DW mg conversion but modified to be for DW=(Cw/0.366)^(1/1.037)
# returns metabolic requirements (mgC) of krill at length over a year for each temperature
#    16.35244=365*(24/535.7)  using carbon conversion of Hofmann & Lascara

# old formula based on Dry Weight was
# return(16.35244*10^(0.02438*Temp-0.1838)*DW^(0.8918-0.0109*Temp))

# formula simplified and based on carbon weight (mg)
#

MetC_EMC_from_CarbonWT<-function(Temp   # vector of temperature
                  ,Cw     # dry weight (mg)
                  ,Met_p  # parameters for reducing metabolism
                  ,Cw_prop_ideal # carbon body weight as a proportion of ideal carbon body weight
                  ){
  MetC_per_year<-10^(0.02438*Temp+1.029783)*(Cw/0.366)^(0.8599807-0.01051109*Temp)
#  MetC_per_year<-16.35244*10^(0.02438*Temp-0.1838)*(Cw/0.366)^((1/1.037)*(0.8918-0.0109*Temp))
if(Cw_prop_ideal<Met_p$Starving_Prop_Ideal_Cw)  {
    if(Cw_prop_ideal<Met_p$Hibernation_Prop_IdealCw) {
      MetC_per_year<-Met_p$Hibernation_Metabolism_prop*MetC_per_year
    } else {MetC_per_year<-Met_p$Starving_Metabolism_prop*MetC_per_year}
  }
 return (MetC_per_year)
}



###############################################################################
###############################################################################

# 10. source(paste(dir_fns,"Growth_Fn_EMC_DW_from_Length.r",sep=""))  # checked


###############################################################################

# 1. Morphometrics - DW from length (mm)

DW_from_Length<-function(L,LW,Row){
  DW<-(LW[Row,3]*(LW[Row,1]*L^LW[Row,2])+LW[Row,4])# DW from mm  where DW units are determined by LW[Row,3]
  DW[DW<0]<-0  # protection from negative weight
  DW
  }

WW_from_Length<-function(L,LW,Row){
  WW<-(LW[Row,1]*L^LW[Row,2])# WW from mm  where WW units are determined by WW[Row,3]
  WW[WW<0]<-0  # protection from negative weight
  WW
  }

###############################################################################
###############################################################################

# 11. source(paste(dir_fns,"Growth_Fn_EMC_CW_from_Length.r",sep=""))  # checked

# 1. Morphometrics - Carbon weight from length (mm)

CW_from_Length<-function(L,LW,Row){
  CW<-0.366*(LW[Row,3]*(LW[Row,1]*L^LW[Row,2])+LW[Row,4])^1.037# Cw mgC from mm
  CW[CW<0]<-0  # protection from negative weight
  CW
  }


###############################################################################

###############################################################################
###############################################################################

# 12. source(paste(dir_fns,"Growth_Fn_EMC_allocate.r",sep=""))        # checked (not vectorised)



###############################################################################

# Growth function - Energetics Moult-Cycle Function
# Constable & Kawaguchi 2012

# allocation routine between reproduction and growth based on proportional allocation

# version 1
# return investment in growth (iG), reproduction (iR)

#################################################

Growth_EMC_allocate<-function(
                 Ration
                ,iRmax     # maximum investment possible in reproduction
                ,iGmax     # maximum investment possible in growth
                ,PropRb4G  # proportion Reproduction before Growth
                ){

      if (PropRb4G>0 & PropRb4G<1){
           if((iRmax/PropRb4G)<=(iGmax/(1-PropRb4G))){
               if(Ration<(iRmax/PropRb4G)) {
                     iG <- (1-PropRb4G)*Ration
                     } else {
                     if(Ration < (iRmax+iGmax)) iG<-(Ration-iRmax) else iG<-iGmax
                     }
                } else {
                if (Ration<= (iGmax/(1-PropRb4G))) iG<-((1-PropRb4G)*Ration) else iG<-iGmax
                }
          } else {
            if(Ration<=(PropRb4G*iRmax)) {
                  iG<-0
                } else {
                  if(Ration<(iGmax+PropRb4G*iRmax)) {
                    iG<-(Ration-PropRb4G*iRmax)
                    } else {
                    iG<-iGmax
                    }
                 }
          }
      iR<-Ration-iG
      if (iR>iRmax) iR<-iRmax
  return(c(iG,iR))

  } # end EMC allocation function


###############################################################################
###############################################################################

# IGR_Atkinson_MaxFood


###############################################################################


# 5a. Atkinson etal 2006 model (Growth is % increase in their model - hence multiply by L)
IGR_Atkinson_MaxFood<-function(L,Temp,Food,B)((B[1]+B[2]*L+B[3]*L^2+B[4]*Food/(B[5]+Food)+B[6]*Temp+B[7]*Temp^2)*L/1000)
Atkinson_B<-c(6.60,-0.385,0.00259,17.53,0.332,0.595,-0.477)


###############################################################################
###############################################################################

# 14. source(paste(dir_fns,"Growth_Fn_EMC_DWfromC.r",sep=""))



###############################################################################

# Dry Weight (mg) from Carbon mgC (Hofmann & Lascara 2000)

DWfromC<-function(Cwt){(Cwt/0.366)^(1/1.037)}

###############################################################################
###############################################################################

# 15. source(paste(dir_fns,"Growth_Fn_EMC_Length_from_DW.r",sep=""))

# Morphometrics - Length (mm) from DW (mg)

Length_from_DW<-function(DW,LW,Row){
  ((DW-LW[Row,4])/(LW[Row,3]*LW[Row,1]))^(1/LW[Row,2])
  }



###############################################################################

###############################################################################
###############################################################################

# 16. source(paste(dir_fns,"Growth_Fn_VB.r",sep=""))



###############################################################################


# von Bertalanffy growth

Growth.VB<-function(P
                   ,Age  # data for input to the model:  Age
                   ,Time  # start time
                   ,dTime # length of time interval
                   ){
      return (P$Linf*(1-exp((-P$K)*(Age+Time+dTime)-P$T0)))
    }

###############################################################################
###############################################################################

# 17. source(paste(dir_fns,"Growth_Fn_PVB.r",sep=""))



###############################################################################

Growth.PVB<-  # von Bertalanffy growth function - part year
           function(P     # parameters  Linf
                                      #, K
                                      #, T0
                                      #, Fs = fraction of year when growth starts (MUST be greater than 0)
                                      #, g = fraction of year over which growth occurs
                   ,Age     # data for input to the model:  Age
                   ,Time  # start time
                   ,dTime # length of time interval
                   ){

                   FgFromY0<-0  # how much initial time in the year (fraction of year from the beginning of the year) does growth occur
                   if((P$g+P$Fs)>1) FgFromY0<-(P$g+P$Fs)-1

                   Yplus<-floor(Age+Time+dTime)
                   Fy<-(Age+Time+dTime)-Yplus
                   Gt<-Fy*0
                   Gt[Fy>=P$Fs]<-(FgFromY0+(Fy[Fy>=P$Fs]-P$Fs))/P$g
                   Gt[Fy<=FgFromY0]<-Fy[Fy<=FgFromY0]/P$g
                   Gt[Fy>FgFromY0 & Fy<P$Fs]<-FgFromY0/P$g
                   Gt[Gt>1]<-1

    return (P$Linf*(1-exp(-P$K*(Yplus+Gt)-P$T0)))
    }

###############################################################################
###############################################################################

# 18. source(paste(dir_fns,"Growth_Fn_PVB_F.r",sep=""))



###############################################################################

Growth.PVB_F<-  # von Bertalanffy growth function - part year - Candy & Kawaguchi 2006
           function(P     # parameters  Linf, K, T0, g(fraction of year), Fs(start time as fraction of year)
                   ,Age     # data for input to the model:  Age
                   ,Time  # start time
                   ,dTime # length of time interval
                   ){

                   FgFromY0<-0 # partial growth increment at beginning of year
                   FendYr<-0 # partial growth increment at end of year
                   if((P$g+P$Fs)>1) {
                       FgFromY0<-(P$g+P$Fs)-1
                       FendYr<-1-P$Fs
                       }

                   Yplus<-floor(Age+Time+dTime)
                   Fy<-(Age+Time+dTime)-Yplus

                   # fraction of growth period

                   Gt<-Fy*0
                   Gt[Fy>=P$Fs]<-(Fy[Fy>=P$Fs]-P$Fs)/P$g
                   Gt[Fy<=FgFromY0]<-(Fy[Fy<=FgFromY0]+FendYr)/P$g
                   Gt[Fy>FgFromY0 & Fy<P$Fs]<-1
                   Gt[Gt>1]<-1

                   # Length at beginning of growth period and growth increment in year
                   Index<-(Yplus>1 & Yplus<length(P$LatAge))
                   L0<-rep(NA,length(Fy))
                   Ldiff<-L0
                   L0[Index]<-P$LatAge[(Yplus[Index]-1)]
                   Ldiff[Index]<-P$LatAge[(Yplus[Index])]-L0[Index]
                   if(FgFromY0>0){
                     L0[Index][Fy[Index]>=P$Fs]<-P$LatAge[Yplus[Index]][Fy[Index]>=P$Fs]
                     Ldiff[Index][Fy[Index]>=P$Fs]<-P$LatAge[(Yplus[Index]+1)][Fy[Index]>=P$Fs]-L0[Index][Fy[Index]>=P$Fs]
                     }

            # cumulative F distribution with 50,50 dof  pf(val,50,50)
    return (L0+pf((2*P$Beta*Gt),50,50)*Ldiff)
    }



###############################################################################
###############################################################################

# 19. source(paste(dir_fns,"Growth_Fn_SVB.r",sep=""))



###############################################################################

Growth.SVB<-  # von Bertalanffy growth function - seasonal - from Siegel 1986 (eqn 4 in Candy & Kawaguchi 2006)
           function(P     # parameters  Linf, K, T0, Theta0, Theta1
                   ,Age     # data for input to the model:  Age
                   ,Time  # start time
                   ,dTime # length of time interval
                   ){

                   A<-Age+Time+dTime

    return (P$Linf*(1-exp(-P$K*((A-P$T0)-P$Theta0/(2*pi)*sin(2*pi*(A-P$Theta1))))))
    }


###############################################################################
###############################################################################

# 20. source(paste(dir_fns,"Reprod_Fn_Eggs_max.r",sep=""))



###############################################################################

# Minimum reproductive output (size of gonads) determined as a proportion of Morris etal 1988
# based on Nicol etal 1995 results for 50mm krill i.e. mean of 1914 eggs for mean length of 50.22 mm

# DW is in grams - need to convert to milligrams before determining eggs

Reprod_EggsMax<-function(Length,LW){
            GonadDW(Length,LW)*25.8/0.6824 # determine number of eggs by converting gonad weight (mg DW) to energy and then dividing
            }                             # by energy content of an egg (see Nicol etal 1995


