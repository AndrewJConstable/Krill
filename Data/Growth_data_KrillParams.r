EMC_p<-list(
   AssimilationEfficiency = 0.8  # 0.4E-2 #0.35E-3  # can use this parameter to adjust for krill density  i.e. Food per m2/Krill per m2  = ~  1/250 (look it up)
                                 # conversion parameters of body DW and gonad DW to be equivalent in units of investment
                                 # from ingested food
  ,Ingest                = list(Max=FALSE            # logical - should ingestion be limited at a maximum food density
                                ,MaxIngest = 1E10) # mgC - max instaneous food (mg)density per m3 (per 1/365 or per day)
  ,DWunits                = 1  # units are in mg DW
  ,Moult                  = list(BodyDwProp = 0.075 # the proportion of DW that is the moult (Nicol etal 1992)
                                ,PropC = 0.238) # proportion of the moult DW that is carbon (Ikeda 1984)
  ,GonadConversion        = 0.551 # proportion converts from GonadDW (mg) to mgC (from Ikeda 1984 for eggs)
  ,GonadDW_to_eggs        = 25.8/0.6824  # data from Nicol et al (1995)
  ,GonadDW_to_WW          = c(2.074666,1.0699)  # mg to mg : a.DW^b   estimated from Morris maxima in LW below
  ,MaturityByAge          = TRUE
  ,Birthday               = 1   # 12 am on 1 January (Note : do not set to 0)
  ,EndOfEggSeason         = 0.5 # pending new recruit model - ensure eggs do not get accounted in wrong year.  those later in year need to wait at least one year
  ,AgeMature              = 1.5
  ,LengthMature           = 0 # 39.28832  #  old estimate 36.173
  ,PropReprodBeforeGrowth = 0.95
  ,MinReprodProp          = 0.1 #0.52 # taken from Nicol etal 1995 as a proportion of the reproduction from Morris etal for the size of krill in Nicol's paper
  ,nMoultsRegress         = 5
  ,nMoultsMature          = 5
  ,DLreprod               = 5  # regression is entered if DeclineDayLength & Daylength< DLreprod and food is limited
  ,CritDayLength          = 2   # regression automatically entered after a moult if the day length is declining and the day length < CritDayLength
  ,Shrink                 = TRUE
  ,ShrinkParams           = list(slope = 1.418E-5
                                 ,int = 1.3E-3
                                   )
  ,Met = list(dt = 1/365  # duration of time step for integrating metabolism (years)
             ,Starving_Prop_Ideal_Cw = 0.99  # proportion of ideal carbon weight when starving considered to begin
             ,Hibernation_Prop_IdealCw = 0.8 # proportion when hibernation considered to begin
             ,Starving_Metabolism_prop = 0.33 #0.33 # proportion of ideal metabolism when starving - i.e. hibernation (Quetin & Ross 1991)
             ,Hibernation_Metabolism_prop = 0.1
             ,Male = list(HiRate = TRUE, Mature_Length = 35.0, Multiplier = 1.1)  # allowing for males to have higher energy expenditure than females - activity and gonads
             ) # proportion of ideal metabolism when hibernating - Meyer
  ,Mortality = list(Nat = 0.8
            ,Starve = list(IndexCrit = (-1.5), Exp = 2, Male_Specific = FALSE)  # index parameters for mortality rate arising from starvation
            ,Larvae = 5.0
            )
  ,LW = matrix(c(
     21.38E-3,2.76,0.195,(-14.0) # MA2   (a,b,d1,d2  in formula from paper) giving DW in mg from length in mm (note original parameters were for g from mm) thus parameters a & d2 were multiplied by 1000
    ,30.20E-3,2.62,0.251,(-13.0) # FA1
    , 9.77E-3,2.98,0.275,(-14.0) # FA4
  ),ncol=4,byrow=TRUE)
  )



# calculations of Gonad DW to WW conversion
# LW<-matrix(c(
#   21.38E-3,2.76,0.195,(-14.0) # MA2   (a,b,d1,d2  in formula from paper) giving DW in mg from length in mm (note original parameters were for g from mm) thus parameters a & d2 were multiplied by 1000
#   ,30.20E-3,2.62,0.251,(-13.0) # FA1
#   , 9.77E-3,2.98,0.275,(-14.0) # FA4
# ),ncol=4,byrow=TRUE)
# 
# tmp_L<-c(35:80)
# tmp_WW<- (LW[3,1]*tmp_L^LW[3,2])-(LW[2,1]*tmp_L^LW[2,2])
# tmp_DW<- (LW[3,3]*(LW[3,1]*tmp_L^LW[3,2])+LW[3,4])-(LW[2,3]*(LW[2,1]*tmp_L^LW[2,2])+LW[2,4])
# plot(log(tmp_DW),log(tmp_WW))
# lines(log(tmp_DW),(1.0699*log(tmp_DW)+0.7298))
# 
# tmp_dw<-c(0:100)/100*500
# plot(tmp_dw,(exp(0.7298)*tmp_dw^1.0699))
# 
# tmp_tst<-(tmp_WW-(2.074666*tmp_DW^1.0699))/tmp_WW
# mean(tmp_tst)
# sd(tmp_tst)
# 
# lm(log(tmp_WW)~log(tmp_DW))
# 
# plot(tmp_DW,(log(tmp_WW)/log(tmp_DW)))
# mean(tmp_WW/tmp_DW)
# sd(tmp_WW/tmp_DW)
# 
