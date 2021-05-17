###########################################################################
###########################################################################
# Growth - Data

###########################################################################
# 1. scenarios

DaysInMonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
CumDaysInMonth<-apply(matrix(c(1:12),nrow=1),2,function(m,DIM){sum(DIM[c(1:m)])},DaysInMonth)
CumPropYrAPeriodStart<-c(0,CumDaysInMonth[c(1:11)]/365)
# 31  59  90 120 151 181 212 243 273 304 334 365

refMonth<-matrix(c(1:12),nrow=1)


# SSRU Elephant Island (EPOC Area 8)
# SSRU South Georgia West (EPOC Area 20)

Area_n<-c(8,20)
Area_Area<-c(1,43080448642.3966)

Area_col<-c("black","red") # plotting colour
Area_lat<-c(-65,-55)
Area_names<-c("EI","SG-W")
Area_lty<-c(1,1)

###########################################################################
# 2. Sea surface temperature (oC)

  SST<-matrix(c( 8, 1, 1.343081858,0.386184037
               , 8, 2, 1.528896831,0.345323414
               , 8, 3, 1.088661308,0.393935477
               , 8, 4, 0.420484764,0.297858027
               , 8, 5,-0.241198962,0.350801156
               , 8, 6,-0.669309516,0.35996916
               , 8, 7,-0.903839814,0.376720732
               , 8, 8,-0.851951408,0.34998953
               , 8, 9,-0.784357636,0.275890476
               , 8,10,-0.596537484,0.191668006
               , 8,11,-0.030453754,0.287945069
               , 8,12, 0.614372895,0.389847948
               ,20, 1, 2.977688308,0.545736361
               ,20, 2, 3.674458423,0.510002269
               ,20, 3, 3.566845192,0.4973703
               ,20, 4, 3.04992636 ,0.466596541
               ,20, 5, 2.15436472 ,0.492427588
               ,20, 6, 1.344083196,0.414378222
               ,20, 7, 0.807527468,0.359595024
               ,20, 8, 0.450114131,0.379167007
               ,20, 9, 0.260947186,0.386936734
               ,20,10, 0.405072008,0.358666728
               ,20,11, 1.013203156,0.525696475
               ,20,12, 1.972289308,0.532591851
               ),ncol=4,byrow=TRUE)


#################################################
# 3. Food in day in the year


# SSRU Elephant Island (EPOC Area 8)
# SSRU South Georgia West new carbon measures


      
##  [8]     

# Chl: cols : 1=SSRU, 2 = month, 3 = mean, 4=SD

# for Area 20, these are actual abundances for the area (mg)).  Area 8 remain at densities mg.m-2

Chl<-matrix(c( 8, 1,0.38491479 ,0.051755452
               , 8, 2,0.36056937 ,0.11383729
               , 8, 3,0.28682625 ,0.059662066
               , 8, 4,0.237112967,0.0663785
               , 8, 5,0.177183111,0.071857447
               , 8, 6,0.0        ,0.0
               , 8, 7,0.0        ,0.0
               , 8, 8,0.138077989,0.010203458
               , 8, 9,0.16895978 ,0.020115181
               , 8,10,0.23142768 ,0.032753194
               , 8,11,0.33732674 ,0.105439134
               , 8,12,0.36952999 ,0.056923504
               ,20, 1,4.598e+16 ,0.051755452
               ,20, 2,5.137e+16 ,0.11383729
               ,20, 3,7.721e+16 ,0.059662066
               ,20, 4,4.155e+16,0.0663785
               ,20, 5,2.771e+15,0.071857447
               ,20, 6,0.000e+00        ,0.0
               ,20, 7,0.000e+00        ,0.0
               ,20, 8,4.120e+15,0.010203458
               ,20, 9,7.118e+15 ,0.020115181
               ,20,10,8.514e+15 ,0.032753194
               ,20,11,2.433e+16 ,0.105439134
               ,20,12,3.399e+16 ,0.056923504
               ),ncol=4,byrow=TRUE)

EnvRecPropYr<-1/12  # for calculating proportion of food available.


# note that the conversion to carbon from chlorophyll is often considered to be 50x.

FoodFromChl<-1  # carbon already calculate 

AreaEnv<-cbind(Chl[,1]
              ,Chl[,2]
              ,DaysInMonth[Chl[,2]]/365
              ,SST[,3]
              ,Chl[,3]*FoodFromChl
              ,CumPropYrAPeriodStart[Chl[,2]]
              ,CumDaysInMonth/365
              )
dimnames(AreaEnv)[[2]]<-c("Area","Month","dt","Temp","Food","CalPropYrStart","CalPropYrEnd") # $Env dimnames col 1 = "dt" (prop of year), 2 = "Temp" (oC), 3 = "Food" (mg), 4 = calendar prop of year at start of period


