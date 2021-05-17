SizeAtAge_p<-list(
   VB     = list(Linf=60, K=0.45, T0=0)
  ,PVB    = list(Linf=60, K=0.45, T0=0, g = 0.25, Fs = 0.0)
  ,SVB_CK = list(Linf=53.89, K=0.6792, T0=0.43, g = 0.25, Fs = 0.0, Theta0 = 2.519, Theta1 = 0.094) # from Candy & Kawaguchi 2006 estimated
  ,SVB    = list(Linf=61, K=0.4728, T0=0.1418, g = 0.25, Fs = 0.0, Theta0 = 0.9598, Theta1 = 0.0272) # from Siegel 1986 we think
  ,PVB_F  = list(LatAge = c(28.0,41.77,48.95,53.55,57.53,61.06) # length at age where age is determined as the length in the no growth period
                ,Beta = 1.239
                ,g = 0.5808 # calculated from Candy & Kawaguchi eqn 6
                ,Fs = 0.8329)
   )
