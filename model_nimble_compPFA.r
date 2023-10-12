model.nimble <- paste("","

model.nimble <- nimbleCode({

# -------------------------------------------------------------------------
#                     ATLANTIC SALMON LIFE CYLE MODEL
# Etienne RIVOT & Maxime OLMOS
## @ Max Olmos
## Version 3/10/2020  
  
# -------------------------------------------------------------------------
",sep="")

model.nimble <- paste(model.nimble,
"

# -------------------------------------------------------------------------
# Guide to principal variable names (critical stages and transitions)
# -------------------------------------------------------------------------

# Index
# -----------------------------

# t = 1:n - years
# r = 1,...,N : stock units
# N = 24
# N.NAC = 6
# N.NEAC = 18


# N1[t,r]								Eggs
# N2[t,r]								Total smolt produced by the spawning deposition N1[t,r]
# N3[t,k,r]								Smolts of age k migrating year t
# N3.tot[t,r]							Total smolts migrating year t (sum over ages k)
# N4[t,r]								PFA
# N5[t,r] / N8[t,r]						Maturing / non maturing PFA

# theta3[t,r]							Survival N3.tot --> N4
# theta4[t,r]							Proportion maturing N4 --> (N5 / N8)

# both theta3 and theta4 are modelled as multivariate random walk in the logit scale

# Maturing PFA --> Spawners 1SW
# -----------------------------
# N5[t,r] --> N6[t,r]					Maturing PFA --> 1SW returns
# C5.NEAC...  C5.NAC...					Sequential catches on maturing fish from NEAC and NAC (see below)
# h5.NEAC...  h5.NEAC...	 			Associated harvest rates
#										Different for NEAC and NAC (see below)

# N6[t,r] --> N7[t,r]					Returns 1SW --> spawners 1SW
# Chw.1SW[t,r]	h.hw.1SW[t,r]			Homewater catches 1SW and harvest rates

# Non maturing PFA --> Spawners 2SW
# -----------------------------
# N8[t,r] --> N9[t,r]					Non maturing PFA --> 2SW returns 
# N8.1[t,r]								Escapement fishery before WG
# N8.2[t,r]								Escapement Greenland fishery
# C8.NEAC...  C8.NAC...					Sequential catches on maturing fish (see below)
#										Different for NEAC and NAC (see below)
# h8.NEAC...  h8.NEAC...	 			Associated harvest rates

# N9[t,r] --> N10[t,r]					Returns 2SW --> spawners 2SW
# Chw.2SW[t,r]	h.hw.2SW[t,r]			Homewater catches 2SW and harvest rates

# Catches, Harvest rates and escapement for sequential fisheries at sea
# Different for NEAC and NAC
# -----------------------------

# NEAC
# ----

# Maturing 
# C5.NEAC.1[t,r]	h5.NEAC.1[t,r]		Faroes 1SW mature

# Non maturing
# C8.NEAC.1[t,r]	h8.NEAC.1[t,r]		Faroes fishery 1SW non mature
# C8.2[t,r]			h8.2[t,r]			Greenland fishery on 2SW (common with NAC)
# C8.NEAC.3[t+1,r]	h8.NEAC.3[t,r]		Faroes fishery on 2SW

# NAC
# ---

# Maturing
# C5.NAC.1[t,r]		h5.NAC.1[t]			NFL fishery on 1SW mature (harvest rate homogeneous among SU)
# C5.NAC.2[t,r]		h5.NAC.2[t,r]		LB fishery on 1SW mature
# C5.NAC.3[t,r]		h5.NAC.3[t,r]		Saint-Pierre et Miquelon fishery on 1SW mature # deleted in CYCLE PFA

# Non maturing
# C8.NAC.1[t,r]		h8.NAC.1[t]		 	NFL fishery on 1SW non mature (harvest rate homogeneous among SU)
# C8.2[t,r]			h8.2[t,r]			Greenland on 2SW (common with NEAC)
# C8.NAC.3[t+1,r]	h8.NAC.3[t]			NFL fishery on 2SW (harvest rate homogeneous among SU)
# C8.NAC.4[t,r]		h8.NAC.4[t,r]		LB fishery on 2SW
# C8.NAC.5[t,r]		h8.NAC.5[t,r]		SPM fishery on 2SW


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# FIXED PARAMETERS AND TIGHT INFORMATIVE PRIORS
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# Precision of additional logNormal noise on some demographic transitions
# modeled using logNormal process noise with a variance arbitrarily fixed
# to a very low value corresponding to a coefficient of variation
# CV.dummy = 0.01 (in the Constants)

tau.dummy <- 1/log(CV.dummy*CV.dummy + 1) 

# Precision for LogNormal random noise on the average proportion of smolts ages
# CV.psm fixed to an arbitrary very low value CV.psm = 0.01 (in the Constants)
 
tau.psm <- 1/log(CV.psm*CV.psm + 1)

# Precision for logNormal obs error on Homewater catches
# Arbitrarily set to relatively low value in the data CV.hw = 0.05 (in the Constants)

tau.hw <- 1/log(CV.hw*CV.hw + 1)

# M = Monthly natural mortality rate M
# Considered constant after PFA
# Is applied for all stages between PFA and returns using duration deltat..

M <- E.M

# One can also use a tight informative prior 
# ~LogNormal with E.M and CV.M fixed in the data

# tau.log.M <- 1/log(CV.M*CV.M + 1)	
# E.log.M <- log(E.M) - 0.5/tau.log.M
# M ~ dlnorm(E.log.M,tau.log.M)

# Inter-annual stochasticity
# Eggs N1 --> total Smolts per cohorts N2
# CV.theta1 fixed to an arbitrarily value (in the Constants)
# Default value: CV.Theta1 = 0.4

tau.theta1 <- 1/log(CV.dummy*CV.dummy + 1)

",sep="")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#   Prior on post-smolts survival - theta3 : Smolt N3.tot[t] --> PFA N4[t+1]
#   and probability to mature the first year at sea - theta4 : N4[t] --> N5[t] + N8[t]
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if(Cov==T)
{
  model.nimble <- paste(model.nimble,"

# Time series modeled as MultiNormal with NxN variance-covariance matrix

# Note :
# Variance - covariance matrix is only estimated starting at year 6
# Year 1:5 - estimates are sensitive to prior on number of fish (smolts)
# not generated by the model
# --> Temporal variations modelled through independent random terms
#     Variance-covariance not estimated
# Year 6:n : Start of multivariate random walk


# Post-smolts survival theta3 : N3.tot[t] --> N4[t+1]
# --------------------------------------------------------

	# Year 1
for (r in 1:N)
{
  for (t in 1:5)
  {
    logit.theta3[t,r] ~ dnorm(0,1)
  }
}

# Years in 6:n

# Prior on the variance-covariance matrix
# Wishart prior on the precision (Omega and ddl fixed in the Constants
# Omega = identity(N) and ddl=N
tau.theta3[1:N, 1:N] ~ dwish(omega[1:N, 1:N],N)

# Multivariate random walk in the logit scale
for (t in 5:(n-1))
{
  logit.theta3[t+1, 1:N] ~ dmnorm(logit.theta3[t,1:N], tau.theta3[1:N,1:N])

}

# Back to natural scale
for (r in 1:N)
{
  for (t in 1:n)
  {
    logit(theta3[t,r]) <- logit.theta3[t,r]
  } 
}


# Probability early maturing PFA N4[t] --> N5[t] + N8[t]
# ------------------------------------------------------------

# Year 1
for (r in 1:N)
{
  for (t in 1:5)
  {
    logit.theta4[t,r] ~ dnorm(0,1)
  }
}

# Years in 6:n

# Prior on the variance-covariance matrix
# Wishart prior on the precision (Omega and ddl fixed in the Constants
# Omega = identity(N) and ddl=N
tau.theta4[1:N, 1:N] ~ dwish(omega[1:N, 1:N],N)

# Multivariate random walk in the logit scale
for (t in 5:(n-1))
{
  logit.theta4[t+1, 1:N] ~ dmnorm(logit.theta4[t,1:N], tau.theta4[1:N,1:N])
  
}

  # Back to natural scale
for (r in 1:N)
{
  for (t in 1:n)
  {
    logit(theta4[t,r]) <- logit.theta4[t,r]
  } 
}

",sep="")
}# end if

if(noCov==T)
{
  model.nimble<-paste(model.nimble," 
# Post-smolts survival theta3 : N3.tot[t] --> N4[t+1]
# --------------------------------------------------------
                      
                      # Year 1
                      for (r in 1:N)
                      {
                      for (t in 1:5)
                      {
                      logit.theta3[t,r] ~ dnorm(0,1)
                      }
                      
                      
                      # Years in 6:n


	sigma.theta3[r] ~ dunif(0,1)
	tau.theta3[r] <- 1/pow(sigma.theta3[r],2)
	                    

                      # Multivariate random walk in the logit scale
                      for (t in 5:(n-1))
                      {
                      logit.theta3[t+1, r] ~ dnorm(logit.theta3[t,r], tau.theta3[r])
                                            }
                      
                      # Back to natural scale

                      for (t in 1:n)
                      {
                      logit(theta3[t,r]) <- logit.theta3[t,r]
                      } 

                      }
                      
                      
                      # Probability early maturing PFA N4[t] --> N5[t] + N8[t]
                      # ------------------------------------------------------------
                      
                      # Year 1
                      for (r in 1:N)
                      {
                      for (t in 1:5)
                      {
                      logit.theta4[t,r] ~ dnorm(0,1)
                      }
                      }
                      
                      # Years in 6:n
                      
                      # Prior on the variance-covariance matrix
	for (r in 1:N)
	{
	sigma.theta4[r] ~ dunif(0,1)
	tau.theta4[r] <- 1/pow(sigma.theta4[r],2)
	
                      # Multivariate random walk in the logit scale
                      for (t in 5:(n-1))
                      {
                      logit.theta4[t+1, r] ~ dnorm(logit.theta4[t,r], tau.theta4[r])
                      }
                      
                      # Back to natural scale
                   
                      for (t in 1:n)
                      {
                      logit(theta4[t,r]) <- logit.theta4[t,r]
                      } 
                      }
                      

",sep="")
} # end if

model.nimble<-paste(model.nimble,"	
  
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# POPULATION DYNAMICS
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# N7, N10 = Number of spawners 1SW and 2SW respectively
# defined as = returns (N6 and N9, respect.) - Homewater catches
# ---------------------------------------------------------
	
	# Harvest rate - Homewater fishery
	# n is the last year of data for 2SW
	# then (n-1) is the last year for which returns are updated by data

	for (r in 1:N)
	{   
		for (t in 1:(n-1))
		{ 
		h.hw.1SW[t,r] ~ dbeta(1,2)
		}
		
		for (t in 1:n)
		{ 
		h.hw.2SW[t,r] ~ dbeta(1,2)
		}
    }
	
	# Harvest rate on delayed spawners
	# Null for all SU except r=23 (RU.KW)
	# Only r = 23 has a likelihood for those catches
	for (t in 1:(n-1))
	{ 
		h.hw.1SW.delSp[t,23] ~ dbeta(1,2)
		for (r in 1:22)
		{
		h.hw.1SW.delSp[t,r] <- 0
		}
		for (r in 24:N)
		{
		h.hw.1SW.delSp[t,r] <- 0
		}	
	}
	
	for (t in 1:n)
	{ 
		h.hw.2SW.delSp[t,23] ~ dbeta(1,2)
		for (r in 1:22)
		{
		h.hw.2SW.delSp[t,r] <- 0
		}
		for (r in 24:N)
		{
		h.hw.2SW.delSp[t,r] <- 0
		}	
	}
	
	# Harvest rate for supplementary mortality in freshwater: only Scotland
	# Null for all SUs
	# except for EA_SC and WE_SC
	
		for (t in 1:(n-1))
		{
	  h.hw.sc.mort.1SW[t,12] ~ dbeta(1,2)
	  h.hw.sc.mort.1SW[t,13] ~ dbeta(1,2)
	  for (r in 1:11)
	 {
	   h.hw.sc.mort.1SW[t,r] <- 0
	}	
	 for (r in 14:N)
	 {
	  h.hw.sc.mort.1SW[t,r] <- 0
	 }	
	}
	
		for (t in 1:(n))
		{
	  h.hw.sc.mort.2SW[t,12] ~ dbeta(1,2)
	  h.hw.sc.mort.2SW[t,13] ~ dbeta(1,2)
	  for (r in 1:11)
	 {
	  h.hw.sc.mort.2SW[t,r] <- 0
	}	
	for (r in 14:N)
	{
	   h.hw.sc.mort.2SW[t,r] <- 0
	 }	
	}
	
  
	# Catches and escapement
	# -------------------------------
 	
	for (r in 1:N)
	{
	
	# 1SW
	# --------
	# Spawners (escapement)
	# delspawn1SW is the proportion of fish that returns to homewated but spawns the year after (delayed spawners)
	# must be considered AFTER the HW catches (see Run Reconstruction model NEAC, lines 3470)
	# (1-delspawn1SW)*(N6[t] - HWcatches[t]) = number of fish N6[t] that will spawn year t
	# delspawn1SW*(N6[t-1]-HWcatches[t-1]) = number of returns of the previous year (N6[t-1]) that will spawn year t (delayed)
	# Those delayed spawners are subject to additional catches (only for r = 23 = RU.KW)
	
	# First year
	#N7[1,r] <-  N6[1,r]*(1- h.hw.1SW[1,r]) 
	N7[1,r] <- scot.pct.mort[1,r]*N6[1,r]*(1- h.hw.1SW[1,r])*(1-h.hw.sc.mort.1SW[1,r])
    
	for (t in 2:(n-1))
    { 
   # N7[t,r] <- 	N6[t,r]*(1- h.hw.1SW[t,r])*(1 - prop.delSp.1SW[t,r]) +  N6[t-1,r]*(1- h.hw.1SW[t-1,r])*prop.delSp.1SW[t-1,r]*(1-h.hw.1SW.delSp[t,r])
	 N7[t,r] <-scot.pct.mort[t,r]*(N6[t,r]*(1- h.hw.1SW[t,r])*(1-h.hw.sc.mort.1SW[t,r])*(1 - prop.delSp.1SW[t,r]) +  N6[t-1,r]*(1- h.hw.1SW[t-1,r])*prop.delSp.1SW[t-1,r]*(1-h.hw.1SW.delSp[t,r]))
	  }
	
	# Homewater catches (likelihood until t = (n-1) only) 
    for (t in 1:(n-1)) 
    {
    Chw.1SW[t,r] <- h.hw.1SW[t,r] * N6[t,r]
	  }
	
	# Catches on delayed spawners
	# Will be 0 for all r except r = 23 = RU.KW
	Chw.1SW.delSp[1,r] <- 1		# Dummy - No delayed spawners defined for the first year
	for (t in 2:(n-1)) 
    {
	Chw.1SW.delSp[t,r] <- N6[t-1,r]*(1- h.hw.1SW[t-1,r])*prop.delSp.1SW[t-1,r]*h.hw.1SW.delSp[t,r]
	}
	
	# Catches on fish from scotland : freshwater mortality + EW catches specific on Scottish fish
	for (t in 1:(n-1))
	{
	  Chw.SC.killed1SW.m[t,r] <-scot.pct.mort[t,r] *N6[t,r]*(1- h.hw.1SW[t,r])*h.hw.sc.mort.1SW[t,r] 	
	}
	
	# 2SW
	# --------
	
	# Spawners (escapement)
	# Idem as 1SW
	# Additional stocking for USA (all 0 except r=6=USA)
    
	# First year
	#N10[1,r] <- N9[1,r]*(1-h.hw.2SW[1,r])
	N10[1,r] <- scot.pct.mort[1,r]*N9[1,r]*(1-h.hw.2SW[1,r])*(1-h.hw.sc.mort.2SW[1,r])
	
	for(t in 2:n)
    {   	
    #N10[t,r] <- N9[t,r]*(1- h.hw.2SW[t,r])*(1 - prop.delSp.2SW[t,r]) + N9[t-1,r]*(1- h.hw.2SW[t-1,r])*prop.delSp.2SW[t-1,r]*(1-h.hw.2SW.delSp[t,r]) + Stocking.2SW[t,r]
	N10[t,r] <- scot.pct.mort[t,r]*(N9[t,r]*(1- h.hw.2SW[t,r])*(1 - prop.delSp.2SW[t,r])*(1-h.hw.sc.mort.2SW[t,r]) + N9[t-1,r]*(1- h.hw.2SW[t-1,r])*prop.delSp.2SW[t-1,r]*(1-h.hw.2SW.delSp[t,r])) + Stocking.2SW[t,r]
	  }
	
	
    # Homewater catches (likelihood until t = n)
    for (t in 1:n)
    {
    Chw.2SW[t,r] <- h.hw.2SW[t,r]*N9[t,r]
    }
    
	# Catches on delayed spawners 2SW
	# Will be 0 for all r except r = 23 = RU.KW
	Chw.2SW.delSp[1,r] <- 1		# Dummy - No delayed spawners defined for the first year
	for (t in 2:n) 
    {
	Chw.2SW.delSp[t,r] <- N9[t-1,r]*(1- h.hw.2SW[t-1,r])*prop.delSp.2SW[t-1,r]*h.hw.2SW.delSp[t,r]
	}
	
   	
	# Catches on fish from scotland : freshwater mortality + EW catches specific on Scottish fish
	for (t in 1:(n))
	{
	 Chw.SC.killed2SW.m[t,r] <- scot.pct.mort[t,r]*N9[t,r]*(1- h.hw.2SW[t,r])*h.hw.sc.mort.2SW[t,r]
	}
	
	}


# N1: Number of eggs from spawners --> N2: total Smolts per cohort
# ---------------------------------------------------------

    for (r in 1:N)
    {  
		for (t in 1:(n-1))
		{
		# N1: nb eggs
		# -----------------
		N1[t,r] <- N7[t,r]*eggs[1,r,t] + N10[t,r]*eggs[2,r,t]
		
		N2[t,r] <- E.theta1*N1[t,r]

		# N1 --> N2 (Smolts)
		# -----------------
		# Mean survival known and fixed to E.theta1 (Constants)
		# CV.theta1 fixed in the data set (Constants)
		
		# Without density dependence
		#log.N2.m[t,r] <- log(a*N1[t,r]) - 0.5/tau.theta1
		#log.N2.m[t,r] <- log(N1[t,r]) - 0.5/tau.theta1
		
		# With density dependence
		# theta1.ddp[t,r] <- a/(1+B[r]*N1[t,r])
		# log.N2.m[t,r] <- log(theta1.ddp[t,r]*N1[t,r]) - 0.5/tau.theta1
		
		#N2[t,r] ~ dlnorm(log.N2.m[t,r],tau.theta1)
		
		# Surv.eggs is for monitoring
		#Surv.eggs[t,r] <- N2[t,r]/N1[t,r]
		}
	}
    

# N3: Smolts distribution by age class (6 age classes)
# ---------------------------------------------------------

# Dirichlet Informative prior
# Information equivalent to the one gained with a sample size = N.Sample.sm
# N.Sample.sm fixed to 100 in the Constants
# To improve computational speed, the Dirichlet is written with Gamma
    
    for (r in 1:N) 
    {
	
	# Proportion of smolts ages
	# k = smolt ages from 1 to nSm=6
	for (k in 1:nSm)
	{
	# +1 is needed to avoid mu.psm = 0 (would eventually crash the gamma()
	mu.psm[k,r] <- p.smolt[k,r]*N.Sample.sm + 1
        
        for (t in 1:n)
        {
		# psm.stoch[t,1:nSm,r] ~ ddirich(mu.psm[1:nSm,r])
        prop_gamma[t,k,r] ~ dgamma(mu.psm[k,r],1)
        psm.stoch[t,k,r] <- prop_gamma[t,k,r]/sum(prop_gamma[t,1:nSm,r])
        }
	}
      
	# N3
	for (t in 1:(n-1))
	{
		for (k in 1:nSm)
		{
		N3[t+1+k,k,r] <- psm.stoch[t,k,r]*N2[t,r]
		}

	# N3 tot : Total smolt migrating each year t
	N3.tot[t,r] <- sum(N3[t,1:nSm,r])
	}  
	
	# Dummy - not used in the model - no effect on the model
	# Just to be sure the array N3 has values in all cells (no cells with NA)
	
	for (s in 1:(nSm-1))
	{
		for (k in 1:s)
		{
		N3[n+1+s,k,r] <- 99
		}
	}	
	
	}

	
# N4 : Smolt --> PFA (N4)
# From N3.tot and survival theta3
# ---------------------------------------------------------

	for (r in 1:N)
	{
		for (t in 1:(n-1))
		{
        log.N4.m[t+1,r] <- log(theta3[t,r]*N3.tot[t,r]) - 0.5/tau.dummy
        N4[t+1,r] ~ dlnorm(log.N4.m[t+1,r], tau.dummy)
		}
	}
	  
# N5 : PFA maturing during the first year at sea
# From N4 (PFA) and theta4 = proba. maturing the first year at sea
# ---------------------------------------------------------
for (r in 1:N)
{
		for (t in 1:n)
		{
		log.N5.m[t,r] <-  log(N4[t,r] * theta4[t,r]) - 0.5/tau.dummy
		N5[t,r] ~ dlnorm(log.N5.m[t,r], tau.dummy)

# N8 : PFA non maturing after the first year at sea
# From N4 (PFA) and (1-theta4)
# ---------------------------------------------------------
        
		log.N8.m[t,r] <- log(N4[t,r] * (1-theta4[t,r])) - 0.5/tau.dummy
		N8[t,r] ~ dlnorm(log.N8.m[t,r], tau.dummy)
		}
	}


# N5 -> N6
# MATURING FISH returning in homewater as 1SW fish
# ---------------------------------------------------------	

# NEAC - Faroes fisheries and returns to N6  
# -----------------------------------------

    # Natural mortality between sequential fisheries
	# PFA -> Faroes m
	theta5.1.NEAC <- exp(-M*deltat.5.1.NEAC)

	# Faroes  -> returns
	theta5.2.NEAC <- exp(-M*(deltat.5.2.NEAC))

	# Exploitation rate - Faroes 1SWm
	# Separate for each SU
	for (r in 1:N.NEAC)
	{   
		for(t in 1:(n-1))
		{ 
		h5.NEAC.1[t,r] ~ dbeta(1,2)
		}
    }
	
	# Catches Faroes 1SWm
	for (t in 1:(n-1))
	{
		for (r in 1:N.NEAC)
		{   
		# Catches
		C5.NEAC.1[t,r] <-  theta5.1.NEAC * h5.NEAC.1[t,r] * N5[t,(r+N.NAC)]
		
		# Proportion to allocate the catches transformed as 
		# parameter for the Dirichlet likelihood
		mu.F1.NEAC.m[t,r] <- (C5.NEAC.1[t,r]/C5.NEAC.1.tot[t])*N.Sample[t]
		}
	
		# Total Catches Faroes 1SW maturing
		C5.NEAC.1.tot[t] <- sum(C5.NEAC.1[t,1:N.NEAC])
	}
	
	#  N6 : Returns 1SW NEAC
    for (r in 1:N.NEAC)
	{
		for(t in 1:(n-1))
		{
		N6[t,r+N.NAC]  <- theta5.2.NEAC * theta5.1.NEAC * (1-h5.NEAC.1[t,r]) * N5[t,r+N.NAC]
		}
	}
    
	
# NAC - LB/NF and SPM fisheries and returns to N6
# ------------------------------------------------

	# Natural mortality between sequential fisheries
	# PFA to Labrador/NFDL 1SWm
	theta5.1.NAC <- exp(-deltat.5.1.NAC*M)
	# Labrador/NFDL -> SPM
	theta5.2.NAC <- exp(-deltat.5.2.NAC*M)

	
	# Prior for harvest rates 

	for (t in 1:(n-1))
	{
		# NFL fishery on 1SW mature (harvest rate homogeneous among SU)
		# h5.NAC.1[t]		
		h5.NAC.1[t] ~ dbeta(1,2)
		
		# LB fishery on 1SW mature
		# h5.NAC.2[t,r]		
		# h homogeneous among all r except separate h for Labrador (r=1)
		h5.NAC.2.other[t] ~ dbeta(1,2)
		h5.NAC.2.lab[t] ~ dbeta(1,2)
		for(r in 2:N.NAC)
		{
		h5.NAC.2[t,r] <- h5.NAC.2.other[t]
		}		
		# Separate h for Labrador (r=1)
		h5.NAC.2[t,1] <- h5.NAC.2.lab[t]
		
		# Saint-Pierre et Miquelon fishery on 1SW mature
		# h5.NAC.3[t,r]		
		# h homogeneous among all r except h=0 for Labrador (r=1)
		#h5.NAC.3.other[t] ~ dbeta(1,2)
		#for(r in 2:N.NAC)
		#{
		#h5.NAC.3[t,r] <- h5.NAC.3.other[t]
		#}
		# Separate h for Labrador (r=1)		
		#h5.NAC.3[t,1] <- 0
	}
	
	# Sequential catches at sea 

	for (t in 1:(n-1)) 
	{	
		for (r in 1:N.NAC)
		{
			# PFA --> NFL fishery on 1SW mature
			C5.NAC.1[t,r] <- h5.NAC.1[t] * theta5.1.NAC * N5[t,r] 

			# LB fishery on 1SW mature
			C5.NAC.2[t,r] <- h5.NAC.2[t,r] * (1-h5.NAC.1[t]) * theta5.1.NAC * N5[t,r]
		}	
	
	# Total catches
	C5.NAC.1.tot[t] <- sum(C5.NAC.1[t,1:N.NAC])
		
	C5.NAC.2.other[t] <- sum(C5.NAC.2[t,2:N.NAC])
	C5.NAC.2.lab[t] <- C5.NAC.2[t,1]
		

	}


	# N6: Returns 1SW NAC
	# Survival after sequential catches + natural mortality theta5.1

	for (r in 1:N.NAC)
	{
		for (t in 1:(n-1))
		{
		N6[t,r] <- theta5.2.NAC * (1-h5.NAC.2[t,r]) * (1-h5.NAC.1[t]) * theta5.1.NAC * N5[t,r]
		}
	}			
	

	
#  NON  MATURING FISH
#  N8 -> N8.1 (escapement before Greenland Fishery)
# ---------------------------------------------------------	

# NEAC 

# Survival rate theta8.1.NEAC for the transition PFA --> Faroes
# followed by Faroes fisheries on non mature fish 1SWnm
# Harvest rate h8.SNEAC.1, variable across r
# ------------------------------------------	
    
    # Survival rate before Faroes fishery
    theta8.1.NEAC <- exp(-M*(deltat.8.1.NEAC))
    
    # Prior on exploitation rate
	# Separate for all regions r
    for (t in 1:(n-1))
    {
		for (r in 1:N.NEAC)
		{
		h8.NEAC.1[t,r] ~ dbeta(1,2)
		}
    }
     
    # Catches Faroes 1SW non maturing
    
    for (t in 1:(n-1))
    {	
		for (r in 1:N.NEAC)
		{
		# Catches
		C8.NEAC.1[t,r] <- h8.NEAC.1[t,r] * theta8.1.NEAC * N8[t,r+N.NAC]
		
		# Proportion to allocate the catches transformed as 
		# parameters for the Dirichlet likelihood
		mu.F1.NEAC.nm[t,r] <- (C8.NEAC.1[t,r]/C8.NEAC.1.tot[t]) * N.Sample[t]
		}
      
    # Total catches
    C8.NEAC.1.tot[t] <- sum(C8.NEAC.1[t,1:N.NEAC])	
    }
		
	# Transition N8 --> N8.1
    for (t in 1:(n-1))
	{  
		for (r in 1:N.NEAC)
		{
		N8.1[t,r+N.NAC] <- theta8.1.NEAC * (1-h8.NEAC.1[t,r]) * N8[t,r+N.NAC]
		}
    }
	
    
# NAC

# Survival theta8.1.NAC for the transition PFA --> LBandNF
# followed by LBandNF fisheries on 1SW non maturing 
# Harvest rate h8.NAC.1, homogeneous across r
# ---------------------------------------------

	# Survival rate
	
	theta8.1.NAC <- exp(-M * deltat.8.1.NAC)

	# Prior on exploitation rate
	# h homogeneous across the 6 regions
	
	for (t in 1:(n-1))
	{
		h8.NAC.1[t] ~ dbeta(1,2)
	}
	
	# Catches 1SW non maturing
	
	for (t in 1:(n-1))
	{	
		for (r in 1:N.NAC)
		{
		C8.NAC.1[t,r] <- h8.NAC.1[t] * theta8.1.NAC * N8[t,r] 
		}
		
		# Total catches
		C8.NAC.1.tot[t] <- sum(C8.NAC.1[t,1:N.NAC])
	}

	# Transitions N8 --> N8.1
	
	for (t in 1:(n-1))
	{	
		for (r in 1:N.NAC)
		{
		N8.1[t,r] <- (1-h8.NAC.1[t]) * theta8.1.NAC	* N8[t,r]	
		}
	}


# N8.1 --> N8.2
# Survival theta8.2[r] for NAC and NEAC
# Followed by Greenland Fisheries on NAC + NEAC stocks
# Harvest rate h8.2
# ---------------------------------------------------------		
	
	# Survival rate before Greenland Fishery - theta8.2[r]
	# calculated with deltat.8.2.NAC and deltat.8.2.NEAC
	
	for (r in 1:N.NAC)
	{
	theta8.2[r] <- 	exp(-M * (deltat.8.2.NAC))
	}
	
	for (r in 1:N.NEAC)
	{
	theta8.2[N.NAC+r] <- exp(-M * (deltat.8.2.NEAC))
	}

	
	# Proportion to allocate the catches between NEAC and NAC complexes
	
	for (t in 1:(n-1))
	{
	  for (k in 1:2)
	  {
	    mu.Gld.comp[t,k] <- (C8.2.comp[t,k]/C8.2.tot[t])*N.Sample[t]
	  }
	  
	  #Total catches
	  C8.2.tot[t] <- sum(C8.2.comp[t,1:2])	
	}

	# Exploitation rates at West Greenland Fishery
	# harvest rates all different for SUs NAC + NEAC
	
	  # NAC
	  for(t in 1:(n-1))
	  { 	
	    for (r in 1:N)
	    {
	    h8.2[t,r] ~ dbeta(1,2)
	    }
	  }
	
	for (t in 1:(n-1))
	{
	  for (r in 1:N.NAC)
	  {
	    mu.Gld.NAC[t,r] <- (C8.2[t,r]/(C8.2.comp[t,2]))*N.Sample[t]
	  }
	  
	  for (r in 1:N.NEAC)
	  {
	    mu.Gld.NEAC[t,r] <- (C8.2[t,(r+N.NAC)]/(C8.2.comp[t,1]))*N.Sample[t]
	  }
	}
	
	
	# Catches Greenland per SU
	for (t in 1:(n-1))
	{	
	  for (r in 1:N)
	  {
	    C8.2[t,r] <- h8.2[t,r] * N8.1[t,r]*theta8.2[r]
	  }
	  
	  # Total catches per complexe
	  C8.2.comp[t,2] <- sum(C8.2[t,1:N.NAC])	
	  C8.2.comp[t,1] <- sum(C8.2[t,(1+N.NAC):N])	
	}
	    # Escapement W Greenland fishery N8.1 --> N8.2
	    for (t in 1:(n-1))
	    {
	      for (r in 1:N)
	      {
	        N8.2[t,r] <-  theta8.2[r] * (1-h8.2[t,r]) * N8.1[t,r]	
	      }
	    }
	    
	
# AFTER Greenland Fisheries
# From N8.2 to returns as 2SW fish N9
# ---------------------------------------------------------	
 
# NEAC
# Faroes 2SW fisheries and return
# ---------------------------------------------------------

    # Survival rate 
	# Greenland -> Faroes 2SW
    theta8.2.1.NEAC <- exp(-M * (deltat.8.2.1.NEAC))
    # Faroes -> Returns
    theta8.2.2.NEAC <- exp(-M * (deltat.8.2.2.NEAC))

   
    # Exploitation rate different for all stock units
	# t = 1 are not used 
    for (r in 1:N.NEAC)
    {
		for(t in 1:n)
		{
		h8.NEAC.3[t,r] ~ dbeta(1,2)
		}
    }
		
    # Catches
	for (r in 1:N.NEAC)
    {
		for(t in 1:(n-1))
		{		 
		C8.NEAC.3[t+1,r] <- h8.NEAC.3[t+1,r] * theta8.2.1.NEAC * N8.2[t,r+N.NAC]
		}
    }
    
	# Total catches
	for(t in 1:(n-1))
    {
    C8.NEAC.3.tot[t+1] <- sum(C8.NEAC.3[t+1,1:N.NEAC])
    }
    
	# Proportion to allocate the catches transformed as 
	# parameters for the Dirichlet likelihood
    # Starts only at year t=2 
	for(t in 1:(n-1))
    {
		for ( r in 1:N.NEAC)
		{
		mu.F2.NEAC[t+1,r] <- (C8.NEAC.3[t+1,r]/C8.NEAC.3.tot[t+1])*N.Sample[t+1]
		}
    }
	
    #  Returns 2SW NEAC
	for (r in 1: N.NEAC) 
    {
		for(t in 1:(n-1))
		{
		N9[t+1,r+N.NAC] <-  theta8.2.2.NEAC * theta8.2.1.NEAC * (1-h8.NEAC.3[t+1,r]) * N8.2[t,r+N.NAC]
		}
    }
    
	
	# Dummy - not used - just to avoid NA in array
	for (r in 1:N.NEAC)
    {
    C8.NEAC.3[1,r] <- 99
	mu.F2.NEAC[1,r] <- C8.NEAC.3[1,r]/C8.NEAC.3.tot[1]
	}
	C8.NEAC.3.tot[1] <- sum(C8.NEAC.3[1,1:N.NEAC])

	
# NAC
# Sequential fisheries LB/NF and SPM
# ----------------------------------------------

	# Greenland -> Labrador/NF 2SW
	theta8.2.1.NAC <- exp(-M * deltat.8.2.1.NAC)
	# Labrador -> SPM 2SW
	theta8.2.2.NAC <- exp(-M * deltat.8.2.2.NAC)
	# SPM -> returns 2SW
	#theta8.2.3.NAC <- exp(-M * deltat.8.2.3.NAC)

	
	# NF on 2SW
	# h8.NAC.3	
	# Homogeneous across the 6 regions
	for (t in 1:n)
	{
	h8.NAC.3[t] ~ dbeta(1,2)
	}

  	# Catches
	for (t in 1:(n-1))
	{
		for (r in 1:N.NAC)
		{
		C8.NAC.3[t+1,r] <- h8.NAC.3[t+1] * theta8.2.1.NAC * N8.2[t,r]
		}
		C8.NAC.3.tot[t+1] <- sum(C8.NAC.3[t+1,1:N.NAC])
	}

	# LB on 2SW
	# h8.NAC.4
	# Homogeneous across the 5 regions except Labrador (r=1)
	# h.NAC.4 for year t=1 not used	
	
	for (t in 1:n)
	{
	h8.NAC.4.other[t] ~ dbeta(1,2)
	h8.NAC.4.lab[t] ~ dbeta(1,2)
		# All regions except Labrador
		for(r in 2:N.NAC)
		{
		h8.NAC.4[t,r] <- h8.NAC.4.other[t]
		}
		# Separate for Lab r =1
		h8.NAC.4[t,1] <- h8.NAC.4.lab[t]
	}
	
	for (t in 1:(n-1))
	{
		for (r in 1:N.NAC)
		{
		C8.NAC.4[t+1,r] <- h8.NAC.4[t+1,r] * (1-h8.NAC.3[t+1]) * theta8.2.1.NAC * N8.2[t,r]
		}
		C8.NAC.4.other[t+1] <- sum(C8.NAC.4[t+1,2:N.NAC])
		C8.NAC.4.lab[t+1] <- C8.NAC.4[t+1,1]
	}
	

	
	# Catches

	#  N9 : Returns 2SW NAC
	for (r in 1:N.NAC)
	{
		for (t in 1:(n-1))
		{
		N9[t+1,r] <-  theta8.2.2.NAC * (1-h8.NAC.4[t+1,r]) * (1-h8.NAC.3[t+1]) * theta8.2.1.NAC * N8.2[t,r]			
		}
	}	

	
	# Dummy - not used in the model - no effect on the model - just to avoid NA in array
	for (r in 1:N.NAC)
	{
		C8.NAC.3[1,r] <- 99
		C8.NAC.4[1,r] <- 99
		C8.NAC.5[1,r] <- 99 
	}
	C8.NAC.3.tot[1] <- sum(C8.NAC.3[1,1:N.NAC])
	C8.NAC.4.other[1] <- sum(C8.NAC.4[1,2:N.NAC])
	C8.NAC.4.lab[1] <- C8.NAC.4[1,1]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# OBSERVATION EQUATIONS 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# Likelihood for returns
# ---------------------------------------------------------	
    
    # 1 SW
    for (t in 1:(n-1))
    {
		for (r in 1:N)
		{
		log.R1SW.m[t,r] ~ dnorm(log(N6[t,r]),log.R1SW.tau[t,r])
		}	
    }
    
    # 2 SW
    for (t in 1:n)
    {		
		for (r in 1:N)
		{
		log.R2SW.m[t,r] ~ dnorm(log(N9[t,r]),log.R2SW.tau[t,r])
		}
    }
 
 
# Likelihood for Homewater catches
# ------------------------------------------------
    # Regular homewater catches
	# --------------------------------------------
	
    # 1 SW (up to t=n-1)
	# ------------------
    for(t in 1:(n-1))
    {
		for (r in 1:N)
		{
		log.hwC1SW.m[t,r] ~ dnorm(log(Chw.1SW[t,r]),tau.hw)
		}
    }
    
    # 2 SW (up to t=n)
	# -----------------
    for(t in 1:n)
    {
		for (r in 1:N)
		{
		log.hwC2SW.m[t,r] ~ dnorm(log(Chw.2SW[t,r]),tau.hw)
		}
    }   

	# Catches on delayed spawners
	# --------------------------------------------
	# Only needed for r = 23 because harvest rates fixed to 0 for all other stock units
	
	# 1SW (up to t=(n-1)
	for (t in 1:(n-1)) 
    {
		log.Chw.1SW.delSp.m[t,23] ~ dnorm(log(Chw.1SW.delSp[t,23]),tau.hw)
	}
	# 2SW (up to t=n) 
	for (t in 1:n) 
    {
		log.Chw.2SW.delSp.m[t,23] ~ dnorm(log(Chw.2SW.delSp[t,23]),tau.hw)
	}
	
	# Catches on scotish returns (natural mortality + cathces)
	# --------------------------------------------
	# Only needed for r = 23 because harvest rates fixed to 0 for all other stock units
	
	# 1SW (up to t=(n-1)
		for (t in 1:(n-1)) 
		{
	  log.SC.killed1SW.m[t,12] ~ dnorm(log(Chw.SC.killed1SW.m[t,12]),tau.hw)
	  log.SC.killed1SW.m[t,13] ~ dnorm(log(Chw.SC.killed1SW.m[t,13]),tau.hw)
	}
	 #2SW (up to t=n) 
	for (t in 1:n) 
	{
	 log.SC.killedMSW.m[t,12] ~ dnorm(log(Chw.SC.killed2SW.m[t,12]),tau.hw)
	log.SC.killedMSW.m[t,13] ~ dnorm(log(Chw.SC.killed2SW.m[t,13]),tau.hw)
	 }
	
	
  
# Likelihood for total Catches at sea
# -----------------------------------------------------------------
		
	# West Greenland NAC + NEAC
	# ----------------------------------------------- 
	
	for (t in 1:(n-1))
	{
		log.CG2.m[t]  ~ dnorm(log(C8.2.tot[t]), log.CG2.tau[t])
	} 
	
    
    # NEAC
	# -----------------------------------------------
    
	for(t in 1:(n-1))
    {	
		# Faroes 1SW maturing
		log.CF1.m.m[t] ~ dnorm(log(C5.NEAC.1.tot[t]),log.CF1.m.tau[t])
	  
		# Faroes 1SW non maturing
		log.CF1.nm.m[t] ~ dnorm(log(C8.NEAC.1.tot[t]),log.CF1.nm.tau[t])
    }
    
    for(t in 2:n)
    {
		# Faroes 2SW
		log.CF2.m[t] ~ dnorm(log(C8.NEAC.3.tot[t]),log.CF2.tau[t])
    }
    	
	
    # NAC
	# -----------------------------------------------   

	# Mature fish
	# Likelihood up to t=(n-1)
	for (t in 1:(n-1))
	{
		log.C1.Nf.3_7.m[t] ~ dnorm(log(C5.NAC.1.tot[t]), log.C1.Nf.3_7.tau[t])
	
		log.C1.LbNf.other_spm.m[t] ~ dnorm(log(C5.NAC.2.other[t]), log.C1.LbNf.other_spm.tau[t])

		log.C1.tot.Lb.m[t] ~ dnorm(log(C5.NAC.2.lab[t]), log.C1.tot.Lb.tau[t])
		
	
		
	}
	
	# Non mature fish
	# Likelihood up to t=(n-1)
	for (t in 1:(n-1))
	{
		log.C1.nm.LbNf.m[t] ~ dnorm(log(C8.NAC.1.tot[t]), log.C1.nm.LbNf.tau[t])
	}
	
	# Likelihood starts at t=2
	# C8.NAC.3.tot[t=1] not defined
	for (t in 2:n)
	{
		log.C2.Nf.3_7.m[t] ~ dnorm(log(C8.NAC.3.tot[t]), log.C2.Nf.3_7.tau[t])
	}
	
	# Likelihood starts at t=2
	# C8.NAC.4.other[t=1] and C8.NAC.4.lab[t=1] not defined
	for (t in 2:n)
	{
		log.C2.LbNf.other_spm.m[t] ~ dnorm(log(C8.NAC.4.other[t]), log.C2.LbNf.other_spm.tau[t])
		log.C2.tot.Lb.m[t] ~ dnorm(log(C8.NAC.4.lab[t]), log.C2.tot.Lb.tau[t])
	}

 
# Likelihood for proportion of each stock unit in the catches
# -----------------------------------------------------------------
        
	for(t in 1:(n-1))
	{
		# 1SWm Faroes
		prop_F1.m[t,1:N.NEAC] ~ ddirch(mu.F1.NEAC.m[t,1:N.NEAC])
      
		# 1SWnm Faroes
		prop_F1.nm[t,1:N.NEAC] ~ ddirch(mu.F1.NEAC.nm[t,1:N.NEAC])
	  
		# West Greenland
		p.NEAC.WG[t,1:2] ~ ddirich(mu.Gld.comp[t,1:2])
		prop_NAC.WG[t,1:N.NAC] ~ ddirich(mu.Gld.NAC[t,1:N.NAC])
		prop_NEAC.WG[t,1:N.NEAC] ~ ddirich(mu.Gld.NEAC[t,1:N.NEAC])
	}	
    
    # 2SW Faroes catches
	for(t in 2:n)
	{
		prop_F2[t,1:N.NEAC] ~ ddirch(mu.F2.NEAC[t,1:N.NEAC])	
	}
    

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# INITIALIZATION OF THE LOOP ON t
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# All number at all stages that are not generated by equations must be initialized in the model 

	# N1    
	# Weakly informative prior for the number of eggs (specific to each region)
	# (the mean mu.N1.pr is calculated in the data block)
	# CV.N1.pr fixed in the Constant (0.8)
	# -------------------------------------------------
    
	tau.N1.pr <- 1/log(CV.N1.pr * CV.N1.pr + 1)

	for (r in 1:N)
	{  
		for(t in 1:(nSm+1))
		{
		log.N1.pr[t,r] <- mu.N1.pr[r] - 0.5/tau.N1.pr
		N1.pr[t,r] ~ dlnorm(log.N1.pr[t,r],tau.N1.pr)
		}
	}	
    
    # N2
	# -------------------------------------------------- 
	
	for (r in 1:N)
	{  
		for(t in 1:(nSm+1))
		{
		# No density dependence
		log.N2.m.pr[t,r] <- log(a*N1.pr[t,r]) - 0.5/tau.theta1
				
		# Density dependence 
		# theta1.ddp.pr[t,r] <- a/(1+B[r]*N1.pr[t,r])
		# log.N2.m.pr[t,r] <- log(theta1.ddp.pr[t,r]*N1.pr[t,r]) - 0.5/tau.theta1
		
		N2.pr[t,r] ~ dlnorm(log.N2.m.pr[t,r],tau.theta1)
		}
	}

    # N3
	# Initialization for all smolts that are not generated by the dynamics
	# The same informative Dirichlet distribution for the proportion of smolts age
	# is used than the one in the model
	# -------------------------------------------------- 
    
	for (r in 1:N)
    { 
		# Allocation of smolt ages for years 64 to 70
		# mu.psm independent of year and defined in the Constants
		for (k in 1:nSm)
		{
			for (t in 1:(nSm+1))
			{
			prop_gamma.pr[t,k,r] ~ dgamma(mu.psm[k,r],1)
			psm.stoch.pr[t,k,r] <- prop_gamma.pr[t,k,r]/sum(prop_gamma.pr[t,1:nSm,r])
			N3.pr[t+k+1,k,r] <- psm.stoch.pr[t,k,r] * N2.pr[t,r]
			}
		}

		# Year 71 corresponds to the 8th line in N3.pr (from 64 to 70)
		# The first year that has to be completed
	  	for(k in 1:nSm)
		{
		N3[1,k,r] <- N3.pr[8,k,r]
		}
	  
		# Filling years 72:77 from N3.pr of years 64 to 70
		for (k in 1:nSm)
		{
			for (kk in k:nSm)
			{
			N3[k+1,kk,r] <- N3.pr[k+8,kk,r]
			}
		}
    }
    
    
    # N4	
	# --------------------------------------------------
	
	for (r in 1:N)
    {
	# Will provide approximately a uniform on [0,1] for theta3
	logit.theta3.pr[r] ~ dnorm(0,0.5)
	logit(theta3.pr[r]) <- logit.theta3.pr[r]
	
	log.N4.m[1,r] <- log(theta3.pr[r]*sum(N3.pr[(2+nSm),1:nSm,r])) - 0.5/tau.dummy
	N4[1,r] ~ dlnorm(log.N4.m[1,r],tau.dummy)
    }

	
    # N9	
	# --------------------------------------------------
	
	for (r in 1:N)
    {
	tau.N9.pr[r] <- pow((max.log.N9[r] - min.log.N9[r])/2,-2)
	mean.log.N9.pr[r] <- (min.log.N9[r] + max.log.N9[r])/2 - 0.5/tau.N9.pr[r]
	log.N9.pr[r] ~ dnorm(mean.log.N9.pr[r], tau.N9.pr[r])
	N9[1,r] <- exp(log.N9.pr[r])
    }

",sep="")
    
model.nimble<-paste(model.nimble,"
}
)"
)

write(model.nimble,file("model.txt"))
