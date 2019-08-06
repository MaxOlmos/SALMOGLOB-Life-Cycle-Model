## SALMOGLOB-Life-Cycle-Model

A life cycle model 

	1. considering spatial covariation of marine life history traits of Atlantic salmon populations 
	For a description of the model please, see : Olmos et al., (2019), Evidence for a spatial coherence in time trends of marine
	life history traits of Atlantic salmon populations in the North Atlantic
	> model.Fish&Fish.txt : is the JAGS code correponding to the baseline configuration of the life cycle model
 
      * density-independent egg-to-smolt survival relationship with a common homogeneous survival rate among SUs, constant over time and modelled with very low inter-annual stochastic variability (F1 configuration)
      
      * no a priori hypotheses about the covariation structure of 13×13 variance-covariance matrices ∑θ3 and ∑θ4 are made (M1 configuration)
      
      * stock assignement data based on genetics are used to assign the origin of the catches to SUs, which allows harvest rates to vary by SU (E1 configuration)

	
	and 
	
	2. investigating the environmental drivers and the demographic mechanisms of the widespread decline in marine survival rates of Atlantic salmon (Salmo salar) over the last four decades
		This model allows 
		a. the investigation of the degree of synchrony in Atlantic salmon post-smolt survival and explicitly quantifies the amount of variance that is captured by trends at various spatial scales. 
		b. the investigation whether the temporal variation in the marine survival can be explained by environmental variation encountered by salmon during the early post-smolt marine phase when salmon use specific transit habitat, or during the later phase of the first year at sea when salmon of different areas aggregate at common feeding areas.
		For a description of the model, please, see : Olmos et al., (20XX), Spatial synchrony in the response of a long range migratory species (Salmo salar) to climate change in the North Atlantic Ocean, currently submitted to Global Change Biology
		> model.enviro.txt : is the nimble code 
 
		* correponding to the  baseline configuration of the life cycle model in Olmos et al 2019 (Fish and Fisheries)
		* quantifying the influence of covariaties in the common space-time domains
