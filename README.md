	PURPOSE
	The SOMNUS_NU_datapull_11Aug26 repository contains code used to create synthetic data to run analyses prespecified on clinicaltrials.gov. 
	This is to ensure that all required data elements and variables will be accounted for in our one-time Northwestern data pull. 
	
	DATA 
	1. Synthetic prescription and clinican-monthly data with prespecified fixed and random effect sizes for a two-part hurdle model with a binomial 
	   and poisson distribution 

	_ANALYSES_
	*PRIMARY*
	1. Two-stage hurdle model testing time x treatement reductions in Z-drug prescribing frequency (binomial) and 5-mg pill counts (poisson)
	2. Other analyses with primary outcome:
	   a. By user type (short- vs. long-term)
	   b. During intervention and post-intervention period
	      i. Both overall, and by user type 
	      
	*SECONDARY*
	1. 