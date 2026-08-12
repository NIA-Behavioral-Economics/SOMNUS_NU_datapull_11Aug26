**PURPOSE**  
The SOMNUS_NU_datapull_11Aug26 repository contains code to create synthetic data to run analyses prespecified on clinicaltrials.gov. 
This is to ensure that all required data elements and variables will be accounted for in a one-time data pull at Northwestern Medicine.

**DATA**  
Synthetic prescription and clinician-monthly visit data with prespecified fixed and random effect sizes for two-part hurdle model with binomial 
and Poisson distributions 

**ANALYSES**  
	*PRIMARY*  
1. Two-stage knotted spline hurdle model testing truncated kmonth x treatment reductions in Z-drug prescribing frequency (binomial) and 5-mg pill counts (Poisson)
2. Other analyses with primary outcome:
   - **a.** During intervention (study months 1 to 18) and post-intervention (study months 19 to 30)
      - **i.** Overall and by user type	 
           
	*SECONDARY*  
	1. Mixed logistic model testing month x treatment change in clinician-monthly CBT-I orders from intervention start to intervention end   
	   a. By user type (short- vs. long-term)  
	2. Two-stage knotted spline hurdle model testing truncated kmonth x treatment reductions in benzodiazepine prescribing frequency (binomial)  
	   and 2-mg pill equivalents (Poisson)  
	3. Mixed logistic model testing truncated kmonth x treatment reductions in discordant Z-drug prescribing (> 35 pills for < 35 days' supply)  
	4. Secondary outcomes also assessed in post-intervention period   
	
	*UNRESTRICTED*  
	1. Hurdle model for primary outcome will be tested using an unrestricted model which includes two-way treatment interactions 
	   (accountable justification x default) and three-way knotted time by treatment interactions (accountable justification x default x kmonth)
	   a. Unrestricted model will be tested separately in intervention and post-intervention periods
	   b. If interaction effects are insignificant, then we will proceed with restricted model (i.e., no treatment interactions) 

**LICENSE**<br>
Schaeffer Center for Health Policy and Economics, University Southern California
	   
**FILES**

  *drug_tables.R*
  
  i. Creates separate tables for each pharmaceutical class (e.g., Z-drugs, benzodiazepines) to derive drug strength for synthetic prescription data 
     a. zdrug_rx (prescription name, strength, route, etc. for Z-drugs)
     b. benzo_rx (same as above for benzodiazepines)
  
  *binomial_visits.R*
  
  i. Creates a flat file (binomial.rds) with aggregated visit counts for every study month per-clinician (444 clinicians * 48 study months = 21,312 rows)
 ii. Variable list:
 		 a. prov_id (clinician ID: 1-444)
 		 b. month (study month: baseline = 1-18, intervention = 19-36, post-intervention = 37-48)
 		 c. total_vsts (total visits per clinician-month: mean = 250, sd = 40)
 		 d. long_vsts (total visits per clinician-month for long-term Z-drug users (> 180 days of Z-drugs in past 365 days): 15% of total visits)
 		 e. short_vsts (total  visits per clinician-month for short-term Z-drug users (<= 180 days of Z-drugs in past 365 days): total_vsts - long_vsts)
 		 f. clinic_id (clinic ID: 1-64)
 		 g. Tx1 (accountable justification: 50% of clinics randomized)
 		 h. Tx2 (default: 50% of clinics randomized)
 		 i. cpres (clinic prescribing (high vs. low): 50% of clinic randomized to 'high' prescribing)
 		 j. mnth (month centered: baseline = -17-0, intervention = 1-18, post-intervention = 19-30)
 		 k. post (study period: baseline = 0, intervention = 1, post-intervention = 2)
 		 l. kmnthTx (intervention study month for knotted spline model: baseline = 0, intervention = 1-18, post-intervention = NA)
 		 m. kmnthFu (post-intervention study month for knotted spline model: baseline = 0, intervention = 0, post-intervention = 19-30)
 		 n. xb (fixed effects linear predictor for binomial Z-drug prescription outcome on log scale)
 		 o. xb_cbti (fixed effects linear predictor for binomial CBT-I outcome on log scale)
 		 p. clinic_re_logi (random clinic intercept for binomial outcomes on log scale: total variance = pi^2/3)
 		 q. provider_re_logi (random clinician intercept for binomial outcomes on log scale: total variable = pi^2/3)
 		 r. ETA_logi_short (linear predictor for binomial Z-drug prescription outcome on log scale including fixed and random effects for long-term users)
 		 s. ETA_logi_long (linear predictor for binomial Z-drug prescription outcome on log scale including fixed and random effects for short-term users)
 		 t. ETA_cbti_short (linear predictor for binomial CBT-I order outcome on log scale including fixed and random effects for short-term users)
 		 u. ETA_cbti_long (linear predictor for binomial CBT-I order outcome on log scale including fixed and random effects for long-term users)
 		 v. ETA_benzo (linear predictor for binomial benzodiazepine prescription outcome on log scale including fixed and random effects)
 		 w. zdrug_short_vsts (randomly generated number of visits for short-term Z-drug users based on probabilities from ETA_logi_short)
 		 x. zdrug_long_vsts (randomly generated number of visits for long-term Z-drug users based on probabilities from ETA_logi_long)
 		 y. benzo_vsts (randomly generated number of visits for benzodiazepine prescriptions based on probabilities from ETA_benzo)
 		 z. cbti_short_vsts (randomly generated number of visits for CBT-I orders for short-term users based on probabilities from ETA_cbti_short)
 	  aa. cbti_long_vsts (randomly generated number of visits for CBT-I orders for long-term users based on probabilities from ETA_cbti_long)
 	  ab. start_date (start date for study month)
 	  ac. end_date (end date for study month)
 	  
 	  *poisson_visits.R*
 	  
 	  i. Combines visit, clinician, and clinic data from binomial.rds with poisson prescription outcomes to create analytic datasets
 	     a. total.rds (analytic dataset for primary Z-drug outcome for all visits, n = 5,321,768)
 	        1. prov_id (clinician ID)
 	        2. clinic_id (clinic ID)
 	        3. Tx1 (accountable justification intervention)
 	        4. Tx2 (default intervention)
 	        5. post (study period)
 	        6. mnth (month centered)
 	        7. kmnthTx (intervention study month)
 	        8. kmnthFu (post-intervention study month)
 	        9. cpresc (high vs. low clinic prescribing)
 	       10. quantity (NA if rx = 0 (i.e., no prescription), otherwise pill quantity randomly generated using same methodology for binomial outcome above)
 	       11. dose (prescription info. merged from zdrug_rx, NA if rx = 0, otherwise destringed drug strength) 
 	       12. rx (Z-drug prescription: 0 = no, 1 = yes)
 	       13. cbti (CBT-I order where rx = 1: 0 = no, 1 = yes)
 	       14. rx_start_date (start date for prescription randomly generated from study month dates)
 	       15. rx_end_date (end date for prescription: start_date + days_supply)
 	       16. days_supply (randomly generated, assumes a minimum of 5 and a 0.30 decrease over time post-intervention for treatment groups)
 	       17. name (prescription name: NA if rx = 0, otherwise populated) 
 	       18. pills (5 mg Z-drug pill equivalents: (dose/5) x quantity, NA if rx = 0, else populated)
 	       19. discordant (Z-drug prescription not within guidelines: > 35 pills for 5-week days' supply, NA if rx = 0, otherwise populated)
 	     b. zdrug_long.rds (n = 798,285)
 	        1. same variables as total.rds on subset of patients that are long-term users (n = 798,285)
 	     c. zdrug_short.rds (4,523,483)
 	        1. same variables as total.rds on subset of patients that are short-term users (n = 4,523,483)
 	     d. benzo.rds (n = 5,321,768)
 	        1. same variables as total.rds with addition of:
 	           i. convert (diazepam conversion factor)
 	          ii. DME (diazepam milligram equivalents: dose*quantity*convert)
 	         iii. pills in this dataset correspond to 2 mg pill equivalents (DME/2)
 	     e. cbti_total.rds (444 clinicians x 30 study months (excludes baseline) = 13,320)
 	        1. prov_id (clinician ID)
 	        2. clinic_id (clinic ID)
 	        3. Tx1 (accountable justification intervention)
 	        4. Tx2 (default intervention)
 	        5. post (study period)
 	        6. mnth (month centered: 1-30)
 	        6. cpresc (clinic prescribing)
 	     f. cbti_short.rds (n = 13,320)
 	        1. same variables as cbti_total on short-term Z-drug users
 	     g. cbti_long (n = 13,320)
 	        1. same variables as cbti_total on long-term Z-drug users
 	     h. pills_long.rds (subset of visits where Z-drug was prescribed for Poisson model for long-term users (n = 50,419))
 	     i. pills_short.rds (subset of visits where Z-drug was prescribed for Poisson model for short-term users (n = 121,199))
 	     j. pills_total.rds (subset of visits where Z-drug was prescribed for Poisson model for all users (n = 171,618))
 	     k. pills_benzo.rds (subset of visits where benzodiazepine was prescribed (n = 248,218)
 	     
 	     *models.R*
 	     
 	     i. Executes models reported on clinicaltrials.gov
 	    ii. For hurdle models, separate glmer models used to estimate binomial (outcome = rx) and poisson (pills) outcomes for Z-drugs and benzodiazepines 
 	   iii. Glmer Poisson used for CBT-I orders 
 	    iv. Glmer binomial used for discordant Z-drug prescribing     

	