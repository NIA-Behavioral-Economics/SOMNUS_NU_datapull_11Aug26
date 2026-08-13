**PURPOSE**  
The SOMNUS_NU_datapull_11Aug26 repository contains code to create synthetic data to run analyses prespecified on clinicaltrials.gov.
This is to ensure that all required data elements and variables will be accounted for in one-time data pull at Northwestern Medicine.

**DATA**  
Synthetic prescription and clinician-monthly visit data with prespecified fixed and random effect sizes for two-part hurdle model with binomial and Poisson distributions.

**ANALYSES**  
*PRIMARY*  
1. Two-stage knotted spline hurdle model testing truncated kmonth x treatment reductions in Z-drug prescribing frequency (binomial) and 5-mg pill counts (Poisson)
2. Other analyses with primary outcome:

   * During intervention (study months 1 to 18) and post-intervention (study months 19 to 30)

     * Overall and by user type

*SECONDARY*  
1. Mixed logistic model testing month x treatment change in clinician-monthly CBT-I orders from intervention start to intervention end

   * By user type (short- vs. long-term)
2. Two-stage knotted spline hurdle model testing truncated kmonth x treatment reductions in benzodiazepine prescribing frequency (binomial) and 2-mg pill equivalents (Poisson)
3. Mixed logistic model testing truncated kmonth x treatment reductions in discordant Z-drug prescribing (> 35 pills for < 35 days' supply)
4. Secondary outcomes also assessed in post-intervention period

*UNRESTRICTED*  
1. Hurdle model for primary outcome will be tested using an unrestricted model which includes two-way treatment interactions (accountable justification x default) and three-way knotted time by treatment interactions (accountable justification x default x kmonth)

   * Unrestricted model will be tested separately in intervention and post-intervention periods
   * If interaction effects are insignificant, then we will proceed with restricted model (i.e., no treatment interactions)

**LICENSE**  
Schaeffer Center for Health Policy and Economics, University Southern California

**FILES**  
* *drug_tables.R*

  * Creates separate tables for each pharmaceutical class (e.g., Z-drugs, benzodiazepines) to derive drug strength for synthetic prescription data

    * `zdrug_rx` (prescription name, strength, route, etc. for Z-drugs)
    * `benzo_rx` (same as above for benzodiazepines)

* *binomial_visits.R*

  * Creates a flat file with aggregated visit counts for every study month per-clinician (444 clinicians * 48 study months = 21,312 rows)
  * `binomial.rds`

    * Variable list:

      * `prov_id` (clinician ID: 1-444)
      * `month` (study month: baseline = 1-18, intervention = 19-36, post-intervention = 37-48)
      * `total_vsts` (total visits per clinician-month: mean = 250, sd = 40)
      * `long_vsts` (total visits per clinician-month for long-term Z-drug users (> 180 days of Z-drugs in past 365 days): 15% of total visits)
      * `short_vsts` (total visits per clinician-month for short-term Z-drug users (<= 180 days of Z-drugs in past 365 days): total_vsts - long_vsts)
      * `clinic_id` (clinic ID: 1-64)
      * `Tx1` (accountable justification: 50% of clinics randomized)
      * `Tx2` (default: 50% of clinics randomized)
      * `cpres` (clinic prescribing (high vs. low): 50% of clinic randomized to 'high' prescribing)
      * `mnth` (month centered: baseline = -17-0, intervention = 1-18, post-intervention = 19-30)
      * `post` (study period: baseline = 0, intervention = 1, post-intervention = 2)
      * `kmnthTx` (intervention study month for knotted spline model: baseline = 0, intervention = 1-18, post-intervention = NA)
      * `kmnthFu` (post-intervention study month for knotted spline model: baseline = 0, intervention = 0, post-intervention = 19-30)
      * `xb` (fixed effects linear predictor for binomial Z-drug prescription outcome on log scale)
      * `xb_cbti` (fixed effects linear predictor for binomial CBT-I outcome on log scale)
      * `clinic_re_logi` (random clinic intercept for binomial outcomes on log scale: total variance = pi^2/3)
      * `provider_re_logi` (random clinician intercept for binomial outcomes on log scale: total variance = pi^2/3)
      * `ETA_logi_short` (linear predictor for binomial Z-drug prescription outcome on log scale including fixed and random effects for long-term users)
      * `ETA_logi_long` (linear predictor for binomial Z-drug prescription outcome on log scale including fixed and random effects for short-term users)
      * `ETA_cbti_short` (linear predictor for binomial CBT-I order outcome on log scale including fixed and random effects for short-term users)
      * `ETA_cbti_long` (linear predictor for binomial CBT-I order outcome on log scale including fixed and random effects for long-term users)
      * `ETA_benzo` (linear predictor for binomial benzodiazepine prescription outcome on log scale including fixed and random effects)
      * `zdrug_short_vsts` (randomly generated number of clinician-monthly visits with a Z-drug prescription for short-term Z-drug users based on probabilities from `ETA_logi_short`)
      * `zdrug_long_vsts` (randomly generated number of clinician-monthly visits with a Z-drug prescription for long-term Z-drug users based on probabilities from `ETA_logi_long`)
      * `benzo_vsts` (randomly generated number of clinician-monthly visits with a benzodiazepine prescription based on probabilities from `ETA_benzo`)
      * `cbti_short_vsts` (randomly generated number of clinician-monthly visits with a CBT-I order for short-term users based on probabilities from `ETA_cbti_short`)
      * `cbti_long_vsts` (randomly generated number of clinician-monthly visits with a CBT-I order for long-term users based on probabilities from `ETA_cbti_long`)
      * `start_date` (start date for study month)
      * `end_date` (end date for study month)

* *poisson_visits.R*

  * Combines visit, clinician, and clinic data from `binomial.rds` with Poisson prescription outcomes to create analytic datasets
  * `total.rds` (analytic dataset for primary Z-drug outcome for all visits, n = 5,321,768)

    * Variable list:

      * `prov_id` (clinician ID)
      * `clinic_id` (clinic ID)
      * `Tx1` (accountable justification intervention)
      * `Tx2` (default intervention)
      * `post` (study period)
      * `mnth` (month centered)
      * `kmnthTx` (intervention study month)
      * `kmnthFu` (post-intervention study month)
      * `cpresc` (high vs. low clinic prescribing)
      * `quantity` (NA if rx = 0 (i.e., no prescription), otherwise pill quantity randomly generated using same methodology for binomial outcome above)
      * `dose` (prescription info. merged from zdrug_rx, NA if rx = 0, otherwise destringed drug strength)
      * `rx` (Z-drug prescription: 0 = no, 1 = yes)
      * `cbti` (CBT-I order where rx = 1: 0 = no, 1 = yes)
      * `rx_start_date` (start date for prescription randomly generated from study month dates)
      * `rx_end_date` (end date for prescription: start_date + days_supply)
      * `days_supply` (randomly generated, assumes a minimum of 5 and a 0.30 decrease over time post-intervention for treatment groups)
      * `name` (prescription name: NA if rx = 0, otherwise populated)
      * `pills` (5 mg Z-drug pill equivalents: (dose/5) x quantity, NA if rx = 0, else populated)
      * `discordant` (Z-drug prescription not within guidelines: > 35 pills for 5-week days' supply, NA if rx = 0, otherwise populated)
  * `zdrug_long.rds` (n = 798,285)

    * Same variables as total.rds on subset of patients that are long-term users 
  * `zdrug_short.rds` (4,523,483)

    * Same variables as total.rds on subset of patients that are short-term users 
  * `benzo.rds` (n = 5,321,768)

    * Same variables as total.rds with addition of:

      * `convert` (diazepam conversion factor)
      * `DME` (diazepam milligram equivalents: dose x quantity x convert)
      * `pills` in this dataset correspond to 2 mg pill equivalents (DME/2)
  * `cbti_total.rds` (444 clinicians x 30 study months (excludes baseline) = 13,320)

    * Variable list

      * `prov_id` (clinician ID)
      * `clinic_id` (clinic ID)
      * `Tx1` (accountable justification intervention)
      * `Tx2` (default intervention)
      * `post` (study period)
      * `mnth` (month centered: 1-30)
      * `cpresc` (clinic prescribing)
  * `cbti_short.rds` (n = 13,320)

    * Same variables as cbti_total for short-term Z-drug users
  * `cbti_long.rds` (n = 13,320)

    * Same variables as cbti_total for long-term Z-drug users
  * `pills_long.rds` (subset of visits where Z-drug was prescribed for Poisson model for long-term users (n = 50,419))
  * `pills_short.rds` (subset of visits where Z-drug was prescribed for Poisson model for short-term users (n = 121,199))
  * `pills_total.rds` (subset of visits where Z-drug was prescribed for Poisson model for all users (n = 171,618))
  * `pills_benzo.rds` (subset of visits where benzodiazepine was prescribed (n = 248,218))

* *models.R*

  * Executes models reported on clinicaltrials.gov (n = 28)
  * For hurdle models, separate glmer models used to estimate binomial (outcome = `rx`) and Poisson (outcome = `pills`) outcomes for Z-drugs and benzodiazepines
  * Glmer Poisson used for CBT-I orders
  * Glmer binomial used for discordant Z-drug prescribing
