#load libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(truncnorm)
library(sqldf)
library(data.table)
set.seed(123654)

#################################### simulate aggregated visit data given by Ji Young #############################################
n_prov <- 444
n_clinic <- 64
nvsts_clinician <- expand.grid(prov_id = 1:n_prov, month = 1:48)

#total visits per clinician-month (mean = 250 (SD = 40))
nvsts_clinician$total_vsts <- round(rnorm(n = nrow(nvsts_clinician), mean = 250, sd = 40))

#remove negative values 
nvsts_clinician$total_vsts <- ifelse(nvsts_clinician$total_vsts < 0, 0, nvsts_clinician$total_vsts)

#visits for long-term users (assume 15% of visits)
nvsts_clinician$long_vsts <- round(nvsts_clinician$total_vsts * 0.15)

#short term visits 
nvsts_clinician$short_vsts <- nvsts_clinician$total_vsts - nvsts_clinician$long_vsts 
###################################################################################################################################

###################################### create XB grid #############################################################################

#assign each clinician a clinic
prov_df <- data.frame(
  prov_id = 1:n_prov,
  clinic_id = rep(1:n_clinic, length.out = n_prov)
)

#assign clinic variables: tx1 (AJ), tx2 (default), cpresc (clinic prescribing high vs. low)
clinic_df <- data.frame(
  clinic_id = 1:n_clinic,
  Tx1 = rbinom(n_clinic, 1, 0.5),
  Tx2 = rbinom(n_clinic, 1, 0.5),
  cpresc = rbinom(n_clinic, 1, 0.5)
)

#join provider/clinic data with clinic assignments 
clinician_demo <- left_join(prov_df, clinic_df, by = "clinic_id")

#add clinician demos. to dis-aggregated non-zdrug visits 
sample <- left_join(nvsts_clinician, clinician_demo, by = "prov_id")

#time covariates 
#mnth centered 
sample$mnth <- sample$month - 18

#post variable 
sample$post <- ifelse(sample$mnth < 1, 0, 
                      ifelse(sample$mnth > 0 & sample$mnth < 19, 1, 2))

#mnth truncated intervention period 
sample$kmnthTx <- ifelse(sample$mnth < 1|sample$mnth >18, 0, sample$mnth)

#mnth truncated post-intervention period 
sample$kmnthFu <- ifelse(sample$mnth < 19, 0, sample$mnth)
#############################################################################################################################

######################################################## add effects ########################################################
#Model part 1: binomial treatment effects for primary outcome and DME 
#time main effects
bmnth <- log(1)
bkmnthTx <- log(0.99)
bkmnthFu <- log(0.99)

#tx and covariate main effects
bcpresc <- log(1.2)

### will not have subgroups (n patient id) for non Z-drug visits if Ji Young provides aggregated values-need work around
btx1 <- log(0.99)
btx2 <- log(0.99)

#time and treatment two-way interactions
bmnth_tx1 <- log(1)   
bmnth_tx2 <- log(1) 
bkmnthTx_tx1 <- log(0.98)
bkmnthTx_tx2 <- log(0.98)
bkmnthFu_tx1 <- log(0.98)
bkmnthFu_tx2 <- log(0.98)

#two-way treatment interaction 
btx1_tx2 <- log(0.999)

#### remove 3-way interactions for now
#three-way treatment interaction 
bmnth_tx1_tx2     <- log(1)
bkmnthTx_tx1_tx2  <- log(0.99)
bkmnthFu_tx1_tx2  <- log(0.99)

#intercepts
#assume baseline z-drug prescribing rate of 2% for short-term users
beta0_logi_short  <- log(0.02/0.98)
#assume baseline z-drug prescribing rate of 7% for long-term users (0.15 * 0.46 = 0.07)
beta0_logi_long  <- log(0.05/0.95)
#assume baseline benzo Rx rate of 4%
beta0_logi_benzo <- log(0.04/0.96)

#CBT-I coefficients 
bcbti_mnth        <- log(1.0001)
bcbti_Tx1         <- log(1.0)
bcbti_Tx2         <- log(1.0)
bcbti_Tx1_mnth    <- log(1.5)
bcbti_Tx2_mnth    <- log(1.5)
bcbti_cpresc      <- log(1.30)
beta0_cbti_long   <- log(0.5/1)
beta0_cbti_short  <- log(0.3/1)

#ICCs
icc_clinic    <- 0.10
icc_provider  <- 0.20

#total variance for random intercepts
###lambda (mean pills) 
lambda <- 30
tvp <- log(1+(1/lambda))
#logistic 
tvl <- pi^2/3

#random intercepts (binomial)
#clinic
var_clinic_logi <- (icc_clinic/(1-icc_clinic))*tvl
clinic_re_logi <- data.frame(
  clinic_id = unique(sample$clinic_id),
  clinic_re_logi = rnorm(length(unique(sample$clinic_id)),0, sqrt(var_clinic_logi)))

#clinician 
var_provider_logi <- (icc_provider/(1-icc_provider))*tvl
provider_re_logi <- data.frame(
  prov_id = unique(sample$prov_id),
  provider_re_logi = rnorm(length(unique(sample$prov_id)), 0, sqrt(var_provider_logi)))

#ETA primary outcome binomial distribution 
sample$xb <-
  #time main effects
  bmnth     * sample$mnth       +
  bkmnthTx  * sample$kmnthTx    +
  bkmnthFu  * sample$kmnthFu    +
  
  #Tx and covariate main effects
  bcpresc  * sample$cpresc      +
  btx1     * sample$Tx1         +
  btx2     * sample$Tx2         +
  
  #tx1 and tx2 two-way interaction 
  btx1_tx2 * sample$Tx1 * sample$Tx2 +
  
  #tx by time two-way interactions 
  bmnth_tx1     * sample$mnth    * sample$Tx1 +
  bmnth_tx2     * sample$mnth    * sample$Tx2 +
  bkmnthTx_tx1  * sample$kmnthTx * sample$Tx1 +
  bkmnthTx_tx2  * sample$kmnthTx * sample$Tx2 +
  bkmnthFu_tx1  * sample$kmnthFu * sample$Tx1 +
  bkmnthFu_tx2  * sample$kmnthFu * sample$Tx2 +
  
  #time by tx1 and tx2 three-way interactions 
  bmnth_tx1_tx2    * sample$mnth    * sample$Tx1 * sample$Tx2 +
  bkmnthTx_tx1_tx2 * sample$kmnthTx * sample$Tx1 * sample$Tx2 +
  bkmnthFu_tx1_tx2 * sample$kmnthFu * sample$Tx1 * sample$Tx2 

#ETA CBT-I 
sample$xb_cbti <- 
  #time main effect
  bcbti_mnth      * sample$mnth              +
  #Tx main effects 
  bcbti_Tx1       * sample$Tx1               +
  bcbti_Tx2       * sample$Tx2               +
  #covariates 
  bcbti_cpresc    * sample$cpresc            +
  #Tx by time two-way interactions           
  bcbti_Tx1_mnth  * sample$Tx1 * sample$mnth +
  bcbti_Tx2_mnth  * sample$Tx2 * sample$mnth 

#merge with sample
sample <- sample %>%
  left_join(clinic_re_logi, by = "clinic_id") %>%
  left_join(provider_re_logi, by = "prov_id")

#add fixed and random intercepts to ETA 
#short term users, logistic
sample$ETA_logi_short <- beta0_logi_short + sample$xb + sample$clinic_re_logi + sample$provider_re_logi
#long term users, logistic
sample$ETA_logi_long <- beta0_logi_long + sample$xb + sample$clinic_re_logi + sample$provider_re_logi
#CBT-I short-term users 
sample$ETA_cbti_short <- beta0_cbti_short + sample$xb_cbti + sample$clinic_re_logi + sample$provider_re_logi
#CBT-I long-term users 
sample$ETA_cbti_long <- beta0_cbti_long + sample$xb_cbti + sample$clinic_re_logi + sample$provider_re_logi

#benzos (assumes same effect sizes except intercept)
sample$ETA_benzo <- beta0_logi_benzo + sample$xb + sample$clinic_re_logi + sample$provider_re_logi

#expected z-drug counts for short term users 
sample$zdrug_short_vsts <- rbinom(n = nrow(sample), size = sample$short_vsts, prob = plogis(sample$ETA_logi_short))
#long-term users
sample$zdrug_long_vsts <- rbinom(n = nrow(sample), size = sample$long_vsts, prob = plogis(sample$ETA_logi_long))
#benzos 
sample$benzo_vsts <- rbinom(n = nrow(sample), size = sample$total_vsts, prob = plogis(sample$ETA_benzo))
#CBT-I short-term user counts 
sample$cbti_short_vsts <- rbinom(n = nrow(sample), size = sample$zdrug_short_vsts, prob = plogis(sample$ETA_cbti_short))
#CBT-I long-term user counts 
sample$cbti_long_vsts <- rbinom(n = nrow(sample), size = sample$zdrug_long_vsts, prob = plogis(sample$ETA_cbti_long))
###############################################################################################################################

##################################### add study dates #########################################################################
#dataframe with start and end date for each month 
month_df <- data.frame(mnth = -17:30)
start <- as.Date("2023-05-14")
month_df$start_date <- start %m+% months(month_df$mnth + 17)

#shift months 1+ by 36 days
month_df$start_date[month_df$mnth >= 1] <-
  month_df$start_date[month_df$mnth >= 1] + 36

month_df$end_date <- c(
  month_df$start_date[-1] - 1,
  as.Date("2027-06-21")
)

#make month 0 end on 11/14/24
month_df$end_date[month_df$mnth == 0] <- as.Date("2024-11-14")

#add dates to sample
sample <- sqldf("select t.*, l.start_date, l.end_date from sample t
                left join 
                month_df l 
                on t.mnth = l.mnth
                order by prov_id, mnth")
#save data 
saveRDS(sample, file = "/directory/binomial.rds")
