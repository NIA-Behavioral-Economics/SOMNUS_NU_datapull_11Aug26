#############################################################################################################
### RUN SOMNUS MODELS WITH SYNTHETIC DATA FOR NU DATA PULL
#############################################################################################################

# changes: changed data structure, models running separately vs. using covariate in model, 
#and time frame (days vs. months) to match ct.gov 
# questions: discrepancy between protocol and ct.gov regarding cbt-i time frame 
# benzo pill equivalents (not DME)
#### 2-mg diazepam pill equivalents. DME divided by 2-unlike 5-mg z-drug equivalents where strength is divided by 5
# THE WAY JI YOUNG GIVES THE DATA SUBGROUP ANALYSES WILL NOT BE STRATIFIED BY USER TYPE (IN OTHER WORDS MODEL WLL BE OVERALL)
# should discordant encounters also include cbt-i prescrbing 
# restricted model: limited only to primary outcome? Also by user type and other subgroups? 
## do we want to know whether clinician removed an order? Current alert data tells us what happens if they keep an order, but we won't 
## know the cases where there was an order that was removed (unless it shows up as canceled?)
## unless in other words 

# LIBRARIES 
library(lme4)
library(sqldf)
library(glmmTMB)
library(dplyr)

# LOAD DATA 
zdrug_short <- readRDS("/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Data/zdrug_short.rds")
zdrug_long <- readRDS("/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Data/zdrug_long.rds")
total <- readRDS("/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Data/total.rds")
benzo <- readRDS("/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Data/benzo.rds")
  
# OUTPUT RESULTS 
sink("/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Results/SOMNUS_models.txt")

# MODEL 1
## PRIMARY OUTCOME: 5 MG PILL COUNTS
### RESTRICTED 
##### INTERVENTION ONLY 
###### HURDLE MODEL 
####### SEPARATING BINOMIAL AND POISSON OUTCOMES BECAUSE GLMER HAS FEWER CONVERGENCE ISSUES THAN GLMMTMB
primary_outcome <- function(dat, datpills) {
  
  # PART 1
  m1_p1 <- glmer(
    rx ~ 
      # time main effects
      mnth + kmnthTx + 
      # tx + covariate main effects 
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx +
      (1|clinic_id/prov_id),
    family = binomial,
    data = dat %>%
    filter(post != 2),
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
  
  # PART 2
  m1_p2 <- glmer(
    pills ~ 
      # time main effects
      mnth + kmnthTx + 
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx +
      (1|clinic_id/prov_id),
    data = datpills %>%
    filter(post != 2),
    family = poisson, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
  return(list(summary(m1_p1), summary(m1_p2)))
}
# CONVERGED
#m1_long <- primary_outcome(zdrug_long, pills_long)
#m1_long
# CONVERGED
#m1_short <- primary_outcome(zdrug_short, pills_short)
#m1_short
# CONVERGED
#m1_total <- primary_outcome(total, pills_total)
#m1_total

# MODEL 2
## SECONDARY OUTCOME: CBT-I ORDERS
### RESTRICTED 
#### INTERVENTION ONLY 
##### POISSON 
secondary_cbti <- function(dat) {
  m2 <- glmer(
    cbti ~ 
      # time main effect
      mnth + 
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth +
      (1|clinic_id/prov_id),
    data = dat %>%
    filter(post != 2),
    family = poisson, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
 return(summary(m2))
}
# CONVERGED
# m2_long <- secondary_cbti(cbti_long)
# CONVERGED
#m2_long
# CONVERGED
#m2_short <- secondary_cbti(cbti_short)
#m2_short
# CONVERGED
#m2_total <- secondary_cbti(cbti_total)
#m2_total

# MODEL 3
## SECONDARY OUTCOME: 2-MG DIAZEPAM PILL EQUIVALENTS
### RESTRICTED 
#### INTERVENTION ONLY 
###### HURDLE MODEL
m3_p1 <- glmer(
    rx ~ 
      # time main effects
      mnth + kmnthTx + 
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx +
      (1|clinic_id/prov_id),
      data = benzo %>%
      filter(post != 2),
    family = binomial, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
)
# CONVERGED
#summary(m3_p1)

m3_p2 <- glmer(
    pills ~ 
      # time main effects
      mnth + kmnthTx + 
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx +
      (1|clinic_id/prov_id),
      data = benzo %>%
      filter(post != 2 & rx == 1),
    family = poisson, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
)
summary(m3_p2)

# MODEL 4
## SECONDARY OUTCOME: Z-DRUG DISCORDANT PRESCRIBING; GREATER THAN 35 PILLS FOR 35 DAYS' SUPPLY OR LESS
### RESTRICTED 
#### INTERVENTION ONLY WHERE Z-DRUG PRESCRIBED 
###### POISSON
m4 <- glmer(
   discordant ~ 
      # time main effects
      mnth + kmnthTx + 
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx +
      (1|clinic_id/prov_id),
    data = total %>%
    filter(post != 2 & rx == 1),
    family = binomial, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
summary(m4)
  
# MODEL 5
## PRIMARY OUTCOME: 5 MG PILL COUNTS
### RESTRICTED 
##### INTERVENTION & POST-INTERVENTION 
###### HURDLE MODEL 
primary_outcome_post <- function(dat, dat_pills) {
    
# PART 1
m5_p1 <- glmer(
      rx ~ 
        # time main effects
        mnth + kmnthTx + kmnthFu +
        # tx + covariate main effects 
        cpresc + Tx1 + Tx2 +
        # tx + time two-way interactions
        Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx + Tx1:kmnthFu + Tx2:kmnthFu +
        (1|clinic_id/prov_id),
      data = dat,
      family = binomial,
      glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
    )
    
# PART 2
m5_p2 <- glmer(
      pills ~ 
        # time main effects
        mnth + kmnthTx + kmnthFu +
        # tx + covariate main effects
        cpresc + Tx1 + Tx2 +
        # tx + time two-way interactions
        Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx + Tx1:kmnthFu + Tx2:kmnthFu + 
        (1|clinic_id/prov_id),
      data = dat_pills,
      family = poisson, 
      glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
    )
    return(list(summary(m5_p1), summary(m5_p2)))
  }
m5_long <- primary_outcome_post(zdrug_long, pills_long)
m5_long
m5_short <- primary_outcome_post(zdrug_short, pills_short)
m5_short
m5_total <- primary_outcome_post(total, pills_total)
m5_total

# MODEL 6
## SECONDARY OUTCOME: CBT-I ORDERS
### RESTRICTED 
#### INTERVENTION & POST-INTERVENTION 
##### POISSON 
secondary_cbti_post <- function(dat) {
  m6 <- glmer(
    cbti ~ 
      # time main effects
      mnth + 
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + 
      (1|clinic_id/prov_id),
    data = dat,
    family = poisson, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
  return(summary(m6))
}
m6_long <- secondary_cbti_post(cbti_long)
m6_long
m6_short <- secondary_cbti_post(cbti_short)
m6_short
m6_total <- secondary_cbti_post(cbti_total)
m6_total 

# MODEL 7
## SECONDARY OUTCOME: 2-MG DIAZEPAM PILL EQUIVALENTS
### RESTRICTED 
#### INTERVENTION & POST-INTERVENTION ONLY 
###### HURDLE MODEL
m7_p1 <- glmer(
  rx ~ 
    # time main effects
    mnth + kmnthTx + kmnthFu + 
    # tx + covariate main effects
    cpresc + Tx1 + Tx2 +
    # tx + time two-way interactions
    Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx + Tx1:kmnthFu + Tx2:kmnthFu + 
    (1|clinic_id/prov_id),
  data = benzo,
  family = binomial, 
  glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
)
summary(m7_p1)

m7_p2 <- glmer(
  pills ~ 
    # time main effects
    mnth + kmnthTx + kmnthFu + 
    # tx + covariate main effects
    cpresc + Tx1 + Tx2 +
    # tx + time two-way interactions
    Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx + Tx1:kmnthFu + Tx2:kmnthFu + 
    (1|clinic_id/prov_id),
  data = benzo %>%
  filter(rx == 1),
  family = poisson, 
  glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
)
summary(m7_p2)

# MODEL 8
## SECONDARY OUTCOME: Z-DRUG DISCORDANT PRESCRIBING; GREATER THAN 35 PILLS FOR 35 DAYS' SUPPLY OR LESS
### RESTRICTED 
#### INTERVENTION & POST-INTERVENTION WHERE Z-DRUG PRESCRIBED 
###### POISSON
m8 <- glmer(
  discordant ~ 
    # time main effects
    mnth + kmnthTx + kmnthFu + 
    # tx + covariate main effects
    cpresc + Tx1 + Tx2 +
    # tx + time two-way interactions
    Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx + Tx1:kmnthFu + Tx2:kmnthFu + 
    (1|clinic_id/prov_id),
  data = total %>%
  filter(rx == 1),
  family = binomial, 
  glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
)
summary(m8)

# MODEL 9
## PRIMARY OUTCOME: 5 MG PILL COUNTS
### UNRESTRICTED 
##### INTERVENTION ONLY 
###### HURDLE MODEL 

# PART 1
  m9_p1 <- glmer(
    rx ~ 
      # time main effects
      mnth + kmnthTx + 
      # tx + covariate main effects 
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx +
      # tx1 + tx2 two-way interaction
      Tx1:Tx2 +
      # tx1 + tx2 + time three-way interactions
      Tx1:Tx2:mnth + Tx1:Tx2:kmnthTx + 
      (1|clinic_id/prov_id),
    family = binomial,
    data = total %>%
    filter(post != 2),
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
  summary(m9_p1)
  
  # PART 2
  m9_p2 <- glmer(
    pills ~ 
      # time main effects
      mnth + kmnthTx + 
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx +
      # tx1 + tx2 two-way interaction
      Tx1:Tx2 +
      # tx1 + tx2 + time three-way interactions
      Tx1:Tx2:mnth + Tx1:Tx2:kmnthTx + 
      (1|clinic_id/prov_id),
    data = pills_total %>%
    filter(post != 2),
    family = poisson, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
  summary(m9_p2)

# MODEL 10
## PRIMARY OUTCOME: 5 MG PILL COUNTS
### UNRESTRICTED 
##### INTERVENTION & POST-INTERVENTION 
###### HURDLE MODEL 
  
# PART 1
  m10_p1 <- glmer(
    rx ~ 
      # time main effects
      mnth + kmnthTx + kmnthFu +
      # tx + covariate main effects 
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx + Tx1:kmnthFu + Tx2:kmnthFu +
      # tx1 + tx2 two-way interaction
      Tx1:Tx2 +
      # tx1 + tx2 + time three-way interactions
      Tx1:Tx2:mnth + Tx1:Tx2:kmnthTx + Tx1:Tx2:kmnthFu + 
      (1|clinic_id/prov_id),
    data = total,
    family = binomial,
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
summary(m10_p1)
  
# PART 2
  m10_p2 <- glmer(
    pills ~ 
      # time main effects
      mnth + kmnthTx + kmnthFu +
      # tx + covariate main effects
      cpresc + Tx1 + Tx2 +
      # tx + time two-way interactions
      Tx1:mnth + Tx2:mnth + Tx1:kmnthTx + Tx2:kmnthTx + Tx1:kmnthFu + Tx2:kmnthFu + 
      # tx1 + tx2 two-way interaction
      Tx1:Tx2 +
      # tx1 + tx2 + time three-way interactions
      Tx1:Tx2:mnth + Tx1:Tx2:kmnthTx + Tx1:Tx2:kmnthFu + 
      (1|clinic_id/prov_id),
    data = pills_total,
    family = poisson, 
    glmerControl(calc.derivs = FALSE, optimizer = "bobyqa")
  )
summary(m10_p2)
sink()









