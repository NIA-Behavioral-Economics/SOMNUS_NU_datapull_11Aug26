library(tidyr)
library(stringr)
library(extraDistr)
library(dplyr)
library(data.table)
library(sqldf)
library(lubridate)

#LOAD DATA 
sample <- readRDS("/directory/binomial.rds")
zdrug_rx <- readRDS("/directory/zdrug_rx.rds")
benzo_rx <- readRDS("/directory/benzo_rx.rds")

#non-rx visits for each dataframe: short-term patients, long-term patients, and patients who received benzos 
sample$non_zdrug_short_vsts <- sample$short_vsts - sample$zdrug_short_vsts
sample$non_zdrug_long_vsts <- sample$long_vsts - sample$zdrug_long_vsts
sample$non_benzo_vsts <- sample$total_vsts - sample$benzo_vsts

#unfurl prescription counts based on probabilities 
expand_visits <- function(rx_vsts, add_cbti, cbti_var,
                          #coefficients 
                          es1, es2, es3, es4, 
                          es5, es6, es7, es8, 
                          es9, es10, es11, es12,
                          es13, es14, es15, es16, 
                          int, lambda, nrx, rxdat, non_rx_vsts) {
  
rx <- sample %>%
    tidyr::uncount({{ rx_vsts }}, .id = "visit") %>%
    group_by(prov_id, mnth) %>%
    mutate(
      cbti = if (add_cbti & first(mnth) > 0) {
        sample(c(rep(1, first({{ cbti_var }})),
                 rep(0, n() - first({{ cbti_var }}))))
      } else {
        0L
      }
    ) %>%
    ungroup() %>%
    select(
      prov_id,
      clinic_id,
      Tx1,
      Tx2,
      post,
      mnth,
      start_date,
      end_date,
      kmnthTx,
      kmnthFu,
      cpresc,
      cbti
    )

  #effect sizes
  bmnth                <- log(es1)
  bkmnthTx             <- log(es2)
  bkmnthFu             <- log(es3)
  bcpresc              <- log(es4)
  btx1                 <- log(es5)
  btx2                 <- log(es6)
  btx1_tx2             <- log(es7)
  bmnth_tx1            <- log(es8)
  bmnth_tx2            <- log(es9)
  bkmnthTx_tx1         <- log(es10) 
  bkmnthTx_tx2         <- log(es11) 
  bkmnthFu_tx1         <- log(es12)
  bkmnthFu_tx2         <- log(es13)
  bmnth_tx1_tx2        <- log(es14)
  bkmnthTx_tx1_tx2     <- log(es15)
  bkmnthFu_tx1_tx2     <- log(es16)
  beta0                <- log(int)
  
  #ETA poisson distribution (pill quantity)
  rx$xb <-
    #time main effects
    bmnth     * rx$mnth                                   +
    bkmnthTx  * rx$kmnthTx                                +
    bkmnthFu  * rx$kmnthFu                                +
    
    #Tx and covariate main effects
    bcpresc  * rx$cpresc                                  +
    btx1     * rx$Tx1                                     +
    btx2     * rx$Tx2                                     +
    
    #tx1 and tx2 two-way interaction 
    btx1_tx2 * rx$Tx1 * rx$Tx2                            +
    
    #tx by time two-way interactions 
    bmnth_tx1     * rx$mnth    * rx$Tx1                   +
    bmnth_tx2     * rx$mnth    * rx$Tx2                   +
    bkmnthTx_tx1  * rx$kmnthTx * rx$Tx1                   +
    bkmnthTx_tx2  * rx$kmnthTx * rx$Tx2                   +
    bkmnthFu_tx1  * rx$kmnthFu * rx$Tx1                   +
    bkmnthFu_tx2  * rx$kmnthFu * rx$Tx2                   +
    
    #time by tx1 and tx2 three-way interactions 
    bmnth_tx1_tx2    * rx$mnth    * rx$Tx1 * rx$Tx2       +
    bkmnthTx_tx1_tx2 * rx$kmnthTx * rx$Tx1 * rx$Tx2       +
    bkmnthFu_tx1_tx2 * rx$kmnthFu * rx$Tx1 * rx$Tx2 
  
  #random intercepts assuming 0.10 ICC for clinic and 0.20 for clinician 
  #total variance for Poisson distribution
  tvp <- log(1+(1/lambda))
  #clinic
  clinic_var <- (0.10/(1-0.10))*tvp
  clinic_re <- data.frame(clinic_id = unique(rx$clinic_id),
  clinic_re = rnorm(length(unique(rx$clinic_id)), mean = 0, sd = sqrt(clinic_var)))
  
  #clinician 
  prov_var <- (0.20/(1-0.20))*tvp
  prov_re <- data.frame(prov_id = unique(rx$prov_id),
  prov_re = rnorm(length(unique(rx$prov_id)), mean = 0, sd = sqrt(prov_var)))
  
  #merge random intercepts
  rx <- rx %>%
    left_join(clinic_re, by = "clinic_id") %>%
    left_join(prov_re, by = "prov_id")
  
  #ETA
  rx$ETA <- exp(beta0 + rx$xb + rx$clinic_re + rx$prov_re)
  #convert ETA to counts 
  #this simulates a truncated Poisson distribution, where every value must be greater than 0
  rx$quantity <- extraDistr::rtpois(nrow(rx), lambda = rx$ETA, a = 1)
  
  #add drug strength 
  #generate random number from 1 to n based on number of possible drugs 
  rx$rn <- sample(1:nrx, size = nrow(rx), replace = TRUE)
  
  #merge with drug datasets
  setDT(rx)
  setDT(rxdat)
  rx <- rxdat[rx, on = .(rn)]
  
  #drug dose 
  #remove characters from strength for z-drugs and benzos
  rx$dose <- as.numeric(str_extract(rx$STRENGTH, "\\d+\\.?\\d*"))
  
  #add hypothetical days' supply and rx dates
  #simulate start date based on month start date (start_date)
  rx <- rx %>%
    mutate(rx_start_date = start_date + floor(runif(n(), min = 0, 
          max = as.numeric(end_date - start_date) + 1)),
  
  #simulate days supply
  #assume baseline mean of 45 days (SD = 15)
  #assume days' supply monthly decrease of 0.30 days (<1%) for treatment*time interaction effects
  mean_days = 45 - if_else(Tx1 == 1 & Tx2 == 1, 0.30 * pmax(kmnthTx, 0), 0) -
          if_else(Tx1 == 1 & Tx2 == 1, 0.30 * pmax(kmnthFu, 0), 0),
        
  #assume minimum of 5 days' supply
  mean_days = pmax(mean_days, 5),
  
  #simulate days' supply
  days_supply = round(rnorm(n(), mean = mean_days, sd = 15)),
  
  #prevent negative values 
  days_supply = pmax(days_supply, 1),
  
  #create prescription end date
  rx_end_date = rx_start_date + days_supply
  
  #unfurl non-rx visits (zero portion of binomial model)
  nonrx <- sample %>%
    tidyr::uncount({{ non_rx_vsts }}) %>%
    dplyr::select(
      prov_id,
      clinic_id,
      Tx1,
      Tx2,
      post,
      mnth,
      kmnthTx,
      kmnthFu,
      cpresc
    )
  
  #add columns to non-rx data
  nonrx$quantity <- NA_integer_
  nonrx$dose <- NA_integer_
  nonrx$rx <- 0
  nonrx$cbti <- NA_integer_
  nonrx$rx_start_date <- as.Date(NA)
  nonrx$rx_end_date <- as.Date(NA)
  nonrx$days_supply <- NA_integer_
  nonrx$NAME <- NA_character_
  
  #select columns from rx data
  rx <- rx %>%
  select(prov_id, 
         clinic_id, 
         Tx1, 
         Tx2, 
         post, 
         mnth, 
         kmnthTx, 
         kmnthFu, 
         cpresc, 
         rx_start_date, 
         rx_end_date, 
         NAME,
         days_supply, 
         quantity, 
         dose, 
         cbti)
  
  rx$rx <- 1
  
  #append rx and non-rx data
  data <- rbind(nonrx, rx)
  return(data)
}
benzo       <- expand_visits(benzo_vsts, FALSE, NULL, 0.995, 0.995, 0.995, 0.98, 
                             0.95, 0.95, 0.97, 1, 1, 0.98, 0.98, 0.98, 0.98, 1, 
                             0.99, 0.99, 30, 30, 123, benzo_rx, non_benzo_vsts)

zdrug_long  <- expand_visits(zdrug_long_vsts, TRUE, cbti_long_vsts, 0.995, 0.995, 0.995, 0.98, 
                             0.95, 0.95, 0.97, 1, 1, 0.98, 0.98, 0.98, 0.98, 1, 0.97, 
                             0.97, 90, 90, 45, zdrug_rx, non_zdrug_long_vsts)

zdrug_short <- expand_visits(zdrug_short_vsts, TRUE, cbti_short_vsts, 0.995, 0.995, 0.995, 0.98, 
                             0.95, 0.95, 0.97, 1, 1, 0.97, 0.97, 0.97, 0.97, 1, 
                             0.96, 0.96, 30, 30, 45, zdrug_rx, non_zdrug_short_vsts)

#combine short- and long-term users for discordant prescribing analysis 
total <- rbind(zdrug_long, zdrug_short)

#CTB-I datasets overall and by user type
cbti <- function (dat) {
  
  #aggregated CBT-I data 
  outdat <- dat %>%
    filter(post != 0) %>%
    group_by(prov_id, 
             clinic_id,
             Tx1,
             Tx2,
             post,
             mnth,
             cpresc) %>%
    summarise(
      cbti = sum(cbti, na.rm = TRUE),
      .groups = 'drop'
    )
  
  #return datasets
  return(outdat)
}
cbti_short <- cbti(zdrug_short)
cbti_long <- cbti(zdrug_long)
cbti_total <- cbti(total)

# SUBSET DATA FOR POISSON OUTCOME 
pills <- function(dat) {
 out <- dat %>%
      filter(rx == 1) %>%
       mutate(pills = round((dose/5)*quantity))
 return(out)
}
pills_short <- pills(zdrug_short)
pills_long <- pills(zdrug_long)
pills_total <- pills(total)

#guideline discordant prescribing 
total$discordant <- ifelse(total$days_supply <= 35 & total$quantity > 35 & total$rx == 1, 1, 0)
total$discordant <- ifelse(total$rx == 0, NA, total$discordant)

#save data 
saveRDS(total, file = "/directory/total.rds")
saveRDS(zdrug_long, file = "/directory/zdrug_long.rds")
saveRDS(zdrug_short, file = "/directory/zdrug_short.rds")

#DME 
#CONVERSION FACTORS
benzo$convert <- ifelse(grepl("ALP|XAN|CLON|KLON", benzo$NAME), 10,
                    ifelse(grepl("ATIV|LORAZEPAM|LOREE", benzo$NAME), 5,
                      ifelse(grepl("CLORAZE|TRAN", benzo$NAME), 0.67,
                        ifelse(grepl("DIAZEPAM|VALI|DIAST|VALTOCO", benzo$NAME), 1,
                          ifelse(grepl("ESTA", benzo$NAME), 7.5,
                            ifelse(grepl("FLUR|PRO|DALM", benzo$NAME), 0.33,
                              ifelse(grepl("HALC|TRIA", benzo$NAME), 20,
                                ifelse(grepl("OXA|SERA", benzo$NAME), 0.33,
                                  ifelse(grepl("REST|TEM|CLOB|SYMPA", benzo$NAME), 0.5,
                                    ifelse(grepl("CHLORDIAZEPOXIDE", benzo$NAME), 0.4,
                                     NA))))))))))

#DIAZEPAM MILLIGRAM EQUIVALENTS 
benzo$DME <- benzo$dose*benzo$quantity*benzo$convert

#2-MG DIAZEPAM PILL EQUIVALENTS 
benzo$DME <- benzo$dose*benzo$quantity*benzo$convert
pills_benzo <- benzo %>%
       filter(rx == 1 & !is.na(DME)) %>%
       mutate(pills = round((DME/2)))

#save benzo data
saveRDS(benzo, file = "/directory/benzo.rds")






















