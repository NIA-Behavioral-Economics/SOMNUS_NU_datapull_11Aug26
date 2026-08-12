library(sqldf)
library(dplyr)
library(readxl)

#import pharamceutical data 
data <- read_excel("/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Data/drugs.xlsx")

#z drugs 
zdrug <-  sqldf("select name, medication_id, generic_name, strength, form, route, `Jeff Notes`,
                  1 as in_zdrug from drugs 
                    where (pharm_subclass_nm like '%Orexin%' or 
                           pharm_subclass_nm like '%GABA%')
                            AND route != 'Intravenous'
                            AND (`Jeff Notes` IS NULL OR `Jeff Notes` != 'Exclude')")
#add rownumber 
zdrug$rn <- seq_len(nrow(zdrug))
zdrug_rx <- zdrug %>% select(NAME, medication_id, STRENGTH, FORM, ROUTE, rn)

#save data 
saveRDS(zdrug_rx, file = "/directory/zdrug_rx.rds")

#benzos 
benzo_rx <- sqldf("select name, medication_id, generic_name, strength, form, route, 1 as in_benzo from drugs 
                    where pharm_subclass_nm like '%Benzo%' and route in ('Oral', 'Intranasal', 'Sublingual', 'Nasal')
                    AND generic_name is not null 
                    AND (`Jeff Notes` is NULL OR `Jeff Notes` != 'Exclude')")
#add rownumber 
benzo_rx$rn <- seq_len(nrow(benzo_rx))
benzo_rx <- benzo_rx %>% select(NAME, medication_id, STRENGTH, FORM, ROUTE, rn)
saveRDS(benzo_rx, file = "/directory/benzo_rx.rds")


