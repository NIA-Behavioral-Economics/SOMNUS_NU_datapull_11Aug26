library(sqldf)
library(dplyr)

#opioids
opioid <- sqldf("select name, medication_id, generic_name, strength, form, route, 1 as in_opioid from drugs 
                    where pharm_class_nm like '%OPIOID%' or
                      (pharm_class_nm like '%BULK%'
                        and pharm_subclass_nm like '%Opioid%') ")

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
saveRDS(zdrug_rx, file = "/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Data/zdrug_rx.rds")

#benzos 
benzo_rx <- sqldf("select name, medication_id, generic_name, strength, form, route, 1 as in_benzo from drugs 
                    where pharm_subclass_nm like '%Benzo%' and route in ('Oral', 'Intranasal', 'Sublingual', 'Nasal')
                    AND generic_name is not null 
                    AND (`Jeff Notes` is NULL OR `Jeff Notes` != 'Exclude')")
#add rownumber 
benzo_rx$rn <- seq_len(nrow(benzo_rx))
benzo_rx <- benzo_rx %>% select(NAME, medication_id, STRENGTH, FORM, ROUTE, rn)
saveRDS(benzo_rx, file = "/schhome/users/epstewar/SOMNUS/NU_data_pull_synth_data/Data/benzo_rx.rds")

#other 
other <- sqldf("select * from
                (select t.name as NAME, t.medication_id, t.generic_name AS GENERIC_NAME, t.strength AS STRENGTH, t.form AS FORM, t.route AS ROUTE, 
                l.in_opioid, n.in_zdrug, b.in_benzo
                from drugs t
                left join 
                opioid l
                on t.name = l.name
                left join 
                zdrug n
                on t.name = n.name 
                left join 
                benzo b
                on t.name = b.name) where in_opioid is null and in_zdrug is null and in_benzo is null")


