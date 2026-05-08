we noticed that in the project we have multiple background database and we have to assign date.time to each database regardless of whether 
this database should be linked to the foreground database or not (e.g., if the foreground system says "hydro, 2050, SSP5-H", then even with the dynamic 
supply chain, it should only traverse through the premise database `2050, SSP5-H` or maybe 2040/2030 SSP5-H if they're included in the project).  
However, if we dynamically comment out those background project databases (w/o assigning a data.time) if it's not supposed to be used by a foreground DB, 
it always has "KeyError: '*database'" so we tested in `dynElec_v4.1_split_BGyear_hydro_reservoir`, always the same KeyError 
Then, we further tested in `dynElec_v4.2*`, and have only ONE foreground process in the foreground DB, now we can sucessfully link it to the ONLY ONE background 
DB with the date.time assigned. 
In addition, in  `dynElec_v4.2*` we tested a foreground system under 2050, as it's assigned as a 10 years uniform distribution, then for a 2050 process, 
it should start with 2046-01-01, however, when we run timex, it always has a warning "Reference date 2026-01-01 00:00:00 is lower than all provided dates"
so appearantly, with only one foreground activity (a 2050/SSP5 hydro), it can be mapped to ONLY one selected background database_dates (just 2050 SSP background DB), 
but still the output file has the date column from 2026 to 2036 even if all foreground / background only inlc. 2050
Thus, we based on these tests, we assume that: 
1. even though all background database are assigned with date.time,  only the relevant ones (where the technosphere can be traversed) will be used 
2. as date.time output is always 2026 - 2036, so we change the data.time output for each dynLCI output to the correct date.time, e.g., 2046-2056 for a 2050 process 
