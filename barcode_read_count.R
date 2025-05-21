#GET BARCODE READ COUNTS BY SAMPLE

#execute script from run directory like so: Rscript ../rv_temp/barcode_read_count.R

library(tidyverse)
bc <- read_csv('barcode_read_count_qc.csv')
m <- read_csv('userInput/metadata.csv.bak', skip = 10)
df <- m %>% select(project.id, specimen.id, specimen.type, panel, barcode) %>% inner_join(bc, by = 'barcode') %>% arrange(read_count)
write_csv(df, 'sample_barcode_details.csv')
low_reads <- filter(df, read_count <50 & specimen.type == "tumor")
print(low_reads, width=150)