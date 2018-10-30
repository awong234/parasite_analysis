# Setup -------------------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rgdal)

# Clean up metadata  ------------------------------------------------------------

metadata = read.csv(file = 'metadata.csv', stringsAsFactors = F)

metadata %>% head
metadata %>% str

metadata$Date = as.Date(metadata$Date, format = '%m/%d/%Y')

# Remove missing location values

missingLocRows = with(data = metadata, expr = {(is.na(Easting) | is.na(Northing))})

metadata = metadata[!missingLocRows,]

# Verify locations - no large deviations

transects = readOGR(dsn = 'DeerTransects2018/DeerTransects2018.shp', layer = 'DeerTransects2018', stringsAsFactors = F)
transects@coords = transects@coords[,c(1,2)]
attr(x = transects@coords, which = 'dimnames') = list(NULL, c("Easting", "Northing"))

metadata %>% filter(Year == 2018)

plot(transects %>% data.frame)

# Deer fecal locations in 2016 and 2018 can potentially be improved by using the transect start location, and displacement along the transect and off of the transect.
# This would require gathering the transect angles (should be easy from access db), and calculating a vector by which to transform the original starting point.



# Load parasite analysis data and check ------------------------------------------------------------

ff = read.csv(file = 'fecal_flukefinder.csv', stringsAsFactors = F)
mb = read.csv(file = 'fecal_MB.csv', stringsAsFactors = F)
quant = read.csv(file = 'fecal_quant.csv', stringsAsFactors = F)

# remove NA's and add metadata

ff = ff[!{ff$Total_eggs %>% is.na},]

ff_data = left_join(ff, metadata, by = c("PK" = "PK"))

mb_data_dsl = mb %>% select(PK:larvae_per_g_dsl) %>% left_join(y = metadata, by = c("PK" = "PK"))
mb_data_dsl = mb_data_dsl[!is.na(mb_data_dsl$Total_larvae_dsl),]

mb_data_fmagna = mb %>% select(PK, total_eggs_fmagna, eggs_per_g_fmagna) %>% left_join(y = metadata, by = c("PK" = "PK"))
mb_data_fmagna = mb_data_fmagna[!is.na(mb_data_fmagna$total_eggs_fmagna),]

