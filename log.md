# 4-17-2019

Model runs for parasite completed; moving on to hds with covariates.

# January - April

Mostly job applications getting in the way of completing the project. Worked on
it here and there. Lots of difficulty finding inlabru syntax to model parasite
intensity as a function of scat density, which itself is modeled under
hierarchical distance sampling.

# 12-21-2018

Downloaded PRISM data and NYS data.

# 11-19-2018

Error found in sample 542 -- it was not collected on a deer distance transect.
Correct name from 5B1-5 to 5B1

Distance sampling transect points updated to start from start of transect (which
is much more accurate a location given the GPS averaging done while the forest
metrics were taken) along the transect using the parallel distance measurement
and away from the transect using the perpendicular distance measurement.

NOTE: the specific *side* on which the scat was relative to the transect is
unknown -- these data were not collected. However, there are a few arguments in favor of random assignment of sides:

1. The model does not care on which side objects are found -- the distance is the key data
2. The distances are often *centimeters* away from the transect -- at the scale of the Adirondacks, this will not change the outcome appreciably
3. This derived location is **much** more accurate than the GPS locations taken, which are very imprecise due to data entry error and GPS error.

# 11-16-2018

Corrected errors in sample info

* Sample ID 460 corrected from 7B1-4 to 7B1-2
* Sample ID 461 corrected from 7B1-5 to 7B1-4
* Sample ID 462 corrected from 7B1-6 to 7B1-5
* Sample ID 480 corrected from 2B4 to 2B4-1
* Sample ID 519 corrected from 15A4 to 15A4-2
* Sample ID 476 corrected from 6C3-3 to 6C3-4
* Sample ID 413 and 414 do not belong to 10B3-3; corrected to 10B3 ambiguous location.
* Sample ID 478 corrected from 2B4-1 to 2B4

# 11-15-2018

Found some errors in distance sampling data entry

* 9A3-2 Start coordinate wrong
* 1B4-5 start northing 10,000 km off
* 2B1-6 imputed end coordinate

Fixed the following (numbers in PK):

* 32, changed Easting value from 666904 to 566904
* 542, changed Northing value from 429289 to 4929289
* 539, changed Northing value from 4921502 to 4929502
* 211, change Northing value from 4844184 to 4944184
* 607, 608 changed Northing value from 48... to 49...
* 11, changed easting value from 572104 to 562104
* 606, changed northing value from 48... to 49...
* 558, changed northing value from 4849192 to 4849792
* 260, changed northing value from 4965246 to 4865246
* 113, changed 4921... to 4928...
* 488, 489, changed easting from 579... to 578... . Changing Northing rather futile, will await displacement from origin.

STOPPING HERE with 2018 points. There is no need to change these as they will be all adjusted to the way more accurate displacement measure.

# 11-14-2018

Verified that metadata 2018 info is same as entered by me. 

# 11-13-2018

Since last update, few changes made. Minimum distances from scat collection locations to deer transect lines are somewhat large, but a lot were collected off of DS transects. 

Today I completed data entry for the distance sampling data, and I'm working on converting the GPS locations to more accurate locations from transect angles, distance along the transect, and distance perpendicular to the transect.

## Bad locations

Using GIS I identify a few bad locations from the 2018 data. They are the following.

* 464, 7B1-6
* 465, 7B5
* 466, 7B5
* 413, 10B3-3
* 414, 10B3-3
* 483, 2B4-4
* 478, 2B4-1
* 531, 15A2
* 532, 12B1-2
* 533, 12B1
* 467, 7B5-3

From other years, it's readily apparent that the following are bad, as they're outside the park

* 32, 10A3, 2016
* 136

The question is where these data should be altered, and how to propagate corrections. For 2017 data, the best place to alter is the Access database, and propagate here. For 2016 and 2018 data, the best place to alter is the copy in the top-level directory here. 

The locations that belong to distance sampling transects can be easily corrected using the location of the start of the transect, and the displacement to the recorded location. Locations *not* belonging to distance sampling transects that are severely out of place will have to be imputed to the center of the moose transect.


# 10-30-2018

## Activities

### `dataCleaning.R`

### `fecal_results_relevant_copy.xlsx`

There exist duplicates in the specimen ID column. There is no POINT to an ID column if there are duplicates. Therefore a new primary key column has been created within this file. Data are NOT WELL MANAGED IN EXCEL, sheesh.

* Imputed 0's into blank cells in fecal quantitative section using select & find > go to special > blanks > set 0 > ctrl-enter

### `metadata.csv`

* Created a new metadata file from `fecal_results_relevant_copy.xlsx`

# 10-29-2018

## Activities

### `fecal_results_relevant_copy.xlsx`

* Imputed NAs into empty cells, 0's where appropriate.
    * Only imputed 0's for DSL and F. magna.
* Removed 'NES' entries replaced with NA.
* Created 'condition' column where sample condition SHOULD have been recorded, instead of some notes column.
* Imputed metadata for indicated records where only Easting was available. Most of these were unique values for Easting, so there is no concern of wrong imputation. The only record without a unique value for Easting was Specimen ID 395, which matches that of 404. Since occur on the same day in the same general area, they can be treated as interchangeable.