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