# Paxson Lake Humpback Whitefish Abundance and Demographics, 2023

## Abstract

This study will use mark-recapture methodology (two-event Petersen experiment) to estimate the abundance and
length composition of humpback whitefish Coregonus pidschian within Paxson Lake. The marking event will occur
in late May and early June 2023 and will consist primarily of using entanglement nets to capture whitefish. The
recapture event will take place from November 2023 through January 2024 and will consist of sampling subsistence-
harvested fish. ADF&G staff will also capture, sacrifice, and sample an additional 225–300 fish during the recapture
event to estimate sex composition, age composition, length composition, length-at-age, age-at-maturity, and length-
at-maturity. Understanding humpback whitefish abundance and demographics will allow for more informed
management of the subsistence fishery.

## Operational Plan

The Operational Plan for this study may be found here:

https://www.adfg.alaska.gov/FedAidPDFs/ROP.SF.3F.2023.01.pdf

## Repository Structure

Materials specific to the report written in 2024 are compiled in the folder **FDS_2024**.  Within this folder are subfolders:

### R

This folder consists of R scripts to either validate or conduct analysis, as well as some auxiliary scripts to supplement
or explore aspects of existing scripts.  These include:

**1_PLHW_MR_data.R**: This script reads the Mark-Recapture data and validates the Excel-based analysis conducted 
by Corey Schwanke and April Behr.  No analysis outputs are directly generated here, but inferences match those 
produced in Excel.  Additional sections have been inserted to validate numbers in the report text.

**1a_MRcheck2.R**: Additional checks in a fresh script.  This script has no real value and should probably be deleted.

**PLHW_age.R**: This script conducts all analysis specific to Age, Length, and Weight of fish in the spawning sample.
This script is used for analysis, and ALL analysis outputs are generated here.  
Additional sections have been inserted to validate numbers in the report text.

**2a_lengthage_sensitivity.R**: An additional script exploring the possible existence of data points that were 
unusually influential to the fit of the Length ~ Age relationship.  Each data point was removed in turn and the 
top two models were re-fit, and inferences compared.  At present, no data points were censored.

**recapr_prep.R**: Initial work developing a generalized function for processing mark-recapture data.  This will likely
be moved to a different folder.


### R_output

This folder consists of tables and figures directly exported from R script PLHW_age.R.  Table and figure numbers are
given by the numbers of tables and figures from the report.

**Fig3_xx** gives plots of Length ~ Age data with an overlayed envelope corresponding the the model description in **xx**

**Tab8.csv** gives p-values from KS tests of lengths and weights from three samples taken at different times during spawning

**Tab10a.csv** and **Tab10b.csv** give two versions of summary statistics of Lengths for each Age bin 

**Tab11_xx** gives a summary of the parameter estimates from model **xx**

**strat_length_props.csv** gives a table of stratified estimates of length proportions


### flat_data

This folder consists of .csv representations of all input data.  These have been copied directly from input data
provided, with no additional corrections applied.

**Event1.csv** and **Event2.csv** give capture data from the two events in the Mark-Recapture experiment

**spawn_sample.csv** gives Age, Length, and Weight during spawning

**ASL_2021.csv** gives Age, Length, and Weight from the previous project (2021), for the sake of comparison


### posts

This folder exists to house large .Rdata objects specific to full posteriors of Bayesian model output.  This is
currently listed in the `.gitignore` document associated with this repository, and is likely not found on the Github 
page.  However, they may be re-created by running the associated model code (which may be time-consuming!)


### one additional .Rdata file

**sensitivity_outputs.Rdata** contains summarized output from the long-running script in 2a_lengthage_sensitivity.R