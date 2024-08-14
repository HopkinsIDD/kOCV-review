# kOCV-review

This is the repository for the paper titled "Protection from Killed Whole-Cell Cholera Vaccines: A Systematic Review and Meta-Analysis". This systematic review builds upon a previous systematic review by Bi et al. titled "Protection against cholera from killed whole-cell oral cholera vaccines: a systematic review and meta-analysis" published on Lancet Infectious Disease in 2017.  

## Data 

* `ve_meansvar.csv`: this file contains data included in the previous search by Bi et al. 2017. Each estimate in this file represents the pooled estimate over the study's follow-up period without stratification by time since vaccination.
* `FollowupPeriodData.csv`: this file contains data included in the previous search by Bi et al. 2017. The estimates in this file are stratified by time since vaccination.
* `extracted_SR_data_2017_2024.csv`: contains data detected through new search.
* `effectiveness_studies.csv`, `efficacy_studies.csv`: contain meta-data of each included study in the current systematic review.
* `propUnder5_data.csv`: contains age-stratified information on each included study.
* `QA_RCT_resolved.csv`: contains the risk of bias assessment of all the included randomized clinical trials after resolving conflicts between two reviewers.

## Code 
* `OCV_SR_update.Rmd`: code used to conduct all the analyses and generation of tables and figures presented in the manuscript.
* `panelPlot_colored.R`, `utils_OCV_SR.R`: functions used in the `OCV_SR_update.Rmd`.

