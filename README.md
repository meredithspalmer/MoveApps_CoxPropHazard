# Cox Proportional Hazards Survival App

Github repository:

https://github.com/meredithspalmer/MoveApps_CoxPropHazard

## Description

Perform Cox Proportional Hazards survival analyses and optional stratified group comparisons, with automatic fallback to Firth's penalized regression when separation is detected.

These analyses can be performed across an entire dataset, within defined time periods, and/or across user-defined data subsets.

## Documentation

This app implements Cox Proportional Hazards (CoxPH) survival modelling. It produces hazard ratio tables, survival curves, Kaplan-Meier survival curves, cumulative hazard plots, Schoenfeld residual diagnostics, and stratified comparisons of per-group survival for different user-defined attributes. 

**Cox Proportional Hazards Model:** The CoxPH model is a semi-parametric method used to estimate the relationship between one or more covariates and survival time, specifically, the risk (hazard) of experiencing an event (e.g., death) at any point in time given that the individual has survived up to that point. This analysis allows for:

- *Right-censoring*, where the exact time of death is unknown for some individuals because they are still alive at the end of the study period, lost to follow-up (e.g., collar failure), or exit the study period alive for other reasons.
- *Staggered entry* (also called left truncation or delayed entry), where individuals enter the study at different times rather than all starting at the same baseline.

**CoxPH Assumptions**: The Cox PH model assumes that the hazard ratio between any two individuals remains *constant over time* (proportional hazards). The app tests this assumption by generating Schoenfeld residuals and will generate a warning if the assumption is violated for any covariate. Users should interpret results cautiously when this assumption does not hold. 

**Firth's Penalized Cox Model:** When complete separation is detected in the standard CoxPH model (i.e., a covariate perfectly predicts survival or death), the app automatically falls back to Firth's penalized Cox regression (`coxphf`). This method applies a likelihood penalization that produces finite, unbiased estimates in cases where maximum likelihood estimation fails to converge. All outputs clearly indicate whether Firth's or standard Cox was used. Separation detection and Firth's fallback applies both to the main Cox model and to all stratified group comparisons, including pairwise comparisons.

**Stratified Comparisons:** When grouping variables are specified, the app produces stratified Kaplan-Meier survival curves, cumulative hazard plots, and likelihood ratio, Wald, and score test results for each covariate. Pairwise comparisons with Bonferroni correction are computed when more than two groups are present. Continuous covariates (e.g., weight) are detected automatically and handled with Cox model tests and a forest plot only — KM curves and group-level outputs are not produced for continuous predictors, as these require discrete groups.

**Cumulative Hazard Plots:** These plots depict the total accumulated risk (e.g., the expected number of mortalities) the population has experienced up to a specific time.

**Forest Plots:** Visual summaries of hazard ratios and 95% confidence intervals for each covariate, with log-scaled axes to reflect the multiplicative nature of hazard ratios.

**Predicted Survival at Covariate Means:** Users have the option to generate survival outputs estimated for a hypothetical individual with mean values on all covariates, providing a useful population-level summary.

**Schoenfeld Residuals:** Users have the option to generate diagnostic plots are used to assess whether the proportional hazards assumption holds. Random scatter around zero across time indicates the assumption is satisfied; a systematic trend may indicate a violation.

**Mortality plots**: Users have the option to generate diagnostic plots showing monthly mortality across the study period.

Data subsetting:

- Data can be subset up to two times based on specific variables (e.g., only females or only adult females) — this then allows the users to perform comparisons within this smaller subset (e.g., comparing survival of adult females with different collar types).
- Users can also enter additional information defining survival years and animal life stages for comparison across these variables. For life stages, data must contain the column "animal\_birth\_hatch\_year" and a data frame mapping age to life stage must be uploaded.
- Furthermore, users can define a study period, censor data to exclude post-capture mortality events, and specify how to handle missing time-stamp information.

Data pre-processing includes:

- Removing or updating empty or invalid data (according to user specification)
- Flagging and handling marked outliers and test data
- Checking for logical errors (e.g., start date after end date)
- Subsetting data to defined study period
- Subsetting data based on attributes according to user specifications

The app then summarizes the data into a per-individual table and figure containing:

- Duration within the study period (accounting for staggered entry and censoring)
- Survival event indicator (1 = death occurred during study period, 0 = censored or survived)
- Gaps in tracking are also indicated (i.e., if an individual was collared and then recollared at a later date)

Finally, the app fits a Cox Proportional Hazards model on the cleaned dataset, generating:

- A table of hazard ratios with 95% confidence intervals and p-values
- A predicted survival curve at covariate means (with 95% confidence intervals) (optional)
- A cumulative hazard curve at covariate means (optional)
- A forest plot of hazard ratios
- Schoenfeld residual diagnostic plots (optional)

When stratified group comparisons are enabled, the app additionally produces per-covariate:

- Likelihood ratio, Wald, and score test results
- Pairwise hazard ratio comparisons (Bonferroni corrected; >2 groups only)
- Stratified Kaplan-Meier survival curves (categorical covariates only)
- Stratified cumulative hazard plots (categorical covariates only)
- Per-group median survival summary tables (categorical covariates only)
- Annual survival rate bar charts (categorical covariates only; when survival years are defined)

### Application scope

#### Generality of app usability

This app has currently been tested on mammals and birds with N\_individuals > 50, N\_deaths > 1, and study duration >1 year, but should be applicable to any dataset containing mortality data of sufficient sample size (see below).

This app allows for staggered entry during the defined study period and for censored data (individuals that are "lost" from the study, e.g., due to equipment failure, and individuals that survive the study period).

#### Required data properties

**Events**: This app can only be used if mortality information (indication of event along with an associated end date) is captured in one of the following columns: `death_comments`, `deployment_end_comments`, `deployment_end_type`, `mortality_location_filled`.

The app will produce a warning and terminate if none of the individuals in the study experienced a mortality event during the study period or data subset.

**Covariates**: The CoxPH model requires at least one covariate to be specified. Covariates must be present as named columns in the data and must contain at least two non-NA levels or values after cleaning. Covariates with only one unique value (or that are entirely NA) will be dropped automatically with a warning. Continuous (numeric) covariates are supported in the main Cox model and will receive a forest plot; group-level outputs such as KM curves are only produced for categorical covariates.

**Sample size**: The required sample size for a CoxPH analysis depends primarily on the total number of deaths, rather than the total number of individuals. In general, the lower the event (death) rate, the higher the number of required individuals. A commonly used rule of thumb is at least 10 events per covariate.

This app does *not* perform a power analysis prior to performing the survival analyses. However, if fewer than 10 mortality events are detected, the app will generate a warning that the model may have low statistical power, potentially resulting in unreliable estimates and poor predictive power.

Note that a larger sample size is required for stratified group comparisons.

**Subsetting data**: Users can select what subset of individuals to include in the study (only females, only individuals alive within a specific survival year, individuals wearing certain collar types, etc.). Data can be subset based on up to two attributes.

**Comparison types**: This app enables the user to identify up to three grouping variables for stratified comparisons. These grouping classifications are cleaned and standardized for capitalization and white-space errors (e.g., "male" and "Male" will register as the same class) but note that other misspecifications (e.g., "male " vs. "M") will be treated as separate groups. The column name of the group comparison value must be spelled and formatted exactly as found in the MoveBank dataset. 

**Survival years**: Users can define a "survival year" for analyses, different from a calendar year in that this period extends from when an animal is typically born to the end of their first year.

**Life stages**: If the user wants to calculate life stage, the input data must contain the column "animal\_birth\_hatch\_year". Users must also upload auxiliary information (see template: https://github.com/meredithspalmer/MoveApps_Survival/blob/master/animal_birth_hatch_year_table.csv) that links individual animal ages to species-specific life stages.

### Input type

`move2::move2_loc`

### Output type

`move2::move2_loc`

### Artefacts

**Tracking history:** Figure (`tracking_history.png`) detailing the start and end dates of each individual during the tracking period (with gaps in collaring noted). *Optional*

**Mortality plot:** Diagnostic plot (`monthly_mortality.png`) depicting mortality rate per month across the study period. *Optional*

**Cox model hazard ratio table:** Table (`coxph_table.csv`) of hazard ratios, 95% confidence intervals, p-values, and an indicator of whether standard Cox or Firth's penalized Cox was used.

**Survival at means:** Table (`surv_at_means.csv`) and plot (`surv_at_means_plot.png`) of the predicted survival function for a hypothetical individual with mean covariate values, with 95% confidence intervals. *Optional* 

*Cumulative hazard at means:* Plot (`cum_hazard_at_means_plot.png`) depicting total accumulated risk over time for a hypothetical individual at covariate means. *Optional* 

*Schoenfeld residuals:* Plot (`schoenfeld_residuals.png`) for assessing the proportional hazards assumption. One panel is generated per covariate. A note is appended when Firth's model was used, as residuals are approximate in that case. *Optional* 

The following artefacts are produced **per grouping covariate** (where `<cov>` is the covariate name):

*Model test statistics:* Table (`model_tests_<cov>.csv`) with likelihood ratio, Wald, and score test chi-squared statistics, degrees of freedom, and p-values.

**Stratified Cox hazard ratios:** Table (`cox_hr_<cov>.csv`) of hazard ratios, 95% confidence intervals, and p-values from a model for the selected covariate. Indicates whether standard Cox or Firth's model was used.

**Stratified forest plot:** Plot (`forest_<cov>.png`) depicting hazard ratios and 95% confidence intervals on a log scale, with significant estimates highlighted.

**Pairwise comparisons:** Table (pairwise_<cov>.csv) with pairwise Wald p-values and Bonferroni-corrected p-values for each pair of groups. Rows indicate whether standard Cox or Firth's model was used for each pair. *Categorical covariates, >2 groups only*

**Stratified Cox hazard ratios:** Table (`cox_hr_<cov>.csv`) of hazard ratios, 95% confidence intervals, and p-values from a model stratified by the selected covariate.

**Stratified KM survival curves:** Plot (`km_by_<cov>.png`) depicting Kaplan-Meier survival curves by group, with 95% confidence intervals, risk table, and Cox likelihood ratio test result in the subtitle. *Categorical covariates only* 

**Per-group survival summary:** Table (`km_summary_<cov>.csv`) with per-group sample sizes, number of events, median survival time, and estimated survival at one year. *Categorical covariates only*

**Cumulative hazard by group:** Plot (`cumhaz_by_<cov>.png`) depicting stratified cumulative hazard curves; parallel lines support the proportional hazards assumption. *Categoical covariates only* 

*Annual survival rates:* Plot (`annual_surv_by_<cov>.png`) and table (`annual_surv_<cov>.csv`) depicting survival rates per group per survival year. *Categorical covariates only; when survival years are defined*

### Settings

`Start date` and `End date`: Temporal limits of the study period. If left null (default), the analysis will encompass the entire study period. Useful for defining a study year. Unit: `date`.

`Fix empty start times` and `Fix empty end times`: Defines how the app handles NA dates (`deploy_on_timestamp` and `deploy_off_timestamp`) at the beginning and end of the study. Options include:

- Replacing with the first/last recorded timestamp (default)
- Replacing with the current date (for end times only)
- Removing missing data

`Censor capture-related mortality`: If capture-related mortality is a concern, this setting allows users to define a number of days post-capture to exclude from the overall analysis. Default is no censoring. Unit: `days`.

`Perform analysis on subset of data` and `Define subset condition`: If the user wants to perform analyses on a subset of the data, they can define the grouping parameter to split on (subset condition; e.g., "sex") and the group of interest to retain (subset definition; e.g., "f"). Default setting is to include all data. There are up to two subsets users can define.

`Cox covariates 1–3` and `Reference levels 1–3`: The grouping variables to include as covariates in the Cox model and for stratified comparisons. Up to three covariates can be specified simultaneously. For each covariate, users can optionally define a reference level (e.g., "f" for female as reference in a sex comparison). If a specified reference level is not found in the data, the default level ordering will be used. Continuous covariates (e.g., weight) are automatically detected and receive Cox model tests and a forest plot only. 

`Survival year start date`: If comparing across survival years (see above), the user can define the day and month that each 'survival year' begins. The code assumes a year runs 365(6) days. Default is null. Unit: `date`.

`Animal birth/hatch year definitions`: Optional auxiliary file a user can upload if they are running analyses by survival year and wish to classify individual life stage within a specific year. This file maps animal age to user-defined life stages. A template can be found here: https://github.com/meredithspalmer/MoveApps_Survival/blob/master/animal_birth_hatch_year_table.csv. App expects files in `.csv` format.

`Monthly mortality plots`: Create optional plots showcasing monthly mortality across years. Default is no plot. Options are include/exclude.

### Most common errors

Please document and send errors to [mspalmer.zool@gmail.com](mailto:mspalmer.zool@gmail.com).

### Null or error handling

See null handling outlined in *Settings*.

The app will log warnings quantifying reductions in sample size due to data cleaning and censoring, as well as indicating how many null timestamps were updated or removed.

The app will log warnings when data includes fewer than 10 mortality events and will terminate if there are no mortality events.

The app will automatically detect and log a fallback to Firth's penalized Cox model if separation is detected during standard Cox model fitting, stratified covariate fitting, or pairwise comparisons. All outputs will indicate which model was used.

Covariates that are entirely NA or have only a single unique value after cleaning will be dropped automatically with a warning, and the model will proceed with the remaining covariates.

Errors may arise due to how mortality events are recorded (what metadata flags mortality and whether mortality is associated with a specific end date). If mortality is logged in a column not mentioned in this documentation, it will not be recognized.