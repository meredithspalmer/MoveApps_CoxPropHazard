library(move2)
library(survival)
library(survminer)
library(ggplot2)   
library(dplyr)
library(lubridate)
library(stringr)
library(sf)
library(forcats)
library(tidyr)
library(purrr) 
library(viridis)
library(ggpubr)
library(coxphf)
library(broom)
library(ggsurvfit)
library(scales)
library(patchwork)

# logger.fatal(), logger.error(), logger.warn(), logger.info(), logger.debug(), logger.trace()

# Survival Function 
rFunction = function(data,  
                     time_period_start, 
                     time_period_end, 
                     censor_capture_mortality,
                     fix_na_start_times, 
                     fix_na_end_times,  
                     subset_condition_1,
                     subset_condition_define_1, 
                     subset_condition_2,
                     subset_condition_define_2, 
                     cox_covariate_1,
                     cox_covariate_2, 
                     cox_covariate_3, 
                     cox_covariate_1_ref, 
                     cox_covariate_2_ref, 
                     cox_covariate_3_ref, 
                     survival_yr_start,
                     animal_birth_hatch_year_table, 
                     life_table_days, 
                     calc_month_mort,
                     calc_tracking_history, 
                     calc_residuals, 
                     calc_artifacts_at_mean, 
                     add_cis, 
                     zoom_to_plot) {
  
  ## Load auxiliary data ------------------------------------------------------ 
  
  if(!is.null(animal_birth_hatch_year_table)){
    animal_birth_hatch_year_table <- read.csv(getAuxiliaryFilePath("animal_birth_hatch_year_table"))
    logger.info("Auxiliary birth/hatch year data loaded.")
  }
  
  
  ## Cleaning and cropping ----------------------------------------------------
  
  data <- dplyr::filter(data, !sf::st_is_empty(data))       # Exclude empty locations
  data <- mt_filter_unique(data)                            # Exclude marked outliers 
  data <- data %>% filter_track_data(is_test == FALSE)      # Exclude data marked "test"
  
  
  ## Aggregate across multiple deployments (where present) ---
  
  # Extract event-level data 
  events <- data |>
    as_tibble() |>
    dplyr::select(any_of(c(
      "deployment_id",
      "individual_local_identifier",   
      "timestamp")))
  
  # Columns to drop from track data (administrative/redundant metadata)
  exclude_cols <- c("acknowledgements",
                    "citation",
                    "contact_person_name",
                    "deployment_local_identifier",
                    "group_id",
                    "has_quota",
                    "i_am_collaborator",
                    "i_am_owner",
                    "i_can_see_data",
                    "i_have_download_access",
                    "individual_comments",
                    "license_type",
                    "main_location",
                    "principal_investigator_name",
                    "serial_no",
                    "study_number_of_deployments",
                    "study_objective",
                    "study_permission",
                    "study_site",
                    "study_type",
                    "suspend_license_terms",
                    "tag_comments",
                    "tag_failure_comments",
                    "tag_local_identifier",
                    "tag_number_of_deployments",
                    "taxon_canonical_name",
                    "taxon_detail",
                    "taxon_ids",
                    "there_are_data_which_i_cannot_see", 
                    "is_test")
  
  tracks <- mt_track_data(data) |>
    mutate(mortality_location_filled = if_else(
      is.na(mortality_location) | st_is_empty(mortality_location),
      0L, 1L)) |>
    dplyr::select(-any_of(exclude_cols))
  
  # Join track attributes to every event row
  use_deployment_join <- all("deployment_id" %in% names(events),
                             "deployment_id" %in% names(tracks)) &&
    any(!is.na(events$deployment_id)) &&
    any(!is.na(tracks$deployment_id))
  
  if (use_deployment_join) {
    events_with_ind <- events |> left_join(tracks, by = "deployment_id")
    
  } else {
    if (!"individual_local_identifier" %in% names(events)) {
      logger.fatal("Cannot join: neither deployment_id nor individual_local_identifier is available in events")
    }
    logger.info("Joining on individual_local_identifier (deployment_id join not possible)")
    events_with_ind <- events |> left_join(tracks, by = "individual_local_identifier")
  }
  
  events_with_ind <- events_with_ind |>
    relocate(any_of(c("individual_id", "individual_local_identifier",
                      "deployment_id", "timestamp")),
             .before = everything())
  
  # Infer column groups dynamically 
  cols_timestamp_min <- c("timestamp_first_deployed_location", "deploy_on_timestamp")
  cols_timestamp_max <- c("timestamp_last_deployed_location",  "deploy_off_timestamp")
  
  cols_already_handled <- c(
    "individual_id", "individual_local_identifier", "deployment_id",
    "timestamp", "geometry", "mortality_location", "mortality_location_filled",
    cols_timestamp_min, cols_timestamp_max
  )
  
  # Character/factor: collapse unique values
  cols_categorical <- events_with_ind |>
    select(where(~ is.character(.) || is.factor(.))) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  # Numeric and units: mean
  cols_numeric <- events_with_ind |>
    select(where(~ is.numeric(.) || inherits(., "units"))) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  # Logical: any() (that is, TRUE if any event flagged TRUE)
  cols_logical <- events_with_ind |>
    select(where(is.logical)) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  # POSIXct/Date not already in explicit min/max lists: take first non-NA value
  cols_timestamp_first <- events_with_ind |>
    select(where(~ inherits(., "POSIXct") || inherits(., "Date"))) |>
    select(-any_of(c(cols_already_handled, cols_timestamp_min, cols_timestamp_max))) |>
    names()
  
  # sfc geometry columns: take first non-empty geometry per individual
  cols_geometry <- events_with_ind |>
    select(where(~ inherits(., "sfc"))) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  logger.info(paste("Categorical columns:",    paste(cols_categorical,     collapse = ", ")))
  logger.info(paste("Numeric/units columns:",  paste(cols_numeric,         collapse = ", ")))
  logger.info(paste("Logical columns:",        paste(cols_logical,         collapse = ", ")))
  logger.info(paste("Timestamp (first) cols:", paste(cols_timestamp_first, collapse = ", ")))
  logger.info(paste("Geometry (first) cols:",  paste(cols_geometry,        collapse = ", ")))
  
  
  # Summarise to individual level
  na_like <- function(x) x[NA_integer_]
  
  # Helper function: first non-empty geometry or an empty point if all empty/NA
  first_non_empty_geom <- function(x) {
    non_empty <- x[!is.na(x) & !st_is_empty(x)]
    if (length(non_empty) > 0) non_empty[[1]] else st_point()
  }
  
  summary_table <- events_with_ind |>
    group_by(individual_id, individual_local_identifier) |>
    summarise(
      first_timestamp = min(as.Date(timestamp), na.rm = TRUE),
      last_timestamp  = max(as.Date(timestamp), na.rm = TRUE),
      n_locations     = n(),
      
      n_deployments = if ("deployment_id" %in% names(events_with_ind))
        n_distinct(deployment_id[!is.na(deployment_id)]) else 1L,
      
      mortality_location_filled =
        if ("mortality_location_filled" %in% names(events_with_ind))
          as.integer(any(mortality_location_filled >= 1, na.rm = TRUE)) else NA_integer_,
      
      across(any_of(cols_timestamp_min),
             ~ if (all(is.na(.))) na_like(.) else min(., na.rm = TRUE)),
      
      across(any_of(cols_timestamp_max),
             ~ if (all(is.na(.))) na_like(.) else max(., na.rm = TRUE)),
      
      # Deployment-level timestamp constants — take first non-NA value
      across(any_of(cols_timestamp_first),
             ~ if (all(is.na(.))) na_like(.) else first(na.omit(.))),
      
      across(any_of(cols_categorical),
             ~ if (all(is.na(.))) NA_character_
             else str_c(unique(.[!is.na(.)]), collapse = " | ")),
      
      across(any_of(cols_numeric),
             ~ if (all(is.na(.))) na_like(.) else mean(., na.rm = TRUE)),
      
      # Logical — TRUE if any event is TRUE, NA if all NA
      across(any_of(cols_logical),
             ~ if (all(is.na(.))) NA else any(., na.rm = TRUE)),
      
      # Geometry — first non-empty geometry, or empty point
      across(any_of(cols_geometry),
             ~ st_sfc(first_non_empty_geom(.), crs = st_crs(.))),
      
      .groups = "drop") |>
    
    mutate(across(any_of(cols_categorical), ~ if_else(. == "", NA_character_, .)),
           across(any_of(c("deploy_on_timestamp", "deploy_off_timestamp")), as.Date))
  
  # Quick checks 
  unhandled_cols <- setdiff(names(events_with_ind),
                            c(cols_already_handled, cols_categorical, cols_numeric,
                              cols_logical, cols_timestamp_first, cols_geometry,
                              "first_timestamp", "last_timestamp", "n_locations", "n_deployments"))
  
  if (length(unhandled_cols) > 0) {
    logger.warn(paste("Columns not summarised (check types):", paste(unhandled_cols, collapse = ", ")))
  }
  
  
  ## Inject NA columns for optional mortality fields absent from this dataset --- 
  if (!"mortality_date" %in% names(summary_table)) summary_table$mortality_date <- NA_Date_
  if (!"mortality_type" %in% names(summary_table)) summary_table$mortality_type <- NA_character_
  if (!"death_comments" %in% names(summary_table)) summary_table$death_comments <- NA_character_
  if (!"deployment_end_comments" %in% names(summary_table)) summary_table$deployment_end_comments <- NA_character_
  if (!"deployment_end_type" %in% names(summary_table)) summary_table$deployment_end_type <- NA_character_
  

  ## Clean dates ---
  
  # Start times  
  if(fix_na_start_times == "timestamp"){
    summary_table <- summary_table %>% 
      mutate(missing_timestamp_start = is.na(deploy_on_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_on_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_on_timestamp = if_else(
        is.na(deploy_on_timestamp),
        as.Date(first_timestamp),
        deploy_on_timestamp)) %>% 
      dplyr::select(-missing_timestamp_start)
    
    if (n_missing > 0) {
      logger.info(
        sprintf("Warning: Replaced %d missing deploy_on_timestamp value%s with first_timestamp.",
                n_missing, if (n_missing == 1) "" else "s"))
    }
  } 
  
  if (fix_na_start_times == "remove"){
    n_missing <- sum(is.na(summary_table$deploy_on_timestamp))
    summary_table <- summary_table %>% filter(!is.na(deploy_on_timestamp))
    
    if (n_missing > 0) {
      logger.info(sprintf("Warning: Removed %d deploy_on_timestamp value%s that were NA.", n_missing,
                          if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  }
  
  # End times 
  if(fix_na_end_times == "timestamp"){
    summary_table <- summary_table %>%
      mutate(missing_timestamp_end = is.na(deploy_off_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_off_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_off_timestamp = if_else(is.na(deploy_off_timestamp),
                                            as.Date(last_timestamp),
                                            deploy_off_timestamp)) %>% 
      dplyr::select(-missing_timestamp_end)
    
    if (n_missing > 0) {
      logger.info(
        sprintf("Warning: Replaced %d missing deploy_off_timestamp value%s with last_timestamp.", 
                n_missing, if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  } 
  
  if (fix_na_end_times == "systime"){
    summary_table <- summary_table %>%
      mutate(missing_timestamp_end = is.na(deploy_off_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_off_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_off_timestamp = if_else(is.na(deploy_off_timestamp),
                                            Sys.Date(), 
                                            deploy_off_timestamp))%>% 
      dplyr::select(-missing_timestamp_end)
    
    if (n_missing > 0) {
      logger.info(sprintf("Warning: Replaced %d missing deploy_off_timestamp value%s with current date.",
                          n_missing, if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  } 
  
  if (fix_na_end_times == "remove"){
    n_missing <- sum(is.na(summary_table$deploy_off_timestamp))
    summary_table <- summary_table %>% filter(!is.na(deploy_off_timestamp))
    
    if (n_missing > 0) {
      logger.info(
        sprintf("Warning: Removed %d deploy_off_timestamp and/or deploy_on_timestamp value%s that were NA.", 
                n_missing, if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  } 
  
  # Remove data for individuals where "deploy_off_timestamp" occurs before "deploy_on_timestamp" 
  n_original <- nrow(summary_table) 
  summary_table <- summary_table %>%
    filter(deploy_off_timestamp >= deploy_on_timestamp)
  n_removed <- n_original - nrow(summary_table)  
  
  if (n_removed > 0) {
    logger.info(
      sprintf("Warning: Removed %d individual%s where deploy_off_timestamp < deploy_on_timestamp.",
              n_removed, if (n_removed == 1) "" else "s"))
  }
  
  
  ## Crop data to user-defined temporal windows ---
  
  # Removed censored data (mortalities within set period of capture)
  if(censor_capture_mortality > 0){
    n_before <- nrow(summary_table)
    
    summary_table <- summary_table %>%
      mutate(raw_deploy_on_timestamp = deploy_on_timestamp) %>%  
      mutate(censor_cutoff = deploy_on_timestamp + lubridate::days(censor_capture_mortality)) %>%
      mutate(remove_due_to_early_end = !is.na(deploy_off_timestamp) & deploy_off_timestamp <= censor_cutoff) %>%
      filter(!remove_due_to_early_end) %>%
      mutate(deploy_on_timestamp = censor_cutoff) %>%
      select(-censor_cutoff, -remove_due_to_early_end)
    
    n_after  <- nrow(summary_table)
    n_removed <- n_before - n_after
    
    if (n_removed > 0) {
      logger.info(paste0("Warning: Removed ", n_removed, " individual(s) because deploy_off_timestamp occurred within ", censor_capture_mortality, " day(s) after deploy_on_timestamp"))
    } 
  }
  
  
  # Crop to study period of interest 
  
  # Save original deploy_off_time 
  summary_table <- summary_table %>% mutate(raw_deploy_off_timestamp = deploy_off_timestamp) 
  
  # Define window 
  effective_start <- if (is.null(time_period_start)) {
    min(summary_table$deploy_on_timestamp, na.rm = TRUE)
  } else {
    time_period_start
  }
  
  effective_end <- if (is.null(time_period_end)) {
    max(summary_table$deploy_off_timestamp, na.rm = TRUE)
  } else {
    time_period_end
  }
  
  # Run updates 
  if(!is.null(time_period_start) | !is.null(time_period_end)){
    
    # Crop to window 
    n_original <- nrow(summary_table) 
    summary_table <- summary_table %>%
      
      # Determine if the deployment overlaps study window 
      mutate(overlaps_study = deploy_on_timestamp <= effective_end & 
               deploy_off_timestamp  >= effective_start) %>%
      filter(overlaps_study | is.na(overlaps_study)) %>%    
      
      # Crop to window 
      mutate(first_timestamp = pmax(deploy_on_timestamp, effective_start, na.rm = TRUE),
             last_timestamp  = pmin(deploy_off_timestamp, effective_end,   na.rm = TRUE)) %>%
      
      # Clean 
      select(-overlaps_study) 
    
    n_removed <- n_original - nrow(summary_table)
    if (n_removed > 0) {
      logger.info(
        sprintf("Warning: %d record%s did not overlap the user-defined study window and were removed.",
                n_removed, if (n_removed == 1) "" else "s"))
    }
  } 
  
  
  ## Calculate entry time and exit time (for staggered entry) ---
  origin_date <- if (is.null(time_period_start) || is.na(time_period_start)) {
    min(summary_table$deploy_on_timestamp, na.rm = TRUE)
  } else {
    time_period_start
  }
  
  summary_table <- summary_table %>%
    mutate(analysis_entry_date = pmax(deploy_on_timestamp, effective_start, na.rm = TRUE),
           analysis_exit_date  = pmin(deploy_off_timestamp, effective_end, na.rm = TRUE),
           entry_time_days = as.numeric(difftime(analysis_entry_date, origin_date, units = "days")),
           exit_time_days  = as.numeric(difftime(analysis_exit_date,  origin_date, units = "days"))) 
  
  
  ## Calculate mortality indicator --------------------------------------------
  # Event = 1 if observed death, 0 if censored or survived 
  
  # Death comments to flag
  positive_pattern <- "dead|death|died|cod|predation|predator|vehicle|collision|killed|poach|poached|shot|hunt|harvest|harvested|mortality"
  
  # Check which columns are present 
  has_death_comments          <- "death_comments"          %in% names(summary_table)
  has_deployment_end_comments <- "deployment_end_comments" %in% names(summary_table)
  has_deployment_end_type     <- "deployment_end_type"     %in% names(summary_table)
  has_mortality_type          <- "mortality_type"          %in% names(summary_table)
  has_mortality_location      <- "mortality_location_filled" %in% names(summary_table)
  
  # Construct mortality event
  summary_table <- summary_table %>%
    
    # Initialise mortality event
    mutate(mortality_event = NA_real_) %>%
    
    # Identify survivors (individuals who last beyond study)
    mutate(survived_beyond_study = !is.na(raw_deploy_off_timestamp) &
             raw_deploy_off_timestamp > as.Date(effective_end),
           mortality_event = if_else(survived_beyond_study, 0L, mortality_event)) %>%
    
    # Update optional columns to remove ambiguity for survivors
    { if (has_death_comments)
      mutate(., death_comments = if_else(survived_beyond_study, "survived beyond study", death_comments))
      else . } %>%
    { if (has_deployment_end_comments)
      mutate(., deployment_end_comments = if_else(survived_beyond_study, "survived beyond study", deployment_end_comments))
      else . } %>%
    { if (has_deployment_end_type)
      mutate(., deployment_end_type = if_else(survived_beyond_study, "survived beyond study", deployment_end_type))
      else . } %>%
    { if (has_mortality_location)
      mutate(., mortality_location_filled = if_else(survived_beyond_study, 0L, mortality_location_filled))
      else . } %>%
    
    # A. death_comments keywords
    { if (has_death_comments)
      mutate(., mortality_event = case_when(
        str_detect(tolower(death_comments), "\\bnot\\b")      ~ 0L,
        str_detect(tolower(death_comments), positive_pattern) ~ 1L,
        mortality_event == 1L                                  ~ 1L,
        TRUE                                                   ~ mortality_event))
      else . } %>%
    
    # B. deployment_end_comments keywords
    { if (has_deployment_end_comments)
      mutate(., mortality_event = case_when(
        mortality_event == 1L ~ 1L,
        str_detect(tolower(deployment_end_comments), positive_pattern) ~ 1L,
        TRUE ~ mortality_event))
      else . } %>%
    
    # C. mortality_type is filled (non-NA)
    { if (has_mortality_type)
      mutate(., mortality_event = case_when(
        mortality_event == 1L  ~ 1L,
        !is.na(mortality_type) ~ 1L,
        TRUE                   ~ mortality_event))
      else . } %>%
    
    # D. mortality_location_filled
    { if (has_mortality_location)
      mutate(., mortality_event = case_when(
        mortality_location_filled >= 1 ~ 1L,
        mortality_event == 1L          ~ 1L,
        TRUE                           ~ mortality_event))
      else . } %>%
    
    # E. deployment_end_type
    { if (has_deployment_end_type)
      mutate(., mortality_event = case_when(
        mortality_event == 1L ~ 1L,
        str_detect(tolower(deployment_end_type), "\\bdead\\b|\\bdeath\\b")           ~ 1L,
        tolower(deployment_end_type) %in%
          c("removal", "other", "unknown", "survived beyond study")                  ~ 0L,
        is.na(deployment_end_type) & is.na(mortality_event)                          ~ 0L,
        TRUE                                                                         ~ mortality_event))
      else
        mutate(., mortality_event = if_else(is.na(mortality_event), 0L, mortality_event))
    } %>%
    
    # Final censoring: remaining NA → 0 if deploy_off_timestamp is present
    mutate(mortality_event = if_else(
      is.na(mortality_event) & !is.na(deploy_off_timestamp),
      0L,
      mortality_event)) %>%
    
    select(-survived_beyond_study) |>
    relocate(mortality_event, .after = any_of("deployment_end_type"))
  
  # Error out: no deaths
  n_mort_events <- sum(summary_table$mortality_event == 1L, na.rm = TRUE)
  if (n_mort_events == 0) {
    logger.fatal("Cannot run survival analysis: no mortality events detected.")
  }
  
  # Warning: few deaths
  if (n_mort_events <= 10) {
    logger.warn(sprintf("Few (%d) deaths detected across entire dataset. Particularly if data is further subset, model may have low statistical power. This could potentially result in unreliable estimates and poor predictive capacity", n_mort_events))
  }
  
  
  ## Calculate survival years (if selected) -----------------------------------
  
  yearly_survival <- NULL
  
  if (!is.null(survival_yr_start)) {
    
    # Survival year start 
    start_month <- month(survival_yr_start)
    start_day   <- mday(survival_yr_start)
    
    # Handle invalid dates (e.g. Feb 29 in non-leap years)
    safe_make_date <- function(year, month, day) {
      date <- suppressWarnings(make_date(year, month, day))
      if (is.na(date)) {
        ym   <- ymd(sprintf("%d-%02d-01", year, month)) %m-% months(1)
        date <- ceiling_date(ym, "month") - days(1)
      }
      date
    }
    
    # Return survival year and period boundaries for a given date
    get_survival_period <- function(date) {
      if (is.na(date)) return(tibble(survival_year = NA_integer_,
                                     period_start  = NA_Date_,
                                     period_end    = NA_Date_))
      y <- year(date)
      period_start_this_yr <- safe_make_date(y, start_month, start_day)
      
      if (date >= period_start_this_yr) {
        survival_year <- y
        period_start  <- period_start_this_yr
        period_end    <- safe_make_date(y + 1, start_month, start_day) - days(1)
      } else {
        survival_year <- y - 1
        period_start  <- safe_make_date(y - 1, start_month, start_day)
        period_end    <- safe_make_date(y,     start_month, start_day) - days(1)
      }
      tibble(survival_year = survival_year,
             period_start  = period_start,
             period_end    = period_end)
    }
    
    # Vectorized helpers
    get_survival_year <- function(date) get_survival_period(date)$survival_year
    
    # Determine survival year range
    date_range <- events_with_ind %>%
      summarise(min_ts = min(timestamp, na.rm = TRUE),
                max_ts = max(timestamp, na.rm = TRUE))
    
    min_year       <- get_survival_year(date_range$min_ts)
    max_year       <- get_survival_year(date_range$max_ts)
    possible_years <- seq(min_year, max_year, by = 1)
    
    logger.info(sprintf("Survival years found in data: %d to %d (%d years total)",
                        min_year, max_year, length(possible_years)))
    
    # Create all possible survival periods
    possible_periods <- tibble(survival_year = possible_years) %>%
      mutate(period_info = map(survival_year, ~ {
        start <- safe_make_date(.x,     start_month, start_day)
        end   <- safe_make_date(.x + 1, start_month, start_day) - days(1)
        tibble(period_start = start, period_end = end)
      })) %>%
      unnest(period_info)
    
    # Infer which columns to carry into yearly_survival
    has_mortality_date          <- "mortality_date"          %in% names(summary_table)
    has_mortality_type          <- "mortality_type"          %in% names(summary_table)
    has_death_comments          <- "death_comments"          %in% names(summary_table)
    has_deployment_end_comments <- "deployment_end_comments" %in% names(summary_table)
    has_deployment_end_type     <- "deployment_end_type"     %in% names(summary_table)
    
    # Exclude join keys and period-specific columns that are recalculated per year
    cols_exclude_from_carry <- c(
      # Join keys
      "individual_id", "individual_local_identifier",
      # Recalculated per survival year
      "entry_time_days", "exit_time_days",
      "analysis_entry_date", "analysis_exit_date",
      # Raw timestamps kept separately — clipped versions used in yearly output
      "raw_deploy_on_timestamp", "raw_deploy_off_timestamp"
    )
    
    carry_cols <- setdiff(names(summary_table), cols_exclude_from_carry)
    logger.info(paste("Columns carried into yearly_survival:", paste(carry_cols, collapse = ", ")))
    
    # Infer individual-level constants by checking n_distinct per individual
    static_cols <- summary_table %>%
      dplyr::select(individual_id, all_of(carry_cols)) %>%
      group_by(individual_id) %>%
      summarise(across(everything(), n_distinct), .groups = "drop") %>%
      dplyr::select(-individual_id) %>%
      summarise(across(everything(), max)) %>%
      pivot_longer(everything(), names_to = "col", values_to = "n_distinct") %>%
      filter(n_distinct == 1) %>%
      pull(col)
    
    logger.info(paste("Individual-level constant columns (used for crossing):",
                      paste(static_cols, collapse = ", ")))
    
    # Build individual × year table
    yearly_survival <- summary_table %>%
      
      # Cross individual-level constants with all possible survival periods
      dplyr::select(individual_id, individual_local_identifier,
                    any_of(static_cols)) %>%
      distinct() %>%
      crossing(possible_periods) %>%
      
      # Join back all remaining carry columns
      left_join(summary_table %>%
                  dplyr::select(individual_id, individual_local_identifier,
                                all_of(setdiff(carry_cols, static_cols))),
                by = c("individual_id", "individual_local_identifier")) %>%
      
      # Clip monitoring interval to the survival period
      mutate(monitor_start    = pmax(deploy_on_timestamp, period_start, na.rm = TRUE),
             monitor_end      = pmin(deploy_off_timestamp, period_end,  na.rm = TRUE),
             active_in_period = monitor_start <= monitor_end &
               !is.na(monitor_start) & !is.na(monitor_end)) %>%
      filter(active_in_period) %>%
      
      # Inject NA columns for any optional mortality columns that are absent,
      # so that case_when() can reference them without erroring
      { if (!has_mortality_date)          mutate(., mortality_date          = NA_Date_)     else . } %>%
      { if (!has_mortality_type)          mutate(., mortality_type          = NA_character_) else . } %>%
      { if (!has_death_comments)          mutate(., death_comments          = NA_character_) else . } %>%
      { if (!has_deployment_end_comments) mutate(., deployment_end_comments = NA_character_) else . } %>%
      { if (!has_deployment_end_type)     mutate(., deployment_end_type     = NA_character_) else . } %>%
      
      # Mortality flags per survival year
      mutate(
        entry_date = monitor_start,
        exit_date  = monitor_end,
        
        # Did mortality occur inside this survival year?
        # Priority 1: mortality_date present and inside the period
        # Priority 2: mortality_date absent or NA — fall back to deploy_off_timestamp
        died_this_year = case_when(
          mortality_event == 1L &
            !is.na(mortality_date) &
            mortality_date >= period_start &
            mortality_date <= period_end
          ~ TRUE,
          
          mortality_event == 1L &
            is.na(mortality_date) &
            !is.na(deploy_off_timestamp) &
            deploy_off_timestamp >= period_start &
            deploy_off_timestamp <= period_end
          ~ TRUE,
          
          TRUE ~ FALSE),
        
        mortality_event = as.integer(died_this_year),
        censored        = as.integer(!died_this_year),
        
        # Reported mortality date: prefer mortality_date if present, else deploy_off
        mortality_date_reported = case_when(
          died_this_year & !is.na(mortality_date) ~ as.Date(mortality_date),
          died_this_year                           ~ as.Date(deploy_off_timestamp),
          TRUE                                     ~ NA_Date_),
        
        # Carry mortality metadata only when flagging a death
        mortality_type          = if_else(died_this_year, mortality_type,          NA_character_),
        death_comments          = if_else(died_this_year, death_comments,          NA_character_),
        deployment_end_comments = if_else(died_this_year, deployment_end_comments, NA_character_),
        deployment_end_type     = if_else(died_this_year, deployment_end_type,     NA_character_),
        
        days_at_risk = as.integer(exit_date - entry_date) + 1L) %>%
      
      select(-active_in_period, -monitor_start, -monitor_end, -died_this_year) %>%
      
      # Calculate per-year entry/exit time in days (relative to each period_start)
      mutate(analysis_entry_date = entry_date,
             analysis_exit_date  = exit_date,
             entry_time_days     = as.numeric(difftime(entry_date, period_start, units = "days")),
             exit_time_days      = as.numeric(difftime(exit_date,  period_start, units = "days"))) %>%
      
      arrange(individual_id, survival_year)
    
  } else {
    logger.info("Survival years not calculated.")
  }
  
  
  ## Calculate life stages per year (if selected) -----------------------------
  
  # Note: this needs auxiliary file to be loaded (errors earlier in code upon loading) 
  # Note: this needs "survival_yr_start" to be defined 
  if(!is.null(animal_birth_hatch_year_table) && is.null(survival_yr_start)){
    logger.error("Calculating life-stage requires survival years to be defined.")
  } 
  
  if(!is.null(survival_yr_start) && !is.null(animal_birth_hatch_year_table)){
    
    # Confirm data exists 
    if (!"animal_birth_hatch_year" %in% names(yearly_survival)) {
      logger.error("Column 'animal_birth_hatch_year' does not exist in the data frame. Cannot compute life stage.")
    }
    
    # Calculate age and age_class
    yearly_survival <- yearly_survival %>%
      mutate(animal_birth_hatch_year = as.numeric(animal_birth_hatch_year),
             age = survival_year - animal_birth_hatch_year,
             age = as.integer(pmax(0, age)))
    
    # Prepare thresholds from your existing table  
    thresholds <- animal_birth_hatch_year_table %>%
      filter(!is.na(year_at_start)) %>%           
      arrange(year_at_start) %>%
      distinct(year_at_start, animal_life_stage)
    
    # Create a named vector for fast look-up
    stage_lookup <- setNames(thresholds$animal_life_stage,
                             thresholds$year_at_start)
    
    # Pull stage for animals with unknown birth year from lookup table
    na_stage <- animal_birth_hatch_year_table %>%
      filter(is.na(year_at_start)) %>%
      pull(animal_life_stage) %>%
      first(default = "adult")
    
    # Dynamic assignment 
    yearly_survival <- yearly_survival %>%
      mutate(matched_threshold = findInterval(age, thresholds$year_at_start),
             animal_life_stage_new = case_when(
               is.na(age)             ~ na_stage,           
               matched_threshold == 0 ~ "unknown",            
               TRUE                   ~ stage_lookup[matched_threshold]),
             animal_life_stage_new = coalesce(
               animal_life_stage_new,
               na_stage)) %>%
      select(-matched_threshold)
    
    # Log age class summary
    age_class_summary <- yearly_survival %>%
      count(animal_life_stage_new) %>%
      arrange(desc(n))
    
    logger.info(sprintf("Individuals by age class: %s",
                        paste(sprintf("%s (n=%d)", age_class_summary$animal_life_stage_new,
                                      age_class_summary$n),
                              collapse = ", ")))
    
  } else {
    logger.info("Life stages not calculated.")
  }
  
  
  ## Subset based on condition (if selected) ----------------------------------
  
  # Lookup: condition name -> list(col, label, coerce_fn, summary_ok, yearly_ok)
  subset_spec <- list(
    sex = list(col = "sex", label = "sex", coerce = identity, summary_ok = TRUE,  yearly_ok = TRUE),
    attachment_type = list(col = "attachment_type", label = "attachment type", coerce = identity, summary_ok = TRUE,  yearly_ok = TRUE),
    model = list(col = "model", label = "model", coerce = as.integer, summary_ok = TRUE, yearly_ok = TRUE),
    lifestage = list(col = "animal_life_stage_new", label = "life stage", coerce = identity, summary_ok = FALSE, yearly_ok = TRUE),
    survival_year = list(col = "survival_year", label = "survival year", coerce = as.integer, summary_ok = FALSE, yearly_ok = TRUE)
  )
  
  # Helper function 
  apply_subset <- function(condition, define, summary_table, yearly_survival, survival_yr_start) {
    if (is.null(condition)) return(list(summary_table = summary_table, yearly_survival = yearly_survival))
    
    spec <- subset_spec[[condition]]
    if (is.null(spec)) stop(paste("Unknown subset condition:", condition))
    
    logger.info(paste0("Subsetting by ", spec$label, " (", define, ")"))
    value <- spec$coerce(define)
    
    if (is.null(survival_yr_start)) {
      if (!spec$summary_ok) {
        logger.fatal("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      } else {
        summary_table <- summary_table %>% filter(.data[[spec$col]] == value)
      }
    } else {
      if (is.null(yearly_survival)) {          
        logger.fatal("yearly_survival is NULL despite survival_yr_start being set. This is unexpected.")
      }
      if (!spec$yearly_ok) {
        logger.fatal("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      } else {
        yearly_survival <- yearly_survival %>% filter(.data[[spec$col]] == value)
      }
    }
    
    list(summary_table = summary_table, yearly_survival = yearly_survival)
  }
  
  # Subset condition 1 ---
  result <- apply_subset(subset_condition_1, subset_condition_define_1, summary_table,
                         yearly_survival, survival_yr_start)
  summary_table    <- result$summary_table
  yearly_survival  <- result$yearly_survival
  
  if (!is.null(subset_condition_1)) {
    if (is.null(survival_yr_start) && nrow(summary_table) == 0) {
      logger.fatal("There are no individuals meeting the first subsetting condition.")
    } else if (!is.null(survival_yr_start) && nrow(yearly_survival) == 0) {
      logger.fatal("There are no individuals meeting the first subsetting condition.")
    }
  }
  
  # Subset condition 2 ---
  result <- apply_subset(subset_condition_2, subset_condition_define_2, summary_table,
                         yearly_survival, survival_yr_start)
  summary_table   <- result$summary_table
  yearly_survival <- result$yearly_survival
  
  if (!is.null(subset_condition_2)) {
    if (is.null(survival_yr_start) && nrow(summary_table) == 0) {
      logger.fatal("There are no individuals meeting both subsetting conditions.")
    } else if (!is.null(survival_yr_start) && nrow(yearly_survival) == 0) {
      logger.fatal("There are no individuals meeting both subsetting conditions.")
    }
  }
  
  
  ## Clean Cox covariates -----------------------------------------------------
  
  # Note: This sets cox_covariate_X to NULL if no cleaned values are available 
  
  # Comparison covariate 1 --- 
  if(!is.null(cox_covariate_1)) {
    logger.info(paste("Comparing across", cox_covariate_1, "..."))
    
    if(is.null(survival_yr_start)){
      logger.info("... with summary table")
      
      non_na_unique <- unique(na.omit(summary_table[[cox_covariate_1]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", cox_covariate_1, ", is entirely NA; no comparison is possible."))
        } else {
          logger.fatal(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_1, "; no comparison is possible."))
        }
        cox_covariate_1 <- NULL
      }
      
      if (!is.null(cox_covariate_1)) {
        n_original    <- nrow(summary_table)
        summary_table <- summary_table[!is.na(summary_table[[cox_covariate_1]]) &
                                         (is.numeric(summary_table[[cox_covariate_1]]) | 
                                          trimws(as.character(summary_table[[cox_covariate_1]])) != ""),]
        
        unique_values <- sort(unique(summary_table[[cox_covariate_1]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(summary_table)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
      
    } else {
      logger.info("... with yearly survival")
      
      non_na_unique <- unique(na.omit(yearly_survival[[cox_covariate_1]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", cox_covariate_1, ", is entirely NA; no comparison is possible."))
        } else {
          logger.fatal(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_1, "; no comparison is possible."))
        }
        cox_covariate_1 <- NULL
      }
      
      if (!is.null(cox_covariate_1)) {
        n_original      <- nrow(yearly_survival)
        yearly_survival <- yearly_survival[!is.na(yearly_survival[[cox_covariate_1]]) & 
                                             (is.numeric(yearly_survival[[cox_covariate_1]]) | 
                                              trimws(as.character(yearly_survival[[cox_covariate_1]])) != ""), ]
        
        unique_values <- sort(unique(yearly_survival[[cox_covariate_1]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(yearly_survival)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
    }
  }
  
  # Comparison covariate 2 --- 
  if(!is.null(cox_covariate_2)) {
    logger.info(paste("Comparing across", cox_covariate_2, "..."))
    
    if(is.null(survival_yr_start)){
      logger.info("... with summary table")
      
      non_na_unique <- unique(na.omit(summary_table[[cox_covariate_2]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", cox_covariate_2, ", is entirely NA; no comparison is possible."))
        } else {
          logger.fatal(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_2, "; no comparison is possible."))
        }
        cox_covariate_2 <- NULL
      }
      
      if (!is.null(cox_covariate_2)) {
        n_original    <- nrow(summary_table)
        summary_table <- summary_table[!is.na(summary_table[[cox_covariate_2]]) & 
                                         (is.numeric(summary_table[[cox_covariate_2]]) | 
                                            trimws(as.character(summary_table[[cox_covariate_2]])) != ""),]
    
        unique_values <- sort(unique(summary_table[[cox_covariate_2]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(summary_table)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
      
    } else {
      logger.info("... with yearly survival")
      
      non_na_unique <- unique(na.omit(yearly_survival[[cox_covariate_2]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", cox_covariate_2, ", is entirely NA; no comparison is possible."))
        } else {
          logger.info(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_2, "; no comparison is possible."))
        }
        cox_covariate_2 <- NULL
      }
      
      if (!is.null(cox_covariate_2)) {
        n_original      <- nrow(yearly_survival)
        yearly_survival <- yearly_survival[!is.na(yearly_survival[[cox_covariate_2]]) & 
                                             (is.numeric(yearly_survival[[cox_covariate_2]]) | 
                                              trimws(as.character(yearly_survival[[cox_covariate_2]])) != ""), ]
        
        unique_values <- sort(unique(yearly_survival[[cox_covariate_2]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(yearly_survival)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
    }
  }
  
  # Comparison covariate 3 ---  
  if(!is.null(cox_covariate_3)) {
    logger.info(paste("Comparing across", cox_covariate_3, "..."))
    
    if(is.null(survival_yr_start)){
      logger.info("... with summary table")
      
      non_na_unique <- unique(na.omit(summary_table[[cox_covariate_3]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", cox_covariate_3, ", is entirely NA; no comparison is possible."))
        } else {
          logger.fatal(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_3, "; no comparison is possible."))
        }
        cox_covariate_3 <- NULL
      }
      
      if (!is.null(cox_covariate_3)) {
        n_original    <- nrow(summary_table)
        summary_table <- summary_table[!is.na(summary_table[[cox_covariate_3]]) & 
                                         (is.numeric(summary_table[[cox_covariate_3]]) | 
                                          trimws(as.character(summary_table[[cox_covariate_3]])) != ""),]
        
        unique_values <- sort(unique(summary_table[[cox_covariate_3]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(summary_table)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
      
    } else {
      logger.info("... with yearly survival")
      
      non_na_unique <- unique(na.omit(yearly_survival[[cox_covariate_3]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", cox_covariate_3, ", is entirely NA; no comparison is possible."))
        } else {
          logger.fatal(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_3, "; no comparison is possible."))
        }
        cox_covariate_3 <- NULL
      }
      
      if (!is.null(cox_covariate_3)) {
        n_original      <- nrow(yearly_survival)
        yearly_survival <- yearly_survival[!is.na(yearly_survival[[cox_covariate_3]]) & 
                                             (is.numeric(yearly_survival[[cox_covariate_3]]) | 
                                              trimws(as.character(yearly_survival[[cox_covariate_3]])) != ""), ]
        
        unique_values <- sort(unique(yearly_survival[[cox_covariate_3]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(yearly_survival)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
    }
  }
  
  
  ## Basic summaries of data ------------------------------------------------
  
  # Plot each individual's tracking history --- 
  
  if(calc_tracking_history == TRUE) {
    
    # Helper (base R, no extra packages)
    clamp <- function(x, lo, hi) max(lo, min(hi, x))
    
    if (is.null(survival_yr_start)) {
      logger.info("Plotting tracking history using summary table...")
      
      # Create deployment summary
      deployment_summary <- summary_table |>
        mutate(deploy_on  = as.POSIXct(deploy_on_timestamp),
               deploy_off = as.POSIXct(deploy_off_timestamp),
               duration_days = round(as.numeric(difftime(deploy_off, deploy_on, units = "days")), 1)) |>
        
        filter(deploy_off > deploy_on,
               !is.na(deploy_on),
               !is.na(deploy_off)) |>
        
        group_by(individual_id, individual_local_identifier) |>          
        mutate(first_start = min(deploy_on, na.rm = TRUE) ) |>
        ungroup() |>
        
        mutate(individual_label = fct_reorder(
          paste(individual_id, individual_local_identifier, sep = " – "),
          first_start)) |>
        
        arrange(first_start, deploy_on) |>
        mutate(plot_start = deploy_on,
               plot_end   = deploy_off)
      
      # Total location count
      n_locs_total <- summary_table |>
        summarise(total = sum(n_locations, na.rm = TRUE)) |>
        pull(total)
      
      # Gap detection
      deployment_summary_with_gaps <- deployment_summary |>
        group_by(individual_label) |>
        arrange(plot_start) |>
        mutate(prev_end  = lag(plot_end),
               gap_start = prev_end,
               gap_end   = plot_start,
               gap_days  = as.numeric(difftime(gap_end, gap_start, units = "days"))) |>
        filter(gap_days > 30, !is.na(gap_days)) |>
        ungroup()
      
      # Build the plot  
      tracking_history <- ggplot(deployment_summary) +
      geom_segment(aes(x = plot_start, xend = plot_end, 
                         y = individual_label, yend = individual_label),
                     linewidth = 3.2, color = "grey") +
      geom_point(aes(x = plot_start, y = individual_label),
                 color = "#1F77B4", size = 3.5) +
      geom_point(aes(x = plot_end, y = individual_label,
                 color = factor(mortality_event)), size = 3.5) +
      scale_color_manual(
        values = c("0" = "#9467BD", "1" = "red"),
        labels = c("0" = "Censored", "1" = "Mortality"),
        name   = "End event") +
  geom_segment(data = deployment_summary_with_gaps,
               aes(x = gap_start + (gap_end - gap_start)/2,
                   xend = gap_start + (gap_end - gap_start)/2,
                   y = as.numeric(individual_label),
                   yend = as.numeric(individual_label) + 0.45),
               color = "grey50", linewidth = 1.2,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  labs(title = "Individual Collared Periods",
       subtitle = sprintf("%d unique individuals • %d visible deployments • %d locations",
                          n_distinct(deployment_summary$individual_id),
                          nrow(deployment_summary),
                          n_locs_total),
       x = "Time",
       y = "Individual") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8, face = "plain"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey50", margin = margin(b = 10)),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12)) +
  scale_x_datetime(date_breaks = "1 year",
                   date_labels = "%Y",
                   expand = expansion(mult = c(0.01, 0.03)))
      
      # Dynamic dimensions based on data
      n_individuals <- n_distinct(deployment_summary$individual_id)
      time_range_years <- as.numeric(difftime(max(deployment_summary$plot_end), 
                                              min(deployment_summary$plot_start), 
                                              units = "days")) / 365.25
      
      # Height: base + per-individual allowance, clamped to reasonable bounds
      plot_height <- clamp(2 + n_individuals * 0.18, 5, 40)
      
      # Width: base + per-year allowance, clamped to reasonable bounds
      plot_width  <- clamp(4 + time_range_years * 1.2, 7, 24)
      
      # Save plot
      png(appArtifactPath("tracking_history.png"),
          width  = plot_width,
          height = plot_height,
          units  = "in", res = 300)
      print(tracking_history)
      dev.off()
      
    }
    
    if (!is.null(survival_yr_start)) {
      logger.info("Plotting tracking history using yearly survival")
      
      mortality_status <- yearly_survival |>
        group_by(individual_id) |>
        summarise(mortality_event = max(mortality_event, na.rm = TRUE), .groups = "drop")
      
      deployment_summary <- yearly_survival |>
        distinct(individual_id,
                 individual_local_identifier,
                 deploy_on_timestamp,
                 deploy_off_timestamp) |>
        left_join(mortality_status, by = "individual_id") |>   # <-- join here
        mutate(deploy_on  = as.POSIXct(deploy_on_timestamp),
               deploy_off = as.POSIXct(deploy_off_timestamp),
               duration_days = round(as.numeric(difftime(deploy_off, deploy_on, units = "days")), 1)) |>
        filter(deploy_off > deploy_on,
               !is.na(deploy_on),
               !is.na(deploy_off)) |>
        group_by(individual_id) |>
        mutate(first_start = min(deploy_on, na.rm = TRUE),
               individual_label = fct_reorder(paste(individual_id, individual_local_identifier, sep = " – "),
                                              first_start)) |>
        ungroup() |>
        arrange(first_start, deploy_on) |>
        mutate(plot_start = deploy_on,
               plot_end   = deploy_off)
      
      # Total location count 
      n_locs_total <- yearly_survival |>
        distinct(individual_id, n_locations) |>
        summarise(total = sum(n_locations, na.rm = TRUE)) |>
        pull(total)
      
      # Plot 
      tracking_history <- ggplot(deployment_summary) +
        geom_segment(aes(x = plot_start, xend = plot_end,
                         y = individual_label, yend = individual_label),
                     linewidth = 3.2, color = "grey") +
        geom_point(aes(x = plot_start, y = individual_label),
                   color = "#1F77B4", size = 3.5) +
        geom_point(aes(x = plot_end, y = individual_label,
                       color = factor(mortality_event)),
                   size = 3.5) +
        scale_color_manual(
          values = c("0" = "#9467BD", "1" = "red"),
          labels = c("0" = "Censored", "1" = "Mortality"),
          name   = "End event") +
        geom_segment(data = deployment_summary |>
                       group_by(individual_label) |>
                       arrange(plot_start) |>
                       mutate(prev_end  = lag(plot_end),
                              gap_start = prev_end,
                              gap_end   = plot_start,
                              gap_days  = as.numeric(difftime(gap_end, gap_start, units = "days"))) |>
                       filter(gap_days > 30, !is.na(gap_days)),
                     aes(x = gap_start + (gap_end - gap_start)/2,
                         xend = gap_start + (gap_end - gap_start)/2,
                         y = as.numeric(individual_label),
                         yend = as.numeric(individual_label) + 0.45),
                     color = "grey50", linewidth = 1.2,
                     arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
        labs(title    = paste0("Individual Collared Periods"),
             subtitle = sprintf("%d unique individuals • %d visible deployments • %d locations",
                                n_distinct(deployment_summary$individual_id),
                                nrow(deployment_summary),
                                n_locs_total),
             x = "Time",
             y = "Individual") +
        theme_minimal(base_size = 12) +
        theme(axis.text.y        = element_text(size = 8, face = "plain"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor   = element_blank(),
              plot.title         = element_text(face = "bold", size = 14),
              plot.subtitle      = element_text(size = 11, color = "grey50", margin = margin(b = 10)),
              axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title         = element_text(size = 12)) +
        scale_x_datetime(date_breaks = "1 year",
                         date_labels = "%Y",
                         expand      = expansion(mult = c(0.01, 0.03)))
      
      # Dynamic dimensions based on data
      n_individuals <- n_distinct(deployment_summary$individual_id)
      time_range_years <- as.numeric(difftime(max(deployment_summary$plot_end), 
                                              min(deployment_summary$plot_start), 
                                              units = "days")) / 365.25
      
      # Height: base + per-individual allowance, clamped to reasonable bounds
      plot_height <- clamp(2 + n_individuals * 0.18, 5, 40)
      
      # Width: base + per-year allowance, clamped to reasonable bounds
      plot_width  <- clamp(4 + time_range_years * 1.2, 7, 24)
      
      # Save plot
      png(appArtifactPath("tracking_history.png"),
          width  = plot_width,
          height = plot_height,
          units  = "in", res = 300)
      print(tracking_history)
      dev.off()
      
    } 
  }
  
  
  # Calculate monthly mortality ---- 
  if (calc_month_mort == TRUE) {
    
    if (is.null(survival_yr_start)) {
      logger.info("Plotting monthly mortality using summary table...")
      
      min_date <- min(summary_table$deploy_on_timestamp, na.rm = TRUE)
      max_date <- max(summary_table$deploy_off_timestamp, na.rm = TRUE)
      full_years <- seq(year(min_date), year(max_date), by = 1)
      
      mortality_data <- summary_table |>
        filter(mortality_event == 1) |>
        mutate(mort_date   = as.Date(mortality_date),
               deploy_off  = as.Date(deploy_off_timestamp),
               death_date  = coalesce(mort_date, deploy_off),
               death_year  = year(death_date),
               death_month = lubridate::month(death_date, label = TRUE, abbr = TRUE),
               death_month = factor(death_month, levels = month.abb, ordered = TRUE)) |>
        dplyr::select(death_year, death_month)
      
      monthly_morts <- mortality_data %>%
        count(death_year, death_month, name = "n_mortalities") %>%
        complete(death_year  = full_years,
                 death_month = factor(month.abb, levels = month.abb, ordered = TRUE),
                 fill = list(n_mortalities = 0)) %>%
        mutate(death_month_num = as.integer(death_month),
               death_month     = fct_relevel(death_month, month.abb))
      
      # Plot 
      monthly_mort_plot <- ggplot(monthly_morts, aes(x = death_month, y = factor(death_year), 
                                                     fill = factor(n_mortalities))) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_brewer(palette  = "PuBu",
                          na.value = "grey92",
                          name     = "Number of\nmortality events",
                          drop     = FALSE) + 
        scale_x_discrete(position = "top") +
        labs(title    = "Monthly Distribution of Confirmed Mortality Events",
             subtitle = paste0(
               "Total events: ", sum(monthly_morts$n_mortalities, na.rm = TRUE),
               " • Time span: ", format(min_date, "%b %Y"),
               " to ", format(max_date, "%b %Y")),
             x = NULL,
             y = "Year") +
        theme_minimal(base_size = 14) +
        theme(panel.grid        = element_blank(),
              axis.ticks        = element_blank(),
              legend.position   = "right",
              legend.title      = element_text(size = 11),
              legend.text       = element_text(size = 10),
              plot.title        = element_text(face = "bold", hjust = 0, size = 16),
              plot.subtitle     = element_text(hjust = 0, size = 12),
              axis.text.x       = element_text(size = 11, face = "bold"),
              axis.text.y       = element_text(size = 11))
      
      # Save plot  
      png(appArtifactPath("monthly_mortality.png"), 
          width = 10, height = 8, 
          units = "in", res = 300)
      print(monthly_mort_plot)
      dev.off()
    }
    
    if (!is.null(survival_yr_start)) {
      
      logger.info("Plotting monthly mortality using yearly survival")
      
      min_date <- min(yearly_survival$deploy_on_timestamp, na.rm = TRUE)
      max_date <- max(yearly_survival$deploy_off_timestamp, na.rm = TRUE)
      full_years <- seq(year(min_date), year(max_date), by = 1)
      
      mortality_data <- yearly_survival |>
        filter(mortality_event == 1) |>
        mutate(mort_date   = as.Date(mortality_date),
               deploy_off  = as.Date(deploy_off_timestamp),
               death_date  = coalesce(mort_date, deploy_off),
               death_year  = year(death_date),
               death_month = lubridate::month(death_date, label = TRUE, abbr = TRUE),
               death_month = factor(death_month, levels = month.abb, ordered = TRUE)) |>
        dplyr::select(death_year, death_month)
      
      monthly_morts <- mortality_data |>
        count(death_year, death_month, name = "n_mortalities") |>
        complete( death_year  = full_years,                                       
                  death_month = factor(month.abb, levels = month.abb, ordered = TRUE),
                  fill = list(n_mortalities = 0)) |>
        mutate(death_month_num = as.integer(death_month),
               death_month     = fct_relevel(death_month, month.abb))
      
      # Plot 
      monthly_mort_plot <- ggplot(monthly_morts, 
                                  aes(x = death_month, 
                                      y = factor(death_year),
                                      fill = factor(n_mortalities))) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_brewer(palette  = "PuBu",
                          na.value = "grey92",
                          name     = "Number of\nmortality events",
                          drop     = FALSE) + 
        scale_x_discrete(position = "top") +
        labs(title = "Monthly Distribution of Confirmed Mortality Events",
             subtitle = paste0(
               "Total events: ", sum(monthly_morts$n_mortalities, na.rm = TRUE),
               " • Time span: ", format(min_date, "%b %Y"),
               " to ", format(max_date, "%b %Y")),
             x = NULL,
             y = "Year") +
        theme_minimal(base_size = 14) +
        theme(panel.grid = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "right",
              legend.title = element_text(size = 11),
              legend.text = element_text(size = 10),
              plot.title = element_text(face = "bold", hjust = 0, size = 16),
              plot.subtitle = element_text(hjust = 0, size = 12),
              axis.text.x = element_text(size = 11, face = "bold"),
              axis.text.y = element_text(size = 11))
      
      # Save plot  
      png(appArtifactPath("monthly_mortality.png"), 
          width = 10, height = 8, 
          units = "in", res = 300)
      print(monthly_mort_plot)
      dev.off()
    }
  }
  
  
  ## Cox Proportional Hazard Analysis -----------------------------------------
  
  force_categorical_cols <- c("survival_year")
  
  # Select fitting data 
  if (is.null(survival_yr_start)) {
    logger.info("Calculating Cox Proportional Hazards using summary table...")
    fitting_data <- summary_table
    
  } else {
    logger.info("Calculating Cox Proportional Hazards using survival table...")
    fitting_data <- yearly_survival
  }
  
  # Warning for no mortality events
  if (sum(fitting_data$mortality_event) == 0) {
    logger.fatal("There are no mortality events in the chosen subset of data.")
  }
  
  # Collect non-NULL covariates
  covariates <- c(cox_covariate_1, cox_covariate_2, cox_covariate_3)
  ref_levels <- list(cox_covariate_1_ref, cox_covariate_2_ref, cox_covariate_3_ref)
  names(ref_levels) <- c("cox_covariate_1", "cox_covariate_2", "cox_covariate_3")
  valid_idx  <- !sapply(covariates, is.null) & !is.na(covariates)
  covariates <- covariates[valid_idx]
  ref_levels <- ref_levels[valid_idx]
  names(ref_levels) <- covariates

  if (length(covariates) == 0) {
    logger.fatal("No valid covariates remain for Cox model (all were NULL or had only one level). Skipping Cox PH analysis.")
    
  } else {
    
    # Apply reference levels where specified
    for (cov in covariates) {
      ref <- ref_levels[[cov]]
      
      # Coerce forced-categorical columns to factor before releveling
      if (cov %in% force_categorical_cols) {
        fitting_data[[cov]] <- factor(fitting_data[[cov]])
      }
      
      if (!is.null(ref) && ref %in% levels(factor(fitting_data[[cov]]))) {
        fitting_data[[cov]] <- relevel(factor(fitting_data[[cov]]), ref = ref)
      } else if (!is.null(ref)) {
        logger.info(paste0("Reference level '", ref, "' not found in covariate '", cov, "' — using default."))
      }
    }
    
  # Build formula dynamically
    cox_formula <- as.formula(paste("Surv(entry_time_days, exit_time_days, mortality_event) ~",
                                    paste(covariates, collapse = " + ")))
    
  # Fit standard Cox 
  firth_used          <- FALSE
  separation_detected <- FALSE
    
  coxph_fit <- withCallingHandlers(
    coxph(cox_formula, data = fitting_data),
    warning = function(w) {
      if (grepl("Loglik converged before variable", conditionMessage(w))) {
        separation_detected <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
    
  # Fall back to Firth's penalized Cox if separation detected
  if (separation_detected) {
    logger.warn("Separation detected in coxph — falling back to Firth's penalized Cox model.")
    coxphf_fit <- coxphf(cox_formula, data = fitting_data, maxstep = 0.1, maxit = 200)
    firth_used  <- TRUE
    
    coef_names <- names(coxphf_fit$coefficients)
    cox.tab <- data.frame(
      term      = coef_names,
      estimate  = exp(coxphf_fit$coefficients[coef_names]),  
      conf.low  = coxphf_fit$ci.lower[coef_names], # already exponentiated           
      conf.high = coxphf_fit$ci.upper[coef_names],# already exponentiated           
      p.value   = coxphf_fit$prob[coef_names],
      row.names = NULL
    )
  } else {
    cox.tab <- tidy(coxph_fit, exponentiate = TRUE, conf.int = TRUE)
  }
  
  # Flag which model used
  cox.tab$model <- ifelse(firth_used, "Firth's Penalized Cox", "Standard Cox")
    
  # Save
  write.csv(cox.tab, file = appArtifactPath("coxph_table.csv"), row.names = FALSE)
    
  
  ## Results at covariate means -----------------------------------------------
  
  if (isTRUE(calc_artifacts_at_mean)) { 
    
    # Predict survival function for subject with means on all covariates --- 
    if (firth_used) {
      coxph_for_survfit <- coxph(cox_formula, data = fitting_data,
                                 init    = coxphf_fit$coefficients[coef_names],
                                 control = coxph.control(iter.max = 0))
      surv.at.means <- survfit(coxph_for_survfit, data = fitting_data)
    } else {
      surv.at.means <- survfit(coxph_fit, data = fitting_data)
    }
    
    # Save
    surv.at.means.tab <- tidy(surv.at.means)
    write.csv(surv.at.means.tab, file = appArtifactPath("surv_at_means.csv"), row.names = FALSE)
    
    
    # Plot of predicted survival --- 
    optional_layers <- list()
    if (isTRUE(add_cis)) {
      optional_layers <- c(optional_layers, list(add_confidence_interval()))
    }
    
    surv_plot <- ggsurvfit(surv.at.means) +
      optional_layers +
      add_risktable(risktable_stats = c("n.risk", "cum.event"),
                    theme = theme_risktable_default(axis.text.y.size = 9)) +
      labs(title    = "Predicted Survival at Covariate Means",
           subtitle = ifelse(firth_used, "Firth's Penalized Cox Model", "Standard Cox Model"),
           x        = "Days",
           y        = "Survival Probability") +
      scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
      scale_x_continuous(expand = c(0.02, 0),
                         breaks = pretty(range(surv.at.means.tab$time), n = 8)) +
      theme_classic(base_size = 12) +
      theme(plot.title         = element_text(face = "bold", size = 14),
            plot.subtitle      = element_text(size = 12, color = "gray50"),
            axis.text          = element_text(color = "black"),
            axis.title         = element_text(size = 12),
            panel.grid.major.y = element_line(color = "gray90"),
            panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
            legend.position    = "none",
            plot.margin        = margin(10, 10, 10, 10))
    
    # Save 
    png(appArtifactPath("surv_at_means_plot.png"), 
        width = 7, height = 6, 
        units = "in", res = 300)
    print(surv_plot)
    dev.off()
    
    ## Plot cumulative hazard curve --- 
    cum_hazard <- ggsurvfit(surv.at.means, type = "cumhaz") +
      optional_layers +
      labs(title    = "Cumulative Hazard at Covariate Means",
           subtitle = ifelse(firth_used, "Firth's Penalized Cox Model", "Standard Cox Model"),
           x        = "Days",
           y        = "Cumulative Hazard") +
      scale_x_continuous(expand = c(0.02, 0),
                         breaks = pretty(range(surv.at.means.tab$time), n = 8)) +
      theme_classic(base_size = 12) +
      theme(plot.title         = element_text(face = "bold", size = 14),
            plot.subtitle      = element_text(size = 12, color = "gray50"),
            axis.text          = element_text(color = "black"),
            axis.title         = element_text(size = 12),
            panel.grid.major.y = element_line(color = "gray90"),
            panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
            legend.position    = "none",
            plot.margin        = margin(10, 10, 10, 10))
    
    # Save   
    png(appArtifactPath("cum_hazard_at_means_plot.png"), 
        width = 7, height = 6, 
        units = "in", res = 300)
    print(cum_hazard)
    dev.off()
  } 
    
  
  ## Diagnostic residuals -----------------------------------------------------
  
  if (isTRUE(calc_residuals)){
    
    ## Plot Schoenfeld residuals ---
    if (!firth_used) {
      ph_test <- cox.zph(coxph_fit)
      ph_plot <- ggcoxzph(ph_test, point.col = "steelblue", point.size = 1.5, point.alpha = 0.5) 
      
      # Save 
      png(appArtifactPath("schoenfeld_residuals.png"), 
          width  = 8, 
          height = 3 * length(covariates), 
          units  = "in",
          res    = 300)
      print(ph_plot)
      dev.off()
      
    } else {
      # Firth's Cox: coxph_for_survfit is a coxph object frozen at Firth's coefficients
      logger.info("PH test is approximate — coefficients are from Firth's model but SE/residuals from coxph")
      ph_test <- cox.zph(coxph_for_survfit)
      print(ph_test)
      
      ph_plots <- suppressWarnings(ggcoxzph(ph_test, point.col = "steelblue", 
                                            point.size  = 1.5, point.alpha = 0.5))
      
      # Add caption to each plot in the list
      ph_plots <- lapply(ph_plots, function(p) {
        p + 
          labs(caption = "Note: Schoenfeld residuals are approximate — model fit via Firth's penalized Cox.") +
          theme(plot.caption = element_text(color = "grey40", size = 8))
      })
      ph_combined <- wrap_plots(ph_plots, ncol = 1)
      
      # Save
      png(appArtifactPath("schoenfeld_residuals.png"), 
          width  = 8, 
          height = 3 * length(covariates), 
          units  = "in",
          res    = 300)
      print(ph_combined)
      dev.off()
    }
  }
  
  } 
  
  ## Stratified Group Comparisons ---------------------------------------------
  
  ## Prepare covariates --- 
  active_covariates <- list(
    list(var = cox_covariate_1, ref = cox_covariate_1_ref),
    list(var = cox_covariate_2, ref = cox_covariate_2_ref),
    list(var = cox_covariate_3, ref = cox_covariate_3_ref)
  )
  
  active_covariates <- Filter(function(x) !is.null(x$var), active_covariates)
  fitting_data <- if (is.null(survival_yr_start)) summary_table else yearly_survival
  all_forest_tabs <- list()
  
  
  for (cov_item in active_covariates) {
    
    # Reset fitting_data at the start of each iteration
    fitting_data <- if (is.null(survival_yr_start)) summary_table else yearly_survival
    
    cov <- cov_item$var
    ref <- cov_item$ref
    logger.info(paste("Producing group comparison outputs for:", cov))
    
    # Detect whether covariate is continuous
    is_continuous <- (is.numeric(fitting_data[[cov]]) || inherits(fitting_data[[cov]], "units")) &&
      !(cov %in% force_categorical_cols)
    
    # Coerce forced-categorical columns to factor
    if (cov %in% force_categorical_cols) {
      fitting_data[[cov]] <- factor(fitting_data[[cov]])
      if (!is.null(ref) && ref %in% levels(fitting_data[[cov]])) {
        fitting_data[[cov]] <- relevel(fitting_data[[cov]], ref = ref)
      }
    }
    
    # Build formula  
    km_formula <- as.formula(paste("Surv(entry_time_days, exit_time_days, mortality_event) ~", cov))
    
    
    ## For continuous covariates ---
    if (is_continuous) {
      logger.info(paste(cov, "is continuous — running Cox model only."))
      
      # Strip units if present 
      if (inherits(fitting_data[[cov]], "units")) {
        fitting_data[[cov]] <- as.numeric(fitting_data[[cov]])
      }
      
      ## Cox model tests (LR, Wald, Score) ---
      separation_detected <- FALSE
      
      cox_strat_fit <- withCallingHandlers(
        coxph(km_formula, data = fitting_data),
        warning = function(w) {
          if (grepl("Loglik converged before variable", conditionMessage(w))) {
            separation_detected <<- TRUE
            invokeRestart("muffleWarning")
          }
        }
      )
      
      if (separation_detected) {
        logger.warn(sprintf(
          "Separation detected for %s — falling back to Firth's penalized Cox for HR estimates.", cov))
        cox_strat_firth <- coxphf(km_formula, data = fitting_data, maxstep = 0.1, maxit = 200)
      }
      
      cox_summary <- summary(cox_strat_fit)
      
      test_results <- data.frame(
        covariate = cov,
        test      = c("Likelihood Ratio", "Wald", "Score"),
        chi_sq    = round(c(cox_summary$logtest["test"],
                            cox_summary$waldtest["test"],
                            cox_summary$sctest["test"]), 3),
        df        = c(cox_summary$logtest["df"],
                      cox_summary$waldtest["df"],
                      cox_summary$sctest["df"]),
        p_value   = round(c(cox_summary$logtest["pvalue"],
                            cox_summary$waldtest["pvalue"],
                            cox_summary$sctest["pvalue"]), 4),
        note      = ifelse(separation_detected,
                           "Separation detected — HR estimates from Firth's model",
                           "Standard Cox"))
      
      write.csv(test_results,
                file      = appArtifactPath(paste0("model_tests_", cov, ".csv")),
                row.names = FALSE)
      
      ## HR table ---
      if (separation_detected) {
        cox_strat_tab <- data.frame(
          term      = names(cox_strat_firth$coefficients),
          estimate  = exp(cox_strat_firth$coefficients),
          conf.low  = cox_strat_firth$ci.lower,
          conf.high = cox_strat_firth$ci.upper,
          p.value   = cox_strat_firth$prob,
          covariate = cov,
          row.names = NULL
        )
      } else {
        cox_strat_tab <- tidy(cox_strat_fit, exponentiate = TRUE, conf.int = TRUE) %>%
          mutate(covariate = cov)
      }
      
      write.csv(cox_strat_tab,
                file      = appArtifactPath(paste0("cox_hr_", cov, ".csv")),
                row.names = FALSE)
      
      # Collect for combined forest plot
      all_forest_tabs[[cov]] <- cox_strat_tab %>%
        mutate(covariate = cov,
               lr_chisq  = test_results$chi_sq[1],
               lr_p      = test_results$p_value[1])
      
      next
    }
    
    
    ## For categorical covariates ---
    
    # Apply reference level if specified
    if (!is.null(ref) && ref %in% unique(fitting_data[[cov]])) {
      fitting_data[[cov]] <- relevel(factor(fitting_data[[cov]]), ref = ref)
    } else {
      fitting_data[[cov]] <- factor(fitting_data[[cov]])
    }
    
    n_groups <- length(levels(fitting_data[[cov]]))
    pal      <- viridis(n_groups, option = "turbo", end = 0.9)
    
    
    ## Cox model tests (LR, Wald, Score) ---
    separation_detected <- FALSE
    
    cox_strat_fit <- withCallingHandlers(
      coxph(km_formula, data = fitting_data),
      warning = function(w) {
        if (grepl("Loglik converged before variable", conditionMessage(w))) {
          separation_detected <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    )
    
    if (separation_detected) {
      logger.warn(sprintf(
        "Separation detected for %s — falling back to Firth's penalized Cox for HR estimates.", cov))
      cox_strat_firth <- coxphf(km_formula, data = fitting_data, maxstep = 0.1, maxit = 200)
    }
    
    cox_summary <- summary(cox_strat_fit)
    
    test_results <- data.frame(
      covariate = cov,
      test      = c("Likelihood Ratio", "Wald", "Score"),
      chi_sq    = round(c(cox_summary$logtest["test"],
                          cox_summary$waldtest["test"],
                          cox_summary$sctest["test"]), 3),
      df        = c(cox_summary$logtest["df"],
                    cox_summary$waldtest["df"],
                    cox_summary$sctest["df"]),
      p_value   = round(c(cox_summary$logtest["pvalue"],
                          cox_summary$waldtest["pvalue"],
                          cox_summary$sctest["pvalue"]), 4),
      note      = ifelse(separation_detected,
                         "Separation detected — HR estimates from Firth's model",
                         "Standard Cox"))
    
    lr_p <- test_results$p_value[1]
    
    write.csv(test_results,
              file      = appArtifactPath(paste0("model_tests_", cov, ".csv")),
              row.names = FALSE)
    
    
    ## Pairwise comparisons (if >2 groups) ---
    if (n_groups > 2) {
      group_levels <- levels(fitting_data[[cov]])
      pairs        <- combn(group_levels, 2, simplify = FALSE)
      
      pairwise_results <- lapply(pairs, function(pair) {
        sub_data        <- fitting_data[fitting_data[[cov]] %in% pair, ]
        sub_data[[cov]] <- droplevels(sub_data[[cov]])
        
        pair_separation <- FALSE
        
        fit <- withCallingHandlers(
          coxph(km_formula, data = sub_data),
          warning = function(w) {
            if (grepl("Loglik converged before variable|did not converge", conditionMessage(w))) {
              pair_separation <<- TRUE
              invokeRestart("muffleWarning")
            }
          }
        )
        
        if (pair_separation) {
          firth_fit <- tryCatch(
            coxphf(km_formula, data = sub_data, maxstep = 0.1, maxit = 200),
            error = function(e) NULL
          )
          
          if (!is.null(firth_fit)) {
            return(data.frame(
              group_1 = pair[1],
              group_2 = pair[2],
              hr      = round(exp(firth_fit$coefficients[1]), 3),
              p_wald  = round(firth_fit$prob[1], 4),
              method  = "Firth"))
          } else {
            return(data.frame(
              group_1 = pair[1],
              group_2 = pair[2],
              hr      = NA_real_,
              p_wald  = NA_real_,
              method  = "Failed"))
          }
        }
        
        s <- summary(fit)
        data.frame(
          group_1 = pair[1],
          group_2 = pair[2],
          hr      = round(exp(coef(fit)[1]), 3),
          p_wald  = round(s$waldtest["pvalue"], 4),
          method  = "Standard Cox")
      })
      
      pairwise_df              <- do.call(rbind, pairwise_results)
      pairwise_df$p_bonferroni <- round(p.adjust(pairwise_df$p_wald, method = "bonferroni"), 4)
      pairwise_df$significant  <- ifelse(pairwise_df$p_bonferroni < 0.05, "Yes", "No")
      
      write.csv(pairwise_df,
                file      = appArtifactPath(paste0("pairwise_", cov, ".csv")),
                row.names = FALSE)
      
      logger.info(sprintf("Pairwise comparisons for %s saved (%d pairs, Bonferroni corrected).",
                          cov, nrow(pairwise_df)))
    }
    
    
    ## HR table ---
    if (separation_detected) {
      cox_strat_tab <- data.frame(
        term      = names(cox_strat_firth$coefficients),
        estimate  = exp(cox_strat_firth$coefficients),
        conf.low  = cox_strat_firth$ci.lower,
        conf.high = cox_strat_firth$ci.upper,
        p.value   = cox_strat_firth$prob,
        covariate = cov,
        row.names = NULL
      )
    } else {
      cox_strat_tab <- tidy(cox_strat_fit, exponentiate = TRUE, conf.int = TRUE) %>%
        mutate(covariate = cov)
    }
    
    write.csv(cox_strat_tab,
              file      = appArtifactPath(paste0("cox_hr_", cov, ".csv")),
              row.names = FALSE)
    
    
    ## KM survival curves ---
    km_fit <- survfit(km_formula, data = fitting_data)
    
    km_plot <- ggsurvfit(km_fit, linewidth = 0.9) +
      { if (add_cis) add_confidence_interval(alpha = 0.12) } +
      scale_color_manual(values = pal) +
      scale_fill_manual(values  = pal) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
      scale_x_continuous(expand = c(0.02, 0),
                         breaks = pretty(range(tidy(km_fit)$time), n = 8)) +
      labs(title    = paste("Kaplan-Meier Survival by", cov),
           subtitle = sprintf("Cox LR test: chi-sq=%.3f, p=%.4f", test_results$chi_sq[1], lr_p),
           x        = "Days",
           y        = "Survival Probability",
           color    = cov,
           fill     = cov) +
      theme_classic(base_size = 12) +
      theme(plot.title         = element_text(face = "bold", size = 14),
            plot.subtitle      = element_text(size = 10, color = "grey40"),
            legend.position    = "bottom",
            panel.grid.major.y = element_line(color = "grey92"))
    
    # Zoom y-axis to data range
    if (isTRUE(zoom_to_plot)) {
      y_min   <- min(km_fit$surv, na.rm = TRUE)
      y_floor <- floor(y_min * 10) / 10
      km_plot <- km_plot +
        coord_cartesian(ylim = c(y_floor, 1)) +
        annotate("text",
                 x        = -Inf,
                 y        = y_floor,
                 label    = paste0("* y-axis zoomed to [", round(y_floor, 2), ", 1]"),
                 hjust    = -0.05,
                 vjust    = -0.5,
                 size     = 3.2,
                 color    = "gray50",
                 fontface = "italic")
    }
    
    png(appArtifactPath(paste0("km_by_", cov, ".png")),
        width = 8, height = 7, units = "in", res = 300)
    print(km_plot)
    dev.off()
    
    
    ## Per-group median survival summary table ---
    km_summary <- tidy(km_fit) %>%
      group_by(strata) %>%
      summarise(n_risk_start = first(n.risk),
                n_events     = sum(n.event, na.rm = TRUE),
                median_surv  = {
                  cross <- time[estimate <= 0.5]
                  if (length(cross) > 0) first(cross) else NA_real_
                },
                surv_at_1yr  = {
                  t365 <- estimate[time <= 365]
                  if (length(t365) > 0) round(last(t365), 3) else NA_real_
                },
                .groups = "drop") %>%
      mutate(lr_p_value = lr_p, covariate = cov)
    
    write.csv(km_summary,
              file      = appArtifactPath(paste0("km_summary_", cov, ".csv")),
              row.names = FALSE)
    
    
    ## Collect for combined forest plot --- 
    all_forest_tabs[[cov]] <- cox_strat_tab %>%
      mutate(covariate = cov,
             lr_chisq  = test_results$chi_sq[1],
             lr_p      = lr_p)
    
    
    ## Cumulative hazard by group ---
    cumhaz_plot <- ggsurvfit(km_fit, type = "cumhaz", linewidth = 0.9) +
      { if (add_cis) add_confidence_interval(alpha = 0.12) } +
      scale_color_manual(values = pal) +
      scale_fill_manual(values  = pal) +
      scale_x_continuous(expand = c(0.02, 0),
                         breaks = pretty(range(tidy(km_fit)$time), n = 8)) +
      labs(title    = paste("Cumulative Hazard by", cov),
           subtitle = "Parallel lines support proportional hazards assumption",
           x        = "Days",
           y        = "Cumulative Hazard",
           color    = cov,
           fill     = cov) +
      theme_classic(base_size = 12) +
      theme(plot.title      = element_text(face = "bold"),
            plot.subtitle   = element_text(size = 10, color = "grey40"),
            legend.position = "bottom")
    
    # Zoom y-axis to data range
    if (isTRUE(zoom_to_plot)) {
      y_max     <- max(-log(km_fit$surv), na.rm = TRUE)
      y_ceiling <- ceiling(y_max * 10) / 10
      cumhaz_plot <- cumhaz_plot +
        coord_cartesian(ylim = c(0, y_ceiling)) +
        annotate("text",
                 x        = -Inf,
                 y        = y_ceiling,
                 label    = paste0("* y-axis zoomed to [0, ", round(y_ceiling, 2), "]"),
                 hjust    = -0.05,
                 vjust    = 1.5,
                 size     = 3.2,
                 color    = "gray50",
                 fontface = "italic")
    }
    
    png(appArtifactPath(paste0("cumhaz_by_", cov, ".png")),
        width = 7, height = 5, units = "in", res = 300)
    print(cumhaz_plot)
    dev.off()
    
    
    ## Annual survival rate bar chart (yearly mode only) ---
    if (!is.null(survival_yr_start)) {
      
      group_vars <- unique(c("survival_year", cov))
      
      annual_surv <- fitting_data %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(n_animals = n_distinct(individual_id),
                  n_deaths  = sum(mortality_event, na.rm = TRUE),
                  surv_rate = 1 - (n_deaths / n_animals),
                  .groups   = "drop")
      
      # Dynamic sizing based on number of groups and years
      n_groups <- length(unique(annual_surv[[cov]]))
      n_years  <- length(unique(annual_surv$survival_year))
      plot_width  <- max(8, n_years * n_groups * 0.4)
      plot_height <- max(6, 5 + ceiling(n_groups / 4) * 0.6)  # extra height for legend rows
      
      annual_plot <- ggplot(annual_surv, aes(x = factor(survival_year),
                                             y = surv_rate,
                                             fill = .data[[cov]])) +
        geom_col(position = position_dodge(0.85), width = 0.75, alpha = 0.88,
                 color = "white", linewidth = 0.3) +
        geom_text(aes(label = paste0("n=", n_animals),
                      y = pmax(surv_rate - 0.04, 0.04)),   
                  position = position_dodge(0.85),
                  vjust    = 1,
                  hjust    = 0.5,
                  size     = 3.0,
                  color    = "white",
                  fontface = "bold") +
        scale_fill_manual(values = pal) +
        scale_y_continuous(limits = c(0, 1.05),
                           breaks = seq(0, 1, by = 0.25),
                           labels = percent_format(accuracy = 1),
                           expand = c(0, 0)) +
        labs(title    = paste("Annual Survival Rate by", cov),
             subtitle = "Bars show naive annual survival; labels show number of tracked individuals",
             x        = "Survival Year",
             y        = "Survival Rate",
             fill     = cov) +
        theme_classic(base_size = 14) +
        theme(plot.title         = element_text(face = "bold", size = 16),
              plot.subtitle      = element_text(size = 12, color = "gray50",
                                                margin = margin(b = 8)),
              axis.text          = element_text(color = "black", size = 12),
              axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
              axis.title         = element_text(size = 13, face = "bold"),
              panel.grid.major.y = element_line(color = "gray90"),
              panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
              legend.position    = "bottom",
              legend.title       = element_text(face = "bold", size = 12),
              legend.text        = element_text(size = 11),
              legend.key.size    = unit(0.6, "cm"),
              legend.box.spacing = unit(0.4, "cm"),
              plot.margin        = margin(10, 15, 10, 10)) +
        guides(fill = guide_legend(nrow           = ceiling(n_groups / 4),
                                   byrow          = TRUE,
                                   title.position = "top"))
      
      png(appArtifactPath(paste0("annual_surv_by_", cov, ".png")),
          width  = plot_width,
          height = plot_height,
          units  = "in", res = 300)
      print(annual_plot)
      dev.off()
    }
  }
  
  ## Combined forest plot across all covariates --------------------------------
  
  if (length(all_forest_tabs) > 0) {
    
    escaped_covs <- gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1",
                         names(all_forest_tabs))
    
    combined_forest_tab <- bind_rows(all_forest_tabs) %>%
      mutate(
        term_clean  = str_remove(term, regex(paste(escaped_covs, collapse = "|"))),
        term_clean  = if_else(term_clean == "", term, term_clean),
        significant = p.value < 0.05,
        facet_label = sprintf("%s\n(LR chi-sq=%.2f, p=%.4f)", covariate, lr_chisq, lr_p),
        term_clean  = fct_inorder(term_clean)
      )
    
    # Build one plot per covariate
    cov_names   <- names(all_forest_tabs)
    panel_plots <- vector("list", length(cov_names))
    panel_nrows <- integer(length(cov_names))
    
    for (i in seq_along(cov_names)) {
      cov_data <- combined_forest_tab %>% filter(covariate == cov_names[i])
      n_terms  <- nrow(cov_data)
      panel_nrows[i] <- n_terms
      
      # Only show legend on the last panel
      show_legend <- i == length(cov_names)
      
      panel_plots[[i]] <- ggplot(cov_data,
                                 aes(x = estimate, y = term_clean, color = significant)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                       height = 0.25, linewidth = 0.8, color = "grey30") +
        geom_point(size = 3.5) +
        scale_color_manual(values = c("TRUE"  = "#D62728",
                                      "FALSE" = "#1F77B4"),
                           labels = c("TRUE"  = "p < 0.05",
                                      "FALSE" = "p \u2265 0.05"),
                           name   = NULL) +
        scale_x_log10(labels = number_format(accuracy = 0.01)) +
        labs(title = unique(cov_data$facet_label),
             x     = if (i == length(cov_names)) "Hazard Ratio (log scale)" else NULL,
             y     = NULL) +
        theme_classic(base_size = 12) +
        theme(plot.title         = element_text(face = "bold", size = 11, hjust = 0),
              axis.text          = element_text(color = "black"),
              axis.title.x       = element_text(size = 12),
              axis.text.x        = if (i == length(cov_names))
                element_text() else element_blank(),
              axis.ticks.x       = if (i == length(cov_names))
                element_line() else element_blank(),
              panel.grid.major.y = element_line(color = "gray90"),
              panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
              strip.background   = element_rect(fill = "grey95", color = NA),
              legend.position    = if (show_legend) "bottom" else "none",
              legend.justification = "left",
              plot.margin        = margin(5, 10, if (i == length(cov_names)) 5 else 2, 10))
    }
    
    # Stack with patchwork, heights proportional to number of terms
    combined_forest_plot <- wrap_plots(panel_plots, ncol = 1,
                                       heights = panel_nrows) +
      plot_annotation(
        title    = "Hazard Ratios \u2014 All Covariates",
        subtitle = ifelse(firth_used,
                          "Firth's Penalized Cox Model",
                          "Standard Cox Model"),
        theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                         plot.subtitle = element_text(size = 12, color = "gray50"))
      )
    
    # Dynamic height: base per panel + per row allowance
    plot_height <- max(4, length(cov_names) * 1.2 + sum(panel_nrows) * 0.45)
    
    png(appArtifactPath("forest_combined.png"),
        width  = 8,
        height = plot_height,
        units  = "in", res = 300)
    print(combined_forest_plot)
    dev.off()
  }
  
  # Pass original to the next app in the MoveApps workflow
  return(data)
}  
