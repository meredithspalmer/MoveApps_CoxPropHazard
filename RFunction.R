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
                     calc_month_mort) {
  
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
  
  # Extract relevant track-level attributes
  desired_cols <- c("animal_birth_hatch_year", "attachment_type",
                    "capture_method", "capture_timestamp", "death_comments",
                    "deploy_off_timestamp", "deploy_on_timestamp", "deployment_comments",
                    "deployment_end_comments", "deployment_end_type", "deployment_id", 
                    "individual_id", "individual_local_identifier",
                    "individual_number_of_deployments", "is_test", "mortality_location",
                    "model", "mortality_type", "mortality_date", "sex", "tag_id", 
                    "timestamp_first_deployed_location", "timestamp_last_deployed_location")
  
  tracks <- mt_track_data(data) |>
    mutate(mortality_location_filled = if_else(
      is.na(mortality_location) | st_is_empty(mortality_location),
      0L, 1L)) |> 
    dplyr::select(any_of(desired_cols))
  
  # Join track attributes to every event row
  use_deployment_join <- all(c("deployment_id") %in% names(events), 
                             "deployment_id" %in% names(tracks)) &&
    any(!is.na(events$deployment_id)) &&
    any(!is.na(tracks$deployment_id))
  
  if (use_deployment_join) {
    events_with_ind <- events |>
      left_join(tracks, by = "deployment_id")
    
  } else {
    if (!"individual_local_identifier" %in% names(events)) {
      logger.fatal("Cannot join: neither deployment_id nor individual_local_identifier is available in events")
    }
    logger.info("Joining on individual_local_identifier (deployment_id join not possible)")
    events_with_ind <- events |> left_join(tracks, by = "individual_local_identifier")
  }
  
  events_with_ind <- events_with_ind |>
    relocate(any_of(c("individual_id", "individual_local_identifier", "deployment_id", "timestamp")),
             .before = everything())
  
  # Summarize timestamps and location count per individual
  summary_table <- events_with_ind |>
    group_by(individual_id, individual_local_identifier) |>
    summarise(first_timestamp = min(as.Date(timestamp), na.rm = TRUE),
              last_timestamp  = max(as.Date(timestamp), na.rm = TRUE),
              n_locations     = n(),
              n_deployments   = 
                if ("deployment_id" %in% names(events_with_ind)) {
                  n_distinct(deployment_id, na.rm = TRUE)
                } else {
                  1L  
                },
              
              # Time-stamp columns: min / max if present
              timestamp_first_deployed_location = 
                if ("timestamp_first_deployed_location" %in% names(events_with_ind))
                  min(timestamp_first_deployed_location, na.rm = TRUE) else NA,
              
              timestamp_last_deployed_location = 
                if ("timestamp_last_deployed_location" %in% names(events_with_ind))
                  max(timestamp_last_deployed_location, na.rm = TRUE) else NA,
              
              deploy_on_timestamp = 
                if ("deploy_on_timestamp" %in% names(events_with_ind)) {
                  if (all(is.na(deploy_on_timestamp))) as.POSIXct(NA) 
                  else min(deploy_on_timestamp, na.rm = TRUE)
                } else as.POSIXct(NA),
              
              deploy_off_timestamp = if ("deploy_off_timestamp" %in% names(events_with_ind)) {
                if (all(is.na(deploy_off_timestamp))) as.POSIXct(NA)
                else max(deploy_off_timestamp, na.rm = TRUE)
              } else as.POSIXct(NA),
              
              # Mortality location column: 1/0 if filled if present 
              mortality_location_filled = if ("mortality_location_filled" %in% names(events_with_ind))
                as.integer(any(mortality_location_filled >= 1, na.rm = TRUE)) else NA_integer_,
              
              # Categorical columns: collapsed unique if present
              sex = if ("sex" %in% names(events_with_ind))
                str_c(unique(sex[!is.na(sex)]), collapse = " | ") else NA_character_,
              
              mortality_type = if ("mortality_type" %in% names(events_with_ind)) {
                str_c(unique(mortality_type[!is.na(mortality_type)]), collapse = " | ")
              } else NA_character_,
              
              mortality_date = if ("mortality_date" %in% names(events_with_ind)) {
                str_c(unique(mortality_date[!is.na(mortality_date)]), collapse = " | ")
              } else NA_character_,
              
              death_comments = if ("death_comments" %in% names(events_with_ind))
                str_c(unique(death_comments[!is.na(death_comments)]), collapse = " | ") 
              else NA_character_,
              
              deployment_end_comments = if ("deployment_end_comments" %in% names(events_with_ind))
                str_c(unique(deployment_end_comments[!is.na(deployment_end_comments)]), collapse = " | ") 
              else NA_character_,
              
              deployment_end_type = if ("deployment_end_type" %in% names(events_with_ind))
                str_c(unique(deployment_end_type[!is.na(deployment_end_type)]), collapse=" | ") 
              else NA_character_,
              
              animal_life_stage = if ("animal_life_stage" %in% names(events_with_ind))
                str_c(unique(animal_life_stage[!is.na(animal_life_stage)]), collapse = " | ")
              else NA_character_,
              
              model = if ("model" %in% names(events_with_ind)) {
                str_c(unique(model[!is.na(model)]), collapse = " | ")
              } else NA_character_,
              
              attachment_type = if ("attachment_type" %in% names(events_with_ind))
                str_c(unique(attachment_type[!is.na(attachment_type)]), collapse = " | ") 
              else NA_character_,
              
              .groups = "drop") |>
    
    mutate(
      
      # Clean empty strings (fill NA) for columns that exist
      across(any_of(c("death_comments",
                      "mortality_type",
                      "mortality_date",
                      "deployment_end_comments",
                      "deployment_end_type",
                      "animal_birth_hatch_year")),
             ~ if_else(. == "", NA_character_, .)),
      
      # Convert deploy timestamps 
      across(any_of(c("deploy_on_timestamp", "deploy_off_timestamp")), as.Date))
  
  
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
                n_missing, if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
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
    n_missing <- sum(is.na(is.na(summary_table$deploy_off_timestamp)))
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
              n_removed, if (n_removed == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
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
      logger.info(paste0("Warning: Removed ", n_removed, " individual(s) because deploy_off_timestamp occurred within ", censor_capture_mortality, " day(s) after deploy_on_timestamp"),
                  call. = FALSE, immediate. = TRUE)
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
                n_removed, if (n_removed == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
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
  
  # Search in data
  summary_table <- summary_table %>%
    
    # Initialize mortality event
    mutate(mortality_event = NA_real_) %>%
    
    # Identify survivors (individuals who last beyond study)
    mutate(survived_beyond_study = !is.na(raw_deploy_off_timestamp) &
             raw_deploy_off_timestamp > as.Date(effective_end),
           mortality_event = if_else(survived_beyond_study, 0L, mortality_event),
           
           # Update columns to remove ambiguity (e.g., if animal dies after study window)
           death_comments = if ("death_comments" %in% names(.)) {
             if_else(survived_beyond_study, "survived beyond study", death_comments)
           } else death_comments,
           
           deployment_end_comments = if ("deployment_end_comments" %in% names(.)) {
             if_else(survived_beyond_study, "survived beyond study", deployment_end_comments)
           } else deployment_end_comments,
           
           deployment_end_type = if ("deployment_end_type" %in% names(.)) {
             if_else(survived_beyond_study, "survived beyond study", deployment_end_type)
           } else deployment_end_type,
           
           mortality_location_filled = if ("mortality_location_filled" %in% names(.)) {
             if_else(survived_beyond_study, 0L, mortality_location_filled)
           } else mortality_location_filled) %>%
    
    # Search for mortality indicators
    # A. "death_comments" keywords
    mutate(mortality_event = case_when(
      "death_comments" %in% names(.) & str_detect(tolower(death_comments), "\\bnot\\b") ~ 0L,
      "death_comments" %in% names(.) & str_detect(tolower(death_comments), positive_pattern) ~ 1L,
      mortality_event == 1L ~ 1L,
      TRUE ~ mortality_event)) %>%
    
    # B: "deployment_end_comments" contains mortality keywords  
    mutate(mortality_event = case_when(
      mortality_event == 1L ~ 1L,  
      "deployment_end_comments" %in% names(.) &
        str_detect(tolower(deployment_end_comments), positive_pattern) ~ 1L,
      TRUE ~ mortality_event)) %>%
    
    # C: "mortality_type" is filled (non-NA) 
    mutate(mortality_event = case_when(
      mortality_event == 1L ~ 1L,   
      "mortality_type" %in% names(.) &
        !is.na(mortality_type) ~ 1L,
      TRUE ~ mortality_event)) %>%
    
    # D. "mortality_location_filled" is filled 
    mutate(mortality_event = case_when(
      "mortality_location_filled" %in% names(.) &
        mortality_location_filled >= 1 ~ 1L,
      mortality_event == 1L ~ 1L,
      TRUE ~ mortality_event)) %>%
    
    # E. "deployment_end_type" indicates censoring vs. death 
    mutate(mortality_event = case_when(
      mortality_event == 1L ~ 1L,
      
      # Mortality indication
      "deployment_end_type" %in% names(.) &
        str_detect(tolower(deployment_end_type), "\\bdead\\b|\\bdeath\\b") ~ 1L,
      
      # Censoring indication
      "deployment_end_type" %in% names(.) &
        tolower(deployment_end_type) %in% c("removal", "other", "unknown", "survived beyond study") ~ 0L,
      
      # Missing column OR NA value → censored
      (!"deployment_end_type" %in% names(.) | is.na(deployment_end_type)) &
        is.na(mortality_event) ~ 0L,
      TRUE ~ mortality_event)) %>%
    
    # Final censoring: remaining NA → 0 (only if we have "deploy_off_timestamp")
    mutate(mortality_event = if_else(
      is.na(mortality_event) & !is.na(deploy_off_timestamp),
      0L,
      mortality_event)) %>%
    
    # Clean data frame  
    select(-survived_beyond_study) %>%
    relocate(mortality_event, .after = deployment_end_type)
  
  # Error out: No deaths 
  n_mort_events <- sum(summary_table$mortality_event == 1, na.rm = TRUE)
  if (n_mort_events == 0) {
    logger.fatal("Cannot run survival analysis: no mortality events detected.",
                 call. = FALSE, immediate. = TRUE)
  }
  
  # Produce warning: Small proportion of deaths 
  if (n_mort_events <= 10) {
    logger.warn(sprintf("Few (%d) deaths detected across entire dataset. Particularly if data is further subset, model may have low statistical power. This could potentially result in unreliable estimates and poor predictive capacity", n_mort_events), call. = FALSE, immediate. = TRUE)
  }
  
  
  ## Calculate survival years (if selected) -----------------------------------
  
  if(!is.null(survival_yr_start)){
    
    # Extract survival year start date 
    start_month <- month(survival_yr_start)
    start_day   <- mday(survival_yr_start)     
    
    # Function to handle invalid dates (e.g., Feb 29 in non-leap years)
    safe_make_date <- function(year, month, day) {
      date <- suppressWarnings(make_date(year, month, day))
      if (is.na(date)) {
        ym <- ymd(sprintf("%d-%02d-01", year, month)) %m-% months(1)
        date <- ceiling_date(ym, "month") - days(1)
      }
      date
    }
    
    # Function to assign survival year and period boundaries 
    get_survival_period <- function(date) {
      if (is.na(date)) return(tibble(survival_year = NA_integer_, period_start = NA_Date_, 
                                     period_end = NA_Date_))
      y <- year(date)
      period_start_this_year <- safe_make_date(y, start_month, start_day)
      
      if (date >= period_start_this_year) {
        survival_year <- y
        period_start  <- period_start_this_year
        period_end    <- safe_make_date(y + 1, start_month, start_day) - days(1)
      } else {
        survival_year <- y - 1
        period_start  <- safe_make_date(y - 1, start_month, start_day)
        period_end    <- safe_make_date(y, start_month, start_day) - days(1)
      }
      
      tibble(survival_year = survival_year,
             period_start  = period_start,
             period_end    = period_end)
    }
    
    # Vectorized helpers
    get_survival_year  <- function(date) get_survival_period(date)$survival_year
    
    # Determine range of survival years present in the data 
    date_range <- events_with_ind %>%
      summarise(min_ts = min(timestamp, na.rm = TRUE),
                max_ts = max(timestamp, na.rm = TRUE))
    min_year <- get_survival_year(date_range$min_ts)
    max_year <- get_survival_year(date_range$max_ts)
    possible_years <- seq(min_year, max_year, by = 1)
    logger.info(sprintf("Survival years found in data: %d to %d (%d years total)",
                        min_year, max_year, length(possible_years)))
    
    # Create all possible survival periods
    possible_periods <- tibble(survival_year = possible_years) %>%
      mutate(period_info = map(survival_year, ~ {
        start <- safe_make_date(.x, start_month, start_day)
        end   <- safe_make_date(.x + 1, start_month, start_day) - days(1)
        tibble(period_start = start, period_end = end)
      })) %>%
      unnest(period_info)
    
    # Create individual-year rows 
    yearly_survival <- summary_table %>%
      
      # Keep only static/animal-level info for crossing
      dplyr::select(individual_id,
                    individual_local_identifier,
                    any_of(c("sex",
                             "animal_birth_hatch_year",
                             "attachment_type",
                             "model"))) %>% 
      distinct() %>% 
      
      # Cross with all possible survival years
      crossing(possible_periods) %>%
      
      # Bring back deployment & mortality info
      left_join(summary_table %>%
                  dplyr::select(individual_id, 
                                individual_local_identifier, 
                                deploy_on_timestamp, 
                                deploy_off_timestamp, 
                                any_of(c("mortality_date", 
                                         "mortality_type", 
                                         "death_comments", 
                                         "deployment_end_comments", 
                                         "deployment_end_type", 
                                         "mortality_event", 
                                         "first_timestamp", "last_timestamp", 
                                         "n_locations", "n_deployments"))), 
                by = c("individual_id", "individual_local_identifier")) %>%
      
      # Clip monitoring interval to the survival period
      mutate(monitor_start = pmax(deploy_on_timestamp,  period_start, na.rm = TRUE),
             monitor_end   = pmin(deploy_off_timestamp, period_end,   na.rm = TRUE),
             
             # Keep only periods where animal was monitored
             active_in_period = monitor_start <= monitor_end & !is.na(monitor_start) & 
               !is.na(monitor_end)) %>%
      filter(active_in_period) %>%
      
      # Final entry / exit dates for this animal-year
      mutate(entry_date = monitor_start,
             exit_date  = monitor_end,
             
             # Did mortality occur **inside** this survival year?
             died_this_year = case_when(
               
               # Priority 1: mortality_date exists and is inside the period
               mortality_event == 1L &
                 !is.na(mortality_date) &
                 mortality_date >= period_start &
                 mortality_date <= period_end
               ~ TRUE,
               
               # Priority 2: mortality_date is NA, but mortality_event == 1 AND deploy_off inside period
               mortality_event == 1L &
                 is.na(mortality_date) &
                 !is.na(deploy_off_timestamp) &
                 deploy_off_timestamp >= period_start &
                 deploy_off_timestamp <= period_end
               ~ TRUE,
               
               # Otherwise: no death this year
               TRUE ~ FALSE),
             
             # Final flags
             mortality_event = as.integer(died_this_year),
             censored        = as.integer(!died_this_year),
             
             # Reported mortality date: prefer original mortality_date, fall back to deploy_off
             mortality_date_reported = case_when(
               died_this_year & !is.na(mortality_date) ~ as.Date(mortality_date),
               died_this_year &  is.na(mortality_date) ~ as.Date(deploy_off_timestamp),
               TRUE                                    ~ NA_Date_),
             
             # Carry mortality metadata only when we flag a death
             mortality_type          = if_else(died_this_year, mortality_type,          NA_character_),
             death_comments          = if_else(died_this_year, death_comments,          NA_character_),
             deployment_end_comments = if_else(died_this_year, deployment_end_comments, NA_character_),
             deployment_end_type     = if_else(died_this_year, deployment_end_type,     NA_character_),
             
             # Days monitored in this survival year
             days_at_risk = as.integer(exit_date - entry_date) + 1L) %>%
      
      # Final cleaning
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
      mutate(age       = survival_year - animal_birth_hatch_year,
             age       = as.integer(pmax(0, age)))
    
    # Prepare thresholds from your existing table  
    thresholds <- animal_birth_hatch_year_table %>%
      filter(!is.na(year_at_start)) %>%           
      arrange(year_at_start) %>%
      distinct(year_at_start, animal_life_stage)
    
    # Create a named vector for fast look-up
    stage_lookup <- setNames(thresholds$animal_life_stage,
                             thresholds$year_at_start)
    
    # Dynamic assignment 
    yearly_survival <- yearly_survival %>%
      mutate(matched_threshold = findInterval(age, thresholds$year_at_start),
             animal_life_stage_new = case_when(
               is.na(age)                                ~ "unknown",           
               matched_threshold == 0                    ~ "unknown",            
               TRUE                                      ~ stage_lookup[matched_threshold]),
             animal_life_stage_new = coalesce(
               animal_life_stage_new,
               animal_birth_hatch_year_table %>%
                 filter(is.na(year_at_start)) %>%
                 pull(animal_life_stage) %>%
                 first(default = "adult"))) %>% 
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
    
    warning(paste0("Subsetting by ", spec$label, " (", define, ")"))
    value <- spec$coerce(define)
    
    if (is.null(survival_yr_start)) {
      if (!spec$summary_ok) {
        warning("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      } else {
        summary_table <- summary_table %>% filter(.data[[spec$col]] == value)
      }
    } else {
      if (!spec$yearly_ok) {
        warning("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      } else {
        yearly_survival <- yearly_survival %>% filter(.data[[spec$col]] == value)
      }
    }
    
    list(summary_table = summary_table, yearly_survival = yearly_survival)
  }
  
  # SUBSET CONDITION 1 ---
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
  
  # SUBSET CONDITION 2 ---
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
    warning(paste("Comparing across", cox_covariate_1, "..."))
    
    if(is.null(survival_yr_start)){
      warning("... with summary table")
      
      non_na_unique <- unique(na.omit(summary_table[[cox_covariate_1]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          warning(paste0("Warning: The grouping variable, ", cox_covariate_1, ", is entirely NA; no comparison is possible."))
        } else {
          warning(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_1, "; no comparison is possible."))
        }
        cox_covariate_1 <- NULL
      }
      
      if (!is.null(cox_covariate_1)) {
        n_original    <- nrow(summary_table)
        summary_table <- summary_table[!is.na(summary_table[[cox_covariate_1]]) & 
                                         trimws(summary_table[[cox_covariate_1]]) != "", ]
        
        unique_values <- sort(unique(summary_table[[cox_covariate_1]]))
        n_values      <- length(unique_values)
        
        warning(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(summary_table)
        if (n_lost > 0) {
          warning(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
      
    } else {
      warning("... with yearly survival")
      
      non_na_unique <- unique(na.omit(yearly_survival[[cox_covariate_1]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          warning(paste0("Warning: The grouping variable, ", cox_covariate_1, ", is entirely NA; no comparison is possible."))
        } else {
          warning(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_1, "; no comparison is possible."))
        }
        cox_covariate_1 <- NULL
      }
      
      if (!is.null(cox_covariate_1)) {
        n_original      <- nrow(yearly_survival)
        yearly_survival <- yearly_survival[!is.na(yearly_survival[[cox_covariate_1]]) & 
                                             trimws(yearly_survival[[cox_covariate_1]]) != "", ]
        
        unique_values <- sort(unique(yearly_survival[[cox_covariate_1]]))
        n_values      <- length(unique_values)
        
        warning(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(yearly_survival)
        if (n_lost > 0) {
          warning(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
    }
  }
  
  # Comparison covariate 2 --- 
  if(!is.null(cox_covariate_2)) {
    warning(paste("Comparing across", cox_covariate_2, "..."))
    
    if(is.null(survival_yr_start)){
      warning("... with summary table")
      
      non_na_unique <- unique(na.omit(summary_table[[cox_covariate_2]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          warning(paste0("Warning: The grouping variable, ", cox_covariate_2, ", is entirely NA; no comparison is possible."))
        } else {
          warning(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_2, "; no comparison is possible."))
        }
        cox_covariate_2 <- NULL
      }
      
      if (!is.null(cox_covariate_2)) {
        n_original    <- nrow(summary_table)
        summary_table <- summary_table[!is.na(summary_table[[cox_covariate_2]]) & 
                                         trimws(summary_table[[cox_covariate_2]]) != "", ]
        
        unique_values <- sort(unique(summary_table[[cox_covariate_2]]))
        n_values      <- length(unique_values)
        
        warning(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(summary_table)
        if (n_lost > 0) {
          warning(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
      
    } else {
      warning("... with yearly survival")
      
      non_na_unique <- unique(na.omit(yearly_survival[[cox_covariate_2]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          warning(paste0("Warning: The grouping variable, ", cox_covariate_2, ", is entirely NA; no comparison is possible."))
        } else {
          warning(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_2, "; no comparison is possible."))
        }
        cox_covariate_2 <- NULL
      }
      
      if (!is.null(cox_covariate_2)) {
        n_original      <- nrow(yearly_survival)
        yearly_survival <- yearly_survival[!is.na(yearly_survival[[cox_covariate_2]]) & 
                                             trimws(yearly_survival[[cox_covariate_2]]) != "", ]
        
        unique_values <- sort(unique(yearly_survival[[cox_covariate_2]]))
        n_values      <- length(unique_values)
        
        warning(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(yearly_survival)
        if (n_lost > 0) {
          warning(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
    }
  }
  
  # Comparison covariate 3 ---  
  if(!is.null(cox_covariate_3)) {
    warning(paste("Comparing across", cox_covariate_3, "..."))
    
    if(is.null(survival_yr_start)){
      warning("... with summary table")
      
      non_na_unique <- unique(na.omit(summary_table[[cox_covariate_3]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          warning(paste0("Warning: The grouping variable, ", cox_covariate_3, ", is entirely NA; no comparison is possible."))
        } else {
          warning(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_3, "; no comparison is possible."))
        }
        cox_covariate_3 <- NULL
      }
      
      if (!is.null(cox_covariate_3)) {
        n_original    <- nrow(summary_table)
        summary_table <- summary_table[!is.na(summary_table[[cox_covariate_3]]) & 
                                         trimws(summary_table[[cox_covariate_3]]) != "", ]
        
        unique_values <- sort(unique(summary_table[[cox_covariate_3]]))
        n_values      <- length(unique_values)
        
        warning(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(summary_table)
        if (n_lost > 0) {
          warning(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
      
    } else {
      warning("... with yearly survival")
      
      non_na_unique <- unique(na.omit(yearly_survival[[cox_covariate_3]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          warning(paste0("Warning: The grouping variable, ", cox_covariate_3, ", is entirely NA; no comparison is possible."))
        } else {
          warning(paste0("Warning: There is only one non-NA comparison covariate in ", cox_covariate_3, "; no comparison is possible."))
        }
        cox_covariate_3 <- NULL
      }
      
      if (!is.null(cox_covariate_3)) {
        n_original      <- nrow(yearly_survival)
        yearly_survival <- yearly_survival[!is.na(yearly_survival[[cox_covariate_3]]) & 
                                             trimws(yearly_survival[[cox_covariate_3]]) != "", ]
        
        unique_values <- sort(unique(yearly_survival[[cox_covariate_3]]))
        n_values      <- length(unique_values)
        
        warning(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(yearly_survival)
        if (n_lost > 0) {
          warning(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
    }
  }
  
  
  ## Basic summaries of data ------------------------------------------------
  
  # Plot each individual's tracking history --- 
  
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
      geom_point(aes(x = plot_end, y = individual_label),
                 color = "#9467BD", size = 3.5) +
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
    
    # Save plot  
    png(appArtifactPath("tracking_history.png"), 
        width = 10, height = 8,
        units = "in", res = 300)
    print(tracking_history)
    dev.off()
    
    
    # Calculate monthly mortality
    if (calc_month_mort == TRUE) {
      
      logger.info("Plotting monthly mortality using summary table...")
      
      min_date <- min(summary_table$deploy_on_timestamp, na.rm = TRUE)
      max_date <- max(summary_table$deploy_off_timestamp, na.rm = TRUE)
      full_years <- seq(year(min_date), year(max_date), by = 1)
      
      mortality_data <- summary_table %>%
        filter(mortality_event == 1) %>%
        mutate(death_date   = as.Date(deploy_off_timestamp),
               death_year   = year(death_date),
               death_month  = month(death_date, label = TRUE, abbr = TRUE),
               death_month  = factor(death_month, levels = month.abb, ordered = TRUE)) %>%
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
        scale_fill_viridis_d(option    = "magma",
                             direction = -1,
                             na.value  = "grey92",
                             name      = "Number of\nmortality events",
                             drop      = FALSE) +
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
              plot.title        = element_text(face = "bold", hjust = 0.5, size = 16),
              plot.subtitle     = element_text(hjust = 0.5, size = 12),
              axis.text.x       = element_text(size = 11, face = "bold"),
              axis.text.y       = element_text(size = 11))
      
      # Save plot  
      png(appArtifactPath("monthly_mortality.png"), 
          width = 10, height = 8, 
          units = "in", res = 300)
      print(monthly_mort_plot)
      dev.off()
      
    } else{
      logger.info("Not plotting monthly mortality")
    } 
  }
  
  if (!is.null(survival_yr_start)) {
    logger.info("Plotting tracking history using yearly survival")
    
    deployment_summary <- yearly_survival |>
      distinct(individual_id,
               individual_local_identifier,
               deploy_on_timestamp,
               deploy_off_timestamp) |>
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
      geom_point(aes(x = plot_end, y = individual_label),
                 color = "#9467BD", size = 3.5) +
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
    
    # Save plot  
    png(appArtifactPath("tracking_history.png"), 
        width = 10, height = 8, 
        units = "in", res = 300)
    print(tracking_history)
    dev.off()
    
    
    # Calculate monthly mortality 
    if (calc_month_mort == TRUE) {
      
      logger.info("Plotting monthly mortality using yearly survival")
      
      min_date <- min(yearly_survival$deploy_on_timestamp, na.rm = TRUE)
      max_date <- max(yearly_survival$deploy_off_timestamp, na.rm = TRUE)
      full_years <- seq(year(min_date), year(max_date), by = 1)
      
      mortality_data <- yearly_survival |>
        filter(mortality_event == 1 | died_this_year == TRUE) |>
        mutate(mort_date   = as.Date(mortality_date),
               deploy_off  = as.Date(deploy_off_timestamp),
               death_date  = coalesce(mort_date, deploy_off),
               death_year  = year(death_date),
               death_month = month(death_date, label = TRUE, abbr = TRUE),
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
        scale_fill_viridis_d(option = "magma",
                             direction = -1,
                             na.value = "grey92",
                             name = "Number of\nmortality events",
                             drop = FALSE) +
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
              plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              axis.text.x = element_text(size = 11, face = "bold"),
              axis.text.y = element_text(size = 11))
      
      # Save plot  
      png(appArtifactPath("monthly_mortality.png"), 
          width = 10, height = 8, 
          units = "in", res = 300)
      print(monthly_mort_plot)
      dev.off()
      
    } else {
      logger.info("Not plotting monthly mortality.")
    }
  } 
  
  
  ## Cox Proportional Hazard Analysis -----------------------------------------
  
  if (is.null(survival_yr_start)) {
    
    logger.info("Calculating Cox Proportional Hazards using summary table...")
    
    # Warning for no mortality events
    if (sum(summary_table$mortality_event) == 0) {
      logger.fatal("There are no mortality events in the chosen subset of data.")
    }
    
    ## Fit model --- 
    fitting_data <- summary_table
    
    # Collect non-NULL covariates
    covariates <- c(cox_covariate_1, cox_covariate_2, cox_covariate_3)
    ref_levels <- list(cox_covariate_1_ref, cox_covariate_2_ref, cox_covariate_3_ref)
    names(ref_levels) <- c(cox_covariate_1, cox_covariate_2, cox_covariate_3)
    covariates <- covariates[!is.null(covariates) & !is.na(covariates)]
    
    # Apply reference levels where specified
    for (cov in covariates) {
      ref <- ref_levels[[cov]]
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
    firth_used         <- FALSE
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
    
    # Apply Firth's if separation warning was triggered
    if (separation_detected) {
      
      logger.warn("Separation detected in coxph — falling back to Firth's penalized Cox model.")
      coxphf_fit <- coxphf(cox_formula, data = fitting_data,
                           maxstep = 0.1,   
                           maxit   = 200) 
      firth_used  <- TRUE
      
      # For Firth's, align vectors by coefficient names 
      coef_names <- names(coxphf_fit$coefficients)
      cox.tab <- data.frame(term      = coef_names,
                            estimate  = exp(coxphf_fit$coefficients[coef_names]),
                            conf.low  = exp(coxphf_fit$ci.lower[coef_names]),
                            conf.high = exp(coxphf_fit$ci.upper[coef_names]),
                            p.value   = coxphf_fit$prob[coef_names],
                            row.names = NULL)
      
    } else {
      cox.tab <- tidy(coxph_fit, exponentiate = TRUE, conf.int = TRUE)
    }
    
    # Flag which model used
    cox.tab$model <- ifelse(firth_used, "Firth's Penalized Cox", "Standard Cox")
    
    # Save
    write.csv(cox.tab, file = appArtifactPath("coxph_table.csv"), row.names = FALSE)
    
    # Predict survival function for subject with means on all covariates
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
  }
  
  if (!is.null(survival_yr_start)) {
    
    logger.info("Calculating KM using survival table...")
    
    # Warning for no mortality events
    if (sum(yearly_survival$mortality_event) == 0) {
      logger.fatal("There are no mortality events in the chosen subset of data.")
    }
    
    ## Fit model --
    fitting_data <- yearly_survival
    
    # Collect non-NULL covariates
    covariates <- c(cox_covariate_1, cox_covariate_2, cox_covariate_3)
    ref_levels <- list(cox_covariate_1_ref, cox_covariate_2_ref, cox_covariate_3_ref)
    names(ref_levels) <- c(cox_covariate_1, cox_covariate_2, cox_covariate_3)
    
    covariates <- covariates[!is.null(covariates) & !is.na(covariates)]
    
    # Apply reference levels where specified
    for (cov in covariates) {
      ref <- ref_levels[[cov]]
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
    
    # Apply Firth's if separation warning was triggered
    if (separation_detected) {
      
      logger.warn("Separation detected in coxph — falling back to Firth's penalized Cox model.")
      coxphf_fit <- coxphf(cox_formula, data = fitting_data,
                           maxstep = 0.1,
                           maxit   = 200)
      firth_used <- TRUE
      
      # For Firth's, align vectors by coefficient names 
      coef_names <- names(coxphf_fit$coefficients)
      
      cox.tab <- data.frame(term      = coef_names,
                            estimate  = exp(coxphf_fit$coefficients[coef_names]),
                            conf.low  = exp(coxphf_fit$ci.lower[coef_names]),
                            conf.high = exp(coxphf_fit$ci.upper[coef_names]),
                            p.value   = coxphf_fit$prob[coef_names],
                            row.names = NULL)
      
    } else {
      cox.tab <- tidy(coxph_fit, exponentiate = TRUE, conf.int = TRUE)
    }
    
    # Flag which model used
    cox.tab$model <- ifelse(firth_used, "Firth's Penalized Cox", "Standard Cox")
    
    # Save
    write.csv(cox.tab, file = appArtifactPath("coxph_table.csv"), row.names = FALSE)
    
    # Predict survival function for subject with means on all covariates
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
  }
  
  
  ## Visualizations of results at covariate means -----------------------------
  
  # Plot of predicted survival --- 
  surv_plot <- ggsurvfit(surv.at.means) +
    add_confidence_interval() +
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
    theme(plot.title       = element_text(face = "bold", size = 14),
          plot.subtitle    = element_text(size = 10, color = "grey40"),
          axis.title       = element_text(face = "bold"),
          legend.position  = "none",
          panel.grid.major.y = element_line(color = "grey92"))
  
  # Save 
  png(appArtifactPath("surv_at_means_plot.png"), 
      width = 7, height = 6, 
      units = "in", res = 300)
  print(surv_plot)
  dev.off()
  
  
  ## Plot cumulative hazard curve --- 
  cum_hazard <- ggsurvfit(surv.at.means, type = "cumhaz") +
    add_confidence_interval() +
    labs(title = "Cumulative Hazard at Covariate Means", x = "Days", y = "Cumulative Hazard") +
    scale_x_continuous(breaks = pretty(range(surv.at.means.tab$time), n = 8)) +
    theme_classic(base_size = 12)
  
  # Save   
  png(appArtifactPath("cum_hazard_at_means_plot.png"), 
      width = 7, height = 6, 
      units = "in", res = 300)
  print(cum_hazard)
  dev.off()
  
  
  ## Plot Schoenfeld residuals ---
  if (!firth_used) {
    ph_test <- cox.zph(coxph_fit)
    print(ph_test)
    
    ph_plot <- ggcoxzph(ph_test,
                        point.col = "steelblue",
                        point.size = 1.5,
                        point.alpha = 0.5) 
    
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
    
    ph_plots <- suppressWarnings(
      ggcoxzph(ph_test,
               point.col   = "steelblue",
               point.size  = 1.5,
               point.alpha = 0.5)
    )
    
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
  
  
  ## Plot forest plot of hazard ratios --- 
  forest_plot <- ggplot(cox.tab, aes(x = estimate, y = term)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_log10() +  # HRs are naturally log-scaled
    labs(title    = "Hazard Ratios with 95% Confidence Intervals",
         subtitle = ifelse(firth_used, "Firth's Penalized Cox Model", "Standard Cox Model"),
         x        = "Hazard Ratio (log scale)",
         y        = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  # Save   
  png(appArtifactPath("forest_plot.png"), 
      width = 7, height = 5,
      units = "in", res = 300)
  print(forest_plot)
  dev.off()
  
  
  ## Stratified Group Comparisons ---------------------------------------------
  
  ## Prepare covariates --- 
  active_covariates <- list(
    list(var = cox_covariate_1, ref = cox_covariate_1_ref),
    list(var = cox_covariate_2, ref = cox_covariate_2_ref),
    list(var = cox_covariate_3, ref = cox_covariate_3_ref)
  )
  active_covariates <- Filter(function(x) !is.null(x$var), active_covariates)
  
  fitting_data <- if (is.null(survival_yr_start)) summary_table else yearly_survival
  
  for (cov_item in active_covariates) {
    cov <- cov_item$var
    ref <- cov_item$ref
    logger.info(paste("Producing group comparison outputs for:", cov))
    
    # Apply reference level if specified
    if (!is.null(ref) && ref %in% unique(fitting_data[[cov]])) {
      fitting_data[[cov]] <- relevel(factor(fitting_data[[cov]]), ref = ref)
    } else {
      fitting_data[[cov]] <- factor(fitting_data[[cov]])
    }
    
    n_groups <- length(levels(fitting_data[[cov]]))
    pal      <- viridis(n_groups, option = "turbo", end = 0.9)
    
    km_formula <- as.formula(paste("Surv(entry_time_days, exit_time_days, mortality_event) ~", cov))
    
    
    ## Cox model tests (LR, Wald, Score) --- 
    cox_strat_fit <- coxph(km_formula, data = fitting_data)
    cox_summary   <- summary(cox_strat_fit)
    
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
                          cox_summary$sctest["pvalue"]), 4))
    
    # LR p-value for plot subtitles
    lr_p <- test_results$p_value[1]
    
    # Save 
    write.csv(test_results,
              file      = appArtifactPath(paste0("model_tests_", cov, ".csv")),
              row.names = FALSE)

    
    ## Pairwise comparisons (if >2 groups) ---
    if (n_groups > 2) {
      group_levels     <- levels(fitting_data[[cov]])
      pairs            <- combn(group_levels, 2, simplify = FALSE)
      
      pairwise_results <- lapply(pairs, function(pair) {
        sub_data        <- fitting_data[fitting_data[[cov]] %in% pair, ]
        sub_data[[cov]] <- droplevels(sub_data[[cov]])
        fit             <- coxph(km_formula, data = sub_data)
        s               <- summary(fit)
        data.frame(group_1 = pair[1],
                   group_2 = pair[2],
                   hr      = round(exp(coef(fit)[1]), 3),
                   p_wald  = round(s$waldtest["pvalue"], 4))
      })
      
      pairwise_df              <- do.call(rbind, pairwise_results)
      pairwise_df$p_bonferroni <- round(p.adjust(pairwise_df$p_wald, method = "bonferroni"), 4)
      pairwise_df$significant  <- ifelse(pairwise_df$p_bonferroni < 0.05, "Yes", "No")
      
      # Save 
      write.csv(pairwise_df,
                file      = appArtifactPath(paste0("pairwise_", cov, ".csv")),
                row.names = FALSE)
      
      logger.info(sprintf("Pairwise comparisons for %s saved (%d pairs, Bonferroni corrected).",
                      cov, nrow(pairwise_df)))
    }
    
    
    ## HR table for specified covariate ---
    cox_strat_tab <- tidy(cox_strat_fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(covariate = cov)
    
    # Save 
    write.csv(cox_strat_tab,
              file      = appArtifactPath(paste0("cox_hr_", cov, ".csv")),
              row.names = FALSE)
    
    
    ## KM survival curves stratified by covariate ---------------------------
    km_fit <- survfit(km_formula, data = fitting_data)
    
    km_plot <- ggsurvfit(km_fit, linewidth = 0.9) +
      add_confidence_interval(alpha = 0.12) +
      add_risktable(risktable_stats = c("n.risk", "cum.event"),
                    theme = theme_risktable_default(axis.text.y.size = 9)) +
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
            legend.position    = "top",
            panel.grid.major.y = element_line(color = "grey92"))
    
    # Save 
    png(appArtifactPath(paste0("km_by_", cov, ".png")),
        width = 8, height = 7, 
        units = "in", res = 300)
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
    
    # Save
    write.csv(km_summary,
              file      = appArtifactPath(paste0("km_summary_", cov, ".csv")),
              row.names = FALSE)
    
    
    ## Forest plot (per-group HRs) ---
    forest_strat <- ggplot(cox_strat_tab, aes(x = estimate, y = term)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                     height = 0.25, linewidth = 0.8, color = "grey30") +
      geom_point(aes(color = p.value < 0.05), size = 3.5) +
      scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "#1F77B4"),
                         labels = c("TRUE" = "p < 0.05",  "FALSE" = "p ≥ 0.05"),
                         name   = NULL) +
      scale_x_log10(labels = number_format(accuracy = 0.01)) +
      labs(title    = paste("Hazard Ratios —", cov),
           subtitle = sprintf("Cox LR test: chi-sq=%.3f, p=%.4f", test_results$chi_sq[1], lr_p),
           x        = "Hazard Ratio (log scale)",
           y        = NULL) +
      theme_classic(base_size = 12) +
      theme(plot.title      = element_text(face = "bold"),
            plot.subtitle   = element_text(size = 10, color = "grey40"),
            legend.position = "top")
    
    # Save   
    png(appArtifactPath(paste0("forest_", cov, ".png")), 
        width = 7, height = 4,
        units = "in", res = 300)
    print(forest_strat)
    dev.off()
    
    
    ## Cumulative hazard by group -------------------------------------------
    cumhaz_plot <- ggsurvfit(km_fit, type = "cumhaz", linewidth = 0.9) +
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
            legend.position = "top")
    
    # Save
    png(appArtifactPath(paste0("cumhaz_by_", cov, ".png")),
        width = 7, height = 5, 
        units = "in", res = 300)
    print(cumhaz_plot)
    dev.off()
    
    
    ## Log cumulative hazard (log-log) plot: visual PH assumption check -----
    logger.info("If parallel lines = proportional hazards assumption holds")
    
    km_tidy <- tidy(km_fit) %>%
      filter(estimate > 0, estimate < 1) %>%
      mutate(log_time     = log(time),
             log_log_surv = log(-log(estimate)),
             group        = strata)
    
    loglog_plot <- ggplot(km_tidy, aes(x = log_time, y = log_log_surv, color = group)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = pal, name = cov) +
      labs(title    = paste("Log-Log Survival Plot —", cov),
           subtitle = "Parallel lines support the proportional hazards assumption",
           x        = "log(Time)",
           y        = "log(-log(Survival))") +
      theme_classic(base_size = 12) +
      theme(plot.title      = element_text(face = "bold"),
            plot.subtitle   = element_text(size = 10, color = "grey40"),
            legend.position = "top")
    
    # Save   
    png(appArtifactPath(paste0("loglog_", cov, ".png")), 
        width = 7, height = 5,
        units = "in", res = 300)
    print(loglog_plot)
    dev.off()
    
    
    ## Annual survival rate bar chart (yearly mode only) --------------------
    if (!is.null(survival_yr_start)) {
      
      annual_surv <- fitting_data %>%
        group_by(survival_year, .data[[cov]]) %>%
        summarise(n_animals = n_distinct(individual_id),
                  n_deaths  = sum(mortality_event, na.rm = TRUE),
                  surv_rate = 1 - (n_deaths / n_animals),
                  .groups   = "drop")
      
      annual_plot <- ggplot(annual_surv, aes(x = factor(survival_year), 
                                             y = surv_rate, fill = .data[[cov]])) +
        geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
        geom_text(aes(label = paste0("n=", n_animals)),
                  position = position_dodge(0.8),
                  vjust = -0.4, size = 3, color = "grey30") +
        scale_fill_manual(values = pal) +
        scale_y_continuous(limits = c(0, 1.05),
                           labels = percent_format(accuracy = 1)) +
        labs(title = paste("Annual Survival Rate by", cov),
             x     = "Survival Year",
             y     = "Survival Rate",
             fill  = cov) +
        theme_classic(base_size = 12) +
        theme( plot.title      = element_text(face = "bold"),
               legend.position = "top",
               axis.text.x     = element_text(angle = 45, hjust = 1))
      
      # Save   
      png(appArtifactPath(paste0("annual_surv_by_", cov, ".png")), 
          width = 8, height = 5,
          units = "in", res = 300)
      print(annual_plot)
      dev.off()
      
      write.csv(annual_surv,
                file      = appArtifactPath(paste0("annual_surv_", cov, ".csv")),
                row.names = FALSE)
    }
  }
  
  # Pass original to the next app in the MoveApps workflow
  return(data)
} 