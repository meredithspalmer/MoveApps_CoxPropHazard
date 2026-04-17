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
                     group_comparison_individual,
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
    
    # C. "mortality_location_filled" is filled 
    mutate(mortality_event = case_when(
      "mortality_location_filled" %in% names(.) &
        mortality_location_filled >= 1 ~ 1L,
      mortality_event == 1L ~ 1L,
      TRUE ~ mortality_event)) %>%
    
    # D. "deployment_end_type" indicates censoring vs. death 
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
  
  # SUBSET CONDITION 1 ---
  
  if (!is.null(subset_condition_1) && subset_condition_1 == "sex") {
    logger.info(paste0("Subsetting by sex (", subset_condition_define_1, ")"))
    
    if (is.null(survival_yr_start)) {
      summary_table <- summary_table %>% filter(sex == subset_condition_define_1)
    } else {
      yearly_survival <- yearly_survival %>% filter(sex == subset_condition_define_1)
    }
  }
  
  
  if (!is.null(subset_condition_1) && subset_condition_1 == "attachment_type") {
    logger.info(paste0("Subsetting by attachment type (", subset_condition_define_1, ")"))
    
    if (is.null(survival_yr_start)) {
      summary_table <- summary_table %>% filter(attachment_type == subset_condition_define_1)
    } else {
      yearly_survival <- yearly_survival %>% filter(attachment_type == subset_condition_define_1)
    }
  }
  
  if (!is.null(subset_condition_1) && subset_condition_1 == "model") {
    logger.info(paste0("Subsetting by model (", subset_condition_define_1, ")"))
    
    if (is.null(survival_yr_start)) {
      summary_table <- summary_table %>% filter(model == subset_condition_define_1) 
    } else {
      yearly_survival <- yearly_survival %>% filter(model == as.integer(subset_condition_define_1))
    }
  }
  
  if (!is.null(subset_condition_1) && subset_condition_1 == "lifestage") {
    logger.info(paste0("Subsetting by life stage (", subset_condition_define_1, ")"))
    
    if (is.null(survival_yr_start)) {
      logger.warn("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      
    } else {
      yearly_survival <- yearly_survival %>% filter(animal_life_stage_new == subset_condition_define_1)
    }
  }
  
  if (!is.null(subset_condition_1) && subset_condition_1 == "survival_year") {
    logger.info(paste0("Subsetting by survival year (", subset_condition_define_1, ")"))
    
    if (is.null(survival_yr_start)) {
      logger.warn("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      
    } else {
      yearly_survival <- yearly_survival %>% filter(survival_year == as.integer(subset_condition_define_1))
    }
  }
  
  if(!is.null(subset_condition_1)){
    if(is.null(survival_yr_start) && nrow(summary_table) == 0){
      logger.fatal("There are no individuals meeting the first subsetting condition.")
    } else if(!is.null(survival_yr_start) && nrow(yearly_survival) == 0){
      logger.fatal("There are no individuals meeting the first subsetting condition.")
    }
  }
  
  
  # SUBSET CONDITION 2 ---
  
  if (!is.null(subset_condition_2) && subset_condition_2 == "sex") {
    logger.info(paste0("Subsetting by sex (", subset_condition_define_2, ")"))
    
    if (is.null(survival_yr_start)) {
      summary_table <- summary_table %>% filter(sex == subset_condition_define_2)
    } else {
      yearly_survival <- yearly_survival %>% filter(sex == subset_condition_define_2)
    }
  } 
  
  if (!is.null(subset_condition_2) && subset_condition_2 == "attachment_type") {
    logger.info(paste0("Subsetting by attachment type (", subset_condition_define_2, ")"))
    
    if (is.null(survival_yr_start)) {
      summary_table <- summary_table %>% filter(attachment_type == subset_condition_define_2)
    } else {
      yearly_survival <- yearly_survival %>% filter(attachment_type == subset_condition_define_2)
    }
  }
  
  if (!is.null(subset_condition_2) && subset_condition_2 == "model") {
    logger.info(paste0("Subsetting by model (", subset_condition_define_2, ")"))
    
    if (is.null(survival_yr_start)) {
      summary_table <- summary_table %>% filter(model == subset_condition_define_2) 
    } else {
      yearly_survival <- yearly_survival %>% filter(model == subset_condition_define_2)
    }
  }
  
  if (!is.null(subset_condition_2) && subset_condition_2 == "lifestage") {
    logger.info(paste0("Subsetting by life stage (", subset_condition_define_2, ")"))
    
    if (is.null(survival_yr_start)) {
      logger.warn("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      
    } else {
      yearly_survival <- yearly_survival %>% filter(animal_life_stage_new == subset_condition_define_2)
    }
  }
  
  if (!is.null(subset_condition_2) && subset_condition_2 == "survival_year") {
    logger.info(paste0("Subsetting by survival year (", subset_condition_define_2, ")"))
    
    if (is.null(survival_yr_start)) {
      logger.warn("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      
    } else {
      yearly_survival <- yearly_survival %>% filter(survival_year == as.integer(subset_condition_define_2))
    }
  }
  
  if(!is.null(subset_condition_2)){
    if(is.null(survival_yr_start) && nrow(summary_table) == 0){
      logger.fatal("There are no individuals meeting both subsetting conditions.")
    } else if(!is.null(survival_yr_start) && nrow(yearly_survival) == 0){
      logger.fatal("There are no individuals meeting both subsetting conditions.")
    }
  }
  
  
  ## Clean user-defined grouping attributes (if selected) ---------------------
  
  if(!is.null(group_comparison_individual) && group_comparison_individual == "sex"){
    logger.info("Grouping by sex...")
    
    if(is.null(survival_yr_start)){
      logger.info("... with summary table")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(summary_table$sex))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # Remove NAs 
      n_original <- nrow(summary_table)
      summary_table <- summary_table[!is.na(summary_table$sex),] 
      
      # Get unique sexes after cleaning
      unique_sexes <- sort(unique(summary_table$sex))
      n_sexes <- length(unique_sexes)
      
      logger.info(sprintf("%d sexes detected after cleaning: %s", n_sexes, paste(unique_sexes, collapse = ", ")),
                  call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(summary_table)
      if (n_lost > 0) {
        logger.warn(sprintf("%d individuals with NA sex removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
      }
      
    } else {
      logger.info("... with yearly survival")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(yearly_survival$sex))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # remove NAs 
      n_original      <- nrow(yearly_survival)
      yearly_survival <- yearly_survival[!is.na(yearly_survival$sex),] 
      
      # Get unique sexes after cleaning
      unique_sexes <- sort(unique(yearly_survival$sex))
      n_sexes      <- length(unique_sexes)
      
      logger.info(sprintf("%d sexes detected after cleaning: %s", n_sexes, paste(unique_sexes, collapse = ", ")),
                  call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(yearly_survival)
      if (n_lost > 0) {
        logger.info(sprintf("%d individuals with NA sex removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
      }
    }
  }
  
  if(!is.null(group_comparison_individual) && group_comparison_individual == "lifestage"){
    logger.info("Grouping by life stage...")
    
    if(is.null(survival_yr_start)){ 
      logger.error("This comparison only makes sense if survival years are defined.")
      
    } else {
      logger.info("... with yearly survival")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(yearly_survival$animal_life_stage_new))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # Remove NAs 
      n_original      <- nrow(yearly_survival)
      yearly_survival <- yearly_survival[!is.na(yearly_survival$animal_life_stage_new),] 
      
      # Get unique life-stages after cleaning
      unique_stages  <- sort(unique(yearly_survival$animal_life_stage_new))
      n_life_stages  <- length(unique_stages)
      
      logger.info(sprintf("%d life-stages detected after cleaning: %s", n_life_stages, 
                          paste(unique_stages, collapse = ", ")),
                  call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(yearly_survival)
      if (n_lost > 0) {
        logger.info(sprintf("%d individuals with NA life-stage removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
      }
    }
  }
  
  if(!is.null(group_comparison_individual) && group_comparison_individual == "model"){
    logger.info("Group by model...")    
    
    if(is.null(survival_yr_start)){
      logger.info("... with summary table")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(summary_table$model))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # Clean data, remove NAs 
      n_original <- nrow(summary_table)
      summary_table <- summary_table %>%
        filter(!is.na(model)) %>%
        mutate(model = str_trim(model),
               model = str_replace_all(model, "\\s+", ""),
               model = str_extract(model, "^[^|]+"),
               model = str_replace(model, "–", "-"))
      
      # Get unique conditions after cleaning
      unique_conditions <- sort(unique(summary_table$model))
      n_conditions <- length(unique_conditions)
      
      logger.info(
        sprintf("%d models detected after cleaning: %s", n_conditions, paste(unique_conditions, collapse = ", ")),
        call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(summary_table)
      if (n_lost > 0) {
        logger.info(sprintf("%d individuals with NA model removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
      }
      
    } else {
      logger.info("... with yearly survival")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(yearly_survival$model))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # Clean data, remove NAs 
      n_original      <- nrow(yearly_survival)
      yearly_survival <- yearly_survival %>%
        filter(!is.na(model)) %>%
        mutate(model = str_trim(model),
               model = str_replace_all(model, "\\s+", ""),
               model = str_extract(model, "^[^|]+"),
               model = str_replace(model, "–", "-"))
      
      # Get unique conditions after cleaning
      unique_conditions <- sort(unique(yearly_survival$model))
      n_conditions      <- length(unique_conditions)
      
      logger.info(sprintf("%d models detected after cleaning: %s", n_conditions, 
                          paste(unique_conditions, collapse = ", ")),
                  call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(yearly_survival)
      if (n_lost > 0) {
        logger.info(sprintf("%d individuals with NA model removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
      }
    } 
  }
  
  if(!is.null(group_comparison_individual) && group_comparison_individual == "attachment_type"){
    logger.info("Grouping by attachment...")
    
    if(is.null(survival_yr_start)){
      logger.info("... with summary table")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(summary_table$attachment_type))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # Clean data, remove NAs 
      n_original    <- nrow(summary_table)
      summary_table <- summary_table %>%
        filter(!is.na(attachment_type)) %>%
        mutate(attachment_type = str_trim(attachment_type),
               attachment_type = str_replace_all(attachment_type, "\\s+", ""),
               attachment_type = str_extract(attachment_type, "^[^|]+"),
               attachment_type = str_replace(attachment_type, "–", "-"))
      
      # Get unique conditions after cleaning
      unique_conditions <- sort(unique(summary_table$attachment_type))
      n_conditions <- length(unique_conditions)
      
      logger.info(sprintf("%d attachment types detected after cleaning: %s", n_conditions,
                          paste(unique_conditions, collapse = ", ")),
                  call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(summary_table)
      if (n_lost > 0) {
        logger.info(sprintf("%d individuals with NA attachment types removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
      }
      
    } else {
      logger.info("... with summary table")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(yearly_survival$attachment_type))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # Clean data, remove NAs 
      n_original      <- nrow(yearly_survival)
      yearly_survival <- yearly_survival %>%
        filter(!is.na(attachment_type)) %>%
        mutate(attachment_type = str_trim(attachment_type),
               attachment_type = str_replace_all(attachment_type, "\\s+", ""),
               attachment_type = str_extract(attachment_type, "^[^|]+"),
               attachment_type = str_replace(attachment_type, "–", "-"))
      
      # Get unique conditions after cleaning
      unique_conditions <- sort(unique(yearly_survival$attachment_type))
      n_conditions      <- length(unique_conditions)
      
      logger.info(sprintf("%d attachment types detected after cleaning: %s", n_conditions, 
                          paste(unique_conditions, collapse = ", ")),
                  call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(yearly_survival)
      if (n_lost > 0) {
        logger.info(sprintf("%d individuals with NA attachment type removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
      }
    } 
  }
  
  if(!is.null(group_comparison_individual) && group_comparison_individual == "survival_year"){
    logger.info("Grouping by survival year...")
    
    if(is.null(survival_yr_start)){ 
      logger.error("This comparison only makes sense if survival years are defined.")
      
    } else {
      logger.info("... with yearly survival")
      
      # Ensure different conditions are present
      non_na_unique <- unique(na.omit(yearly_survival$survival_year))
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.warn("Warning: The grouping variable is entirely NA; no comparison is possible.")
        } else {
          logger.warn("Warning: There is only one non-NA grouping variable; no comparison is possible.")
        }
      }
      
      # Remove NAs 
      n_original      <- nrow(yearly_survival)
      yearly_survival <- yearly_survival[!is.na(yearly_survival$survival_year),] 
      
      # Get unique life-stages after cleaning
      unique_surv_yrs  <- sort(unique(yearly_survival$survival_year))
      n_surv_yrs       <- length(unique_surv_yrs)
      
      logger.info(sprintf("%d survival years detected after cleaning: %s", n_surv_yrs, 
                          paste(unique_surv_yrs, collapse = ", ")),
                  call. = FALSE, immediate. = TRUE)
      
      n_lost <- n_original - nrow(yearly_survival)
      if (n_lost > 0) {
        logger.info(sprintf("%d individuals with NA survival year removed from study.", n_lost),
                    call. = FALSE, immediate. = TRUE)
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
    png(appArtifactPath("tracking_history.png"), width = 10, height = 8, units = "in", res = 300)
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
      png(appArtifactPath("monthly_mortality.png"), width = 10, height = 8, units = "in", res = 300)
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
    png(appArtifactPath("tracking_history.png"), width = 10, height = 8, units = "in", res = 300)
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
      png(appArtifactPath("monthly_mortality.png"), width = 10, height = 8, units = "in", res = 300)
      print(monthly_mort_plot)
      dev.off()
      
    } else {
      logger.info("Not plotting monthly mortality.")
    }
  } 
  
  
  ## Survival Analysis: No comparisons ----------------------------------------
  
  if (is.null(survival_yr_start)) {
    
    logger.info("Calculating KM using summary table...")
    
    # Warning for not mortality 
    if(sum(summary_table$mortality_event) > 0){
      logger.warn("There are no mortality events in the chosen subset of data.")
    }
    
    # Fit model 
    fitting_data <- summary_table
    km_fit <- survfit(Surv(entry_time_days, exit_time_days, mortality_event) ~ 1, 
                      data = summary_table)
    
    # Data for life table 
    lt.length.out <- ceiling(max(summary_table$exit_time_days)/life_table_days)
    times <- round(seq(min(summary_table$entry_time_days), max(summary_table$exit_time_days), 
                       length.out = lt.length.out))
    
  } 
  
  if(!is.null(survival_yr_start)) {
    
    logger.info("Calculating KM using survival table...")
    
    # Warning for no mortality 
    if(sum(yearly_survival$mortality_event) == 0){
      logger.warn("There are no mortality events in the chosen subset of data.")
    }
    
    # Fit model 
    fitting_data <- yearly_survival
    km_fit <- survfit(Surv(days_at_risk, mortality_event) ~ 1, data = yearly_survival)
    
    # Data for life table 
    max_days <- max(yearly_survival$days_at_risk, na.rm = TRUE)
    lt.length.out <- ceiling(max_days / life_table_days) + 1   
    times <- round(seq(0, max_days, length.out = lt.length.out))
  }
  
  # Generate life table ---
  s <- summary(km_fit, times = times)
  life_table <- data.frame(time_days        = s$time,
                           n_risk           = s$n.risk,
                           n_event          = s$n.event,
                           survival_prob    = s$surv,
                           std_err          = s$std.err,
                           lower_95         = s$lower,
                           upper_95         = s$upper)
  
  # Save 
  write.csv(life_table, file = appArtifactPath("life_table.csv"), row.names = F)
  
  # Plot KM curve ---
  n.ind <- nrow(fitting_data)
  n.events <- nrow(fitting_data[fitting_data$mortality_event == 1,])
  n.days <- as.numeric(summary(km_fit)$table["median"])
  med <- survminer::surv_median(km_fit)
  
  km_curve <- ggsurvplot(
    km_fit,
    data = fitting_data,
    title = "Kaplan-Meier Survival Curve",
    subtitle = paste0("N = ", n.ind, ", Events = ", n.events, ", Median Survival = ", 
                      med$median, " days"),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    risk.table = FALSE,
    conf.int = TRUE,
    censor.shape = "|",
    censor.size = 3,
    legend = "none",
    pval = FALSE,
    surv.median.line = "hv",        
    palette = c("#E69F00", "#56B4E9"),
    ggtheme = theme_classic(base_size = 12) + 
      theme(plot.title         = element_text(face = "bold", size = 14), 
            plot.subtitle      = element_text(size = 12, color = "gray50"),
            axis.text          = element_text(color = "black"),
            panel.grid.major.y = element_line(color = "gray90"), 
            panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
            line               = element_line(linewidth = 0.1),
            plot.margin        = margin(10, 10, 10, 10)))
  
  # Save plot  
  png(appArtifactPath("km_survival_curve.png"), width = 10, height = 8, units = "in", res = 300)
  print(km_curve)
  dev.off()
  
  
  # Plot cumulative hazard curve --- 
  cum_hazard <- ggsurvplot(
    km_fit,
    data = fitting_data,
    fun = "cumhaz",
    conf.int = TRUE,
    risk.table = FALSE,
    cumevents = FALSE,                 
    pval = FALSE,                     
    xlab = "Time (days)",
    ylab = "Cumulative Hazard",
    title = "Cumulative Hazard",
    subtitle = paste0("N = ", n.ind, ", Events = ", n.events), 
    palette = c("#E69F00", "#56B4E9"),
    legend = "none",
    ggtheme = theme_classic(base_size = 12) + 
      theme(plot.title   = element_text(face = "bold", size = 14), 
            plot.subtitle = element_text(size = 12, color = "gray50"),
            axis.text    = element_text(color = "black"),
            panel.grid.major.y = element_line(color = "gray90"), 
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            plot.margin  = margin(10, 10, 10, 10)))
  
  # Save plot 
  png(appArtifactPath("cumulative_hazard_plot.png"), width = 10, height = 8, units = "in", res = 300)
  print(cum_hazard)
  dev.off()
  
  
  ## Survival Analysis: Group comparisons -------------------------------------
  
  # Dynamically update comparison variables 
  grouping_labels <- c("model"                         = "Tag Model",
                       "sex"                           = "Sex",
                       "attachment"                    = "Tag Attachment", 
                       "survival_year"                 = "Survival Year",
                       "animal_life_stage_new"         = "Lifestage")
  
  # For yearly_summary data 
  if (!is.null(group_comparison_individual) && !is.null(survival_yr_start)){
    logger.info(paste0("Comparing ", group_comparison_individual, " using yearly survival data"))
    
    # Update headings if necessary
    if(group_comparison_individual == "lifestage"){group_comparison_individual <- "animal_life_stage_new"}
    
    # Fit survival object 
    surv_formula <- as.formula(paste("Surv(days_at_risk, mortality_event) ~", group_comparison_individual))
    km_fit_comp <- survfit(surv_formula, data = yearly_survival) 
    
    # Log-Rank test --- 
    test <- survdiff(surv_formula, data=yearly_survival)
    
    # Extract components
    groups       <- names(test$n)
    n_total      <- sum(test$n)
    events_total <- sum(test$obs)
    chisq_val    <- round(test$chisq, 2)
    df_val       <- length(test$n) - 1
    p_val        <- 1 - pchisq(test$chisq, df_val)
    p_formatted  <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
    
    grouping_var <- grouping_labels[[group_comparison_individual]]
    if (is.na(grouping_var)) {
      logger.error("Unknown grouping variable: ", group_comparison_individual)
    }
    
    # Summary table 
    per_group <- tibble(!!grouping_var   := sub(".*=", "", groups),
                        `N`                      = test$n,
                        `Events`                 = test$obs,
                        `Expected events`        = round(test$exp, 2),
                        `O/E ratio`              = round(test$obs / test$exp, 2)) %>%
      mutate(`N (events)` = sprintf("%d (%d)", N, Events), .keep = "unused")
    
    summary_row <- tibble(!!grouping_var          := "Overall",
                          `N (events)`             = sprintf("%d (%d)", n_total, events_total),
                          `Chisq (log-rank)`       = chisq_val,
                          `df`                     = df_val,
                          `p-value`                = p_formatted)
    
    logrank_table <- bind_rows(per_group, summary_row)
    logrank_table <- logrank_table %>%
      dplyr::select(any_of(grouping_var), `N (events)`, `Expected events`, `O/E ratio`,
                    `Chisq (log-rank)`, df, `p-value`)
    
    # Save
    write.csv(logrank_table, file = appArtifactPath("logrank_table_statistics.csv"), row.names = F)
    
    
    # KM comparison plots ---
    km_fit_comp <- surv_fit(surv_formula, data = yearly_survival)  
    
    # Check if any groups have N == 1 and remove 
    filter_singleton_strata_refit <- function(fit, data) {
      
      grouping_var <- all.vars(formula(fit))[3]    
      n_per_group <- fit$n
      group_names <- sub(".*=", "", names(fit$strata) %||% "Overall")
      
      singletons <- n_per_group == 1
      removed <- group_names[singletons]   
      
      if (any(singletons)) {
        logger.info(paste("Removed the following singleton group(s) (N=1):\n",
                          paste(" •", removed, collapse = "\n ")), call. = FALSE)
        
        # Filter data using clean group names
        data_clean <- data[!data[[grouping_var]] %in% removed, , drop = FALSE]
        
        # Refit 
        km_clean <- surv_fit(formula(fit), data = data_clean)
        return(km_clean)
      }
      return(fit)
    }
    
    km_fit_clean <- filter_singleton_strata_refit(km_fit_comp, yearly_survival)
    
    title_text <- paste0("Kaplan-Meier Survival Curve: ", grouping_var)
    group_rows <- logrank_table %>% filter(!!sym(grouping_var) != "Overall") 
    group_rows_clean <- group_rows %>% filter(!str_detect(`N (events)`, "^1\\s*\\("))
    removed_groups <- group_rows %>% filter(str_detect(`N (events)`, "^1\\s*\\(")) %>%
      pull(1)
    
    if (length(removed_groups) > 0) {
      logger.info(paste0("The following group(s) had N=1 and were removed from the table:\n • ",
                         paste(removed_groups, collapse = "\n • ")), call. = FALSE)
    }
    
    if (nrow(group_rows_clean == 1)){
      logger.fatal("There is only one group left; unable to perform comparisons.")
      
    } else {
      subtitle_text <- paste0(paste(sprintf("N_%s: %s", 
                                            group_rows_clean[[grouping_var]],
                                            group_rows_clean$`N (events)`),
                                    collapse = ", "),
                              "\nP-value: ", 
                              logrank_table$`p-value`[logrank_table[[grouping_var]] == "Overall"])
      
      old_strata_names <- names(km_fit_clean$strata)
      new_strata_names <- sub(".*=", "", old_strata_names)
      names(km_fit_clean$strata) <- new_strata_names
      
      n_groups <- length(names(km_fit_clean$strata))
      my_palette <- viridis(n_groups, option = "turbo")
      
      # Plot 
      km_comp_curve <- ggsurvplot(km_fit_clean,
                                  data = yearly_survival,
                                  title = title_text,
                                  subtitle = subtitle_text,
                                  conf.int = TRUE,
                                  risk.table = FALSE,
                                  palette = my_palette, 
                                  xlab = "Days at risk",
                                  ylab = "Survival probability",
                                  legend.title = grouping_var,
                                  legend = "bottom",
                                  legend.labs = levels(summary_table[[group_comparison_individual]]),
                                  censor.shape = "|",
                                  censor.size = 4,
                                  font.main = c(14, "bold", "black"),
                                  font.x = 12, font.y = 12, font.tickslab = 11, 
                                  ggtheme = theme_classic(base_size = 12) +
                                    theme(plot.title   = element_text(face = "bold", size = 14),
                                          plot.subtitle = element_text(size = 12, color = "gray50"),
                                          axis.text    = element_text(color = "black"),
                                          panel.grid.major.y = element_line(color = "gray90"),
                                          panel.border = element_rect(color="black", fill=NA, linewidth=0.5),
                                          plot.margin  = margin(10, 10, 10, 10)))
      
      km_comp_curve$plot <- km_comp_curve$plot + 
        guides(color = guide_legend(nrow = ifelse(n_groups <= 4, 1, 2), byrow = TRUE,
                                    title.position = "top"))
      
      # Save plot 
      png(appArtifactPath("km_comparison_curves.png"), width = 10, height = 8, units = "in", res = 300)
      print(km_comp_curve)
      dev.off()
      
      
      ## Cumulative hazard comparison plots ---
      
      # Prepare statistics for subtitle 
      n_per_group      <- km_fit_clean$n
      sum_fit          <- summary(km_fit_clean)
      events_per_group <- tapply(sum_fit$n.event, sum_fit$strata, sum, na.rm = TRUE)
      
      # Clean strata names  
      clean_strata <- gsub("^.*=", "", names(km_fit_clean$strata))
      
      # Create one-line subtitle for groups: "N(Group): N (events)"
      subtitle_parts <- mapply(function(group, n, ev) {
        sprintf("N_%s: %d(%d)", group, n, ev)
      },
      clean_strata, n_per_group, events_per_group, SIMPLIFY = FALSE)
      
      groups_line <- paste(subtitle_parts, collapse = ", ")
      test <- surv_pvalue(km_fit_clean, data = yearly_survival)
      pval_text <- sprintf("P-value: %.3f", test$pval)
      subtitle_text <- paste0(groups_line, "\n", pval_text)
      
      # Plot 
      cum_hazard_comp <- ggsurvplot(km_fit_clean,
                                    data         = yearly_survival,
                                    fun          = "cumhaz",
                                    conf.int     = TRUE,
                                    censor.shape = "|",
                                    censor.size  = 4,
                                    title        = paste0("Cumulative Hazard by ", grouping_var),
                                    subtitle     = subtitle_text,   
                                    xlab         = "Days at risk",
                                    ylab         = "Cumulative Hazard",
                                    legend       = "bottom",
                                    legend.title = grouping_var,
                                    legend.labs  = levels(summary_table[[group_comparison_individual]]),
                                    palette      = my_palette,
                                    risk.table   = FALSE,
                                    cumevents    = FALSE,
                                    font.main    = c(14, "bold", "black"),
                                    font.x       = 12, font.y = 12, font.tickslab = 11,
                                    ggtheme = theme_classic(base_size = 12) +
                                      theme(plot.title    = element_text(face = "bold", size = 14),
                                            plot.subtitle = element_text(size = 12, color = "gray50"),
                                            axis.text     = element_text(color = "black"),
                                            panel.grid.major.y = element_line(color = "gray90"),
                                            panel.border  = element_rect(color = "black", fill=NA, linewidth=0.5),
                                            plot.margin   = margin(10, 10, 10, 10),
                                            legend.position = "bottom",
                                            legend.title  = element_text(size = 11),
                                            legend.text   = element_text(size = 10)))
      
      cum_hazard_comp$plot <- cum_hazard_comp$plot + 
        guides(color = guide_legend(nrow = ifelse(n_groups <= 4, 1, 2), byrow = TRUE,
                                    title.position = "top"))
      
      # Save plot  
      png(appArtifactPath("cum_hazard_comparison_plot.png"), width = 10, height = 8, units = "in", res = 300)
      print(cum_hazard_comp)
      dev.off()
    }
  } 
  
  # For summary table data 
  if(!is.null(group_comparison_individual) && is.null(survival_yr_start)) {
    logger.info(paste0("Comparing ", group_comparison_individual, " using summary table data"))
    
    # Fit survival object ---
    summary_table$time_at_risk <- summary_table$exit_time_days - summary_table$entry_time_days
    formula_str <- paste("Surv(time_at_risk, mortality_event) ~", group_comparison_individual)
    surv_formula <- as.formula(formula_str)
    km_fit_comp <- survfit(surv_formula, data = summary_table)
    
    
    # Log-Rank test --- 
    test <- survdiff(surv_formula, data=summary_table)
    
    # Extract components
    groups <- names(test$n)
    n_total <- sum(test$n)
    events_total <- sum(test$obs)
    chisq_val <- round(test$chisq, 2)
    df_val <- length(test$n) - 1
    p_val <- 1 - pchisq(test$chisq, df_val)
    p_formatted <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
    
    # Dynamically update comparison variables 
    grouping_labels <- c("model"                         = "Tag model",
                         "sex"                           = "Sex",
                         "attachment"                    = "Tag Attachment")
    grouping_var <- grouping_labels[[group_comparison_individual]]
    if (is.na(grouping_var)) {
      logger.error("Unknown grouping variable: ", group_comparison_individual)
    }
    
    # Summary table 
    per_group <- tibble(!!grouping_var   := sub(".*=", "", groups),
                        `N`               = test$n,
                        `Events`          = test$obs,
                        `Expected events` = round(test$exp, 2),
                        `O/E ratio`       = round(test$obs / test$exp, 2)) %>%
      mutate(`N (events)` = sprintf("%d (%d)", N, Events), .keep = "unused")
    
    
    summary_row <- tibble(!!grouping_var          := "Overall",
                          `N (events)`             = sprintf("%d (%d)", n_total, events_total),
                          `Chisq (log-rank)`       = chisq_val,
                          `df`                     = df_val,
                          `p-value`                = p_formatted)
    
    logrank_table <- bind_rows(per_group, summary_row)
    logrank_table <- logrank_table %>%
      dplyr::select(any_of(grouping_var), 
                    `N (events)`,
                    `Expected events`,
                    `O/E ratio`,
                    `Chisq (log-rank)`,
                    df,
                    `p-value`)
    
    # Save
    write.csv(logrank_table, file = appArtifactPath("logrank_table_statistics.csv"), row.names = F)
    
    
    # KM comparison plots ---
    km_fit_comp <- surv_fit(surv_formula, data = summary_table)   
    
    # Check if any groups have N == 1 and remove 
    filter_singleton_strata_refit <- function(fit, data) {
      grouping_var <- all.vars(formula(fit))[3]   
      n_per_group <- fit$n
      group_names <- sub(".*=", "", names(fit$strata) %||% "Overall")
      
      # Identify groups where N == 1 
      singletons <- n_per_group == 1
      removed <- group_names[singletons]    
      
      if (any(singletons)) {
        logger.info(paste("Removed the following singleton group(s) (N=1):\n",
                          paste(" •", removed, collapse = "\n ")), call. = FALSE)
        
        # Filter data using clean group names
        data_clean <- data[!data[[grouping_var]] %in% removed, , drop = FALSE]
        
        # Refit 
        km_clean <- surv_fit(formula(fit), data = data_clean)
        
        return(km_clean)
      }
      return(fit)
    }
    
    km_fit_clean <- filter_singleton_strata_refit(km_fit_comp, summary_table)
    
    # Titles and formatting 
    title_text <- paste0("Kaplan-Meier Survival Curve: ", grouping_var)
    
    group_rows <- logrank_table %>% 
      filter(!!sym(grouping_var) != "Overall") 
    
    # Remove any rows where N = 1 
    group_rows_clean <- group_rows %>%
      filter(!str_detect(`N (events)`, "^1\\s*\\("))
    
    removed_groups <- group_rows %>%
      filter(str_detect(`N (events)`, "^1\\s*\\(")) %>%
      pull(1)
    
    if (length(removed_groups) > 0) {
      logger.info(paste0("The following group(s) had N=1 and were removed from the table:\n • ",
                         paste(removed_groups, collapse = "\n • ")), call. = FALSE)
    }
    
    if (nrow(group_rows_clean == 1)){
      logger.fatal("There is only one group left; unable to perform comparisons.")
      
    } else { 
      subtitle_text <- paste0(paste(sprintf("N_%s: %s", 
                                            group_rows_clean[[grouping_var]],
                                            group_rows_clean$`N (events)`),
                                    collapse = ", "),
                              "\nP-value: ", 
                              logrank_table$`p-value`[logrank_table[[grouping_var]] == "Overall"])
      
      old_strata_names <- names(km_fit_clean$strata)
      new_strata_names <- sub(".*=", "", old_strata_names)
      names(km_fit_clean$strata) <- new_strata_names
      
      n_groups <- length(names(km_fit_clean$strata))
      my_palette <- viridis(n_groups, option = "turbo")
      
      # Plot 
      km_comp_curve <- ggsurvplot(km_fit_clean,
                                  data = summary_table,
                                  title = title_text,
                                  subtitle = subtitle_text,
                                  conf.int = TRUE,
                                  risk.table = FALSE,
                                  palette = my_palette, 
                                  xlab = "Days at risk",
                                  ylab = "Survival probability",
                                  legend.title = grouping_var,
                                  legend = "bottom",
                                  legend.labs = levels(summary_table[[group_comparison_individual]]),
                                  censor.shape = "|",
                                  censor.size = 4,
                                  font.main = c(14, "bold", "black"),
                                  font.x = 12, font.y = 12, font.tickslab = 11, 
                                  ggtheme = theme_classic(base_size = 12) +
                                    theme(
                                      plot.title   = element_text(face = "bold", size = 14),
                                      plot.subtitle = element_text(size = 12, color = "gray50"),
                                      axis.text    = element_text(color = "black"),
                                      panel.grid.major.y = element_line(color = "gray90"),
                                      panel.border = element_rect(color="black", fill=NA, linewidth=0.5),
                                      plot.margin  = margin(10, 10, 10, 10)))
      
      km_comp_curve$plot <- km_comp_curve$plot + 
        guides(color = guide_legend(nrow = ifelse(n_groups <= 4, 1, 2), byrow = TRUE,
                                    title.position = "top"))
      
      # Save plot 
      png(appArtifactPath("km_comparison_curves.png"), width = 10, height = 8, units = "in", res = 300)
      print(km_comp_curve)
      dev.off()
      
      
      ## Cumulative hazard comparison plots ---
      
      # Prepare statistics for subtitle 
      n_per_group      <- km_fit_clean$n
      sum_fit          <- summary(km_fit_clean)
      events_per_group <- tapply(sum_fit$n.event, sum_fit$strata, sum, na.rm = TRUE)
      
      # Clean strata names  
      clean_strata <- gsub("^.*=", "", names(km_fit_clean$strata))
      
      # Create one-line subtitle for groups: "N(Group): N (events)"
      subtitle_parts <- mapply(function(group, n, ev) {
        sprintf("N_%s: %d(%d)", group, n, ev)
      },
      clean_strata, n_per_group, events_per_group, SIMPLIFY = FALSE)
      
      groups_line <- paste(subtitle_parts, collapse = ", ")
      test <- surv_pvalue(km_fit_clean, data = summary_table)
      pval_text <- sprintf("P-value: %.3f", test$pval)
      subtitle_text <- paste0(groups_line, "\n", pval_text)
      
      # Plot 
      cum_hazard_comp <- ggsurvplot(km_fit_clean,
                                    data = summary_table,
                                    fun          = "cumhaz",
                                    conf.int     = TRUE,
                                    censor.shape = "|",
                                    censor.size  = 4,
                                    title        = paste0("Cumulative Hazard by ", grouping_var),
                                    subtitle     = subtitle_text,   
                                    xlab         = "Days at risk",
                                    ylab         = "Cumulative Hazard",
                                    legend       = "bottom",
                                    legend.title = grouping_var,
                                    legend.labs  = levels(summary_table[[group_comparison_individual]]),
                                    palette      = my_palette,
                                    risk.table   = FALSE,
                                    cumevents    = FALSE,
                                    font.main    = c(14, "bold", "black"),
                                    font.x       = 12, font.y = 12, font.tickslab = 11,
                                    ggtheme = theme_classic(base_size = 12) +
                                      theme(plot.title    = element_text(face = "bold", size = 14),
                                            plot.subtitle = element_text(size = 12, color = "gray50"),
                                            axis.text     = element_text(color = "black"),
                                            panel.grid.major.y = element_line(color = "gray90"),
                                            panel.border  = element_rect(color = "black", fill=NA, linewidth=0.5),
                                            plot.margin   = margin(10, 10, 10, 10),
                                            legend.position = "bottom",
                                            legend.title  = element_text(size = 11),
                                            legend.text   = element_text(size = 10)))
      
      cum_hazard_comp$plot <- cum_hazard_comp$plot + 
        guides(color = guide_legend(nrow = ifelse(n_groups <= 4, 1, 2), byrow = TRUE,
                                    title.position = "top"))
      
      # Save plot 
      png(appArtifactPath("cum_hazard_comparison_plot.png"), width = 10, height = 8, units = "in", res = 300)
      print(cum_hazard_comp)
      dev.off()
    }
  }
  
  # Pass original to the next app in the MoveApps workflow
  return(data)
} 