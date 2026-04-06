#===============================================================================================#
#===============================================================================================#
#   Catalogue of climatic sensitivity patterns for temperate tree species in Central Europe  ####
#                                     Analysis part of script
#===============================================================================================#
#===============================================================================================#

# Necessary packages
library(terra) # raster work
library(mapview) # map visualization
library(dplR) # tree ring processing
library(lubridate) # package for DOY adding
library(dendroTools) # package for tree-ring daily correlations
library(reshape2) # data structure changing
library(dplyr) # selecting
library(tidyverse)

# Set of working directory
setwd("C:/Users/Jiri/Desktop/Climate_signal_change/")

#==================================================================================#
####                  Loading necessary data                   ####
#==================================================================================#

# Loading of list of sites (New list with 950 sites is filtered, so each site have 70 or more years in 2010)
Loc_list<-readRDS("Data/Loc_list2.rds") |> 
  vect(geom = c("LON", "LAT"), crs="EPSG:4326")

# Loading of daily EOBS climatic grids
EOBS_daily_T_AOI <- rast("Data/Clima_grids/EOBS_daily_T_AOI.tif")
EOBS_daily_P_AOI <- rast("Data/Clima_grids/EOBS_daily_P_AOI.tif")

# Selection of layers since 1.1.1961 to 31.12.2010
start_date <- as.Date("1950-01-01")
end_date   <- as.Date("2010-12-31")
EOBS_daily_T_AOI <- EOBS_daily_T_AOI[[which(time(EOBS_daily_T_AOI) >= start_date & time(EOBS_daily_T_AOI) <= end_date)]]
EOBS_daily_P_AOI <- EOBS_daily_P_AOI[[which(time(EOBS_daily_P_AOI) >= start_date & time(EOBS_daily_P_AOI) <= end_date)]]

#==================================================================================#
####                  1. Calculation of daily climatic signal                   ####
#==================================================================================#

# Loading of PER and DOY selection functions
source("Skripts/Accompanying/Needed_functions.R")

# Empty data-sets for tree-ring data and daily correlations
Daily_CAT<- data.frame(SPECIES = NA, SITE = NA, CLIM = NA, DOY = NA, PP = NA, NN = NA, PN = NA, NP = NA)
Daily_cor<- data.frame(SPECIES=NA, SITE=NA, DOY=NA, TEMP_EARLY = NA, TEMP_LATE = NA, PREC_EARLY = NA, PREC_LATE = NA)

# Trick for switching among folders
Databases<- unique(Loc_list$DATABASE)

start<- now()

for (folder in c(1:length(Databases))) {
  
  Database<- Databases[folder]
  
  Set<- Loc_list[Loc_list$DATABASE == Database,]
  
  for (i in 1:c(nrow(Set))) {
    
    print(paste("Database =", Database, "Site =",i, "Done =", (values(Set[i,"NUM"]) |>  as.numeric()/nrow(Loc_list))*100, "%"))
    
    #===================================================================#
    #####     Loading of tree-ring data and chronology building     #####
    #===================================================================#
    
    serie<- read.rwl(paste("C:/Users/Jiri/Desktop/Climate_signal_change/Data/TRW_data/", Database, "/", unlist(values(Set[i,"RWL"])), sep = ""), format = "tucson")
    
    chronology <- detrend(serie, method = "Spline", nyrs = 30) |> chron(biweight = T, prewhiten = T)
    chronology<- subset(chronology, subset=rownames(chronology)<=2010 & rownames(chronology)>=1950)
    chronology$std<- NULL
    
    #===================================================================#
    #####        Extraction of daily climatic data from grids       #####
    #                       Detrend of daily data                       #
    #===================================================================#
    
    Temperature_daily<- terra::extract(EOBS_daily_T_AOI, Set[i,], ID=FALSE) |> t() |> as.data.frame()
    T_det<- Clima_detrend(Temperature_daily, "temp")
    
    Precipitation_daily<- terra::extract(EOBS_daily_P_AOI, Set[i,], ID=FALSE) |> t() |> as.data.frame()
    P_det<- Clima_detrend(Precipitation_daily, "prec")
    
    #===================================================================#
    #####              Calculation of daily correlations           #####
    #===================================================================#
    
    # For temperature
    TEMP_EARLY <- daily_response(chronology, env_data = T_det, previous_year = T, 
                               subset_years = c(1951, 1980), tidy_env_data = T, method = "cor", 
                               aggregate_function = "mean", reference_window = "middle", 
                               remove_insignificant = F, lower_limit = 10, upper_limit = 90,
                               row_names_subset = TRUE) |> 
                  pluck("calculations") |> 
                  melt() |> 
                  pivot_wider(names_from = Var1, values_from = value) |> 
                  select(-Var2)
      

    TEMP_LATE <- daily_response(chronology, env_data = T_det, previous_year = T, 
                                  subset_years = c(1981, 2010), tidy_env_data = T, method = "cor", 
                                  aggregate_function = "mean", reference_window = "middle", 
                                  remove_insignificant = F, lower_limit = 10, upper_limit = 90,
                                  row_names_subset = TRUE) |> 
                  pluck("calculations") |> 
                  melt() |> 
                  pivot_wider(names_from = Var1, values_from = value) |> 
                  select(-Var2)
      
    
      TEMP_CLASS <- map2_df(TEMP_EARLY, TEMP_LATE, ~{
        case_when(
          is.na(.x) | is.na(.y) ~ NA_character_,
          .x > 0 & .y > 0       ~ "PP",
          .x < 0 & .y < 0       ~ "NN",
          .x > 0 & .y < 0       ~ "PN",
          .x < 0 & .y > 0       ~ "NP",
          TRUE                 ~ NA_character_  # covers zeros
        )
      })
      
      
      cols <- colnames(TEMP_CLASS)
      
      TEMP_percent <- TEMP_CLASS |>
        rowwise() |>
        mutate(
          n_valid = sum(!is.na(c_across(all_of(cols)))),
          
          PP = 100 * sum(c_across(all_of(cols)) == "PP", na.rm = TRUE) / n_valid,
          NN = 100 * sum(c_across(all_of(cols)) == "NN", na.rm = TRUE) / n_valid,
          PN = 100 * sum(c_across(all_of(cols)) == "PN", na.rm = TRUE) / n_valid,
          NP = 100 * sum(c_across(all_of(cols)) == "NP", na.rm = TRUE) / n_valid
        ) |>
        ungroup() |>
        select(PP, NN, PN, NP) |> 
        add_column(CLIM = "TEMP", DOY = c(1:726))
    
      
      # For Precipitation
      PREC_EARLY <- daily_response(chronology, env_data = P_det, previous_year = T, 
                                   subset_years = c(1951, 1980), tidy_env_data = T, method = "cor", 
                                   aggregate_function = "mean", reference_window = "middle", 
                                   remove_insignificant = F, lower_limit = 10, upper_limit = 90,
                                   row_names_subset = TRUE) |> 
        pluck("calculations") |> 
        melt() |> 
        pivot_wider(names_from = Var1, values_from = value) |> 
        select(-Var2)
      
      
      
      PREC_LATE <- daily_response(chronology, env_data = P_det, previous_year = T, 
                                  subset_years = c(1981, 2010), tidy_env_data = T, method = "cor", 
                                  aggregate_function = "mean", reference_window = "middle", 
                                  remove_insignificant = F, lower_limit = 10, upper_limit = 90,
                                  row_names_subset = TRUE) |> 
        pluck("calculations") |> 
        melt() |> 
        pivot_wider(names_from = Var1, values_from = value) |> 
        select(-Var2)
      
      
      PREC_CLASS <- map2_df(PREC_EARLY, PREC_LATE, ~{
        case_when(
          is.na(.x) | is.na(.y) ~ NA_character_,
          .x > 0 & .y > 0       ~ "PP",
          .x < 0 & .y < 0       ~ "NN",
          .x > 0 & .y < 0       ~ "PN",
          .x < 0 & .y > 0       ~ "NP",
          TRUE                 ~ NA_character_  # covers zeros
        )
      })
      
      
      cols <- colnames(PREC_CLASS)
      
      PREC_percent <- PREC_CLASS |>
        rowwise() |>
        mutate(
          n_valid = sum(!is.na(c_across(all_of(cols)))),
          
          PP = 100 * sum(c_across(all_of(cols)) == "PP", na.rm = TRUE) / n_valid,
          NN = 100 * sum(c_across(all_of(cols)) == "NN", na.rm = TRUE) / n_valid,
          PN = 100 * sum(c_across(all_of(cols)) == "PN", na.rm = TRUE) / n_valid,
          NP = 100 * sum(c_across(all_of(cols)) == "NP", na.rm = TRUE) / n_valid
        ) |>
        ungroup() |>
        select(PP, NN, PN, NP) |> 
        add_column(CLIM = "PREC", DOY = c(1:726))
      
      
      # Combination together
      COR_clim<- rbind(TEMP_percent, PREC_percent)
      COR_clim$SITE<- unlist(values(Set[i, "SITE"]))
      COR_clim$SPECIES<- unlist(values(Set[i, "SPECIES2"]))

      # To combine all sites together
      Daily_CAT<- rbind(Daily_CAT, COR_clim)
      
      
      Daily_means<- data.frame(
        SITE = unlist(values(Set[i, "SITE"])),
        SPECIES = unlist(values(Set[i, "SPECIES2"])),
        DOY = c(1:726),
        TEMP_EARLY = rowMeans(TEMP_EARLY, na.rm = T),
        TEMP_LATE = rowMeans(TEMP_LATE, na.rm = T),
        PREC_EARLY = rowMeans(PREC_EARLY, na.rm = T),
        PREC_LATE = rowMeans(PREC_LATE, na.rm = T))
      
      Daily_cor<- rbind(Daily_cor, Daily_means)
  }
}

Daily_CAT<- Daily_CAT[-1,]
Daily_cor<- Daily_cor[-1,]

end<- now()

saveRDS(Daily_CAT, "Data/Temporary_data/Daily_CAT.rds")
saveRDS(Daily_cor, "Data/Temporary_data/Daily_cor.rds")
