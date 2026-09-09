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
library(tidyr)

# Set of working directory
setwd("C:/Users/Jiri/Desktop/Climate_signal_change/")

Loc_list<-readRDS("Data/Loc_list2.rds")

cor_change_categories <- function(early, late, threshold = 0.27) {
  
  early_mat <- as.matrix(early)
  late_mat  <- as.matrix(late)
  
  # Difference: late - early
  change <- late_mat - early_mat
  
  # Direction of change
  category <- matrix(
    NA_character_,
    nrow = nrow(early_mat),
    ncol = ncol(early_mat)
  )
  
  # Strong changes
  category[early_mat > 0 & late_mat > 0 & abs(change) > threshold] <- "PP"
  category[early_mat < 0 & late_mat < 0 & abs(change) > threshold] <- "NN"
  category[early_mat > 0 & late_mat < 0 & abs(change) > threshold] <- "PN"
  category[early_mat < 0 & late_mat > 0 & abs(change) > threshold] <- "NP"
  
  # Weak changes
  category[!is.na(change) & abs(change) <= threshold] <- "WEAK"
  
  # Results
  result <- tibble(
    DOY = seq_len(nrow(early_mat)),
    
    EARLY = rowMeans(early_mat, na.rm = TRUE),
    LATE  = rowMeans(late_mat, na.rm = TRUE),
    CHANGE = LATE - EARLY,
    
    PP = rowSums(category == "PP", na.rm = TRUE),
    NN = rowSums(category == "NN", na.rm = TRUE),
    PN = rowSums(category == "PN", na.rm = TRUE),
    NP = rowSums(category == "NP", na.rm = TRUE),
    WEAK = rowSums(category == "WEAK", na.rm = TRUE)
    )
  
  return(result)
}

# Empty data-sets for tree-ring data and daily correlations
Daily_cor<- data.frame(SPECIES2=NA, SITE=NA, CLIM = NA, DOY=NA, EARLY = NA, LATE = NA, CHANGE = NA, 
                       PP = NA, NN = NA, PN = NA, NP = NA, WEAK = NA)

for (i in c(1:nrow(Loc_list))) {
  
  print(paste("Site = ", i, "Done =", i/950*100))
  
  
  #=============================#
  #         Temperature         #
  #=============================#
  
  TEMP_EARLY<- readRDS(paste("C:/Users/Jiri/Desktop/Climate_signal_change/Data/Interim_data/Site_cor/", 
                             Loc_list[i, "SITE"], "_TEMP_EARLY.rds", sep = ""))

  TEMP_LATE<- readRDS(paste("C:/Users/Jiri/Desktop/Climate_signal_change/Data/Interim_data/Site_cor/", 
                             Loc_list[i, "SITE"], "_TEMP_LATE.rds", sep = ""))
  
  TEMP_change <- cor_change_categories(
    TEMP_EARLY,
    TEMP_LATE,
    threshold = 0.27
  )
  
  TEMP_change$SITE<- Loc_list[i, "SITE"]
  TEMP_change$SPECIES2<- Loc_list[i, "SPECIES2"]
  TEMP_change$CLIM<- "TEMP"
  
  #=============================#
  #         Precipitation       #
  #=============================#
  
  PREC_EARLY<- readRDS(paste("C:/Users/Jiri/Desktop/Climate_signal_change/Data/Interim_data/Site_cor/", 
                             Loc_list[i, "SITE"], "_PREC_EARLY.rds", sep = ""))

  PREC_LATE<- readRDS(paste("C:/Users/Jiri/Desktop/Climate_signal_change/Data/Interim_data/Site_cor/", 
                            Loc_list[i, "SITE"], "_PREC_LATE.rds", sep = ""))
  
  PREC_change <- cor_change_categories(
    PREC_EARLY,
    PREC_LATE,
    threshold = 0.27
  )
  
  PREC_change$SITE<- Loc_list[i, "SITE"]
  PREC_change$SPECIES2<- Loc_list[i, "SPECIES2"]
  PREC_change$CLIM<- "PREC"
  
  Daily_cor<- rbind(Daily_cor, TEMP_change, PREC_change)
  
}

Daily_cor<- Daily_cor[-1,]
Daily_cor<- Daily_cor |> filter(DOY>=121 & DOY<=639)

saveRDS(Daily_cor, "Data/Interim_data/Daily_cor_smart.rds")
