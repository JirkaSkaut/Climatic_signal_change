
# Data operations
library(purrr)
library(tidyverse)
library(terra)
library(mgcv)
library(spdep)
library(car)
library(spatialreg)
library(Hmisc)

# plotting of figures
library(ggplot2)
library(ggh4x)
library(ggtext)
library(ggnewscale)
library(ggpubr)
library(viridis)
library(rworldmap)
library(tidyterra)

setwd("C:/Users/Jiri/Desktop/Climate_signal_change/")

#==============================================#
###      1. Explaining factors          ####
#==============================================#

# Setting area of interest
AOI <- matrix(
  c(6, 45,
    6, 55,
    27, 55,
    27, 45,
    6, 45), ncol = 2, byrow = TRUE) |>
  vect(type = "polygons", crs = "EPSG:4326")

# Loading main list of sites
Loc_list<-readRDS("Data/Loc_list.rds")

# Making spatial layer
Loc_sp<-readRDS("Data/Loc_list.rds") |> 
  vect(geom = c("LON", "LAT"), crs="EPSG:4326")

plet(AOI) |> points(Loc_sp)

#============#
#### CWB ####
#===========#

# (Climatic water balance) raster
CWB_windows<- rast("Data/Clima_grids/CWB_time_windows.tif")
Loc_list$CWB<- terra::extract(mean(CWB_windows), Loc_sp, ID = F)[,1]

#================#
#### d_CWB ####
#===============#

# Climatic water balance difference weighted by Annual precipitation
PREC<- rast("Data/Clima_grids/EOBS_monthly_sum_P.tif") |> 
  crop(AOI) |> 
  tapp("years", "sum") |> 
  mean()

d_CWB<- (CWB_windows$`1980–2010` - CWB_windows$`1950–1980`)/PREC
Loc_list$d_CWB<- terra::extract(d_CWB, Loc_sp, ID = FALSE)[,1]

#============#
#### CCI ####
#===========#

# Conrad's continentality index
CCI<- rast("Data/Clima_grids/EOBS_yearly_CCI.tif") |> 
  mean()
Loc_list$CCI<- terra::extract(CCI, Loc_sp, ID = FALSE)[,1]

#=================================================#
###      2. Loading daily calculated data      ####
#================================================#

Loc_list<-readRDS("Data/Loc_list.rds")

Daily_cor<- readRDS("Data/Interim_data/Daily_cor_smart.rds") |> 
  left_join(Loc_list[,c("SITE","CWB", "d_CWB", "AGE", "CEC", "CCI", "LAT", "LON")], by = "SITE")

species_colors <- c(
  "ABAL" = "#A6CEE3",
  "PCAB" = "#377EB8",
  "PISY" = "#FFEA61",
  "FASY" = "#4DAF4A", 
  "QUsp" = "#a7e065")

var_colors <- c(
  "AGE" = "#7fc97f",
  "d_CWB" = "#6699CC",
  "CWB" = "#f7b448",
  "CCI" = "#E73F74")

cat_colors <- c(
  "WEAK" = "#9b9d9e",
  "PP" = "#5ba3f5",
  "NN" = "#4788d1",
  "NP" = "#F0E442",
  "PN" = "#E69F00")

#==============================================#
###       3. Drivers percentage         ####
#==============================================#

CHANGE_SITE <- Daily_cor |>
  group_by(SPECIES2, SITE, CLIM) |>
  summarise(
    across(c(PP, NN, PN, NP, WEAK),
           \(x) mean(x, na.rm = TRUE)),
    .groups = "drop") |> 
  left_join(Loc_list[,c("SITE", "CWB", "d_CWB", "AGE", "CCI","LAT", "LON")], by = "SITE")

CHANGE_SITE$CODE <- paste(CHANGE_SITE$SPECIES2, CHANGE_SITE$CLIM, sep = "_")
CODE<- as.data.frame(unique(CHANGE_SITE$CODE)); colnames(CODE)<- "CODE"

VAR_models<- data.frame(VAR = NA, Model = NA, CODE = NA, Z_val = NA, P_val = NA)
P_models<- data.frame(CODE = NA, Model = NA, Wald = NA, W_p = NA, AIC = NA, AIC_lm = NA, 
                      pR2 = NA, RMSE = NA)

for (i in c(1:nrow(CODE))) {
  
  print(i)
  
  # Selection of species and climate var
  Set<- CHANGE_SITE[CHANGE_SITE$CODE == CODE[i, "CODE"],]
  
  # Calculating neighbours
  coords <- as.matrix(Set[, c("LON", "LAT")])
  knb <- knn2nb(knearneigh(coords, k = 3))
  lw <- nb2listw(knb, style = "W")
  
  PP_model <- errorsarlm(PP ~ CWB + d_CWB + AGE  + CCI, listw = lw, data = Set) |> 
    summary()
  
  PP_p<- data.frame(CODE = CODE[i, "CODE"], 
                    Model = "PP",
                    Wald = PP_model[["Wald1"]][["statistic"]][["Wald statistic"]], 
                    W_p = PP_model[["Wald1"]][["p.value"]][["Wald statistic"]], 
                    AIC = AIC(PP_model), 
                    AIC_lm = PP_model[["AIC_lm.model"]],
                    pR2 = cor(as.data.frame(predict(PP_model))[,1],
                                   Set$PP)^2,
                    RMSE = (sqrt(mean((Set$PP - as.data.frame(predict(PP_model))[,1])^2))) / sd(Set$PP)
                    )
  
  
  PP_var<- PP_model|>
    coef() |>
    as.data.frame() |>
    rownames_to_column("VAR") |>
    select(VAR, `z value`, `Pr(>|z|)`) |>
    filter(VAR != "(Intercept)") |>
    rename(
      Z_val = `z value`,
      P_val = `Pr(>|z|)`) |>
    mutate(
      Model = "PP",
      CODE = CODE[i, "CODE"])
  
  NN_model <- errorsarlm(NN ~ CWB + d_CWB + AGE + CCI, listw = lw, data = Set) |> 
    summary()
  
  NN_p<- data.frame(CODE = CODE[i, "CODE"], 
                    Model = "NN",
                    Wald = NN_model[["Wald1"]][["statistic"]][["Wald statistic"]], 
                    W_p = NN_model[["Wald1"]][["p.value"]][["Wald statistic"]], 
                    AIC = AIC(NN_model), 
                    AIC_lm = NN_model[["AIC_lm.model"]],
                    pR2 = cor(as.data.frame(predict(NN_model))[,1],
                                   Set$NN)^2,
                    RMSE = (sqrt(mean((Set$NN - as.data.frame(predict(NN_model))[,1])^2))) / sd(Set$NN)
                    )
  
  NN_var<- NN_model|>
    summary() |>
    coef() |>
    as.data.frame() |>
    rownames_to_column("VAR") |>
    select(VAR, `z value`, `Pr(>|z|)`) |>
    filter(VAR != "(Intercept)") |>
    rename(
      Z_val = `z value`,
      P_val = `Pr(>|z|)`) |>
    mutate(
      Model = "NN",
      CODE = CODE[i, "CODE"])
  
  PN_model <- errorsarlm(PN ~ CWB + d_CWB + AGE + CCI, listw = lw, data = Set) |> 
    summary()
  
  PN_p<- data.frame(CODE = CODE[i, "CODE"], 
                    Model = "PN",
                    Wald = PN_model[["Wald1"]][["statistic"]][["Wald statistic"]], 
                    W_p = PN_model[["Wald1"]][["p.value"]][["Wald statistic"]], 
                    AIC = AIC(PN_model), 
                    AIC_lm = PN_model[["AIC_lm.model"]],
                    pR2 = cor(as.data.frame(predict(PN_model))[,1],
                                   Set$PN)^2,
                    RMSE = (sqrt(mean((Set$PN - as.data.frame(predict(PN_model))[,1])^2))) / sd(Set$PN)
                    )
  
  PN_var<- PN_model|>
    summary() |>
    coef() |>
    as.data.frame() |>
    rownames_to_column("VAR") |>
    select(VAR, `z value`, `Pr(>|z|)`) |>
    filter(VAR != "(Intercept)") |>
    rename(
      Z_val = `z value`,
      P_val = `Pr(>|z|)`) |>
    mutate(
      Model = "PN",
      CODE = CODE[i, "CODE"])
  
  NP_model <- errorsarlm(NP ~ CWB + d_CWB + AGE + CCI, listw = lw, data = Set) |> 
    summary()
  
  NP_p<- data.frame(CODE = CODE[i, "CODE"], 
                    Model = "NP",
                    Wald = NP_model[["Wald1"]][["statistic"]][["Wald statistic"]], 
                    W_p = NP_model[["Wald1"]][["p.value"]][["Wald statistic"]], 
                    AIC = AIC(NP_model), 
                    AIC_lm = NP_model[["AIC_lm.model"]],
                    pR2 = cor(as.data.frame(predict(NP_model))[,"fit"],
                                   Set$NP)^2,
                    RMSE = (sqrt(mean((Set$NP - as.data.frame(predict(NP_model))[,1])^2))) / sd(Set$NP)
                    )
  
  NP_var<- NP_model|>
    summary() |>
    coef() |>
    as.data.frame() |>
    rownames_to_column("VAR") |>
    select(VAR, `z value`, `Pr(>|z|)`) |>
    filter(VAR != "(Intercept)") |>
    rename(
      Z_val = `z value`,
      P_val = `Pr(>|z|)`) |>
    mutate(
      Model = "NP",
      CODE = CODE[i, "CODE"])
  
  VAR_models<- rbind(VAR_models, PP_var, NN_var, PN_var, NP_var)
  P_models<- rbind(P_models, PP_p, NN_p, PN_p, NP_p)
  
}

VAR_models<- VAR_models[-1,]
P_models<- P_models[-1,]

write.table(P_models, "Graphs/Table_S1.txt", sep = "\t", dec = ",", row.names = F)

#==============================================#
###    4. Drivers daily   ####
#==============================================#

DOYs<- data.frame(DOY = c(121:639))

VARs_residuals<- data.frame(CODE = NA, SITE = NA, DOY = NA, CHANGE = NA, CAT = NA, Changes_res = NA)
VARs_spat_models<- data.frame(VAR = NA, DOY = NA, CODE = NA, Z_val = NA, P_val = NA)

for (i in c(1:nrow(CODE))) {
  
  print(i)
  
  # Selecting CODE
  Set<- Daily_cor[Daily_cor$CODE == CODE[i, "CODE"],]
  
  for (j in c(1:nrow(DOYs))) {
    
    print(paste("CODE:", i, "DOY: ", DOYs[j, "DOY"]))
    
    # Selecting DOY
    Set_DOY<- subset(Set, DOY == DOYs$DOY[j])
    
    # Detrending from R2 size
    model1<- lm(abs(CHANGE) ~ abs(LATE), data = Set_DOY)
    
    # Getting residuals
    Set_DOY$Changes_res<- residuals(model1)
    Set_DOY$Changes_res<- ifelse (Set_DOY$CHANGE > 0, sqrt(Set_DOY$Changes_res ^ 2), -1* sqrt(Set_DOY$Changes_res ^ 2))
    
    VARs_residuals<- rbind(VARs_residuals, Set_DOY[,c("CODE", "SITE", "DOY", "CHANGE", "CAT", "Changes_res")])
    
    # Spatial regression
    
    # Calculating neighbours
    coords <- as.matrix(Set_DOY[, c("LON", "LAT")])
    knb <- knn2nb(knearneigh(coords, k = 3))
    lw <- nb2listw(knb, style = "W")
    
    model2<- errorsarlm(Changes_res ~ CWB + d_CWB + AGE + CCI, listw = lw, data = Set_DOY) |>
      summary() |>
      coef() |>
      as.data.frame() |>
      rownames_to_column("VAR") |>
      select(VAR, `z value`, `Pr(>|z|)`) |>
      filter(VAR != "(Intercept)") |>
      rename(
        Z_val = `z value`,
        P_val = `Pr(>|z|)`) |>
      mutate(
        DOY = as.numeric(DOYs$DOY[j]),
        CODE = CODE[i, "CODE"])
    
    # Saving results
    VARs_spat_models <- rbind(VARs_spat_models, model2)
    
  }
}

VARs_residuals<- VARs_residuals[-1,]
VARs_spat_models<- VARs_spat_models[-1,]

saveRDS(VARs_residuals, "Data/Interim_data/VARs_residuals2.rds")
saveRDS(VARs_spat_models, "Data/Interim_data/VARs_spat_models2.rds")

#===========================================================================================#
###                          Fig. 1 (Map and climatic niches)                           ####
#===========================================================================================#

#=====================================================================#
####                   Map of study sites                        #####
#=====================================================================#

# Loading of list of sites (New list with 950 sites is filtered, so each site have 70 or more years in 2010)
Loc_sp<-readRDS("Data/Loc_list2.rds") |> 
  vect(geom = c("LON", "LAT"), crs="EPSG:4326")
Loc_sp$TREE_TYPE_order<- factor(Loc_sp$TREE_TYPE, levels = c("Conifer", "Broadleaf"))

# Loading the CWB raster
CWB_mean<- rast("Data/Clima_grids/CWB_time_windows.tif") |> mean()

# Creating the polygon of area of interest
AOI <- matrix(
  c(6, 45,
    6, 55,
    27, 55,
    27, 45,
    6, 45), ncol = 2, byrow = TRUE) |>
  vect(type = "polygons", crs = "EPSG:4326")

# Downloading and cropping state boundaries
Countries<- crop(vect(countriesLow), AOI)

# Loading digital elevation model, croping and hillshade calculating
DEM<- rast("C:/Users/Jiri/Desktop/Science/FG_layers/World/DEM of Europe/data/commonData/Data0/euro_dem/w001001.adf") |> crop(AOI)
Slope<- terrain(DEM, "slope", unit = "radians")
Aspect<- terrain(DEM, "aspect", unit = "radians")
Hillshade<- shade(Slope, Aspect, angle = 45, direction = 315, normalize = T)

Map<- ggplot()+
  
  geom_spatvector(data = AOI, fill = "#8bcaf7", color = "black", linewidth = 0.7) +
  
  new_scale_fill()+
  
  geom_spatraster(data = Hillshade, interpolate =T)+
  scale_fill_gradientn(colours = grey(0:100 / 100), na.value = "transparent")+
  guides(fill = "none")+
  
  new_scale_fill()+
  
  geom_spatraster(data = CWB_mean,interpolate =T, alpha = 0.7)+
  scale_fill_viridis(option="magma", direction=-1, na.value=NA)+
  labs(fill = "CWB (mm)")+
  theme(legend.title=element_text(size=14))+
  
  new_scale_fill()+
  geom_spatvector(data = Countries, fill = NA, color = "black", linewidth = 0.7) +
  
  new_scale_fill()+
  geom_spatvector(data = Loc_sp, aes(fill=SPECIES2, shape=TREE_TYPE), size=4)+
  
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = species_colors)+
  
  labs(x= "Longitude (°)", y= "Latitude (°)", shape = "")+
  guides(fill = "none")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 18))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-5))+
  theme(legend.key.width= unit(12, 'mm'))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

#Map
#ggsave("Graphs/Graphical_abstract/Map.tiff", height = 200, width = 250, units = "mm", dpi = 300)

#====================================================================================#
####                   Climatic niche of each species                        #####
#====================================================================================#

# Loading aggregated monthly rasters
EOBS_monthly_mean_T<- rast("Data/Clima_grids/EOBS_monthly_mean_T.tif")
EOBS_monthly_sum_P<- rast("Data/Clima_grids/EOBS_monthly_sum_P.tif")

# Aggregation of monthly rasters to yearly and general layer
EOBS_mean_T<- tapp(EOBS_monthly_mean_T, "years", mean) |> mean()
EOBS_mean_P<- tapp(EOBS_monthly_sum_P, "years", sum) |> mean()

# Loading shapefiles of specie's areal
ABAL<- vect("Data/Species_areals/ABAL_areal.shp")
PISY<- vect("Data/Species_areals/PISY_areal.shp")
PCAB<- vect("Data/Species_areals/PCAB_areal.shp")
FASY<- vect("Data/Species_areals/FASY_areal.shp")
QUsp<- vect("Data/Species_areals/QUsp_areal.shp")

SPECIES<- list(ABAL, PISY, PCAB, FASY, QUsp)
SPECIES2<- c("ABAL", "PISY", "PCAB", "FASY", "QUsp")

Clima_niches<- data.frame(SPECIES2 = NA, TEMP = NA, PREC = NA)

for (i in c(1:length(SPECIES))) {
  
  Areal_whole<- SPECIES[[i]]
  
  #Areal<- crop(Areal_whole, AOI)
  Areal<- Areal_whole
  
  Mean_t_areal<- crop(EOBS_mean_T, Areal, mask = T)
  Mean_p_areal<- crop(EOBS_mean_P, Areal, mask = T)
  
  Coords<- as.data.frame(crds(Mean_t_areal)) |>  vect(geom = c("x", "y"), crs="EPSG:4326")
  
  Coords$TEMP<- terra::extract(Mean_t_areal, Coords)[2]
  Coords$PREC<- terra::extract(Mean_p_areal, Coords)[2]
  Coords$SPECIES2<- SPECIES2[i]
  
  Clima_niches<- rbind(Clima_niches, as.data.frame(Coords)[,c(3, 1, 2)])
  
}

Clima_niches<- Clima_niches[-1,]
Clima_niches$TREE_TYPE<- ifelse(Clima_niches$SPECIES %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")
Clima_niches$TREE_TYPE_order<- factor(Clima_niches$TREE_TYPE, levels = c("Conifer", "Broadleaf"))

# Extracting climatic coordinates for sites
Loc_sp$TEMP<- terra::extract(EOBS_mean_T, Loc_sp)[2] |> unlist() |>  as.numeric()
Loc_sp$PREC<- terra::extract(EOBS_mean_P, Loc_sp)[2] |> unlist() |>  as.numeric()
Loc_sp$TREE_TYPE_order<- factor(Loc_sp$TREE_TYPE, levels = c("Conifer", "Broadleaf"))

ggplot()+
  
  geom_bin_2d(data = Clima_niches, aes(x = TEMP, y = PREC))+
  scale_fill_gradient(low = "#fcfc95", high = "#0814c4", na.value = NA)+
  new_scale_fill()+
  
  geom_point(data = Loc_sp, aes(x = TEMP, y = PREC, fill=SPECIES2, shape=TREE_TYPE), size = 3)+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(limits = c(-6, 18), breaks = c(-5, 0, 5, 10, 15))+
  scale_y_continuous(limits = c(100, 3400))+
  
  facet_nested(~TREE_TYPE_order+SPECIES2)+
  labs(x = "MAT (°C)", y = "MAP (mm)")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "none", legend.margin=margin(t=-5))+
  theme(legend.key.width= unit(8, 'mm'))

ggsave("Graphs/Graphical_abstract/Niches.tiff", height = 100, width = 320, units = "mm", dpi = 300)

#====================================================================================#
####                   Gradients                        #####
#====================================================================================#

Variables<- Loc_list[, c("SITE", "SPECIES2", "CWB", "d_CWB", "AGE", "CCI")] |> 
  pivot_longer(values_to = "VAL", names_to = "VAR", 3:6)

ggplot(Variables)+
  geom_histogram(aes(x = VAL, fill = SPECIES2), bins = 20, color = "black")+
  
  facet_wrap(~VAR, scales = "free", ncol = 2, 
             labeller = labeller(VAR = c("AGE" = "AGE (years)", 
                                         "CCI" = "CCI (index)",
                                         "CWB" = "CWB (mm)",
                                         "d_CWB" = "∆CWB (%)")))+
  
  scale_fill_manual(values = species_colors)+
  labs(x= "", y= "Number of sites", fill = "")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(l=-50))+
  theme(panel.grid = element_blank())

ggsave("Graphs/Graphical_abstract/Factors coverage.tiff", height = 230, width = 150, units = "mm", dpi = 300)

#===================================================================#
###                   Fig. 2 (Mean correlations)                 ####
#===================================================================#

Sig_change_cor<- Daily_cor |> 
  filter(abs(CHANGE) > 0.27) |>
  pivot_longer(values_to = "COR", names_to = "T_WIN", cols = c(5,6))
Sig_change_cor$T_WIN_f<- ifelse(Sig_change_cor$T_WIN == "EARLY", "EARLY 1951–1980", "LATE 1981–2010")
Sig_change_cor$T_WIN_f<- factor(Sig_change_cor$T_WIN_f, c("EARLY 1951–1980", "LATE 1981–2010"))
Sig_change_cor$CLIM2<- ifelse(Sig_change_cor$CLIM=="TEMP", "Temperature", "Precipitation")
Sig_change_cor$TREE_TYPE<- ifelse(Sig_change_cor$SPECIES2 %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")

Sig_change_cor<- Sig_change_cor |>
  group_by(CLIM2, DOY, SPECIES2, TREE_TYPE, T_WIN_f) |>
  summarise(
    mean_COR = mean(COR, na.rm=T),
    sd_COR = sd(COR, na.rm=T),
    .groups="drop")

breaks<- c(121, 182, 244, 305, 366, 425, 486, 547, 609)+15 # plus 15 to shift it to middle of month
labels <- c("May", "Jul", "Sep", "Nov", "JAN", "MAR", "MAY", "JUL", "SEP")

Mean_corr_gg<- ggplot(Sig_change_cor)+
  
  geom_line(aes(x = DOY, y = mean_COR, color = T_WIN_f))+
  
  geom_ribbon(aes(x = DOY, ymin = mean_COR - sd_COR, ymax = mean_COR + sd_COR, fill = T_WIN_f), 
              alpha = 0.2)+
  
  scale_color_manual(values = c("#40a843", "#b0680b"))+
  scale_fill_manual(values = c("#40a843", "#b0680b"))+
  
  facet_nested(TREE_TYPE~SPECIES2~CLIM2)+
  
  geom_vline(xintercept = 172, color = "#888a8c")+
  geom_vline(xintercept = 365)+
  geom_vline(xintercept = 537, color = "#888a8c")+
  
  geom_hline(yintercept = 0, linetype = "dotted")+
  
  labs(x = "", y = "Mean correlation", fill = "", color = "")+
  scale_x_continuous(breaks = breaks, labels = labels)+
  scale_y_continuous(limits = c(-0.7, 0.7))+
  
  theme(strip.text.y = element_blank())+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black", angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-20))+
  theme(legend.key.width= unit(10, 'mm'))

#ggsave("Graphs/Fig. 2 (Mean correlations).tiff", height = 200, width = 250, units = "mm", dpi = 300)

#====================================================================================#
####                   Mean change                        #####
#====================================================================================#

Mean_change <- Daily_cor |>
  filter(abs(CHANGE) > 0.27) |>
  group_by(SPECIES2, CLIM, DOY) |>
  summarise(
    Mean_CHANGE = mean(CHANGE, na.rm = TRUE),
    SD_CHANGE   = sd(CHANGE, na.rm = TRUE),
    .groups = "drop")

Mean_change$CLIM2<- ifelse(Mean_change$CLIM=="TEMP", "Temperature", "Precipitation")
Mean_change$TREE_TYPE<- ifelse(Mean_change$SPECIES2 %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")

breaks<- c(121, 182, 244, 305, 366, 425, 486, 547, 609)+15 # plus 15 to shift it to middle of month
labels <- c("May", "Jul", "Sep", "Nov", "JAN", "MAR", "MAY", "JUL", "SEP")

Mean_change_gg<- ggplot()+
  
  geom_area(data = Mean_change, aes(x = DOY, y = Mean_CHANGE, fill = CLIM2),
            alpha = 0.6, color = "black")+
  
  facet_nested(TREE_TYPE~SPECIES2~"Correlation change")+
  geom_vline(xintercept = 172, color = "#888a8c")+
  geom_vline(xintercept = 365)+
  geom_vline(xintercept = 537, color = "#888a8c")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  
  scale_x_continuous(breaks = breaks, labels = labels)+
  scale_color_manual(values = c("#3965f7", "#c72e35"))+
  scale_fill_manual(values = c("#3965f7", "#c72e35"))+
  
  #theme(strip.text.y = element_blank())+
  labs(x = "", y = "Mean correlation difference ∆COR", fill = "") +
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black", angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-20))+
  theme(legend.key.width= unit(10, 'mm'))+
  theme(plot.margin = margin(t=5, r=5, b=5, l=0))

ggarrange(Mean_corr_gg, Mean_change_gg, ncol = 2, widths = c(1.5, 1), align = "h",
          common.legend = F, legend = "bottom", 
          labels = c("A","B"), font.label = list(size = 20))+
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("Graphs/Fig. 2 (Correlation change).tiff", height = 250, width = 350, units = "mm", dpi = 300)

#===================================================================#
###             Fig. 3 (Year_profiles)             ####
#===================================================================#


#====================================================================================#
####                   Year profiles                        #####
#====================================================================================#

# Calculating percentages in each category for DOYs
Profiles_year <- Daily_cor |>
  group_by(SPECIES2, DOY, CLIM) |>
  summarise(
    PP   = mean(PP,   na.rm = TRUE),
    NN   = mean(NN,   na.rm = TRUE),
    PN   = mean(PN,   na.rm = TRUE),
    NP   = mean(NP,   na.rm = TRUE),
    WEAK = mean(WEAK, na.rm = TRUE),
    .groups = "drop") |>
  mutate(
    TOTAL = PP + NN + PN + NP + WEAK,
    
    PP   = PP   / TOTAL * 100,
    NN   = NN   / TOTAL * 100,
    PN   = PN   / TOTAL * 100,
    NP   = NP   / TOTAL * 100,
    WEAK = WEAK / TOTAL * 100) |>
  pivot_longer(
    cols = c(PP, NN, PN, NP, WEAK),
    names_to = "CAT",
    values_to = "COUNT") |>
  select(SPECIES2, DOY, CLIM, CAT, COUNT)

Profiles_year$CLIM2<- ifelse(Profiles_year$CLIM=="TEMP", "Temperature", "Precipitation")
Profiles_year$TREE_TYPE<- ifelse(Profiles_year$SPECIES2 %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")
Profiles_year$CAT2 <- factor(Profiles_year$CAT, levels = c("WEAK","PP", "NN", "NP", "PN"))

breaks<- c(121, 182, 244, 305, 366, 425, 486, 547, 609)+15 # plus 15 to shift it to middle of month
labels <- c("May", "Jul", "Sep", "Nov", "JAN", "MAR", "MAY", "JUL", "SEP")

ggplot()+
  
  geom_area(data = Profiles_year, aes(x = DOY, y = COUNT, fill = CAT2, color = CAT2), 
            linewidth = 0.2, color = "black")+
  scale_fill_manual(values = cat_colors,     
                    labels = c(PP = "PP", NN = "NN", NP = "NP", PN = "PN", WEAK = "Stationary"))+
  
  geom_vline(xintercept = 172, color = "#888a8c")+
  geom_vline(xintercept = 365)+
  geom_vline(xintercept = 537, color = "#888a8c")+
  geom_hline(yintercept = 50, linetype = "dotted")+
  
  facet_nested(TREE_TYPE~SPECIES2~CLIM2)+
  
  coord_cartesian(ylim = c(0, 100))+
  guides(color = "none")+
  labs(x= "", y= "Mean percentage of correlations (%)", fill = "Shift direction: ")+

  scale_x_continuous(breaks = breaks, labels = labels)+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black", angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-20))+
  theme(plot.margin = margin(t=5, r=5, b=5, l=2))

ggsave("Graphs/Fig. 3 (Percantage profiles).tiff", height = 250, width = 250, units = "mm", dpi = 300)

#===================================================================#
###             Fig. 4 (Drivers)             ####
#===================================================================#

VAR_models<- VAR_models |> separate(CODE, into = c("SPECIES2", "CLIM"), sep = "_")
VAR_models$CLIM2<- ifelse(VAR_models$CLIM=="TEMP", "Temperature", "Precipitation")
VAR_models$TREE_TYPE<- ifelse(VAR_models$SPECIES2 %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")
VAR_models$Model_f <- factor(VAR_models$Model, levels = c("PP", "NN", "NP", "PN"))
VAR_models$Model2<- ifelse(VAR_models$Model %in% c("PP", "NN"), "Persistent", "Reversed")
VAR_models$Model2<- factor(VAR_models$Model2, levels = c("Persistent", "Reversed"))

ggplot() +
  geom_bar(data = VAR_models,
           aes(x = Model_f, y = Z_val, fill = VAR, alpha = P_val < 0.05, color = P_val < 0.05), 
           stat = "identity", position = position_dodge()) +
  
  geom_vline(xintercept = 1.5, linetype = "dotted", color = "#888a8c")+
  
  scale_fill_manual(values = var_colors, name = "Variable: ",
                    labels = c("AGE" = "AGE", 
                               "CCI" = "CCI", 
                               "CWB" = "CWB",
                               "d_CWB" ="∆CWB")) +
  scale_color_manual(values = c("#cdcfd1", "black"))+
  scale_alpha_manual(
    values = c("TRUE" = 1, "FALSE" = 0.2),
    labels = c("FALSE" = "No", "TRUE" = "Yes"),
    name = "Significant: ") +
  
  geom_hline(yintercept = 0) +
  
  facet_nested(TREE_TYPE ~ SPECIES2 ~ CLIM2 + Model2,scales = "free_x") +
  labs(x = "", y = "Z value") +
  
  guides(alpha = guide_legend(order = 1, nrow = 1),
         fill  = guide_legend(order = 2, nrow = 1),
         color = "none") +
  
  theme(
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 15),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5),
    legend.box.margin = margin(t = -10),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5))

ggsave("Graphs/Fig. 4 (Drivers).tiff", height = 230, width = 230, units = "mm", dpi = 300)

#===================================================================#
###             Fig. 5 (Year_profiles_vars)             ####
#===================================================================#

# Preparing driving factors for ggplot
VARs_spat_models<- readRDS("Data/Interim_data/VARs_spat_models2.rds")
VARs_spat_models$CLIM2<- ifelse(VARs_spat_models$CLIM=="TEMP", "Temperature", "Precipitation")
VARs_spat_models$TREE_TYPE<- ifelse(VARs_spat_models$SPECIES2 %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")

breaks<- c(121, 182, 244, 305, 366, 425, 486, 547, 609)+15 # plus 15 to shift it to middle of month
labels <- c("May", "Jul", "Sep", "Nov", "JAN", "MAR", "MAY", "JUL", "SEP")

ggplot() +
  
  geom_bar(data = VARs_spat_models, aes(x = DOY, y = Z_val, fill = VAR, alpha = P_val < 0.05),
           stat = "identity", position = "stack") +
  
  scale_fill_manual(values = var_colors, name = "Variable: ",
                    labels = c("AGE" = "AGE", 
                               "CCI" = "CCI", 
                               "CWB" = "CWB",
                               "d_CWB" ="∆CWB")) +
  
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3),
                     labels = c("FALSE" = "No", "TRUE" = "Yes"),
                     name = "Significant: ") +
  
  facet_nested(TREE_TYPE ~ SPECIES2 ~ CLIM2) +
  
  geom_vline(xintercept = 172, color = "#888a8c") +
  geom_vline(xintercept = 365) +
  geom_vline(xintercept = 537, color = "#888a8c") +
  
  scale_x_continuous(breaks = breaks, labels = labels) +
  
  labs(x = "", y = "Z value") +
  
  guides(alpha = guide_legend(order = 1, nrow = 1),
         fill  = guide_legend(order = 2, nrow = 1)) +
  
  theme(
    strip.placement = "outside",
    axis.ticks = element_line(color = "black"),
    strip.text = element_text(size = 15),
    axis.text.x = element_text(size = 15, color = "black",
                               angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -5),
    legend.box.margin = margin(t = -10),
    legend.key.width = unit(10, "mm"),
    
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5))

ggsave("Graphs/Fig. 5 (Year_profiles_vars).jpg", height = 250, width = 250, units = "mm", dpi = 300)

#===================================================================#
###             Fig. S1 (TEMP_COR)             ####
#===================================================================#

# Adjustments for ggplot
Temp_cor<- Daily_cor[Daily_cor$CLIM == "TEMP",] |> 
  pivot_longer(values_to = "COR", names_to = "T_WIN", cols = c(5,6))
Temp_cor$T_WIN_f<- ifelse(Temp_cor$T_WIN == "EARLY", "EARLY 1951–1980", "LATE 1981–2010")
Temp_cor$T_WIN_f<- factor(Temp_cor$T_WIN_f, c("EARLY 1951–1980", "LATE 1981–2010"))
Temp_cor$CLIM2<- ifelse(Temp_cor$CLIM=="TEMP", "Temperature", "Precipitation")
Temp_cor$TREE_TYPE<- ifelse(Temp_cor$SPECIES %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")

breaks<- c(121, 182, 244, 305, 366, 425, 486, 547, 609)+15 # plus 15 to shift it to middle of month
labels <- c("May", "Jul", "Sep", "Nov", "JAN", "MAR", "MAY", "JUL", "SEP")

ggplot(Temp_cor)+
  geom_line(aes(x = DOY, y = COR, group = SITE, color = CWB))+
  scale_color_viridis(option="magma", direction=-1, name = "CWB (mm)")+
  facet_nested(TREE_TYPE~SPECIES~T_WIN_f)+
  
  geom_vline(xintercept = 172, color = "#888a8c")+
  geom_vline(xintercept = 365)+
  geom_vline(xintercept = 537, color = "#888a8c")+
  
  geom_hline(yintercept = 0, linetype = "dotted")+
  
  labs(x = "", y = "Correlation")+
  scale_x_continuous(breaks = breaks, labels = labels)+
  scale_y_continuous(limits = c(-0.7, 0.7))+
  
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black", angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-20))+
  theme(legend.key.width= unit(10, 'mm'))

ggsave("Graphs/Fig. S1 (TEMP_COR).tiff", height = 250, width = 300, units = "mm", dpi = 300)

#===================================================================#
###             Fig. S2 (PREC_COR)             ####
#===================================================================#

# Adjustments for ggplot
Prec_cor<- Daily_cor[Daily_cor$CLIM == "PREC",] |> 
  pivot_longer(values_to = "COR", names_to = "T_WIN", cols = c(5,6))
Prec_cor$T_WIN_f<- ifelse(Prec_cor$T_WIN == "EARLY", "EARLY 1951–1980", "LATE 1981–2010")
Prec_cor$T_WIN_f<- factor(Prec_cor$T_WIN_f, c("EARLY 1951–1980", "LATE 1981–2010"))
Prec_cor$CLIM2<- ifelse(Prec_cor$CLIM=="TEMP", "Temperature", "Precipitation")
Prec_cor$TREE_TYPE<- ifelse(Prec_cor$SPECIES %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")

ggplot(Prec_cor)+
  geom_line(aes(x = DOY, y = COR, group = SITE, color = CWB))+
  scale_color_viridis(option="magma", direction=-1, name = "CWB (mm)")+
  facet_nested(TREE_TYPE~SPECIES~T_WIN_f)+
  
  geom_vline(xintercept = 172, color = "#888a8c")+
  geom_vline(xintercept = 365)+
  geom_vline(xintercept = 537, color = "#888a8c")+
  
  geom_hline(yintercept = 0, linetype = "dotted")+
  
  labs(x = "", y = "Correlation")+
  scale_x_continuous(breaks = breaks, labels = labels)+
  scale_y_continuous(limits = c(-0.7, 0.7))+
  
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black", angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-20))+
  theme(legend.key.width= unit(10, 'mm'))

ggsave("Graphs/Fig. S2 (PREC_COR).tiff", height = 250, width = 300, units = "mm", dpi = 300)

#===============================================================#
###                   Fig. S3 (Percentage bars)                ####
#===============================================================#

Change_bars <- Daily_cor |>
  group_by(SPECIES2, CLIM) |>
  summarise(
    PP   = mean(PP,   na.rm = TRUE),
    NN   = mean(NN,   na.rm = TRUE),
    PN   = mean(PN,   na.rm = TRUE),
    NP   = mean(NP,   na.rm = TRUE),
    WEAK = mean(WEAK, na.rm = TRUE),
    .groups = "drop") |>
  mutate(
    TOTAL = PP + NN + PN + NP + WEAK,
    
    PP   = PP   / TOTAL * 100,
    NN   = NN   / TOTAL * 100,
    PN   = PN   / TOTAL * 100,
    NP   = NP   / TOTAL * 100,
    WEAK = WEAK / TOTAL * 100) |>
  pivot_longer(
    cols = c(PP, NN, PN, NP, WEAK),
    names_to = "CAT",
    values_to = "COUNT") |>
  select(SPECIES2, CLIM, CAT, COUNT)

Change_bars$CLIM2<- ifelse(Change_bars$CLIM=="TEMP", "Temperature", "Precipitation")
Change_bars$TREE_TYPE<- ifelse(Change_bars$SPECIES2 %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")
Change_bars$CAT2 <- factor(Change_bars$CAT, levels = c("PP", "NN", "NP", "PN"))
Change_bars<- Change_bars[Change_bars$CAT != "WEAK",]

ggplot(Change_bars, aes(x = "", y= COUNT, fill = CAT2))+
  geom_bar(stat="identity", color = "black")+
  facet_nested(CLIM2~TREE_TYPE+SPECIES2, scales = "free", space='free')+
  scale_fill_manual(values = cat_colors)+
  
  geom_text(aes(label = round(COUNT, 0)),
            position = position_stack(vjust = 0.5), size = 4, color = "black")+
  
  labs(x= "", y= "Mean percentage (%)", fill = "Shift direction: ")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-20))+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

ggsave("Graphs/Fig. 3 (Change_bars).tiff", height = 170, width = 250, units = "mm", dpi = 300)

#==========================================================================================================#
###                                 Fig. S4 (Percentage maps)                                  ####
#==========================================================================================================#

Sites<- pivot_longer(CHANGE_SITE, values_to = "PER", names_to = "CAT", cols = 4:7)
Sites<- Sites |>  vect(geom = c("LON", "LAT"), crs="EPSG:4326")
Sites$TREE_TYPE<- ifelse(Sites$SPECIES2 %in% c("PISY", "PCAB", "ABAL"), "Conifer", "Broadleaf")
Sites$CAT_f = factor(Sites$CAT, levels=c('PP','NN','PN','NP'))

ggplot()+
  
  geom_spatvector(data = AOI, fill = "#8bcaf7", color = "black", linewidth = 0.7) +
  
  new_scale_fill()+
  
  theme(legend.title=element_text(size=14))+
  
  new_scale_fill()+
  geom_spatvector(data = Countries, fill = "grey", color = "black") +
  
  new_scale_fill()+
  
  geom_spatvector(data = Sites, aes(fill=PER, shape=TREE_TYPE), size=3)+
  scale_fill_viridis(option="viridis", direction=-1, na.value=NA)+
  
  scale_shape_manual(values = c(21, 24))+
  
  facet_wrap(~CAT_f)+
  labs(x= "Longitude (°)", y= "Latitude (°)", fill = "Percentage (%)")+
  guides(shape = "none")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 12))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-5))+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

ggsave("Graphs/Fig. S4 (Percentage maps).tiff", height = 230, width = 280, units = "mm", dpi = 300)

#===============================================================#
###                   Fig. S5 (Factors maps)                ####
#===============================================================#

# Downloading and cropping state boundaries
Countries<- crop(vect(countriesLow), AOI)

#====================================================================================#
####                   CWB map                        #####
#====================================================================================#
CWB_map<- ggplot()+
  
  geom_spatvector(data = AOI, fill = "#8bcaf7", color = "black", linewidth = 0.7) +
  
  new_scale_fill()+
  
  geom_spatraster(data = mean(CWB_windows),interpolate =T)+
  scale_fill_viridis(option="magma", direction=-1, na.value=NA)+
  labs(fill = "CWB (mm)")+
  theme(legend.title=element_text(size=14))+
  
  new_scale_fill()+
  geom_spatvector(data = Countries, fill = NA, color = "black") +
  
  facet_wrap(~"Climatic water balance")+
  labs(x= "Longitude (°)", y= "Latitude (°)", shape = "")+
  guides(fill = "none")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 12))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-5))+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

#====================================================================================#
####                   d CWB map                        #####
#====================================================================================#

d_CWB_map<- ggplot()+
  
  geom_spatvector(data = AOI, fill = "#8bcaf7", color = "black", linewidth = 0.7) +
  
  new_scale_fill()+
  
  geom_spatraster(data = d_CWB,interpolate =T)+
  scale_fill_gradient2(low = "#db1304", mid = "#fcfab3", high = "#266ad1", midpoint = 0, na.value=NA)+
  labs(fill = "Relative ∆CWB")+
  theme(legend.title=element_text(size=14))+
  
  new_scale_fill()+
  geom_spatvector(data = Countries, fill = NA, color = "black") +
  
  facet_wrap(~"Relative ∆ of climatic water balance")+
  labs(x= "Longitude (°)", y= "Latitude (°)", shape = "")+
  guides(fill = "none")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 12))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-5))+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

#====================================================================================#
####                   CCI map                        #####
#====================================================================================#

CCI_map<- ggplot()+
  
  geom_spatvector(data = AOI, fill = "#8bcaf7", color = "black", linewidth = 0.7) +
  
  new_scale_fill()+
  
  geom_spatraster(
    data = CCI, interpolate = T) +
  
  scale_fill_gradient2(low = "#343aed", mid = "#fcfab3", high = "#ebe70e", 
                       midpoint = 30, na.value=NA)+
  
  labs(fill = "CCI")+
  theme(legend.title=element_text(size=14))+
  
  new_scale_fill()+
  geom_spatvector(data = Countries, fill = NA, color = "black") +
  
  facet_wrap(~"Conrad's continentality index")+
  labs(x= "Longitude (°)", y= "Latitude (°)", shape = "")+
  guides(fill = "none")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 12))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-5))+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

#====================================================================================#
####                   AGE map                        #####
#====================================================================================#

AGE_map<- ggplot()+
  
  geom_spatvector(data = AOI, fill = "#8bcaf7", color = "black", linewidth = 0.7) +
  
  new_scale_fill()+
  
  theme(legend.title=element_text(size=14))+
  
  new_scale_fill()+
  geom_spatvector(data = Countries, fill = "grey", color = "black") +
  
  new_scale_fill()+
  
  geom_spatvector(data = Loc_sp, aes(fill=AGE, shape=TREE_TYPE), size=4)+
  
  scale_shape_manual(values = c(21, 24))+
  scale_fill_gradient2(low = "#edac58", mid = "#fcfab3", high = "#0c9103", 
                       midpoint = 150, na.value=NA)+
  
  facet_wrap(~"Mean age of trees")+
  labs(x= "Longitude (°)", y= "Latitude (°)", fill = "Age (years)")+
  guides(shape = "none")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 12))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t=-5))+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

ggarrange(CWB_map, d_CWB_map, CCI_map, AGE_map, ncol = 2, nrow = 2,
          labels = c("A", "B", "C", "D"), font.label = list(size = 25))+
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("Graphs/Fig. S5 (Factors maps).tiff", height = 250, width = 250, units = "mm", dpi = 300)

#==========================================================================================================#
###                                 Fig. S6 (Tree-ring chronologies)                                  ####
#==========================================================================================================#

Tree_data<-readRDS("Data/Interim_data/Tree_data.rds")

Tree_data$YEAR<- as.numeric(rownames(Tree_data))
Tree_gg<- pivot_longer(Tree_data, names_to = "SITE", values_to = "VAL", 1:ncol(Tree_data)-1)
Tree_gg<- left_join(Tree_gg, Loc_list[,c(4, 6, 7, 16)], by = "SITE")
Tree_gg$TREE_TYPE_order<- factor(Tree_gg$TREE_TYPE, levels = c("Conifer", "Broadleaf"))
Tree_gg<-Tree_gg[!(Tree_gg$SITE %in% c("P810104PCAB")),]

ggplot(Tree_gg, aes(x=YEAR, y=VAL, group=SITE, colour = CWB))+
  geom_line(linewidth = 0.2)+
  
  geom_vline(xintercept = 1980)+
  
  facet_nested_wrap(~TREE_TYPE_order+SPECIES2, axes = "all", nrow = 2)+
  scale_color_viridis(option="magma", direction=-1, na.value=NA)+
  
  labs(x ="", y= "Tree-ring index", color = "CWB (mm)")+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(strip.text  = element_text(size = 15, color = "black"))+
  theme(axis.text.x = element_text(size = 15, color = "black"))+
  theme(axis.text.y = element_text(size = 15, color = "black"))+
  theme(panel.spacing = unit(1.5, "lines"))+
  theme(axis.title.x = element_text(size = 15, color = "black"))+
  theme(axis.title.y = element_text(size = 15, color = "black"))+
  theme(legend.title = element_text(size = 15, color = "black"))+
  theme(legend.position = "bottom")+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.text = element_text(size = 12))

ggsave("Graphs/Fig. S6 (Tree_rings).tiff", height = 180, width = 280, units = "mm", dpi = 300)

#===============================================================#
###                   Fig. S7 (Breakpoints)                ####
#===============================================================#

Clima_breakpoints<- readRDS("Data/Interim_data/Clima_breakpoints.rds")

BP_plot <- Clima_breakpoints |>
  filter(Breakpoint_year >= 1955,
         Breakpoint_year <= 2005) |>
  mutate(
    Month = factor(month.abb[Month], levels = month.abb))
BP_plot$CLIM2<- ifelse(BP_plot$CLIMA=="TEMP", "Temperature", "Precipitation")

ggplot(BP_plot[BP_plot$CLIM2 == "Temperature",]) +
  geom_histogram(aes(x = Breakpoint_year),
                 binwidth = 1,
                 boundary = 1950.5,
                 color = "black",
                 fill = "#d95f02") +
  
  scale_x_continuous(
    breaks = seq(1950, 2010, 10),
    limits = c(1950, 2010)) +
  labs(x ="", y= "Number of sites")+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(strip.text  = element_text(size = 15, color = "black"))+
  theme(axis.text.x = element_text(size = 15, color = "black"))+
  theme(axis.text.y = element_text(size = 15, color = "black"))+
  theme(panel.spacing = unit(1.5, "lines"))+
  theme(axis.title.x = element_text(size = 15, color = "black"))+
  theme(axis.title.y = element_text(size = 15, color = "black"))+
  theme(legend.title = element_text(size = 15, color = "black"))+
  theme(legend.position = "bottom")+
  theme(legend.key.width= unit(8, 'mm'))+
  theme(legend.text = element_text(size = 12))

ggsave("Graphs/Fig. S7 (Breakpoints).tiff", height = 100, width = 150, units = "mm", dpi = 300)

#===============================================================#
###                   Fig. S8 (Multicollinearity)                ####
#===============================================================#

Intercor <- rcorr(as.matrix(Loc_list[, c("CWB", "d_CWB", "AGE", "CCI")]))

Cor_mat <- data.frame(Intercor$r,
                      VARs = c("CWB", "d_CWB", "AGE", "CCI")) |>
  pivot_longer(-VARs,
               names_to = "VARs2",
               values_to = "COR")

P_mat <- data.frame(Intercor$P,
                    VARs = c("CWB", "d_CWB", "AGE", "CCI")) |>
  pivot_longer(-VARs,
               names_to = "VARs2",
               values_to = "P_val")


Cor_mat$P_val<- P_mat$P_val
Cor_mat<- Cor_mat[-c(1:4, 6:8, 11, 12, 16),]

order<- c("CWB", "d_CWB", "AGE", "CCI")

COR_gg<- ggplot()+ 
  geom_tile(data = Cor_mat, aes(factor(VARs, level=order), factor(VARs2, level=order), fill= COR), color="black")+
  scale_fill_gradient2(low = "#d7191c", mid= "#ffffbf", high = "#2c7bb6", name="Correlation")+
  
  geom_text(data = Cor_mat, aes(VARs, VARs2, label = round(COR, digits = 2)), color = "black", size = 5) +
  
  scale_x_discrete(labels = c("AGE" = "AGE", 
                              "CCI" = "CCI",
                              "d_CWB" ="∆CWB")) +
  
  scale_y_discrete(labels = c("AGE" = "AGE",
                              "CWB" = "CWB",
                              "d_CWB" ="∆CWB")) +
  labs(x= "", y= "")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "none", legend.margin=margin(t=-25))+
  theme(legend.title = element_blank())+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

VIF<- data.frame(VIF = vif(lm(mean_RW ~ AGE + d_CWB + CWB + CCI, data = Loc_list)),
                 VAR = c("AGE", "d_CWB", "CWB", "CCI"))

VIF_gg<- ggplot()+geom_col(data = VIF, aes(x = VAR, y = VIF, fill = VAR), 
                           color = "black", linewidth = 0.2, fill = var_colors)+
  
  scale_x_discrete(labels = c("AGE" = "AGE", 
                              "CCI" = "CCI", 
                              "CWB" = "CWB",
                              "d_CWB" ="∆CWB")) +
  
  labs(x= "", y= "Variance inflation factor")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15, color="black"))+
  theme(axis.text.y = element_text(size = 15, color="black"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "none", legend.margin=margin(t=-25))+
  theme(legend.title = element_blank())+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

ggarrange(COR_gg, VIF_gg, ncol = 2, widths = c(1.3, 1), align = "h",
          labels = c("A","B"), font.label = list(size = 20))+
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("Graphs/Fig. S8 (Multicolinearity).tiff", height = 120, width = 220, units = "mm", dpi = 300)
