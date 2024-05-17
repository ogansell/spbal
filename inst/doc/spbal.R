## ----include = FALSE----------------------------------------------------------
library(bookdown)
knitr::opts_chunk$set( 
  fig.dim = c(8, 6),
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(spbal)
library(ggplot2)
library(gridExtra)

## ----srs1---------------------------------------------------------------------
# Create a random sample of 20 with a seed of 99 from a population of 100.
rand_samp <- spbal::SRS(seed = 99, total_rows = 100, sample_size = 20)
rand_samp

## ----BASex1a------------------------------------------------------------------
# Use the North Carolina shapefile supplied in the sf R package.
shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
shp_gates <- shp_file[shp_file$NAME == "Gates",]
shp_gates
# Vertically aligned master sample bounding box.
bb <- spbal::BoundingBox(shapefile = shp_gates)

## ----BASex2a,warning=FALSE, message=FALSE-------------------------------------
set.seed(511)
n_samples <- 20
result <- spbal::BAS(shapefile = shp_gates,
                     n = n_samples,
                     boundingbox = bb)
BAS20 <- result$sample

## ----BASex3a------------------------------------------------------------------
BAS20

## ----BASex1c------------------------------------------------------------------
# 1. Plot using ggplot
gates <- sf::st_as_sf(shp_gates, coords = c("longitude", "latitude"))

ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_text(data = BAS20,
            size = 3,
            aes(label = spbalSeqID,
                vjust = -0.5,
                geometry = geometry,
                x = after_stat(x),
                y = after_stat(y), color = "ID"),
            stat = "sf_coordinates") +
  geom_point(data = BAS20,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Samples"),
             stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red", "black"),
                     name = "Legend") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "Spatially balanced sample using a equal probability BAS design.",
       subtitle = "Total of 20 survey sites from Gates, North Carolina",
       caption = "spbal Equal probability BAS sample.")

## ----BASex4a,message=FALSE, warning=FALSE-------------------------------------
n_samples <- 50
result2 <- spbal::BAS(shapefile = shp_gates,
                      n = n_samples,
                      boundingbox = bb,
                      seeds = result$seed)
BAS50 <- result2$sample
BAS50

# Check, first n_samples points in both samples must be the same.
all.equal(BAS20$geometry, BAS50$geometry[1:20])

## ----BASex5a, warning=FALSE, message=FALSE------------------------------------
## Convert foreign object to an sf object for ggplot.
gates <- sf::st_as_sf(shp_gates, coords = c("longitude", "latitude"))
ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_text(data = BAS50,
            size = 2,
            aes(label = BAS50$spbalSeqID,
                vjust = -0.7,
                geometry = geometry,
                x = after_stat(x),
                y = after_stat(y), color = "BAS Order"),
            stat = "sf_coordinates") +
  geom_point(data = BAS50,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour = "BAS n = 20 + 30"),
                 stat = "sf_coordinates")+
  geom_point(data = BAS20,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour = "BAS n = 20"),
                 stat = "sf_coordinates"
  ) +
 scale_color_manual(values = c("red","blue","black"),
                    name = "Legend") +  
 theme_bw() +
 labs(x = "Latitude",
      y = "Longitude",
      title = "BAS samples from the Gates county.")


## ----BASex6a, results='hide', message=FALSE-----------------------------------
strata <- c("Gates", "Northampton", "Hertford")
n_strata <- c("Gates" = 20, "Northampton" = 20, "Hertford" = 10)
shp_subset <- shp_file[shp_file[["NAME"]] %in% strata,]

## ----BASex7a, results='hide', message=FALSE-----------------------------------
bb_strata <- spbal::BoundingBox(shapefile = shp_subset)

## ----BASex8a, results='hide', message=FALSE-----------------------------------
set.seed(511)
result3 <- spbal::BAS(shapefile = shp_subset,
                      n = n_strata,
                      boundingbox = bb_strata,
                      stratum = "NAME")
BASMaster <- result3$sample
gates_samp <- BASMaster[BASMaster[["NAME"]] %in% "Gates",]

## ----BASex9a------------------------------------------------------------------
gates_samp 

## ----BASex10a-----------------------------------------------------------------
strat <- sf::st_as_sf(shp_subset, coords = c("longitude", "latitude"))

gates_samp       <- BASMaster[BASMaster[["NAME"]] %in% "Gates",]
northampton_samp <- BASMaster[BASMaster[["NAME"]] %in% "Northampton",]
hertford_samp    <- BASMaster[BASMaster[["NAME"]] %in% "Hertford",]

ggplot() +
  geom_sf() +
  geom_sf(data = strat, size = 4, shape = 23) +
  geom_point(data = gates_samp,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Gates"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = northampton_samp,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Northampton"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = hertford_samp,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Hertford"),
             stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red", "blue", "green"),
                     name = "Legend") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "Stratified samples from BAS master samples.",
       subtitle = "50 sites from North Carolina, U.S.A, 20 from Gates, \n20 from Northampton and 10 from Hertford.",
       caption = "spbal Stratified BAS sample.")


## ----BASex11f-----------------------------------------------------------------
n_strata <- c("Gates" = 20, "Northampton" = 20, "Hertford" = 20)
result4 <- spbal::BAS(shapefile = shp_subset,
                      n = n_strata,
                      boundingbox = bb_strata,
                      seeds = result3$seed,
                      stratum = "NAME")
BASMaster <- result4$sample
gates_samp2 <- BASMaster[BASMaster[["NAME"]] %in% "Gates",]
gates_samp2

## ----BASex12a-----------------------------------------------------------------
# Ensure gates_samp is equal to the first 10 sites in gates_samp2. Must return TRUE.
all.equal(gates_samp$geometry[1:20], gates_samp2$geometry[1:20])

## ----BASex13a-----------------------------------------------------------------
gates_samp4       <- BASMaster[BASMaster[["NAME"]] %in% "Gates",]
northampton_samp4 <- BASMaster[BASMaster[["NAME"]] %in% "Northampton",]
hertford_samp4    <- BASMaster[BASMaster[["NAME"]] %in% "Hertford",]

ggplot() +
  geom_sf() +
  geom_sf(data = strat, size = 4, shape = 23) +
  geom_point(data = gates_samp4,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Gates"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = northampton_samp4,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Northampton"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = hertford_samp4,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Hertford"),
             stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red", "blue", "green"),
                     name = "Legend") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "Stratified samples from BAS master samples.",
       subtitle = "60 sites from North Carolina, U.S.A, 20 from Gates, \n20 from Northampton and 20 from Hertford.",
       caption = "spbal Stratified BAS samples.")

## ----BASex14g-----------------------------------------------------------------
set.seed(511)

n_panels <- c(20, 20, 20)
result5 <- spbal::BAS(shapefile = shp_gates,
                      panels = n_panels,
                      boundingbox = bb)
BASpanel <- result5$sample
BASpanel

## ----BASex15a-----------------------------------------------------------------
panel_2 <- spbal::getPanel(BASpanel, 2)
panel_2 <- panel_2$sample
panel_2

## ----BASex16a, warning=FALSE--------------------------------------------------
# to extract the sample points associated with a specific panelid, we can use the following code:
panel_1 <- BASpanel[BASpanel$panel_id == 1,]
panel_2 <- BASpanel[BASpanel$panel_id == 2,]
panel_3 <- BASpanel[BASpanel$panel_id == 3,]

# or use the spbal::getPanel() function.
#panel_1 <- spbal::getPanel(BASpanel, 1)
#panel_2 <- spbal::getPanel(BASpanel, 2)
#panel_3 <- spbal::getPanel(BASpanel, 3)

ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_point(data = panel_1,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-1 (n = 20)"),
                 stat = "sf_coordinates"
  ) +
  geom_point(data = panel_2,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-2 (n = 20)"),
                 stat = "sf_coordinates"
  ) +
  geom_point(data = panel_3,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-3 (n = 20)"),
                 stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red", "blue", "black"),
                                     name = "Panels") + 
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "Panel Design from BAS Master Sample",
       subtitle = "Total of 60 survey sites from Gates, North Carolina, U.S.A",
       caption = "spbal Non-overlapping Panel Design.")

## ----BASex17a, results='hide', message=FALSE----------------------------------
set.seed(511)
n_panels <- c(20, 20, 20)
n_panel_overlap <- c(0, 5, 5)
result6 <- spbal::BAS(shapefile = shp_gates,
                      panels = n_panels,
                      panel_overlap = n_panel_overlap,
                      boundingbox = bb)
BASpanel <- result6$sample

## ----BASex18a-----------------------------------------------------------------
panel_2 <- spbal::getPanel(BASpanel, 2)
panel_2 <- panel_2$sample
panel_2[1:5,]

## ----BASex19a-----------------------------------------------------------------
panel_1 <- spbal::getPanel(BASpanel, 1)
panel_2 <- spbal::getPanel(BASpanel, 2)
panel_3 <- spbal::getPanel(BASpanel, 3)

ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_jitter(width = 10, height = 10) +
  geom_point(data = sf::st_jitter(panel_1$sample, 0.02),
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-1"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = sf::st_jitter(panel_2$sample, 0.02),
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-2"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = panel_3$sample,
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-3"),
             stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red", "blue", "green"),
                     name = "Overlapped Panels") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "Spatially balanced sample using a overlapping panel \ndesign, of three panels, each containing 20 survey sites.",
       subtitle = "Total of 50 survey sites from Gates, North Carolina. Panel Overlap = (0, 5, 5)",
       caption = "spbal Overlapping Panel Design.")

## ----HFex1a, message=FALSE, warning=FALSE, results='hide'---------------------
set.seed(511)
result6 <- spbal::HaltonFrame(shapefile = shp_gates, 
                              J = c(3, 2),
                              boundingbox = bb)
Frame <- result6$hf.pts.shp
Grid <- result6$hg.pts.shp

## ----HFex1b-------------------------------------------------------------------
# Grid - Halton grid over Gates county.
ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_text(data = Grid,
            size = 3,
            aes(label = ID,
                vjust = -0.5,
                geometry = x,
                x = after_stat(x),
                y = after_stat(y), color = "ID"),
            stat = "sf_coordinates") +
  geom_point(data = Grid,
             aes(geometry = x,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Samples"),
             stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red", "black"),
                     name = "Legend") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "A Halton grid, B = 2^3 * 3^2, over Gates, North Carolina, U.S.A.",
       subtitle = "A total of 432 points.",
       caption = "spbal Halton Grid example.")

## ----HFex1c-------------------------------------------------------------------
# Frame - Halton frame over Gates county.
ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_text(data = Frame,
            size = 3,
            aes(label = ID,
                vjust = -0.5,
                geometry = x,
                x = after_stat(x),
                y = after_stat(y), color = "ID"),
            stat = "sf_coordinates") +
  geom_point(data = Frame,
             aes(geometry = x,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Samples"),
             stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red", "black"),
                     name = "Legend") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "A Halton frame, B = 2^3 * 3^2, over Gates, North Carolina, U.S.A.",
       subtitle = "Showing sample points within the Gates shapefile.",
       caption = "spbal Halton Frame example.")

## ----HFex1d, results='hide', message=FALSE, warning=FALSE---------------------
set.seed(511)
result7 <- spbal::HaltonFrame(shapefile = shp_gates,
                              J = c(8, 5),
                              boundingbox = bb)
Frame <- result7$hf.pts.shp

## ----HFex1e, message=FALSE, warning=FALSE-------------------------------------
n_samples <- 25
FrameSample <-getSample(shapefile = Frame, 
                        n = n_samples)
FrameSample <- FrameSample$sample
FrameSample[1:10, c("x", "spbalSeqID")]

## ----HFex1f, warning=FALSE----------------------------------------------------
ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_text(data = FrameSample, 
            size = 3,
            aes(label = spbalSeqID,
                vjust = -0.5,
                geometry = x,
                x = after_stat(x),
                y = after_stat(y), color = "Frame Order"),
                stat = "sf_coordinates") +
  geom_point(data = FrameSample, 
             aes(geometry = x,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Samples"),
                 stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("black", "red"),
                     name = "Legend") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "First 25 sites from the Gates Halton Frame.")

## ----HFex1g, warning=FALSE, message=FALSE-------------------------------------
n_samples <- 20
FrameSample <-getSample(shapefile = Frame, 
                        n = n_samples, 
                        randomStart = TRUE)
FrameSample <- FrameSample$sample
FrameSample[1:10, c("x", "spbalSeqID")]

## ----HFex3a-------------------------------------------------------------------
set.seed(511)

# Three panels, of 20 samples each.
panels <- c(20, 20, 20)

# second panel overlaps first panel by 5, and third panel 
# overlaps second panel by 5.
panel_overlap <- c(0, 5, 5)

# generate the sample.
samp <- spbal::HaltonFrame(J = c(4, 3),
                           boundingbox = bb,
                           panels = panels,
                           panel_overlap = panel_overlap,
                           shapefile = shp_gates)

# get halton frame data from our sample.
samp3 <- samp$hf.pts.shp
samp3

panelid <- 1
olPanel_1 <- spbal::getPanel(samp3, panelid)

panelid <- 2
olPanel_2 <- spbal::getPanel(samp3, panelid)

panelid <- 3
olPanel_3 <- spbal::getPanel(samp3, panelid)

# Plot using ggplot2
ggplot() +
  geom_sf() +
  geom_sf(data = gates, size = 4, shape = 23) +
  geom_jitter(width = 10, height = 10) +
  geom_point(data = sf::st_jitter(olPanel_1$sample, 0.02),
             aes(geometry = x,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-1"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = sf::st_jitter(olPanel_2$sample, 0.02),
             aes(geometry = x,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-2"),
             stat = "sf_coordinates"
  ) +
  geom_point(data = olPanel_3$sample,
             aes(geometry = x,
                 x = after_stat(x),
                 y = after_stat(y), colour= "Panel-3"),
             stat = "sf_coordinates"
  ) +
  scale_color_manual(values = c("red",
                                "blue",
                                "green"),
                     name = "Overlapped Panels") +  # Define color scale and legend title
  theme_bw() +
  labs(x = "Latitude",
       y = "Longitude",
       title = "Halton Frame sample using a overlapping panel \ndesign, of three panels, each containing 20 survey sites.",
       subtitle = "Total of 50 survey sites from Gates, North Carolina, U.S.A. Panel Overlap = (0, 5, 5)",
       caption = "spbal Halton Frame Overlapping Panel Design.")

## ----HFex4a, results='hide', message=FALSE------------------------------------
strata <- c("Gates", "Northampton", "Hertford")
n_strata <- c("Gates" = 20, "Northampton" = 30, "Hertford" = 10)
shp_subset <- shp_file[shp_file[["NAME"]] %in% strata,]

## ----HFex5a, results='hide', message=FALSE------------------------------------
bb_strata <- spbal::BoundingBox(shapefile = shp_subset)

## ----HFex6a, results='hide', message=FALSE, warning=FALSE---------------------
set.seed(511)
result9 <- spbal::HaltonFrame(shapefile = shp_subset,
                              N = n_strata,
                              J = c(8, 5),
                              boundingbox = bb_strata,
                              stratum = "NAME")
Frame <- result9$hf.pts.shp

## ----HFex7a, warning=FALSE----------------------------------------------------
n_samples <- 10
hertford_samp <- spbal::getSample(Frame, 
                                  n = n_samples, 
                                  strata = "Hertford", 
                                  stratum = "NAME")
hertford_samp <- hertford_samp$sample
hertford_samp[1:10, c("NAME", "spbalSeqID", "x")]

## ----HFex8a, warning=FALSE, message=FALSE, results='hide'---------------------
set.seed(511)
panels <- c(20, 20, 20)
n_panel_overlap <- c(0, 5, 5)
result10 <- spbal::HaltonFrame(shapefile = shp_gates,
                               panels = panels,
                               panel_overlap = n_panel_overlap,
                               J = c(8, 5),
                               boundingbox = bb)
HaltonFramePanel <- result10$hf.pts.shp

## ----HFex9a, warning=FALSE, message=FALSE-------------------------------------
panelid <- 1
SitesPanel_1 <- spbal::getPanel(HaltonFramePanel, panelid)
SitesPanel_1 <- SitesPanel_1$sample
SitesPanel_1[1:10, c("x", "spbalSeqID", "panel_id")]

## ----HIPex1a------------------------------------------------------------------
# set random seed
base::set.seed(511)

# define HIP parameters.
pop <- matrix(stats::runif(5000*2), nrow = 5000, ncol = 2)
n <- 20
its <- 7

# Convert the population matrix to an sf point object.
sf_points <- sf::st_as_sf(data.frame(pop), coords = c("X1", "X2"))
dim(sf::st_coordinates(sf_points))

# generate HIP sample.
result <- spbal::HIP(population = sf_points, 
                     n = n, 
                     iterations =  its)

# HaltonIndex
HaltonIndex <- result$HaltonIndex
# verify all spread equally, should be TRUE.
(length(unique(table(HaltonIndex))) == 1)

# Population Sample
HIPsample <- result$sample
HIPsample

## ----HIPex2a------------------------------------------------------------------
HIPoverSample <- result$overSample
HIPoverSample[1:10, c("geometry", "spbalSeqID")]

OverSampleSize <- dim(HIPoverSample)[1]
OverSampleSize

## ----HIPex3a------------------------------------------------------------------
# compare the HIP sample and oversample, they will be the same.
all.equal(HIPsample$geometry[1:n], HIPoverSample$geometry[1:n])

