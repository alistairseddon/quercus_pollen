library("here")
library("tidyverse")
library("signal")
library("mgcv")
source(here("Rscripts", "functions.R"))
?sgolayfilt
##
data_to_prep <- tibble(.filename = c("0035.D", "0036.D", "0037.D"))

# Get Chromatogram Data
chromatograms <- data_to_prep %>% 
  mutate(chromatogram  = map(.x = .filename, .f = import_chromatogram_manual)) %>% 
  unnest(cols = c(chromatogram)) 

# Do SG filter - need to look into this to understand parameter settings
chromatograms <- chromatograms %>%
   group_by(pollen.code) %>% 
  mutate(SG = map(.x = intensity, .f = sgolayfilt ))
 
que_1 <- plot_chromatogram_gg(.x = chromatograms, facet = TRUE)

# Do baseline correction.
# Multiple options to explore. Here I recommend a gam

for_baseline_correction <- chromatograms %>% 
  #unnest(cols = c(time, intensity, SG)) %>% 
  dplyr::filter(pollen.code =="E839")

baseline_corrected <- chromatograms %>% 
  mutate(baseline_est = map2(.x =time, .y = SG, .f = baseline.gam )) %>% 
  unnest(cols =c(time, intensity, SG, baseline_est)) %>% 
  mutate(intensity_bc = SG - baseline_est)

que_1 <- plot_chromatogram_gg(.x = baseline_corrected, facet = TRUE, baseline = TRUE)
que_1
# Make a plot
que_1_bc <- plot_bc_chromatogram_gg(.x = baseline_corrected, facet = TRUE)
que_1_bc






