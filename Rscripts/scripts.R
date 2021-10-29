library("here")
library("tidyverse")
library("signal")
library("mgcv")
library("MALDIquant")


source(here("Rscripts", "functions.R"))

##
data_to_prep <- tibble(.filename = c("0035.D", "0036.D", "0037.D"))
data_log <- read_csv2(file = here("data", "GC-MS_log.csv" ))


#############################
# Get Chromatogram Data
chromatograms <- data_log %>% 
  mutate(chromatogram  = map(.x = Code, .f = import_chromatogram_manual)) %>% 
  unnest(cols = c(chromatogram)) 

que_1 <- plot_chromatogram_gg(.x = chromatograms, facet = FALSE)
que_1 <- plot_chromatogram_gg(.x = chromatograms, facet = TRUE)



#############################
# add in the maldi objects as column in the tibble
chromatograms <- chromatograms %>% 
  mutate(maldi = pmap(.l = list(.x = time,
                                .y = intensity,
                                GC.code = GC.code,
                                pollen.code = pollen.code), 
                      .f = parse_chrom_MALDI))

# check working
str(chromatograms$maldi[[3]])

#############################
# use maldiquant workflow to smooth spectra and work on baseline correction procedures
chrom_maldi <- as.list(chromatograms$maldi)
chrom_maldi <- transformIntensity(chrom_maldi,method="sqrt")

# chrom_maldi_sg <- smoothIntensity(chrom_maldi, method="SavitzkyGolay", halfWindowSize=10)

baseline <- map(.x = chrom_maldi, .f = estimateBaseline, method="SNIP", iterations=100)
plot(chrom_maldi[[2]])
lines(baseline[[2]], col="red", lwd=2)

baseline_removed <- map(.x = chrom_maldi, .f = removeBaseline, method="SNIP", iterations=100)
plot(baseline_removed[[3]])
# check if this should be included ?calibrateIntensity

# aligning chromatograms- does this do much?
aligned_chromatograms <- alignSpectra(baseline_removed, halfWindowSize=20, 
                        SNR=2, tolerance=0.002, warpingMethod="lowess")

plot(aligned_chromatograms[[3]])





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




# # standardise by the number of pollen grains in the sample
# chromatograms <- chromatograms %>% 
#   unnest(cols = c(intensity)) %>% 
#   mutate(intensity =intensity/Grains ) %>% 
#   nest(intensity = intensity)

# Because of the variable number of grains we will need to standardise against a peak in the sample





