library("here")
library("tidyverse")
library("signal")
library("mgcv")
library("MALDIquant")
library("vegan")


source(here("Rscripts", "functions.R"))

##
data_to_prep <- tibble(.filename = c("0035.D", "0036.D", "0037.D"))
data_log <- read_csv2(file = here("data", "GC-MS_log.csv" ))


#############################
# Get Chromatogram Data
chromatograms <- data_log %>% 
  mutate(chromatogram  = map(.x = Code, .f = import_chromatogram_manual)) %>% 
  unnest(cols = c(chromatogram)) %>% 
  dplyr::filter(Grains > 250) %>% 
  group_by(Tree)

que_1 <- plot_chromatogram_gg(.x = chromatograms, facet = FALSE)
que_2 <- plot_chromatogram_gg(.x = chromatograms, facet = TRUE)

#############################
# Standardise accoring to n grains

chromatograms <- chromatograms %>% mutate(intensity = map2(.x = intensity, .y = Grains, 
                                          .f = function(.x, .y ){
                                            .x/.y
                                          }))
que_3 <- plot_chromatogram_gg(.x = chromatograms, facet = TRUE)
#############################
# add in the maldi objects as column in the tibble
chromatograms <- chromatograms %>% 
  mutate(maldi = pmap(.l = list(.x = time,
                                .y = intensity,
                                GC.code = GC.code,
                                pollen.code = pollen.code,
                                Tree = Tree), 
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
sample_codes <- factor(sapply(aligned_chromatograms, function(x)metaData(x)$Tree))

#### average spectra
avgSpectra <- averageMassSpectra(aligned_chromatograms, labels=sample_codes,  method="mean")
# peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=20, SNR=2)
# peaks <- binPeaks(peaks, tolerance=0.002)
# peaks <- filterPeaks(peaks, minFrequency=0.25)
# featureMatrix <- intensityMatrix(peaks, avgSpectra)
# rownames(featureMatrix) <- names(avgSpectra)
avgChrom <- tibble(species = names(avgSpectra), time = map(avgSpectra, "mass"), intensity = map(avgSpectra, "intensity") )
plot_avg_chromatogram_gg(avgChrom)
plot_avg_chromatogram_gg(avgChrom, facet= TRUE)

#### raw samples for PCA plot
peaks <- detectPeaks(aligned_chromatograms, method="MAD", halfWindowSize=20, SNR=2)
peaks <- binPeaks(peaks, tolerance=0.002)
peaks <- filterPeaks(peaks, minFrequency=0.25)
featureMatrix <- intensityMatrix(peaks, aligned_chromatograms)

devtools::install_github("gavinsimpson/ggvegan")

library(ggvegan)

pca <- vegan::rda(featureMatrix)
str(summary(pca))

pc1_importance <- round(summary(pca)$cont$importance[2,1]*100, 1)
pc2_importance <- round(summary(pca)$cont$importance[2,2]*100, 1)

autoplot(pca, layers = c("species", "sites", "biplot"), arrows= TRUE)
pca_fort <- fortify(pca, display = "sites") %>% 
  mutate(taxa =  chromatograms$Tree)

ggplot(pca_fort, aes(x = PC1, y = PC2, colour = taxa, group = taxa)) +
  #scale_color_viridis_d() +
  coord_equal() +
  theme_bw() +
  geom_hline(linetype = "dashed", yintercept = 0,  col = "darkgrey") +
  geom_vline(linetype = "dashed", xintercept = 0, col = "darkgrey") +
  geom_point(show.legend = FALSE, size= 3) +
  xlab(paste0("PC1 (", pc1_importance, "%)")) +
  ylab(paste0("PC2 (", pc2_importance, "%)")) 


pca_fort_species <- fortify(pca, display = "species") 

PC1 <- arrange(pca_fort_species, desc(abs(PC1)))[1:12,]
PC2 <- arrange(pca_fort_species, desc(abs(PC2)))[1:12,]














par(mfrow = c(3,1))
plot(avgSpectra[[1]])
plot(avgSpectra[[2]])
plot(avgSpectra[[3]])







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





