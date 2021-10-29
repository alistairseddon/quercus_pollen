import_chromatogram_manual <- function(.filename){
  
  path <- here("data", "GC-MS_raw_files", .filename, "CHROMTAB.CSV")
  # imports the actual data
  raw_chromatogram <- suppressWarnings(read.csv(path, 
                                                skip =2, 
                                                row.names= NULL, 
                                                col.names = c("Time", "Intensity") ) )
  # NB. not sure why we are getting the warning here but it doesn't seem to matter
  
  # imports the meta data and extracts the two sample codes
  meta_data <- read.csv(path, 
                        nrows =1, 
                        row.names= NULL)
  
  sample_code_GC <- meta_data$File
  sample_code_pollen <- meta_data$Sample
  
  chromatogram <- tibble(GC.code = sample_code_GC,
                         pollen.code = sample_code_pollen,
                          time = list(raw_chromatogram$Time), 
                         intensity = list(raw_chromatogram$Intensity),
                         ) 
  chromatogram
  
}




###########################
# Plot a chromatogram using ggplot
plot_chromatogram_gg <-function(.x = chromatograms, facet= FALSE) {
  # .x is a chromatograms tibble with nested time and intensity columns
  p <- .x %>% 
    unnest(cols = c(time, intensity, SG)) %>% 
    ggplot(aes(x = time, y = intensity, groups = pollen.code)) +
    geom_line(aes(col = pollen.code)) +
    geom_line(aes(x = time, y = SG), col = "black") +
    xlab("Time (minutes)") +
    ylab("Intensity (m/z)") +
    theme_bw()

    if(facet == TRUE) {
      p <- p +  
        facet_wrap(vars(pollen.code), dir = "v") +
        theme(legend.position = "none")
    }
 p 
}


######
baseline.gam <- function(.x = for_baseline_correction$time, .y = for_baseline_correction$SG){
  SG <- .y
  time <- .x
  bs_fit <- bam(SG ~ s(time), family = Gamma()) 
  # plot(bs_fit, resid = TRUE)
  predict(bs_fit, type = "response")
}

