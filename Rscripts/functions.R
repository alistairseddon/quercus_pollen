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
# Plot a chromatogram using ggplot using a MALDIquant object
plot_chromatogram_gg <-function(.x) {
  
  data.frame(Time = mass(.x), Intensity = intensity(.x))    %>% 
    ggplot(aes( x = Time, y = Intensity)) +
    geom_line() +
    xlab("Time (minutes)") +
    ylab("Intensity (m/z)") +
    theme_bw()
}


######


