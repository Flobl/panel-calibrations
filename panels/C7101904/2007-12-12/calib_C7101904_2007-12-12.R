# Initial calibration of Spectralon disk
# Panel C7101904
# ASD Sphere (2017 rental)

# Load libraries
library(googledrive)
library(tidyverse)
library(spectrolab)

# Panel ID
panel_id <- 'C7101904'
calibration_date <- '2007-12-12'
  
# default working directory
default_wd <- getwd()

# set new working directory
setwd(paste0(default_wd, '/', 'panels/', panel_id, '/', calibration_date))

# Sign in as Caboscience and download initial calibration file
file_name <- 'c7101904_specchioformat.txt'
drive_download(file_name, overwrite = T)

# Read the raw calibration data
calib_raw <- read.table(file_name, header = T) %>% 
  as_tibble() %>% 
  rename(refl = c7101904)
  
# Plot the original data
calib_raw_plot <- ggplot(calib_raw, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Reflectance') +
  theme_bw()
calib_raw_plot

# Zoom around 400 nm
calib_raw_plot +
  coord_cartesian(xlim = c(250, 500))

# Plot sensor matching regions
calib_raw_plot +
  coord_cartesian(xlim = c(850, 875))

# Convert data to wide format
calib_raw_wide <- calib_raw %>% 
  spread(key = wvl, value = refl) %>% 
  mutate(SN = panel_id) %>% 
  as.data.frame()
calib_raw_spectra <- as.spectra(calib_raw_wide, name_idx = ncol(calib_raw_wide))                    
plot(calib_raw_spectra)

 

# Resample to 1-nm and smooth
calib_smooth_0.1 <- smooth(calib_raw_spectra, spar = 0.1)
calib_smooth_0.2 <- smooth(calib_raw_spectra, spar = 0.2)
calib_smooth_0.3 <- smooth(calib_raw_spectra, spar = 0.3)
calib_smooth_0.4 <- smooth(calib_raw_spectra, spar = 0.4)
calib_smooth_0.5 <- smooth(calib_raw_spectra, spar = 0.5)
calib_smooth_0.6 <- smooth(calib_raw_spectra, spar = 0.6)
calib_smooth_0.7 <- smooth(calib_raw_spectra, spar = 0.7)
calib_smooth_0.8 <- smooth(calib_raw_spectra, spar = 0.8)

# Combine smoothed spectra
calib_comb <- combine(calib_smooth_0.1,
                      calib_smooth_0.2) %>% 
  combine(calib_smooth_0.3) %>% 
  combine(calib_smooth_0.4) %>% 
  combine(calib_smooth_0.5) %>%
  combine(calib_smooth_0.6) %>%
  combine(calib_smooth_0.7) %>% 
  combine(calib_smooth_0.8)

# Convert back to data frame and semi-wide format
calib_df <- as.data.frame(calib_comb) %>% 
  as_data_frame() %>% 
  select(-sample_name) %>% 
  mutate(smooth = seq(0.1, 0.8, by = 0.1)) %>% 
  gather(key = wvl, value = refl, -smooth) %>% 
  mutate(wvl = as.integer(wvl)) %>% 
  dplyr::arrange(smooth, wvl) %>% 
  spread(smooth, value = refl) %>% 
  left_join(calib_raw, by = 'wvl') %>%
  rename(refl_raw = refl)
  

# Plot match data and smoothers
calib_smooth_plot <- ggplot(calib_df, aes(x = wvl, y = refl_raw)) +
  geom_line(colour = 'black') +
  geom_line(aes(y = `0.4`), colour = 'orange') +
  geom_line(aes(y = `0.5`), colour = 'red') +
  geom_line(aes(y = `0.6`), colour = 'green') +
  geom_line(aes(y = `0.7`), colour = 'blue') +
  ylab('Reflectance') +
  xlab('Wavelength (nm)')
calib_smooth_plot

# Zoom in
calib_smooth_plot +
  coord_cartesian(xlim = c(250, 500))

# 0.5 seems best
calib_smooth_plot_0.5 <- ggplot(calib_df, aes(x = wvl, y = refl_raw)) +
  geom_line(colour = 'black') +
  geom_line(aes(y = `0.5`), colour = 'orange') +
  ylab('Reflectance') +
  xlab('Wavelength (nm)')
calib_smooth_plot_0.5
ggsave('calib_match_smooth.png', width = 6, height = 4.5, dpi = 600)

# Create a matrix to store smoothed reflectance curve
calib_smooth_mat <- matrix(c(calib_df$wvl, signif(calib_df$`0.5`, 4)), byrow = F, ncol = 2, nrow = nrow(calib_df))
smooth_file_name <- paste0(panel_id, '-', calibration_date, '.calib')

# save a new calib file to google drive
write.table(file = smooth_file_name,
            calib_smooth_mat,
            quote = F,
            row.names = F,
            col.names = F,
            sep = '\t')

# Upload to Google Drive folder
parent_path <- '/CABO/DATA/SPECTROSCOPY/PANELS/'
path_name <- paste0(parent_path, panel_id, '/', calibration_date, '/')
drive_upload(smooth_file_name, path = path_name)
drive_upload('calib_match_smooth.png', path = path_name)
