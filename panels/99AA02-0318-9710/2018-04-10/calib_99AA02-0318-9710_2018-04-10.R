# Initial calibration of Spectralon disk
# Panel 99AA02-0318-9710
# Laboratory reference disk (used to calibrate other disks)

# Load libraries
library(googledrive)
library(tidyverse)
library(spectrolab)

# Panel ID
panel_id <- '99AA02-0318-9710'
calibration_date <- '2018-04-10'
  
# default working directory
default_wd <- getwd()

# set new working directory
setwd(paste0(default_wd, '/', 'panels/', panel_id, '/', calibration_date))

# Sign in as Caboscience and download initial calibration file
file_name <- paste0('SRS-99-010-', panel_id, '.txt')
drive_download(file_name, overwrite = T)

# Read the raw calibration data
calib_raw <- read.table(file_name) %>% 
  as_tibble() %>% 
  rename(wvl = V1, refl = V2)
  
# remove six really high values >.997
calib_raw2 <- filter(calib_raw, refl < 0.997)

# subset of data with refl > 0.997
calib_sub <- filter(calib_raw, refl >= 0.997)

# Plot the original data
calib_raw_plot <- ggplot(calib_raw, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Reflectance') +
  theme_bw()
calib_raw_plot

# Zoom around 860 nm
calib_raw_plot +
  coord_cartesian(xlim = c(750, 1000))

# Plot sensor matching regions
calib_raw_plot +
  coord_cartesian(xlim = c(850, 875))

# Convert data to wide format
calib_raw_wide <- calib_raw2 %>% 
  spread(key = wvl, value = refl) %>% 
  mutate(SN = paste0('SRS-99-010-', panel_id)) %>% 
  as.data.frame()
calib_raw_spectra <- as.spectra(calib_raw_wide, name_idx = ncol(calib_raw_wide))                    
plot(calib_raw_spectra)

# Match sensors
calib_match <- match_sensors(calib_raw_spectra, splice_at = 860,
                             interpolate_wvl = 100)
plot(calib_match)

# Convert to data frame
calib_match_df <- as.data.frame(calib_match) %>% 
  as_data_frame() %>% 
  select(-sample_name) %>% 
  gather(key = wvl, value = refl) %>% 
  mutate(wvl = as.integer(wvl)) %>% 
  bind_rows(calib_sub) %>% 
  rename(refl_match = refl) %>% 
  arrange(wvl)
 

# Resample to 1-nm and smooth
wvls <- 250:2500
calib_smooth_0.1 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.1)
calib_smooth_0.2 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.2)
calib_smooth_0.3 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.3)
calib_smooth_0.4 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.4)
calib_smooth_0.5 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.5)
calib_smooth_0.6 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.6)
calib_smooth_0.7 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.7)
calib_smooth_0.8 <- spectrolab::resample(calib_match, new_wvls = wvls, spar = 0.8)

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
  rename(refl_raw = refl) %>%
  left_join(calib_match_df, by = 'wvl')
  

# Plot match data and smoothers
calib_smooth_plot <- ggplot(calib_df, aes(x = wvl, y = refl_raw)) +
  geom_line(colour = 'black') +
  geom_line(aes(y = refl_match), colour = 'blue') +
  geom_line(aes(y = `0.5`), colour = 'red') +
  geom_line(aes(y = `0.4`), colour = 'green') +
  geom_line(aes(y = `0.3`), colour = 'orange') +
  ylab('Reflectance') +
  xlab('Wavelength (nm)')
calib_smooth_plot

# Zoom in
calib_smooth_plot +
  coord_cartesian(xlim = c(250, 500))

# 0.4 seems best
calib_smooth_plot_0.4 <- ggplot(calib_df, aes(x = wvl, y = refl_raw)) +
  geom_line(colour = 'black') +
  geom_line(aes(y = refl_match), colour = 'blue') +
  geom_line(aes(y = `0.4`), colour = 'orange') +
  ylab('Reflectance') +
  xlab('Wavelength (nm)')
calib_smooth_plot_0.4
ggsave('calib_match_smooth.png', width = 6, height = 4.5, dpi = 600)

# Create a matrix to store smoothed reflectance curve
calib_smooth_mat <- matrix(c(calib_df$wvl, signif(calib_df$`0.4`, 4)), byrow = F, ncol = 2, nrow = nrow(calib_df))
smooth_file_name <- paste0('SRS-99-010-', panel_id, '-', calibration_date, '.calib')

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
