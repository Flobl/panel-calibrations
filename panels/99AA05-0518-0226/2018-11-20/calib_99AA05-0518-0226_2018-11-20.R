# Second calibration of Spectralon disk
# Panel 99AA05-0518-0226
# Sphere S/N 402 replacement reference disk
# 2018-11-20: measurement dirty and clean

# Load libraries
library(googledrive)
library(tidyverse)
library(spectrolab)
library(googlesheets)

# Panel ID
panel_id <- '99AA05-0518-0226'
calibration_date <- '2018-11-20'

# default prefix for working directory
prefix_wd <- '/Users/etienne/Documents/ProjetsRecherche/CABO/panel-calibrations/panels/'

# set new working directory
setwd(paste0(prefix_wd, panel_id, '/', calibration_date))
getwd()

# Get Antoine's metadata (created by Etienne)
gs_auth() # Log on with LEFO
metadata <- gs_url('https://docs.google.com/spreadsheets/d/1ubYPl-CWu0ioNDR5CSZmf9rhLLAzQMB3FDlIRhSEkzI/edit#gid=0')
meta <- gs_read_csv(metadata) %>% 
  arrange(prefix, spectrum_number)
  meta
  
# get file names
meta <- meta %>%
  mutate(filename = paste0(prefix, '_', spectrum_number, '.sig'))

# get working folder in which spectra files are stored, in Google Drive
working_folder <- '2018-11-13-Calibration SPC_99AA05-0518-0226'

# download all files from working folder to local directory (LEFO account)
files <- drive_ls(working_folder) %>% 
  filter(name %in% meta$filename)
mypath <- paste0(getwd(), '/spectra')
for (i in 1:nrow(files)) {
  drive_download(files[i,], path = paste0(mypath,'/',files$name[i]), overwrite = T )
}


# function to extract the wavelengths + counts from .sig file
get_sig_data <- function(filename) {
  con <- file(filename, open = 'r')
  sig_all <- readLines(con)
  data_line <- which(sig_all == "data= ")
  close(con)
  con2 <- file(filename, open = 'r')
  readLines(con2, n = data_line)
  spec_data <- readLines(con2)
  close(con2)
  spec_data <- gsub(',', '.', spec_data)
  spec_data2 <- strsplit(spec_data, '  ')
  spec_data_num <- lapply(spec_data2, function(x) as.numeric(x))
  spec_data_mat <- do.call(rbind, spec_data_num)
  colnames(spec_data_mat) <- c('wvl', 'ref', 'tar', 'refl')
  return(as.data.frame(spec_data_mat) )
}

# apply function to each file
filenames <- list.files(mypath)
paths <- paste0(mypath,'/',filenames)
spec_data_all <- lapply(paths, get_sig_data) ; names(spec_data_all) <- filenames


# convert list to long data frame with filename column
for (i in 1:length(spec_data_all)) spec_data_all[[i]]$filename <- names(spec_data_all)[i]
spec_data_all_df <- bind_rows(spec_data_all) %>% as_data_frame() %>% 
  mutate(refl = refl / 100)

# merge with metadata
meta_sub <- meta %>% 
  select(filename, spectrum_number, type, group, spec_num, cleaning) %>% 
  mutate(spec_fac = factor(spec_num))
spec_all <- meta_sub %>% 
  left_join(spec_data_all_df, by = 'filename')

# check detector overlap
all_wvls <- filter(spec_all, filename == 'gr111318_0000.sig')$wvl
plot(all_wvls, type = 'l')

# find detector overlap regions
all_wvls[512:524] # remove 513:523
rem1 <- all_wvls[513:523]
all_wvls[768:774] # remove 769:773
rem2 <- all_wvls[769:773]
rem <- c(rem1, rem2)

# remove those wavelengths
spec_all <- spec_all %>% 
  filter(!wvl %in% rem)

# Plot the spectra, before cleaning
spec_dirty <- filter(spec_all, cleaning == 'before')
spectra_plot_before <- ggplot(spec_dirty, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = spec_fac, linetype = type, group = filename)) +
  facet_wrap(~ group) +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance') +
  coord_cartesian(ylim = c(.90, 1.10))
spectra_plot_before
ggsave('spectra_before_cleaning.pdf', plot = spectra_plot_before,
       width = 6, height = 5)

# Take average of all reps
spec_dirty_avg <- spec_dirty %>% 
  filter(type == 'target') %>% 
  group_by(group, wvl) %>% 
  summarise_at(vars(refl), mean) %>% 
  group_by(wvl) %>% 
  summarise_at(vars(refl), mean)

# Plot average dirty spectra
spec_dirty_avg_plot <- ggplot(spec_dirty_avg, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance')
spec_dirty_avg_plot


# download calibrated reflectance of reference disk 99AA02-0318-9710
ref_panel_id <- '99AA02-0318-9710'
ref_calib <- read.table('/Users/etienne/Documents/ProjetsRecherche/CABO/panel-calibrations/panels/99AA02-0318-9710/2018-04-10/SRS-99-010-99AA02-0318-9710-2018-04-10.calib') %>% 
  as_data_frame() %>% 
  rename(wvl = V1, refl = V2)
ref_calib

# convert to wide format and to spectra
ref_calib_wide <- ref_calib %>% 
  spread(key = wvl, value = refl) %>% 
  mutate(SN = ref_panel_id) %>% 
  as.data.frame()
ref_calib_spectra <- as.spectra(ref_calib_wide, name_idx = ncol(ref_calib_wide))                    
plot(ref_calib_spectra)

# resample to same wvls as target spectra
wvls <- filter(spec_dirty_avg, wvl <= 2500)$wvl
ref_calib_spectra_resamp <- resample(ref_calib_spectra, new_wvls = wvls)
plot(ref_calib_spectra_resamp)

# convert to long
ref_calib_spectra_resamp_long <- ref_calib_spectra_resamp %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  select(-sample_name) %>% 
  gather(key = wvl, value = refl) %>% 
  mutate(wvl = as.double(wvl))

# download previous calibrated reflectance of target panel
tar_panel_id <- panel_id
tar_calib <- read.table('/Users/etienne/Documents/ProjetsRecherche/CABO/panel-calibrations/panels/99AA05-0518-0226/2018-06-27/SRS-99-010-99AA05-0518-0226-2018-06-27.calib') %>% 
  as_data_frame() %>% 
  rename(wvl = V1, refl = V2)
tar_calib

# Plot the spectra, after cleaning
spec_clean <- filter(spec_all, cleaning == 'after')
spectra_plot_after <- ggplot(spec_clean, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = spec_fac, linetype = type, group = filename)) +
  facet_wrap(~ group) +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance') +
  coord_cartesian(ylim = c(.90, 1.10))
spectra_plot_after
ggsave('spectra_after_cleaning.pdf', plot = spectra_plot_after,
       width = 8, height = 5)

# Take average pf all reps
spec_clean_avg <- spec_clean %>% 
  filter(type == 'target') %>% 
  group_by(group, wvl) %>% 
  summarise_at(vars(refl), mean) %>% 
  group_by(wvl) %>% 
  summarise_at(vars(refl), mean)

# Plot average clean spectra
spec_clean_avg_plot <- ggplot(spec_clean_avg, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance')
spec_clean_avg_plot

# combine the different spectra
spec_clean_df <- spec_clean_avg %>% 
  rename(rel_refl_clean = refl)
spec_dirty_df <- spec_dirty_avg %>% 
  rename(rel_refl_dirty = refl)
tar_calib_df <- tar_calib %>% 
  rename(abs_refl_tar = refl)
ref_calib_df <- ref_calib_spectra_resamp_long %>% 
  rename(abs_refl_ref = refl)
spec_comb <- spec_clean_df %>% 
  left_join(spec_dirty_df) %>% 
  left_join(ref_calib_df) %>% 
  na.omit()

# add absolute reflectances
spec_comb <- spec_comb %>% 
  mutate(abs_ref_clean = rel_refl_clean * abs_refl_ref,
         abs_ref_dirty = rel_refl_dirty * abs_refl_ref) %>% 
  select(-rel_refl_clean, -rel_refl_dirty, -abs_refl_ref) %>% 
  rename(clean = abs_ref_clean, dirty = abs_ref_dirty)

# Gather df
spec_comb_long <- spec_comb %>%
  gather(key = calib, value = refl, -wvl)

# plot them
spec_comb_plot <- ggplot(spec_comb_long, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = calib)) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance')
spec_comb_plot


# Linear interpolation and smooth
new_wav <- 339:2499
interpol <- function(x) {
  x <- na.omit(x)
  inter <- approx(x$wvl, x$refl, new_wav)
  inter_smoothed <- smooth.spline(inter$y, spar = 0.5)$y
  tmp <- data_frame(wvl = inter$x, refl = inter$y, refl_smooth = inter_smoothed)
  return(tmp)
}

spec_interpol <- spec_comb_long %>%
  group_by(calib) %>% 
  do(interpol(.))

# Plot the clean and dirty spectra
spec_interpol_long <- spec_interpol %>% 
  rename(raw = refl, smooth = refl_smooth) %>% 
  gather(key = type, value = refl, -calib, -wvl) %>% 
  ungroup()

# Plot it
spec_interpol_plot <- ggplot(spec_interpol_long, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = type)) +
  facet_wrap(~ calib) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance')
spec_interpol_plot
ggsave('spec_interpol_plot.png', plot = spec_interpol_plot, width = 7, height = 3, dpi = 600)


# Plot all (including initial)
tar_calib <- tar_calib %>% 
  mutate(calib = '2018-06-27') %>% 
  select(calib, wvl, refl)
spec_three <- spec_interpol_long %>% 
  filter(type == 'smooth') %>% 
  select(-type) %>% 
  bind_rows(tar_calib) %>% 
  mutate(Calibration = calib,
         Calibration = replace(Calibration, Calibration == 'clean', '2018-11-20 (cleaned)'),
         Calibration = replace(Calibration, Calibration == 'dirty', '2018-11-13 (dirty)'))

# Plot it
spec_all_plot <- ggplot(spec_three, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = Calibration)) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance')
spec_all_plot
ggsave('99AA05-0518-0226_spec_all_plot.png', plot = spec_all_plot, width = 6, height = 3, dpi = 600)


# Create a matrix to store smoothed reflectance curve
calib_df <- spec_three %>% 
  filter(calib == 'clean')
calib_smooth_mat <- matrix(c(calib_df$wvl, signif(calib_df$refl, 4)), byrow = F, ncol = 2, nrow = nrow(calib_df))
smooth_file_name <- paste0('SRS-99-010-', tar_panel_id, '-', calibration_date, '.calib')

# save a new calib file to google drive
write.table(file = smooth_file_name,
            calib_smooth_mat,
            quote = F,
            row.names = F,
            col.names = F,
            sep = '\t')

# Upload to Google Drive folder
parent_path <- '/CABO/DATA/SPECTROSCOPY/PANELS/'
path_name <- paste0(parent_path, tar_panel_id, '/', calibration_date, '/')
gs_auth(new_user = T) # cabo-science
drive_upload(smooth_file_name, path = path_name)
drive_upload('99AA05-0518-0226_spec_all_plot.png', path = path_name)
drive_upload('spec_interpol_plot.png', path = path_name)
