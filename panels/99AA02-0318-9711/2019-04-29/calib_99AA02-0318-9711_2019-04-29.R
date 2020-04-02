##########################################
#### Calibration files for spectralon ####
##########################################

library(googledrive)
library(tidyverse)
library(spectrolab)
library(googlesheets4)

# Panel ID
panel_id <- "99AA02-0318-9711"
calibration_date <- "2019-04-29" ## this is the date post-cleaning measurements where made
precalibration_date <- "2019-04-23" # this is the date pre-cleaning measurements where made
panel <- "9711"
prefix_wd <- paste0("./panels/", panel_id, "/",calibration_date)

prefix_before <- paste0(precalibration_date,"_before_cleaning/gr042319") 
prefix_after <- paste0(calibration_date,"_after_cleaning/gr042919")

### Currently Google Drive download does not want to work
### Instead: Download .sig files manually and add them to the before/after cleaning folders 

dir.create(paste0(prefix_wd, "/spectra/",precalibration_date,"_before_cleaning"))
dir.create(paste0(prefix_wd, "/spectra/",calibration_date,"_after_cleaning"))

# download all files from working folder to local directory (paste path to google drive folder)
 # files <- drive_ls("https://drive.google.com/drive/u/1/folders/1Yk4Q5J96QcUtG0nMT34iCbz0jXQlX-Y3") %>%
 #   filter(name %in% meta$filename)
 # 
 # for (i in 1:nrow(files)){
 #   drive_download(files[i,], path = paste0(prefix_wd, "/spectra/", files$name[i]))
 #   }

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
filenames <- list.files(paste0(prefix_wd, "/spectra"),recursive = T,pattern = ".sig")
paths <- paste0(prefix_wd, "/spectra/",filenames)
filenames
spec_data_all <- lapply(paths, get_sig_data) 
names(spec_data_all) <- filenames

# convert list to long data frame with filename column
for (i in 1:length(spec_data_all)) spec_data_all[[i]]$filename <- names(spec_data_all)[i]
spec_data_all_df <- bind_rows(spec_data_all) %>% as_tibble() %>% 
  mutate(refl = refl / 100)

# Check readme.txt files in folders with .sig files 
# make sure that order of measurements (before + after cleaning) is the same
# if not: re-create file
meta <- read_csv(paste0("./panels/order_measurements/order_measurements.csv"))

# make file names before
meta_before <- meta %>%
  mutate(filename = paste0(prefix_before, '_', spectrum_number, '.sig'))

# make file names after
meta_after <- meta %>%
  mutate(filename = paste0(prefix_after, '_', spectrum_number, '.sig'))

meta_x <- rbind(meta_before,meta_after)

# merge spectra with metadata
meta_sub <- meta_x %>% 
  select(filename, spectrum_number, type, group, spec_num,rep) %>% 
  mutate(spec_cont= factor(spec_num))

table(meta_sub$filename %in% spec_data_all_df$filename)

spec_all <- meta_sub %>% 
  left_join(spec_data_all_df, by = 'filename')%>% 
  mutate(cleaning=sapply(strsplit(spec_data_all_df$filename,"_"),"[",2))

# check detector overlap
all_wvls <- filter(spec_all, filename == spec_all$filename[1])$wvl
plot(all_wvls, type = 'l')
abline(v=c(513,524))

# find detector overlap regions
all_wvls[512:525] # remove 513:524
rem1 <- all_wvls[513:523]

all_wvls[768:775] # remove 769:774
rem2 <- all_wvls[769:773]

rem <- c(rem1, rem2)

# remove those wavelengths
spec_all <- spec_all %>% 
  filter(!wvl %in% rem)

#### Test if wvls are increasing tester si différences entre les longueurs d'ondes tjrs>0 pcq on vx pas retourner en arrère (loop)
test_wvls <- filter(spec_all, filename == spec_all$filename[1])$wvl
all(diff(test_wvls) >= 0)
test_wvls[which(diff(test_wvls) < 0)] #check for repeated/smaller wavelengts

if (!(all(diff(test_wvls) >= 0))){### we want this to be TRUE
  print("Careful! Duplicated wvls still present")
} else ("Great! Duplicated wvls removed")

# Plot the spectra, before cleaning
spec_dirty <- filter(spec_all, cleaning == 'before', !(type=="stray"))
spectra_plot_before <- ggplot(spec_dirty, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = spec_cont, linetype = type, group = filename)) +
  facet_wrap(~ group) +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance') +
  coord_cartesian(ylim = c(.90, 1.10))
spectra_plot_before

ggsave(paste0(prefix_wd,'/spectra_before_cleaning.pdf'), plot = spectra_plot_before,
       width = 8, height = 5)

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

# read calibrated reflectance of reference disk 99AA02-0318-9710
ref_panel_id <- '99AA02-0318-9710'

ref_calib <- read.table(paste0("./panels/",ref_panel_id, "/2018-04-10/SRS-99-010-",
                               ref_panel_id,"-2018-04-10.calib")) %>%
  as_tibble() %>%
  rename(wvl = V1, refl = V2)

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

# load previous calibrated reflectance of target panel
tar_ls <- list.files(paste0("./panels/", panel_id, "/"),
                     pattern = ".calib",recursive = T, full.names = T)

tar_ls <- tar_ls[2:6]

tar_calib <- read.table(tar_ls[1]) %>% 
  as_tibble() %>% 
  rename(wvl = V1, refl = V2)
tar_calib

tar_calib2 <- read.table(tar_ls[2]) %>% 
  as_tibble() %>% 
  rename(wvl = V1, refl = V2)

tar_calib3 <- read.table(tar_ls[3]) %>% 
  as_tibble() %>% 
  rename(wvl = V1, refl = V2)

tar_calib4 <- read.table(tar_ls[4]) %>% 
  as_tibble() %>% 
  rename(wvl = V1, refl = V2)

tar_calib5 <- read.table(tar_ls[5]) %>% 
  as_tibble() %>% 
  rename(wvl = V1, refl = V2)

# Plot the spectra, after cleaning
spec_clean <- filter(spec_all, cleaning == 'after',!(type=="stray")) #we don't want to plot stray light
spectra_plot_after <- ggplot(spec_clean, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = spec_cont, linetype = type, group = filename)) +
  facet_wrap(~ group) +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance') +
  coord_cartesian(ylim = c(.90, 1.10))
spectra_plot_after
ggsave(paste0(prefix_wd,'/spectra_after_cleaning.pdf'), plot = spectra_plot_after,
       width = 8, height = 5)

# Take average 
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

# Calculate average straylight spectra dirty
stray_dirty <- filter(spec_all, cleaning == 'before', type=="stray")

# Take average
stray_dirty_avg <- stray_dirty %>% 
  group_by(group, wvl) %>% 
  summarise_at(vars(refl), mean) %>% 
  group_by(wvl) %>% 
  summarise_at(vars(refl), mean)

# Plot average dirty stray spectra
ggplot(stray_dirty_avg, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance')

# Calculate average straylight spectra clean
stray_clean <- filter(spec_all, cleaning == 'after', type=="stray")

# Take average
stray_clean_avg <- stray_clean %>% 
  group_by(group, wvl) %>% 
  summarise_at(vars(refl), mean) %>% 
  group_by(wvl) %>% 
  summarise_at(vars(refl), mean)

# Plot average stray clean spectra
ggplot(stray_clean_avg, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance')

# Average clean reference
ref_clean_avg <- spec_clean %>% 
  filter(type == 'ref') %>% 
  group_by(group, wvl) %>% 
  summarise_at(vars(refl), mean) %>% 
  group_by(wvl) %>% 
  summarise_at(vars(refl), mean)

# Plot average clean spectra ref
ggplot(ref_clean_avg, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance')

# Average dirty reference
ref_dirty_avg <- spec_dirty %>% 
  filter(type == 'ref') %>% 
  group_by(group, wvl) %>% 
  summarise_at(vars(refl), mean) %>% 
  group_by(wvl) %>% 
  summarise_at(vars(refl), mean)

# Plot average dirty spectra reference
ggplot(ref_dirty_avg, aes(x = wvl, y = refl)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance')

##################################
# Combine the different spectra
spec_clean_df <- spec_clean_avg %>% 
  rename(rel_refl_clean = refl)
spec_dirty_df <- spec_dirty_avg %>% 
  rename(rel_refl_dirty = refl)
ref_clean_df <- ref_clean_avg %>% 
  rename(ref_clean = refl)
ref_dirty_df <- ref_dirty_avg %>% 
  rename(ref_dirty = refl)
ref_calib_df <- ref_calib_spectra_resamp_long %>% 
  rename(abs_refl_ref = refl)
stray_clean_df <- stray_clean_avg %>% 
  rename(stray_clean = refl)
stray_dirty_df <- stray_dirty_avg %>% 
  rename(stray_dirty = refl)

tar_calib_df <- tar_calib %>% 
  rename(abs_refl_tar = refl)
tar_calib2_df <- tar_calib2 %>% 
  rename(abs_refl_tar2 = refl)
tar_calib3_df <- tar_calib3 %>% 
  rename(abs_refl_tar3 = refl)
tar_calib4_df <- tar_calib4 %>% 
  rename(abs_refl_tar4 = refl) 
tar_calib5_df <- tar_calib5 %>% 
  rename(abs_refl_tar5 = refl) 

spec_comb <- spec_clean_df %>% 
  left_join(spec_dirty_df) %>% 
  left_join(ref_clean_df) %>% 
  left_join(ref_dirty_df) %>% 
  left_join(stray_clean_df) %>% 
  left_join(stray_dirty_df) %>% 
  left_join(ref_calib_df) %>% 
  left_join(stray_clean_df) %>% 
  na.omit()

 # add absolute reflectances
spec_comb <- spec_comb %>% 
  mutate(abs_refl_clean = (rel_refl_clean-stray_clean)/(ref_clean-stray_clean) * abs_refl_ref,
         abs_refl_dirty = (rel_refl_dirty-stray_dirty)/(ref_dirty-stray_dirty) * abs_refl_ref) %>% 
  select(-rel_refl_clean, -rel_refl_dirty, -abs_refl_ref) %>% 
  rename(clean = abs_refl_clean, dirty = abs_refl_dirty)

ggplot(spec_comb, aes(x = wvl, y = clean)) +
  geom_line() +
  xlab('Wavelength (nm)') 

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
tail(spec_comb_long$wvl) ### round up 
head(spec_comb)
new_wav <- 339:2498
interpol <- function(x) {
  x <- na.omit(x)
  inter <- approx(x$wvl, x$refl, new_wav)
  inter_smoothed <- smooth.spline(inter$y, spar = 0.5)$y
  tmp <- tibble(wvl = inter$x, refl = inter$y, refl_smooth = inter_smoothed)
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
spec_interpol_plot <- ggplot(subset(spec_interpol_long,calib=="clean" |calib=="dirty"), aes(x = wvl, y = refl)) +
  geom_line(aes(colour = type)) +
  facet_wrap(~ calib) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance')
spec_interpol_plot

ggsave(paste0(prefix_wd,'/spec_interpol_plot.png'), 
       plot = spec_interpol_plot,
       width = 7, height = 3, dpi=600)

# Plot all (including initial)
tar_ls
tar_calib <- tar_calib %>% 
  mutate(calib = '2018-07-19') %>% 
  select(calib, wvl, refl)
tar_calib2 <- tar_calib2 %>% 
  mutate(calib = '2018-11-13') %>% 
  select(calib, wvl, refl)
tar_calib3 <- tar_calib3 %>% 
  mutate(calib = '2019-05-27') %>% 
  select(calib, wvl, refl)
tar_calib4 <- tar_calib4 %>% 
  mutate(calib = '2019-06-25') %>% 
  select(calib, wvl, refl)
tar_calib5 <- tar_calib5 %>% 
  mutate(calib = '2019-07-22') %>% 
  select(calib, wvl, refl)

spec_three <- spec_interpol_long %>%  ##here: operation needed: uptade clean and dirty date
  filter(type == 'smooth') %>% 
  filter(calib == 'clean'| calib == 'dirty') %>% 
  select(-type) %>% 
  bind_rows(tar_calib, tar_calib2, tar_calib3, tar_calib4, tar_calib5) %>% 
  mutate(Calibration = calib,
         Calibration = replace(Calibration, Calibration == 'clean', '2019-07-22 (cleaned)'),
         Calibration = replace(Calibration, Calibration == 'dirty', '2019-07-19 (dirty)'))

# Plot it toutes les anciennes calibrations nettoyées + sale/propre de la calibration en cours
spec_all_plot <- ggplot(spec_three, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = Calibration)) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance')
spec_all_plot
ggsave(paste0(prefix_wd,'/spec_all_plot.png'), plot = spec_all_plot, 
       width = 6, height = 3, dpi = 600)

# Create a matrix to store smoothed reflectance curve
calib_df <- spec_three %>% 
  filter(calib == 'clean')
calib_smooth_mat <- matrix(c(calib_df$wvl, signif(calib_df$refl, 4)), byrow = F, ncol = 2, nrow = nrow(calib_df))
smooth_file_name <- paste0(prefix_wd,"/SRS-99-010-", panel_id, '-', calibration_date, '.calib')

# save a new calib file to google drive
write.table(file = smooth_file_name,
            calib_smooth_mat,
            quote = F,
            row.names = F,
            col.names = F,
            sep = '\t')

### Upload manually

# # Upload to Google Drive folder
# parent_path <- '/CABO/DATA/SPECTROSCOPY/PANELS/'
# path_name <- paste0(parent_path, "/",panel_id, '/', calibration_date, '/')
# gs_auth() # cabo-science
# 
# drive_upload(smooth_file_name, path = path_name)


