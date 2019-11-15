##########################################
#### Calibration files for spectralon ####
##########################################

library(googledrive)
library(tidyverse)
library(spectrolab)
library(googlesheets)
source("./R_scripts/CABO_processing_functions.R")

# Panel ID
panel_id <- "99AA02-0318-9711"
calibration_date <- "2019-06-25"
panel <- "9711"

# Read File Info
meta <- read_csv("./Processing_chain/order_measurements/order_measurements.csv") %>% 
  arrange(prefix, spectrum_number)
  meta
  
# get file names
meta <- meta %>%
  mutate(filename = paste0(prefix, '_', spectrum_number, '.sig'))


# download all files from working folder to local directory (paste path to google drive folder)
files <- drive_ls("https://drive.google.com/drive/folders/18emY__sRwL0-JeIannvnSJH3t6YSzsZa") %>% 
  filter(name %in% meta$filename)

# mypath <- paste0(getwd(), '/spectra')
dir.create(paste0("./Processing_chain/",panel,"_panel/R_input/",calibration_date,"_spec"))

for (i in 1:nrow(files)){
  drive_download(files[i,], path = paste0("./processing_chain/",panel,"_panel/R_input/",
                                          calibration_date,"_spec/", files$name[i]))
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
filenames <- list.files(paste0("./Processing_chain/",
                               panel,"_panel/R_input/",
                               calibration_date,"_spec"))
paths <- paste0("./Processing_chain/",panel,"_panel/R_input/",calibration_date,
                "_spec/",filenames)

spec_data_all <- lapply(paths, get_sig_data) 
names(spec_data_all) <- filenames

# convert list to long data frame with filename column
for (i in 1:length(spec_data_all)) spec_data_all[[i]]$filename <- names(spec_data_all)[i]
spec_data_all_df <- bind_rows(spec_data_all) %>% as_data_frame() %>% 
  mutate(refl = refl / 100)

# merge with metadata
meta_sub <- meta %>% 
  select(filename, spectrum_number, type, group, spec_num, cleaning,rep) %>% 
  mutate(spec_fac = factor(spec_num))
spec_all <- meta_sub %>% 
  left_join(spec_data_all_df, by = 'filename')

# check detector overlap
all_wvls <- filter(spec_all, filename == 'gr062519_0000.sig')$wvl
plot(all_wvls, type = 'l')
abline(v=524)

# find detector overlap regions
all_wvls[512:525] # remove 513:524
rem1 <- all_wvls[513:524]

all_wvls[765:780] # remove 769:773
rem2 <- all_wvls[769:774]

# 
# all_wvls[600] # remove 769:773
# rem2 <- all_wvls[769:774]
# rem <- c(rem1, rem2)

rem <- c(rem1, rem2)

# remove those wavelengths
spec_all <- spec_all %>% 
  filter(!wvl %in% rem)

#### Test if wvls are increasing 
test_wvls <- filter(spec_all, filename == 'gr062519_0000.sig')$wvl
all(diff(test_wvls) >= 0)
test_wvls[which(diff(test_wvls) < 0)]


# download calibrated reflectance of reference disk 99AA02-0318-9710
ref_panel_id <- '99AA02-0318-9710'
panel_file <-drive_ls("https://drive.google.com/drive/folders/1MHQBy6P2qY67li4y59stmPhHaHnSxTEx") %>% 
  filter(grepl(".calib$", name))

for (i in 1:nrow(panel_file)){
  drive_download(panel_file[i,], 
                 path = paste0("./Processing_chain/",panel,"_panel/docu/",calibration_date,"/",
                               panel_file$name[i]))
}

ref_calib <- read.table(paste0("./Processing_chain/",panel,"_panel/docu/",calibration_date,"/", panel_file$name[1])) %>%
  as_data_frame() %>%
  rename(wvl = V1, refl = V2)

### or 
# ref_calib <- read.table("./Processing_chain/lab_panel/SRS-99-010-99AA02-0318-9710-2018-04-10.calib") %>%
#   as_data_frame() %>%
#   rename(wvl = V1, refl = V2)
ref_calib

# convert to wide format and to spectra
ref_calib_wide <- ref_calib %>%
  spread(key = wvl, value = refl) %>%
  mutate(SN = ref_panel_id) %>%
  as.data.frame()
ref_calib_spectra <- as.spectra(ref_calib_wide, name_idx = ncol(ref_calib_wide))
plot(ref_calib_spectra)

# resample to same wvls as target spectra
wvls <- filter(spec_all[spec_all$spectrum_number=="0000",], wvl <= 2500)$wvl
ref_calib_spectra_resamp <- resample(ref_calib_spectra, new_wvls = wvls)
plot(ref_calib_spectra_resamp)

# convert to long
ref_calib_spectra_resamp_long <- ref_calib_spectra_resamp %>%
  as.data.frame() %>%
  as_tibble() %>%
  select(-sample_name) %>%
  gather(key = wvl, value = refl) %>%
  mutate(wvl = as.double(wvl))

## make sure funcion get_refl_trans_LL() is sources
panels_calibs <- ref_calib_spectra_resamp_long

out <- get_refl_trans_LL(spec_all)

outi <- out %>%
  left_join(spec_all)%>%
  select(-refl,-ref,-tar)

outi$rep <- as.factor(outi$rep)

# Plot the spectra, after cleaning
spec_clean <- filter(outi, cleaning == 'after')

spectra_plot_after <- ggplot(spec_clean, aes(x = wvl, y = value)) +
  geom_line(aes(colour = rep, linetype = type, group = filename)) +
  facet_wrap(~ group) +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance') +
  coord_cartesian(ylim = c(.90, 1.10))
spectra_plot_after
ggsave(paste0("./Processing_chain/", panel,"_panel/",
              calibration_date,"/spectra_after_cleaning.pdf"), plot = spectra_plot_after,
       width = 8, height = 5)

# Take average of all reps
spec_clean_avg <- spec_clean %>% 
  filter(type == 'target') %>% 
  group_by(group, wvl) %>% 
  summarise_at(vars(value), mean) %>% 
  group_by(wvl) %>% 
  summarise_at(vars(value), mean)

# Plot average clean spectra
spec_clean_avg_plot <- ggplot(spec_clean_avg, aes(x = wvl, y = value)) +
  geom_line() +
  xlab('Wavelength (nm)') +
  ylab('Relative reflectance')
spec_clean_avg_plot

# Linear interpolation and smooth
tail(spec_clean_avg$wvl) ### round up 
new_wav <- 339:2498
interpol <- function(x) {
  x <- na.omit(x)
  inter <- approx(x$wvl, x$value, new_wav)
  inter_smoothed <- smooth.spline(inter$y, spar = 0.6)$y
  tmp <- data_frame(wvl = inter$x, refl = inter$y, refl_smooth = inter_smoothed)
  return(tmp)
}

spec_interpol <- interpol(spec_clean_avg)

# Plot the clean and dirty spectra
spec_interpol_long <- spec_interpol %>% 
  rename(raw = refl, smooth = refl_smooth) %>% 
  gather(key = type, value = refl, -wvl) %>% 
  ungroup()

# Plot it
spec_interpol_plot <- ggplot(spec_interpol_long, aes(x = wvl, y = refl)) +
  geom_line(aes(colour = type)) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance')
spec_interpol_plot

ggsave(paste0("./Processing_chain/", panel,"_panel/",
              calibration_date,"/spec_interpol_plot.png"), 
       plot = spec_interpol_plot, width = 4.5, height = 3, dpi = 600)

# Create a matrix to store smoothed reflectance curve
calib_df <- spec_interpol_long %>% 
  filter(type == "smooth") %>% 
  select(-type) 
calib_smooth_mat <- matrix(c(calib_df$wvl, signif(calib_df$refl, 4)), 
                           byrow = F, ncol = 2, nrow = nrow(calib_df))
smooth_file_name <- paste0("./Processing_chain/", panel,"_panel/",
       calibration_date,"/SRS-99-010-", panel_id, '-', calibration_date, '.calib')

paste0("./Processing_chain/", panel,"_panel/",
       calibration_date,"/spec_interpol_plot.png")


# save a new calib file to google drive
write.table(file = smooth_file_name,
            calib_smooth_mat,
            quote = F,
            row.names = F,
            col.names = F,
            sep = '\t')
