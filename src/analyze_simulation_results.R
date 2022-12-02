library(tidyverse)

data_folder <- "E:/Dropbox/ling/research/bracing_paper/data/sim1_12_1"

contacts_df <- tibble()
excitation_df <- tibble()
failures_df <- tibble()

contacts_headers <- c(
  "sim", "front", "mid", "back", "lat_left", "lat_right", "coronal"
)
excitation_headers <- c(
  "sim", "GGP", "GGM", "GGA", "STY", "GH", "MH", "HG", "TRANS", "VERT", "SL", "IL"
)
failures_headers <- c(
  "GGP", "GGM", "GGA", "STY", "GH", "MH", "HG", "TRANS", "VERT", "SL", "IL"
)

for (f in list.files(data_folder)) {
 full_path <- file.path(data_folder, f)
 data <- read_csv(full_path, col_names=FALSE)

 if (str_detect(f, "contacts")) {
   contacts_df <- rbind(contacts_df, data)
 } else if (str_detect(f, "failed")) {
   failures_df <- rbind(failures_df, data)
 } else if (str_detect(f, "excitations")) {
   excitation_df <- rbind(excitation_df, data)
 }
}

colnames(contacts_df) <- contacts_headers
colnames(excitation_df) <- excitation_headers
colnames(failures_df) <- failures_headers

success_df <- inner_join(contacts_df, excitation_df, by = "sim")
contacts <- success_df %>%
  filter(
    front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0
  )

bracing <- contacts %>%
  filter(lat_left > 0 & lat_right > 0)

mean_activations <- bracing %>%
  select(GGP, GGM, GGA, STY, GH, MH, HG, TRANS, VERT, SL, IL) %>%
  colMeans()
