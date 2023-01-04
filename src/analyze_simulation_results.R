library(tidyverse)
library(lme4)

# Change this as necessary
setwd("E:/git_repos/bracing_simulations/")

five_data_folder <- "data/sim3_12_19"
ten_data_folder <- "data/sim5_12_22"

# Simulation results are split over three types of files.
# * contacts files record tongue/palate contact at the end of each
#   successful simulation
# * excitation files record the muscle excitation at the end of each successful
#   simulation
# * failure files record the muscle excitation at the end of each failed
#   simulation
# Each worker (a computer process that runs simulations in parallel with other
# workers) has its own version of each of these files.
contacts_df <- tibble()
excitation_df <- tibble()
failures_df <- tibble()

# These files don't include headers, so we define our own
contacts_headers <- c(
  "sim", "front", "mid", "back", "lat_left", "lat_right", "coronal", "condition"
)
excitation_headers <- c(
  "sim", "GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL", "IL",
  "JAW.OPEN", "JAW.CLOSE", "condition"
)
failures_headers <- c(
  "GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL", "IL",
  "JAW.OPEN", "JAW.CLOSE", "condition"
)

# Go through every file in the 5mm folder and plug it into the appropriate dataframe
for (f in list.files(five_data_folder)) {
  full_path <- file.path(five_data_folder, f)
  data <- read_csv(full_path, col_names=FALSE,col_types = cols())
  data$condition <- '5mm'
  if (str_detect(f, "contacts")) {
    contacts_df <- rbind(contacts_df, data)
  } else if (str_detect(f, "failed")) {
    failures_df <- rbind(failures_df, data)
  } else if (str_detect(f, "excitations")) {
    excitation_df <- rbind(excitation_df, data)
  }
}

# Ditto for the 10mm folder
for (f in list.files(ten_data_folder)) {
  full_path <- file.path(ten_data_folder, f)
  data <- read_csv(full_path, col_names=FALSE,col_types = cols())
  data$condition <- '10mm'
  if (str_detect(f, "contacts")) {
    contacts_df <- rbind(contacts_df, data)
  } else if (str_detect(f, "failed")) {
    failures_df <- rbind(failures_df, data)
  } else if (str_detect(f, "excitations")) {
    excitation_df <- rbind(excitation_df, data)
  }
}

# Add our headers to the DF
colnames(contacts_df) <- contacts_headers
colnames(excitation_df) <- excitation_headers
colnames(failures_df) <- failures_headers

# Combine the contacts and excitation dfs to get one containing all the
# successful simulation results
success_df <- inner_join(contacts_df, excitation_df, by = c("sim", "condition"))

# Keep only simulations that resulted in tongue/palate contact
contacts <- success_df %>%
  filter(front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0)

# Keep only simulations that resulted in bilateral bracing
bracing <- contacts %>%
  filter(lat_left > 0 & lat_right > 0)

# Mean activation of each muscle in bracing outcomes
mean_activations <- bracing %>%
  select(condition, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL) %>%
  group_by(condition) %>%
  summarise(across(c(GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL), mean, na.rm=TRUE))

# Mean number of lateral contacts in bracing outcomes
mean_contacts <- bracing %>%
  select(condition, front, mid, back, lat_left, lat_right, coronal) %>%
  group_by(condition) %>%
  summarise(across(c(front, mid, back, lat_left, lat_right, coronal), mean, na.rm=TRUE))

# Make a df in a format more amenable to plotting
pivot_data <- mean_activations %>%
  pivot_longer(-condition, names_to="muscle", values_to="mean") %>%
  mutate(type=ifelse(muscle %in% c("GGP", "GGM", "MH", "VERT", "SL"),
                     'Bracing Agonist', 'Bracing Antagonist'))

# Print some useful statistics
print("Successful sims")
print(success_df %>% group_by(condition) %>% count())

print("Contact outcomes")
print(contacts %>% group_by(condition) %>% count())

print("Bracing outcomes")
print(bracing %>% group_by(condition) %>% count())

print("Mean activations")
print(mean_activations)

print("Mean contacts")
print(mean_contacts)

# Compare activation by agonist/antagonist category
cat_data <- bracing %>%
  select(condition, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL)

cat_pivot <- cat_data %>%
  pivot_longer(-condition, names_to="muscle", values_to="activation") %>%
  mutate(type=ifelse(muscle %in% c("GGP", "GGM", "MH", "VERT", "SL"), 'agonist', 'antagonist'))

print("Mean activation across conditions by agonist/antagonist")
cat_pivot %>% group_by(condition, type) %>% summarize(mean=mean(activation))

# Fit a logistic regression model
# Tidy up the data a bit
activation_df_stats <- success_df %>%
  mutate(
    contact = front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0, 'contact', 'no contact',
    bracing = lat_left > 0 & lat_right > 0, 'bracing', 'no bracing'
  ) %>%
  select(condition, contact, bracing, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL)
activation_df_stats$condition <- factor(activation_df_stats$condition, levels=c("5mm", "10mm"))

# Fit the model
m <- glm(bracing ~ (GGP + GGM + GGA + STY + MH + HG + TRANS + VERT + SL + IL) * condition,
         data = activation_df_stats, family = 'binomial')
summary(m)

# Check if there were 10mm successes that weren't 5mm successes
# Are there contacts in 10mm that aren't in 5mm
five_contacts <- contacts %>% filter(condition == '5mm')
ten_contacts <- contacts %>% filter(condition == '10mm')

print(nrow(inner_join(five_contacts, ten_contacts,
                      by=c("GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL","IL"))))
print(nrow(ten_contacts))

# Are there bracing outcomes in 10mm that weren't in 5mm
five_bracing <- bracing %>% filter(condition == '5mm')
ten_bracing <- bracing %>% filter(condition == '10mm')

print(nrow(inner_join(five_bracing, ten_bracing,
                      by=c("GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL","IL"))))
print(nrow(ten_bracing))

# Visualization

# Histogram of activation level counts. This isn't used in the paper
activation_df <- success_df %>%
  mutate(
    contact = ifelse(front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0,
                     'contact', 'no contact'),
    bracing = ifelse(lat_left > 0 & lat_right > 0, 'bracing', 'no bracing'),
    group = str_c(contact, bracing, sep = ' & ')
  ) %>%
  select(condition, contact, bracing, group, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL)

activation_df %>%
  filter(bracing == 'bracing') %>%
  mutate(condition = factor(condition, levels = c("5mm", "10mm"))) %>%
  gather(key, value, -condition, -contact, -bracing, -group) %>%
  ggplot(aes(x=value, fill=condition)) +
  geom_bar(aes(y=..prop..), position=position_dodge(preserve = "single")) +
  facet_wrap(~key, nrow=2) +
  xlab("Activation level") +
  ylab("Proportion of bracing outcomes")
ggsave("images/activation_histogram.png")

# Barplot of mean levels
pivot_data %>%
  mutate(condition = factor(condition, levels = c("5mm", "10mm"))) %>%
  ggplot(aes(x=factor(muscle,
                      levels = c("HG", "TRANS", "STY", "GGA", "IL", "VERT", "GGM", "MH", "SL", "GGP")),
             y=mean, fill=condition)) +
  geom_bar(position=position_dodge(preserve = "single"), stat='identity') +
  xlab("Muscle") +
  ylab("Mean activation level for bilateral bracing outcomes") +
  scale_fill_discrete(name = "Jaw Opening") +
  facet_grid(~ type, scale='free') +
  theme(text = element_text(size=20))
ggsave("images/activation_barplot.png")

