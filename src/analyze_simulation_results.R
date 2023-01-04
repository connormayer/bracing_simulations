library(tidyverse)
library(lme4)

#setwd("C:/Users/conno/git_repos/bracing_simulations")
setwd("E:/git_repos/bracing_simulations/")


five_data_folder <- "data/sim3_12_19"
ten_data_folder <- "data/sim5_12_22"

contacts_df <- tibble()
excitation_df <- tibble()
failures_df <- tibble()

contacts_headers <- c(
  "sim", "front", "mid", "back", "lat_left", "lat_right", "coronal", "condition"
)
excitation_headers <- c(
  "sim", "GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL", "IL", "JAW.OPEN", "JAW.CLOSE", "condition"
)
failures_headers <- c(
  "GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL", "IL", "JAW.OPEN", "JAW.CLOSE", "condition"
)

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

colnames(contacts_df) <- contacts_headers
colnames(excitation_df) <- excitation_headers
colnames(failures_df) <- failures_headers

success_df <- inner_join(contacts_df, excitation_df, by = c("sim", "condition"))

contacts <- success_df %>%
  filter(front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0)

bracing <- contacts %>%
  filter(lat_left > 0 & lat_right > 0)

mean_activations <- bracing %>%
  select(condition, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL) %>%
  group_by(condition) %>%
  summarise(across(c(GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL), mean, na.rm=TRUE))

mean_contacts <- bracing %>%
  select(condition, front, mid, back, lat_left, lat_right, coronal) %>%
  group_by(condition) %>%
  summarise(across(c(front, mid, back, lat_left, lat_right, coronal), mean, na.rm=TRUE))

sd_activations <- bracing %>%
  select(condition, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL) %>%
  group_by(condition) %>%
  summarise(across(c(GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL), sd, na.rm=TRUE))

pivot_mean <- mean_activations %>%
  pivot_longer(-condition, names_to="muscle", values_to="mean")
pivot_sd <- sd_activations %>%
  pivot_longer(-condition, names_to="muscle", values_to="sd")

pivot_data <- inner_join(pivot_mean, pivot_sd, by=c("condition", "muscle"))
pivot_data <- pivot_data %>%
  mutate(type=ifelse(muscle %in% c("GGP", "GGM", "MH", "VERT", "SL"), 'Bracing Agonist', 'Bracing Antagonist'))

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

# Compare activation by bracer/non-bracer category
cat_data <- bracing %>%
  select(condition, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL)

cat_pivot <- cat_data %>%
  pivot_longer(-condition, names_to="muscle", values_to="activation") %>%
  mutate(type=ifelse(muscle %in% c("GGP", "GGM", "MH", "VERT", "SL"), 'agonist', 'antagonist'))

print("Mean activation across conditions by agonist/antagonist")
cat_pivot %>% group_by(condition, type) %>% summarize(mean=mean(activation))

# Stats
activation_df_stats <- success_df %>%
  mutate(
    contact = front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0, 'contact', 'no contact',
    bracing = lat_left > 0 & lat_right > 0, 'bracing', 'no bracing'
  ) %>%
  select(condition, contact, bracing, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL) 

activation_df_stats$condition <- factor(activation_df_stats$condition, levels=c("5mm", "10mm"))
m <- glm(bracing ~ (GGP + GGM + GGA + STY + MH + HG + TRANS + VERT + SL + IL) * condition, data = activation_df_stats, family = 'binomial')
summary(m)

# Check if there were 10mm successes that weren't 5mm successes
five_contacts <- contacts %>% filter(condition == '5mm')
ten_contacts <- contacts %>% filter(condition == '10mm')

print(nrow(inner_join(five_contacts, ten_contacts, by=c("GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL","IL"))))
print(nrow(ten_contacts))

five_bracing <- bracing %>% filter(condition == '5mm')
ten_bracing <- bracing %>% filter(condition == '10mm')

print(nrow(inner_join(five_bracing, ten_bracing, by=c("GGP", "GGM", "GGA", "STY", "MH", "HG", "TRANS", "VERT", "SL","IL"))))
print(nrow(ten_bracing))

# Check mean activation level across both conditions
mean(as.matrix(five_bracing[,9:18]))
mean(as.matrix(ten_bracing[,9:18]))

# Visualization

# Histogram of activation level counts

activation_df <- success_df %>%
  mutate(
    contact = ifelse(front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0, 'contact', 'no contact'),
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

