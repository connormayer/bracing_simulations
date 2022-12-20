library(tidyverse)
library(lme4)
library(ggbiplot)

setwd("E:/git_repos/bracing_simulations/")

five_data_folder <- "E:/git_repos/bracing_simulations/data/sim3_12_19"
ten_data_folder <- "E:/git_repos/bracing_simulations/data/sim4_12_20"

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

sd_activations <- bracing %>%
  select(condition, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL) %>%
  group_by(condition) %>%
  summarise(across(c(GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL), sd, na.rm=TRUE))

pivot_mean <- mean_activations %>%
  pivot_longer(-condition, names_to="muscle", values_to="mean")
pivot_sd <- sd_activations %>%
  pivot_longer(-condition, names_to="muscle", values_to="sd")

pivot_data <- inner_join(pivot_mean, pivot_sd, by=c("condition", "muscle"))

print("Successful sims")
print(success_df %>% group_by(condition) %>% count())

print("Contact outcomes")
print(contacts %>% group_by(condition) %>% count())

print("Bracing outcomes")
print(bracing %>% group_by(condition) %>% count())

print("Mean activations")
print(mean_activations)

# Visualization

# Dimensionality-reduced plot
activation_df <- success_df %>%
  mutate(
    contact = ifelse(front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0, 'contact', 'no contact'),
    bracing = ifelse(lat_left > 0 & lat_right > 0, 'bracing', 'no bracing'),
    group = str_c(contact, bracing, sep = ' & ')
  ) %>%
  select(condition, contact, bracing, group, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL)

activation_five <- activation_df %>% filter(condition == '5mm')
activation_ten <- activation_df %>% filter(condition == '10mm')

pc_five <- activation_five %>%
  select(-condition, -group, -contact, -bracing) %>%
  prcomp(center=TRUE, scale. = TRUE)

pc_ten <- activation_ten %>%
  select(-condition, -group, -contact, -bracing) %>%
  prcomp(center=TRUE, scale. = TRUE)

plot_five <- ggbiplot(pc_five, obs.scale = 1, var.scale = 1,
                      groups = activation_five$group) +
  labs(title="5mm jaw aperture")
plot_five$layers <- c(plot_five$layers, plot_five$layers[[1]])
ggsave("images/5mm_pca.png", plot=plot_five)

plot_ten <- ggbiplot(pc_ten, obs.scale = 1, var.scale = 1,
                     groups = activation_ten$group) +
  labs(title="10mm jaw aperture")
plot_ten$layers <- c(plot_ten$layers, plot_ten$layers[[1]])
ggsave("images/10mm_pca.png", plot=plot_ten)

activation_five <- activation_df %>% filter(condition == '5mm' & contact == 'contact')
activation_ten <- activation_df %>% filter(condition == '10mm' & contact == 'contact')

pc_five <- activation_five %>%
  select(-condition, -group, -contact, -bracing) %>%
  prcomp(center=TRUE, scale. = TRUE)

pc_ten <- activation_ten %>%
  select(-condition, -group, -contact, -bracing, -HG) %>%
  prcomp(center=TRUE, scale. = TRUE)

plot_five <- ggbiplot(pc_five, obs.scale = 1, var.scale = 1,
                      groups = activation_five$group) +
  labs(title="5mm jaw aperture")
plot_five$layers <- c(plot_five$layers, plot_five$layers[[1]])
ggsave("images/5mm_pca_contact.png", plot=plot_five)

plot_ten <- ggbiplot(pc_ten, obs.scale = 1, var.scale = 1,
                     groups = activation_ten$group) +
  labs(title="10mm jaw aperture")
plot_ten$layers <- c(plot_ten$layers, plot_ten$layers[[1]])
ggsave("images/10mm_pca_contact.png", plot=plot_ten)

# Histogram of activation level counts

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
  ggplot(aes(x=muscle, y=mean, fill=condition)) +
  geom_bar(position=position_dodge(preserve = "single"), stat='identity') +
  xlab("Muscle") +
  ylab("Mean activation level")
ggsave("images/activation_barplot.png")

# Stats?
activation_df_stats <- success_df %>%
  mutate(
    contact = front > 0 | mid > 0 | back > 0 | lat_left > 0 | lat_right > 0 | coronal > 0, 'contact', 'no contact',
    bracing = lat_left > 0 & lat_right > 0, 'bracing', 'no bracing'
  ) %>%
  select(condition, contact, bracing, GGP, GGM, GGA, STY, MH, HG, TRANS, VERT, SL, IL)

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
