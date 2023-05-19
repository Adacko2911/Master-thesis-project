#---------------------------------------------------
# Summary statistics on empirical data
#---------------------------------------------------

library(tidyverse)
library(ggplot2)
library(foreign)


# Read data
#insert directory to the data (Note the data is not included as it is only available on inquiry)
bipolar_data <- foreign::read.spss("ESM bipolar cleaned data_UU project.sav", to.data.frame = TRUE)
bipolar_data <- bipolar_data %>%
  dplyr::select(patient_id, time,
                bs_diary_5, bs_diary_13, bs_diary_22,
                bs_diary_15, bs_diary_9, bs_diary_10,
                bs_diary_7, bs_diary_11, bs_diary_17,
                bs_diary_8, bs_diary_14, bs_diary_16) %>%
  drop_na() %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id()) %>%
  ungroup() %>%
  arrange(patient_id, time) %>%
  dplyr::select(-time) %>%
  as.matrix()

# Create summary statistics on empirical data to be used in PPCs
true_data <- bipolar_data %>%
  as.data.frame() %>%
  group_by(patient_id) %>%
  mutate(occasion = row_number()) %>%
  gather(variable, value, -patient_id, -occasion)

# Summary statistics on aggregated (group-level) data
group_statistics <- true_data %>%
  mutate(variable = factor(variable, levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                                "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                                "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                                "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                           labels = c("down","dread nrest of day","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "switch and focus",
                                      "extremely well", "full of ideas", "thoughts are racing"))) %>%
  group_by(variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            min_value = min(value, na.rm = TRUE),
            max_value = max(value, na.rm = TRUE),
            iqr_value = IQR(value, na.rm = TRUE),
            prop_zero = mean(value == 0, na.rm = TRUE))

saveRDS(group_statistics,"./example/outputs/group_statistics.rds")

# Summary statistics on patient-level data
patient_statistics <- true_data %>%
  mutate(variable = factor(variable,
                           levels = c("bs_diary_5", "bs_diary_13", "bs_diary_22",
                                      "bs_diary_15", "bs_diary_9", "bs_diary_10",
                                      "bs_diary_7", "bs_diary_11", "bs_diary_17",
                                      "bs_diary_8", "bs_diary_14", "bs_diary_16"),
                           labels = c("down","dread nrest of day","worry",
                                      "inadequate", "tired", "content",
                                      "agitated", "irritated", "switch and focus",
                                      "extremely well", "full of ideas", "thoughts are racing"))
  ) %>%
  group_by(patient_id, variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            min_value = min(value, na.rm = TRUE),
            max_value = max(value, na.rm = TRUE),
            iqr_value = IQR(value, na.rm = TRUE),
            prop_zero = mean(value == 0, na.rm = TRUE))

saveRDS(patient_statistics,"./example/outputs/patient_statistics.rds")
