#==============================================================================#

library(devtools)

# Preferred parallel library
library(tidyverse)
library(mHMMbayes)
library(viridis)
library(simHMM)


# load data
#insert directory to the data (Note the data is not included as it is only available on inquiry)
data <- haven::read_sav("ESM bipolar cleaned data_UU project.sav")

#------------------------------------------------------------------------------#

# All
train_df <- data%>%
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

#------------------------------------------------------------------------------#

# Load output
out_medhmm <- readRDS(paste0("./example/outputs/out_medHMM_m",4,"_12dv_it2000_c1.rds"))
out_mhmm <- readRDS(paste0("./example/outputs/out_mHMM_m",4,"_12dv_it2000_c1.rds"))


#------------------------------------------------------------------------------#
# Figure R.4: validated scales (decoding)

# Load weekly data
#(Note the data is not included as it is only available on inquiry)
bs_asrm <- haven::read_sav("bs_asrm_cleaned_UU.sav") %>%
  arrange(patient_id)
#(Note the data is not included as it is only available on inquiry)
bs_qids <- haven::read_sav("bs_qids_cleaned_UU.sav") %>%
  arrange(patient_id)

# Put together
bs_scales <- full_join(bs_asrm, bs_qids) %>%
  dplyr::select(patient_id, bs_asrm_open_from, bs_asrm_tot, bs_qids_tot) %>%
  rename("open_from" = "bs_asrm_open_from") %>%
  group_by(patient_id) %>%
  mutate(from_time_day = format(as.POSIXct(lag(open_from, 1, default = min(open_from)-(3600*24*7))),format = "%Y-%m-%d"),
         to_time_day = format(as.POSIXct(open_from-1),format = "%Y-%m-%d"),
         patient_id = factor(patient_id)) %>%
  ungroup()



# Get decoding

states_mhmm <- do.call(rbind, lapply(1:length(out_mhmm$sample_path), function(s){
  apply(out_mhmm$sample_path[[s]], 1, function(r) {
    which.max(table(factor(r, levels = 1:4))/length(r))
    }) %>%
    as.data.frame() %>%
    rename("state" = ".") %>%
    mutate(patient_id = s,
           occasion = row_number(),
           model = "mHMM")
}))

states_medhmm <- do.call(rbind, lapply(1:length(out_medhmm$sample_path), function(s){
  apply(out_medhmm$sample_path[[s]], 1, function(r) {
    which.max(table(factor(r, levels = 1:4))/length(r))
  }) %>%
    as.data.frame() %>%
    rename("state" = ".") %>%
    mutate(patient_id = s,
           occasion = row_number(),
           model = "medHMM")
}))


states <- rbind(states_mhmm, states_medhmm)

data_labelled <- left_join(as.data.frame(train_df) %>% group_by(patient_id) %>%
                             mutate(occasion = row_number()), states,by=c("patient_id","occasion")) %>%

  as.data.frame()

# Add week info
data_labelled <- left_join(data_labelled, data %>%
                             drop_na() %>%
                             group_by(patient_id) %>%
                             mutate(patient_id = cur_group_id(),
                                    occasion = row_number()) %>%
                             ungroup() %>%
                             arrange(patient_id) %>%
                             dplyr::select(patient_id, occasion, open_from, bs_diary_open_from, bs_diary_date))


# Validated scales
state_data <- data_labelled %>%
  select(patient_id, state, bs_diary_open_from, occasion, model) %>%
  rename("open_from" = "bs_diary_open_from") %>%
  group_by(patient_id) %>%
  mutate(state = factor(state, levels = 1:4, labels = c("euthymic","manic","mixed","depressive")),
         from_time = lag(open_from, 1, default = min(open_from)-(3600*3)),
         to_time = open_from-1,
         patient_id = factor(patient_id),
         from_occ = occasion-1,
         to_occ = occasion) %>%
  mutate(day = format(as.POSIXct(to_time), format = "%Y-%m-%d")) %>%
  ungroup()

scale_data <- bs_scales %>%
  group_by(patient_id) %>%
  mutate(patient_id = cur_group_id(),
         patient_id = as.factor(patient_id),
         occasion = row_number(),
         day = from_time_day,
         bs_asrm_tot = case_when(is.na(bs_asrm_tot) ~ 99,
                                 !is.na(bs_asrm_tot) ~ bs_asrm_tot),
         bs_qids_tot = case_when(is.na(bs_qids_tot) ~ 99,
                                 !is.na(bs_qids_tot) ~ bs_qids_tot)
  ) %>%
  ungroup() %>%
  select(-occasion, -open_from)

join_data <- full_join(state_data, scale_data) %>%
  tidyr::fill(bs_asrm_tot:bs_qids_tot, .direction = "down") %>%
  gather(variable, value, -patient_id,-state,-open_from,-occasion,-from_time,-to_time,-from_occ,-to_occ,-from_time_day,-to_time_day,-day,-model) %>%
  mutate(value = case_when(value == 99 ~ NA_real_,
                           value != 99 ~ value))

p <- ggplot() +
  geom_line(data = join_data %>%
              group_by(patient_id) %>%
              mutate(variable = factor(variable,
                                       levels = c("bs_asrm_tot","bs_qids_tot"),
                                       labels = c("ASRM","QIDS")),
                     patient_id = factor(patient_id, levels = 1:20, labels = paste0("patient ",1:20))),
            aes(x = from_occ, y = value, colour = variable)) +
  geom_rect(data = state_data %>%
              mutate(patient_id = factor(patient_id, levels = 1:20, labels = paste0("patient ",1:20)),
                     ymin = rep(c(25,32), nrow(state_data)/2), ymax = rep(c(30, 37), nrow(state_data)/2)), aes(xmin = from_occ, xmax = to_occ,
                                                                                                            ymin = ymin, ymax = ymax, fill = state)) +

  scale_fill_manual(values = c("#e49444","#d1615d","#85b6b2","#5778a4")) +
  geom_hline(yintercept = 6, linetype = "dashed") +
  facet_wrap(patient_id~., scales = "free", ncol = 4) +
  theme_minimal() +
  labs(colour = "Questionnaire", fill = "Mood state") +
  theme(legend.position = "bottom") +
  ggtitle(label = "Temporal alignment between mood states and weekly symptom scores") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(label = "Weekly symptom score") +
  xlab(label = "Measurement occasion")

p

ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "./example/plots/decoding_14.pdf", width = 14, height = 11, units = "in", dpi = 300)



# Selecting patients

p <- ggplot() +
  geom_line(data = join_data %>%
              group_by(patient_id) %>%
              mutate(variable = factor(variable,
                                       levels = c("bs_asrm_tot","bs_qids_tot"),
                                       labels = c("ASRM","QIDS")),
                     patient_id = factor(patient_id, levels = 1:20, labels = paste0("patient ",1:20))) %>%
              filter(patient_id %in% paste0("patient ",c(3,4,14))),
            aes(x = from_occ, y = value, colour = variable)) +
  geom_rect(data = state_data %>%
              mutate(patient_id = factor(patient_id, levels = 1:20, labels = paste0("patient ",1:20)),
                     ymin = rep(c(25,32), nrow(state_data)/2), ymax = rep(c(30, 37), nrow(state_data)/2)) %>%
              filter(patient_id %in% paste0("patient ",c(3,4,14))), aes(xmin = from_occ, xmax = to_occ,
                                                                                                               ymin = ymin, ymax = ymax, fill = state)) +

  scale_fill_manual(values = c("#e49444","#d1615d","#85b6b2","#5778a4")) +
  geom_hline(yintercept = 6, linetype = "dashed") +
  facet_wrap(patient_id~., scales = "free", ncol = 4) +
  theme_minimal() +
  labs(colour = "Questionnaire", fill = "Mood state") +
  theme(legend.position = "bottom") +
  ggtitle(label = "Temporal alignment between mood states and weekly symptom scores") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(label = "Weekly symptom score") +
  xlab(label = "Measurement occasion")

p


ggsave(plot = p + theme(text = element_text(size = 14)),
       filename = "./example/plots/decoding_3patients_14.pdf", width = 14, height = 4, units = "in", dpi = 300)
