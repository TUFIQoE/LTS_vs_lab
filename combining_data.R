library(tidyverse)
library(magrittr)

rm(list = ls())

lts_data <- read_csv("phase4_current_results_p1204.csv", col_types = "cnnnnnnnn-n----") 
lts_data %<>%
  filter(n == 7, week < 14) %>%
  select(-n)
lts_data$group = "lts"

lab_data <- read_csv("phase4_lab.csv", col_types = "cnnc")
lab_data %<>%
  left_join(lts_data %>% 
              group_by(week) %>% 
              summarize(mo = first(mo), tu = first(tu), we = first(we), 
                        th = first(th), fr = first(fr), sa = first(sa)), 
            by = "week")
all_data <- rbind(lts_data, lab_data)
all_data <- all_data[!(all_data$external_id %in% c("61", "65", "42", "32", "41", "25", "20")), ]
s1 <- all_data %>% filter(group == "lts") %>%
  group_by(week) %>%
  dplyr::summarise(mean_q2_lts = mean(q2))
s2 <- all_data %>% filter(group != "lts") %>%
  group_by(week) %>%
  dplyr::summarise(mean_q2_lab = mean(q2))
scatter <- left_join(s1, s2)

ggplot(scatter, aes(x = mean_q2_lts, y = mean_q2_lab)) + geom_point()  + ylim(2,5)+ xlim(2,5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")

s3 <- all_data %>% filter(group == "labPl") %>%
  group_by(week) %>%
  dplyr::summarise(mean_q2_labPl = mean(q2))
s4 <- all_data %>% filter(group == "labNor") %>%
  group_by(week) %>%
  dplyr::summarise(mean_q2_labNor = mean(q2))
scatter_labs <- left_join(s3, s4)

ggplot(scatter_labs, aes(x = mean_q2_labPl, y = mean_q2_labNor)) + geom_point()  + ylim(2,5)+ xlim(2,5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")

scatter_Pl <- left_join(s1, s3)
scatter_Nor <- left_join(s1, s4)
scatter_all <- left_join(scatter_Pl, scatter_Nor)

ggplot(scatter_Pl, aes(x = mean_q2_lts, y = mean_q2_labPl)) + geom_point()  + ylim(2,5)+ xlim(2,5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")

ggplot(scatter_Nor, aes(x = mean_q2_lts, y = mean_q2_labNor)) + geom_point()  + ylim(2,5)+ xlim(2,5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")

colors <- c("Polish" = "blue", "Norwegian" = "red")
ggplot(scatter_all) + geom_point(aes(x = mean_q2_lts, y = mean_q2_labPl, color = "Polish")) + 
  geom_point(aes(x = mean_q2_lts, y = mean_q2_labNor, color = "Norwegian"), alpha = 0.5) + ylim(2,5)+ xlim(2,5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + labs(
    x = "MOS from LTS",
    y = "MOS from labs",
    color = "Lab") +
  scale_color_manual(values = colors) + theme_bw()

saveRDS(all_data, file = "all_data_1204.RDS")
write.csv(all_data,"lab_and_lts.csv", row.names = FALSE)
