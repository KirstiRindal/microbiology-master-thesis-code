#Significant OTUs 
#from: https://www.youtube.com/watch?v=IoeSrLXeSos&ab_channel=RiffomonasProject

library(broom)

#test this with the raw count data, with 10% prevalence filter: 
prev_filtered_long_raw_baseline_data <- pivot_longer(raw_baseline_data, 
                                                     cols = -c(p_ID, group, diagnose, sex, age, BMI), 
                                                     names_to = "bacteria_OTU",
                                                     values_to = "abundance") |> 
  #remember to filter subgroup before prevalence filtering! 
  mutate(n_patients = n_distinct(p_ID)) |>
  group_by(bacteria_OTU) |> 
  filter(sum(abundance > 0) > 0.1 * n_patients) |> #236 
  select(-n_patients) |> 
  select(-c(group, sex, age, BMI))

# #test for all four groups at once, with Kruskal-Wallis test:
# sig_OTUs <- prev_filtered_long_raw_baseline_data |> 
#   group_by(bacteria_OTU) |>
#   nest(data = -bacteria_OTU) |> 
#   mutate(test = map(.x = data, ~  kruskal.test(abundance ~ diagnose, data = .x) |> 
#                       tidy())) |>
#   unnest(test) |> 
#   mutate(p.adjust = p.adjust(p.value, method = "BH")) |> 
#   filter(p.adjust < 0.05) |> 
#   select(bacteria_OTU, p.adjust) 
# 
# prev_filtered_long_raw_baseline_data |> 
#   inner_join(sig_OTUs, by = "bacteria_OTU") |> 
#   ggplot(aes(x = abundance+1, y = bacteria_OTU, color = diagnose, fill = diagnose)) + 
#   geom_jitter(aes(shape = diagnose), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3), 
#               size = 1) +
#   #add median and 50% CI for HC and IBS: 
#   stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
#                geom = "pointrange", position = position_dodge(width = 0.8),
#                color = "black", size = 0.2, show.legend = FALSE) +
#   #fix the appearance: 
#   scale_color_manual(name = "Group:",
#                      labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
#                      values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   scale_shape_manual(name = "Group:",
#                      labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
#                      values = c(21, 22, 23, 24)) +
#   scale_fill_manual(name = "Group:",
#                     labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
#                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   labs(x = "Relative abundance", y = "") + 
#   theme(legend.key.height = unit(0.4, "cm"), 
#         legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) + 
#   #log transform the x axis 
#   scale_x_log10()

# #IBS-D vs HC in "absolute abundance"
# long_IBS_D <- prev_filtered_long_raw_baseline_data |> 
#   filter(diagnose %in% c("HC", "IBS-D"))
# 
# sig_OTUs <- long_IBS_D |> 
#   group_by(bacteria_OTU) |>
#   nest(data = -bacteria_OTU) |> 
#   mutate(test = map(.x = data, ~  wilcox.test(abundance ~ diagnose, data = .x) |> 
#                       tidy())) |>
#   unnest(test) |> 
#   mutate(p.adjust = p.adjust(p.value, method = "BH")) |> 
#   filter(p.adjust < 0.05) |> 
#   select(bacteria_OTU, p.adjust) 
# 
# long_IBS_D |> 
#   inner_join(sig_OTUs, by = "bacteria_OTU") |> 
#   ggplot(aes(x = abundance+1, y = bacteria_OTU, color = diagnose, fill = diagnose)) + 
#   geom_jitter(aes(shape = diagnose), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3), 
#               size = 1) +
#   #add median and 50% CI for HC and IBS: 
#   stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
#                geom = "pointrange", position = position_dodge(width = 0.8),
#                color = "black", size = 0.2, show.legend = FALSE) +
#   #fix the appearance: 
#   scale_color_manual(name = "Group:",
#                      labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
#                      values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   scale_shape_manual(name = "Group:",
#                      labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
#                      values = c(21, 22, 23, 24)) +
#   scale_fill_manual(name = "Group:",
#                     labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
#                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   labs(x = "Relative abundance", y = "") + 
#   theme(legend.key.height = unit(0.4, "cm"), 
#         legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) + 
# #log transform the x axis 
# scale_x_log10()
# 
# #more results when the data is not relative abundance. 

#transform the data to relative abundance 
rel_abund <- decostand(raw_baseline_data[ ,-c(1:6)], method = "clr", pseudocount=1)

#change BG_ID in meta_data to character: 
meta_data$p_ID <- as.character(meta_data$p_ID)
  
rel_abund <- rel_abund |> 
  as.data.frame() |>
  #the first column needs to be BG_ID
  rownames_to_column(var = "p_ID") |> 
  #join with meta data to get diagnose info 
  inner_join(meta_data, by = "p_ID") |>  
  select(-c(group, sex, age, BMI, ONT_code))


rel_abund_long <- pivot_longer(rel_abund, 
                          cols = -c(p_ID, diagnose), 
                          names_to = "bacteria_OTU",
                          values_to = "abundance") |> 
  #10% prevalence filter
  mutate(n_patients = n_distinct(p_ID)) |>
  group_by(bacteria_OTU) |> 
  filter(sum(abundance > 0) > 0.1 * n_patients) |> #236 OTUs keept. 
  select(-n_patients)

#test one of the sub groups of IBS vs HC: 

#1: IBS-C vs HC
#filter out HC and IBS-C: 
rel_abund_long_IBS_C <- rel_abund_long |> 
  filter(diagnose %in% c("HC", "IBS-C"))

sig_OTUs <- rel_abund_long_IBS_C |> 
  group_by(bacteria_OTU) |>
  nest(data = -bacteria_OTU) |> 
  mutate(test = map(.x = data, ~  wilcox.test(abundance ~ diagnose, data = .x) |> 
                      tidy())) |>
  unnest(test) |> 
  mutate(p.adjust = p.adjust(p.value, method = "BH")) |> 
  filter(p.adjust < 0.05) |> 
  select(bacteria_OTU, p.adjust) 

IBSC_plot <- rel_abund_long_IBS_C |> 
  inner_join(sig_OTUs, by = "bacteria_OTU") |> 
  ggplot(aes(x = abundance+1, y = bacteria_OTU, color = diagnose, fill = diagnose)) + 
  geom_jitter(aes(shape = diagnose), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3), 
              size = 1, alpha = 0.5) +
  #add median and 50% CI for HC and IBS: 
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", position = position_dodge(width = 0.8),
               color = "black", size = 0.2, show.legend = FALSE) +
  #fix the appearance: 
  scale_color_manual(name = "Group:",
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_shape_manual(name = "Group:",
                     values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x = "Abundance (centered log ratio)", y = "") + 
  theme(legend.key.height = unit(0.4, "cm"), 
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) 
IBSC_plot

#save the plot:
ggsave("IBSC_plot.png", plot = IBSC_plot, width = 18, height = 15, units = "cm")

#2: IBS-D vs HC:
#filter out HC and IBS-D: 
rel_abund_long_IBS_D <- rel_abund_long |> 
  filter(diagnose %in% c("HC", "IBS-D"))

sig_OTUs <- rel_abund_long_IBS_D |> 
  group_by(bacteria_OTU) |>
  nest(data = -bacteria_OTU) |> 
  mutate(test = map(.x = data, ~  wilcox.test(abundance ~ diagnose, data = .x) |> 
                      tidy())) |>
  unnest(test) |> 
  mutate(p.adjust = p.adjust(p.value, method = "BH")) |> 
  filter(p.adjust < 0.05) |> 
  select(bacteria_OTU, p.adjust) 

IBSD_plot <- rel_abund_long_IBS_D |> 
  inner_join(sig_OTUs, by = "bacteria_OTU") |> 
  ggplot(aes(x = abundance+1, y = bacteria_OTU, color = diagnose, fill = diagnose)) + 
  geom_jitter(aes(shape = diagnose), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3), 
              size = 1, alpha = 0.5) +
  #add median and 50% CI for HC and IBS: 
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", position = position_dodge(width = 0.8),
               color = "black", size = 0.2, show.legend = FALSE) +
  #fix the appearance: 
  scale_color_manual(name = "Group:",
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_shape_manual(name = "Group:",
                     values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x = "Abundance (centered log ratio)", y = "") + 
  theme(legend.key.height = unit(0.4, "cm"), 
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) 
IBSD_plot

#save the plot:
ggsave("IBSD_plot.png", plot = IBSD_plot, width = 18, height = 15, units = "cm")

#3: IBS-M vs HC:
#filter out HC and IBS-M: 
rel_abund_long_IBS_M <- rel_abund_long |> 
  filter(diagnose %in% c("HC", "IBS-M"))
  
sig_OTUs <- rel_abund_long_IBS_M |> 
  group_by(bacteria_OTU) |>
  nest(data = -bacteria_OTU) |> 
  mutate(test = map(.x = data, ~  wilcox.test(abundance ~ diagnose, data = .x) |> 
                      tidy())) |>
  unnest(test) |> 
  mutate(p.adjust = p.adjust(p.value, method = "BH")) |> 
  filter(p.adjust < 0.05) |> 
  select(bacteria_OTU, p.adjust) 

IBSM_plot <-rel_abund_long_IBS_M |> 
  inner_join(sig_OTUs, by = "bacteria_OTU") |> 
  ggplot(aes(x = abundance+1, y = bacteria_OTU, color = diagnose, fill = diagnose)) + 
  geom_jitter(aes(shape = diagnose), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3), 
              size = 1, alpha = 0.5) +
  #add median and 50% CI for HC and IBS: 
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", position = position_dodge(width = 0.8),
               color = "black", size = 0.2, show.legend = FALSE) +
  #fix the appearance: 
  scale_color_manual(name = "Group:",
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_shape_manual(name = "Group:",
                     values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x = "Abundance (centered log ratio)", y = NULL) + 
  theme(legend.key.height = unit(0.4, "cm"), 
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) 
IBSM_plot

#save the plot:
ggsave("IBSM_plot.png", plot = IBSM_plot, width = 18, height = 15, units = "cm")

#CHANGE: 
#italicize the OTU names
#remove the underscores from the OTU names

#no clear pattern. 
#very patchy: many of the participants do not have the OTUs. 
#high level of XXX -> IBS symptoms. 
#not one OTU by itself, but in combination. 
#not causation, but association.
#loss of OTU or gain of OTU causing the symptoms? 

#limitaions to this method: 
#does not take into acount that the data is compositional.
#are subject to inflated false discovery rates
