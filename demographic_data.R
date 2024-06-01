#meta data: demographics: number, sex, age, BMI, IBS-SSS. 

# #HC, IBS-C, IBS-M, IBS-D
# diagnose_group_plot <- ggplot(meta_data, aes(x = group, fill = diagnose)) +
#   geom_bar(position = "stack", alpha = 0.5) +
#   labs(title = "",
#        x = "",
#        y = "Number",
#        fill = "") +
#   scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   ylim(0, 60)
# diagnose_group_plot

#count the number in each diagnose: 
meta_data |> 
  group_by(group) |> 
  count() 

#are the differences in number of each group significantly different?
#chi-square test:
chisq.test(number$n) #p-value = 0.03706

#sex:

#Count the number of females in each of the four diagnosis groups: 
meta_data |> 
  group_by(diagnose) |> 
  count(sex)

#are the differences in number of each sex between the groups significant? 
#chi-square test:
chisq.test(table(meta_data$diagnose, meta_data$sex)) #p-value = 0.2042

#I want to make a stacked bar plot that include both the number of individuals in each of the four groups, and their sex: 
meta_data$diagnose_sex <- interaction(meta_data$diagnose, meta_data$sex, sep = " ")

#I want to make a barplot of HC, with one column for female and one for male: 
n_sex_plot <- ggplot(meta_data, aes(fill=diagnose,  x=sex)) + 
  geom_bar(alpha = 0.5) + 
  theme_classic(base_size = 18) +
  facet_wrap(~group) +
  labs(title = NULL, x = NULL, y = "Number", fill = NULL) +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863"))+
  #remove the top from the facet_wrap:
  theme(strip.text = element_blank())
n_sex_plot

ggsave("n_sex_plot.png", plot = n_sex_plot, 
       width = 13, height = 10, units = "cm", dpi = 300)

#age: 
#calculate the median ages for each group
medians_age <- meta_data |>
  group_by(group) |>
  summarise(median_age = median(age, na.rm = TRUE))

#visualize: density plot or boxplot? 

#density plot: by group 
age_density_plot <- ggplot(meta_data, aes(x = age, fill = group, color = group)) +
  geom_density(alpha = 0.5) +
  labs(x = "Age", y = NULL) + 
  scale_fill_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
  #add median lines for each group: 
  geom_vline(data = medians_age, aes(xintercept = median_age, color = group), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) + 
  guides(fill = "none", color = "none") + 
  xlim(0,80)
age_density_plot

# #Box plot: by diagnose 
# age_box_plot <- ggplot(meta_data, aes(x = diagnose, y = age, fill = diagnose)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   #geom_violin(alpha = 0.5) +
#   geom_jitter(width = 0.3, aes(color = diagnose)) +
#   labs(title = "", x = "", y = "Age (years)", fill = "Group:") +
#   scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   guides(fill = "none", color = "none") 
# age_box_plot

#calculate the range and median age for each of the four diagnosis groups: 
meta_data |> 
  group_by(diagnose) |> 
  summarise(min_age = min(age, na.rm = TRUE),
            max_age = max(age, na.rm = TRUE),
            median_age = median(age, na.rm = TRUE))

#are the median age between the IBS group and the HC significantly different?
#t-test:
t.test(meta_data$age ~ meta_data$group) #p-value = 0.606

#BMI: 
#calculate the median ages for each group
medians_BMI <- meta_data |>
  group_by(group) |>
  summarise(median_BMI = median(BMI, na.rm = TRUE))

#Density plot: 
BMI_density_plot <- ggplot(meta_data, aes(x = BMI, fill = group, color = group)) +
  geom_density(alpha = 0.5) +
  labs(x = "BMI", y = NULL, fill = NULL) +
  scale_fill_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
  #add median lines for each group
  geom_vline(data = medians_BMI, aes(xintercept = median_BMI, color = group), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
  guides(color = "none") +
  xlim(10,40)
BMI_density_plot

# #Box plot: by diagnose
# BMI_box_plot <- ggplot(meta_data, aes(x = diagnose, y = BMI, fill = diagnose)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   #geom_violin(alpha = 0.5) +
#   geom_jitter(width = 0.3, aes(color = diagnose)) +  
#   labs(title = "", x = "", y = "BMI", fill = "Group:") +
#   scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +  # Match color scale to fill
#   guides(fill = "none", color = "none")  #hide legend
# BMI_box_plot

#calculate the range and median BMI for each of the four diagnosis groups: 
meta_data |> 
  group_by(diagnose) |> 
  summarise(min_BMI = min(BMI, na.rm = TRUE),
            max_BMI = max(BMI, na.rm = TRUE),
            median_BMI = median(BMI, na.rm = TRUE))

#are the median BMI between the IBS group and the HC significantly different?
#t-test:
t.test(meta_data$BMI ~ meta_data$group) #p-value = 0.606

#combined plots: sex, age and BMI 
age_BMI_plot <-  age_density_plot | BMI_density_plot
age_BMI_plot

ggsave("age_BMI_plot.png", plot = age_BMI_plot, width = 15, height = 10, units = "cm", dpi = 300)

# Calculate the median ages for each group
medians <- IBS_SSS_data |> 
  group_by(group) |> 
  summarise(median_SSS = median(IBS_SSS, na.rm = TRUE))

# # Density plot for IBS-SSS: 
# IBS_SSS_p  <- ggplot(IBS_SSS_data, aes(x = IBS_SSS, fill = group)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "", x = "IBS-SSS score", y = "Density", fill = "Group:") +
#   theme_bw() +
#   scale_fill_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   # Add median lines for each group
#   geom_vline(data = medians, aes(xintercept = median_SSS, color = group), linetype="dashed", size = 0.7) +
#   scale_color_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) + # To ensure that the color of the lines matches the color of the fills.
#   guides(color = FALSE) + # This line removes the color legend
#   xlim(0, 500) +
#   # Add threshold line at x value 175
#   geom_vline(xintercept = treshold, linetype = "solid", color = "black", size = 0.6) 
# IBS_SSS_p

#box plot for IBS-SSS by diagnose: 
IBS_SSS_boxplot <- ggplot(IBS_SSS_data, aes(x = diagnose, y = IBS_SSS, fill = diagnose)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.3, aes(color = diagnose)) +
  labs(x = NULL, y = "IBS-SSS score", fill = "Diagnose:") +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  geom_hline(yintercept = 175, linetype = "dashed", color = "black", linewidth = 0.7) + #threshold. 
  geom_hline(yintercept = 75, linetype = "dashed", color = "black", linewidth = 0.7) + #threshold.
  geom_hline(yintercept = 300, linetype = "dashed", color = "black", linewidth = 0.7) + #threshold.
  guides(fill = "none", color = "none") +
  #I want to add starts to the boxplot to indicate the significance of the pairwise comparisons: 
  geom_signif(comparisons = list(c("HC", "IBS-C")), map_signif_level = TRUE, textsize = 6, vjust = 0.5, step_increase = 0.05) +
  geom_signif(comparisons = list(c("HC", "IBS-M")), map_signif_level = TRUE, textsize = 6, vjust = 0.5, step_increase = 0.05) +
  geom_signif(comparisons = list(c("HC", "IBS-D")), map_signif_level = TRUE, textsize = 6, vjust = 0.5, step_increase = 0.05) +
  geom_signif(comparisons = list(c("IBS-C", "IBS-M")), map_signif_level = TRUE, textsize = 6, vjust = 0.5, step_increase = 0.05) +
  geom_signif(comparisons = list(c("IBS-C", "IBS-D")), map_signif_level = TRUE, textsize = 6, vjust = 0.5, step_increase = 0.05) +
  geom_signif(comparisons = list(c("IBS-M", "IBS-D")), map_signif_level = TRUE, textsize = 6, vjust = 0.5, step_increase = 0.05) 
IBS_SSS_boxplot

ggsave("IBS_SSS_boxplot.png", plot = IBS_SSS_boxplot, width = 18, height = 10, units = "cm", dpi = 300)

#are the IBS-SSS scores between the four groups significantly different?
kruskal.test(IBS_SSS ~ diagnose, data = IBS_SSS_data) 

#pairwise comparison: 
pairwise.wilcox.test(IBS_SSS_data$IBS_SSS, IBS_SSS_data$diagnose, p.adjust.method = "bonferroni")

