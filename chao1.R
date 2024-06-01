#I want to calculate the expected number of species in each individual 
#and group it by diagnose to get the "true" number of species in each group.

#Chao1 from vegan: estimateR or specpool
est_species <- specpool(filtered_comm) #for the whole dataset. 
#est_species_rare <- specpool(min_filt_rarefied_data) 

est_species <- estimateR(filtered_comm) |> 
  t() |> 
  as.data.frame() |>
  rownames_to_column("p_ID") |> 
  #BG_ID needs to be numeric: 
  mutate(BG_ID = as.numeric(p_ID)) |>
  #join with meta data to get diagnose information: 
  left_join(meta_data, by = c("p_ID" = "p_ID")) 

# est_species_2_rare <- estimateR(min_filt_rarefied_data) |> #for individual samples.
#   t() |> 
#   as.data.frame() |>
#   rownames_to_column("p_ID") |> 
#   #BG_ID needs to be numeric: 
#   mutate(p_ID = as.numeric(p_ID)) |>
#   #join with meta data to get diagnose information: 
#   left_join(meta_data, by = c("p_ID" = "p_ID"))

#box plot of observed species richness:
observed_speces_richness <- est_species |> 
  ggplot(aes(x = diagnose, y = S.obs, fill = diagnose)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3) +
  labs(x = NULL, y = "Observed species richness") + 
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  theme(legend.position = "none") 
observed_speces_richness

ggsave("observed_species_richness.png", plot = observed_speces_richness, 
       width = 18, height = 10, units = "cm", dpi = 300)

#box plot of estimated species richness:
estimated_species_richness_boxplot <- est_species |>
  ggplot(aes(x = diagnose, y = S.chao1, fill = diagnose)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3) +
  labs(x = NULL, y = "Chao1 (expected species richness)") +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  theme(legend.position = "none") +
  labs(tag = "A")
estimated_species_richness_boxplot

ggsave("expected_species_richness.png",plot = estimated_species_richness_boxplot,
       width = 18, height = 10, units = "cm", dpi = 300)

#is it a significant difference between the groups?
#Kruskal-Wallis rank sum test:
kruskal.test(S.chao1 ~ diagnose, data = est_species) #p-value = 0.1232

#pairwise Wilcoxon rank sum test:
pairwise.wilcox.test(est_species$S.chao1, est_species$diagnose, p.adjust.method = "BH") #no significant differences between the groups.

# #both observed and expected species richness in one plot: 
# ggplot(est_species_2) +
#   geom_jitter(aes(x = diagnose, y = S.obs, color = diagnose), width = 0.3, alpha = 0.2) +
#   geom_jitter(aes(x = diagnose, y = S.chao1, color = diagnose), width = 0.3) +
#   geom_boxplot(aes(x = diagnose, y = S.chao1, fill = diagnose), alpha = 0.5, outlier.shape = NA) + 
#   labs(x = NULL, y = "Expected and observed number of species") +
#   scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   theme(legend.position = "none")