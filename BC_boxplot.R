#Boxplot presenting BC distances in the four group: 

#calculate Bray-Curtis dissimilarity: 
dist <- avgdist(filtered_comm, sample = 2836, dmethod = "bray") 

#create sample lookup, with join id: 
sample_lookup <- meta_data |> 
  select(p_ID, diagnose) |> 
  #BG_ID must be character
  mutate(p_ID = as.character(p_ID)) 

#are the spread in the groups significantly different from each other? 

#Multivariate homogeneity of groups dispersions (variances)
mod_median <- betadisper(dist, sample_lookup$diagnose, type = "median", bias.adjust = TRUE) 
anova(mod_median)
permutest(mod_median, pairwise = TRUE)

mod.HSD_median <- TukeyHSD(mod_median)
plot(mod.HSD_median)

plot(mod_median)
boxplot(mod_median)

#extract the distances from the betadisper object to make boxplot with ggplot: 
median_data <- data.frame(mod_median$distances) |> 
  rownames_to_column(var = "p_ID") |>
  #join with sample_lookup to get diagnose information:
  inner_join(sample_lookup, by=c("p_ID"="p_ID")) 

#calculate median distances between the samples and the median of the group:
median_dist <- median_data |> 
  group_by(diagnose) |> 
  summarise(mod_median.distances = median(mod_median.distances))

dist_median_boxplot <- ggplot(median_data, aes(x=diagnose, y=mod_median.distances, fill = diagnose)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, alpha = 0.7) +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x="", y="Distance to median") +
  guides(fill="none", color = "none")
dist_median_boxplot

ggsave("dist_median_boxplot.png", plot = dist_median_boxplot, 
       width = 18, height = 10, units = "cm", dpi = 300) 

#significant difference between the median of the groups? 
kruskal.test(mod_median.distances ~ diagnose, data = median_dist)

#pairwise comparison of the median distances:
pairwise.wilcox.test(median_data$mod_median.distances, median_data$diagnose, p.adjust.method = "BH")

#centroid: 
mod_centroid <- betadisper(dist, sample_lookup$diagnose, type = "centroid", bias.adjust = TRUE)
anova(mod_centroid)
permutest(mod_centroid, pairwise = TRUE)

mod.HSD_centroid <- TukeyHSD(mod_centroid)
plot(mod.HSD_centroid)

plot(mod_centroid)
boxplot(mod_centroid)

#extract the distances from the betadisper object to make boxplot with ggplot: 
centroid_data <- data.frame(mod_centroid$distances) |> 
  rownames_to_column(var = "p_ID") |>
  #join with sample_lookup to get diagnose information:
  inner_join(sample_lookup, by=c("p_ID"="p_ID")) 

dist_centroid_boxplot <- ggplot(centroid_data, aes(x=diagnose, y=mod_centroid.distances, fill = diagnose)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, alpha = 0.7) +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x="", y="Distance to centroid") +
  guides(fill="none", color = "none")
dist_centroid_boxplot

ggsave("dist_centroid_boxplot.png", plot = dist_centroid_boxplot, 
       width = 18, height = 10, units = "cm", dpi = 300)

#inspiration: https://www.youtube.com/watch?v=jLVKJ_n6Qd0&ab_channel=RiffomonasProject

dist_data <- dist|> 
  as.matrix() |> 
  as_tibble(rownames="p_ID") |>
  #pivot to long format with column names as BG_ID_2 and values as BC distance:
  pivot_longer(cols=-p_ID, names_to = "p_ID_2", values_to = "distances") |>
  filter(BG_ID < BG_ID_2) |> #remove self-comparisons
  #join with sample_lookup to get diagnose information: 
  inner_join(sample_lookup, by=c("p_ID"="p_ID")) |> 
  relocate(p_ID, diagnose, p_ID_2, distances) |> 
  #join so that we get the diagnose for BG_ID_2 as well:
  inner_join(sample_lookup, by=c("p_ID_2"="p_ID")) |> 
  mutate(diagnose_2 = diagnose.y) |> 
  mutate(diagnose = diagnose.x) |> 
  select(-c(diagnose.x, diagnose.y)) |> 
  relocate(p_ID, diagnose, p_ID_2, diagnose_2, distances) |> 
  #I only want to do inter diagnose comparisons, so I need to filter out the intra diagnose comparisons:
  mutate(comparison = case_when(
    diagnose == "HC" & diagnose_2 == "HC" ~ "HC",
    diagnose == "IBS-D" & diagnose_2 == "IBS-D" ~ "IBS-D",
    diagnose == "IBS-C" & diagnose_2 == "IBS-C" ~ "IBS-C",
    diagnose == "IBS-M" & diagnose_2 == "IBS-M" ~ "IBS-M")) |> 
  #count(comparison) |> #this is not correct..? 
  drop_na()

#make the boxplot: 
BC_boxplot <- ggplot(dist_data, aes(x=comparison, y=distances, fill = comparison)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, alpha = 0.7) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x="", y="Bray-Curtis dissimilarity") +
  guides(fill="none", color = "none") + 
  scale_y_continuous(limits = c(0,1)) 
BC_boxplot

ggsave("BC_boxplot.png", plot = BC_boxplot, 
       width = 18, height = 10, units = "cm", dpi = 300) 

#are the median of the groups significantly different from each other?
dist_data |> 
  group_by(comparison) |> 
  summarize(median = median(distances))

#significance test:
kruskal.test(median ~ comparison, data = medians)

#pairwise comparison of the median distances:
pairwise.wilcox.test(dist_data$distances, dist_data$comparison, p.adjust.method = "BH")
