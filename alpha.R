#Alpha diversity: 

set.seed(19981103)

#compare 1) raw count data, the 2) rarefied data and the 3) minimally filtered rarefied data: 

# raw_baseline_data #no filtering, not rarefied
# rarefied_data <- rrarefy(select(raw_baseline_data, -(p_ID:BMI)), sample = 2836) #no filtering, rarefied 
min_filt_rarefied_data <- filtered_raw_baseline_data |> 
  column_to_rownames(var = "p_ID") |>
  select(-group, -diagnose, -age, -sex, -BMI) |>
  rrarefy(sample = 2836) |> 
  as.data.frame() #minimal filtering, rarefied data.

#species richness: the number of OTUs in each sample. specnumber finds the number of species. 
# richness <- specnumber(comm) |> as.data.frame()
# richness_2 <- specnumber(rarefied_data)
richness_3 <- specnumber(min_filt_rarefied_data) 

#Shannon's index:
#shannon <- diversity(raw_baseline_data[ ,-1:-7])
#shannon_2 <- diversity(rarefied_data)
shannon_3 <- diversity(min_filt_rarefied_data, index = "shannon")

#Pielou's Evenness index: 
#J = H'/ln(S) where H' is Shannon Weiner diversity and S is the total number of species in a sample, across all samples in dataset.
# evenness <- shannon/log(richness) 
# evenness_2 <- shannon_2/log(richness_2)
evenness_3 <- shannon_3/log(richness_3)

#Simpson's index:
simpson_3 <- diversity(min_filt_rarefied_data, index = "simpson")
#ENS = 1/simpson
ENS_3 <- 1/simpson_3 

# #Hill numbers 
# hill <- exp(shannon)
# hill_2 <- exp(shannon_2)
hill_3 <- exp(shannon_3)

#Create alpha diversity dataframe including environmental data
# alpha <- cbind(shannon = shannon, richness = richness, pielou = evenness, hill = hill, chao = chao, meta_data)
# alpha_2 <- cbind(shannon = shannon_2, richness = richness_2, pielou = evenness_2, hill = hill_2, chao = chao_2, meta_data)
alpha_3 <- cbind(shannon = shannon_3, richness = richness_3, pielou = evenness_3, 
                 simpson = simpson_3, hill = hill_3, ENS = ENS_3, meta_data)

#richness boxplot:
richness_boxplot <- ggplot(alpha_3, aes(x = diagnose, y = richness)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, size = 1.5, show.legend = FALSE) +
  scale_color_manual(name = "Diagnose group:", 
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_fill_manual(name = "Diagnose group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  ylab("Species richness") + 
  xlab("") +
  guides(fill = "none", color = "none") +
  ylim(50, 150)
richness_boxplot

#do the four groups have a significant difference in species richness?
kruskal.test(richness ~ diagnose, data = alpha_3) #p-value = 0.6105

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction
pairwise.wilcox.test(alpha_3$richness, alpha_3$diagnose, p.adjust.method = "BH") #p-values > 0.05

#evenness boxplot:
evenness_boxplot <- ggplot(alpha_3, aes(x = diagnose, y = pielou)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, size = 1.5, show.legend = FALSE) +
  scale_color_manual(name = "Diagnose group:", 
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_fill_manual(name = "Diagnose group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  ylab("Pielou's Evenness index") + 
  xlab("") +
  labs(tag = "B") +
  guides(fill = "none", color = "none") +
  ylim(0, 1)
evenness_boxplot

ggsave("evenness_boxplot.png", plot = evenness_boxplot, width = 18, height = 10, units = "cm", dpi = 300)

#do the four groups have a significant difference in Pielou's evenness index?
kruskal.test(pielou ~ diagnose, data = alpha_3) #p-value = 0.6168

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction
pairwise.wilcox.test(alpha_3$pielou, alpha_3$diagnose, p.adjust.method = "BH") #p-values > 0.05

#shannon boxplot:
shannon_boxplot <- ggplot(alpha_3, aes(x = diagnose, y = shannon)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, size = 1.5, show.legend = FALSE) +
  scale_color_manual(name = "Group:",
                     labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_shape_manual(name = "Group:",
                     labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
                     values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Group:",
                    labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  ylab("Shannon's index") +
  xlab("") +
  #remove the legend
  guides(fill = "none", color = "none")
shannon_boxplot

#do the four groups have a significant difference in Shannon index?
kruskal.test(shannon ~ diagnose, data = alpha_3) #p-value = p-value = 0.6332

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
pairwise.wilcox.test(alpha_3$shannon, alpha_3$diagnose, p.adjust.method = "BH") #p-values > 0.05

#Simpson's boxplot: 
simpson_boxplot <- ggplot(alpha_3, aes(x = diagnose, y = simpson)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, size = 1.5, show.legend = FALSE) +
  scale_color_manual(name = "Group:", 
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_fill_manual(name = "Group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  ylab("Simpson's index") + 
  xlab("") +
  labs(tag = "C") +
  #remove the legend
  guides(fill = "none", color = "none") 
simpson_boxplot

#do the four groups have a significant difference in Simpson index?
kruskal.test(simpson ~ diagnose, data = alpha_3) #p-value = 0.6977

#pairwise comparisons using Wilcoxon rank sum test with continuity correction
pairwise.wilcox.test(alpha_3$simpson, alpha_3$diagnose, p.adjust.method = "BH") #p-values > 0.05

ENS_boxplot <- ggplot(alpha_3, aes(x = diagnose, y = ENS)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, size = 1.5, alpha = 0.7, show.legend = FALSE) +
  scale_color_manual(name = "Group:", 
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_fill_manual(name = "Group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  ylab("Effective number of species (Simpson)") + 
  xlab("") +
  #remove the legend
  guides(fill = "none", color = "none")
ENS_boxplot

#do the four groups have a significant difference in Simpson index?
kruskal.test(simpson ~ diagnose, data = alpha_3) #p-value = 0.6977

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
pairwise.wilcox.test(alpha_3$simpson, alpha_3$diagnose, p.adjust.method = "BH") #p-values > 0.05

ggsave("alpha.png", plot = alpha_plots, 
       width = 18, height = 12, units = "cm", dpi = 300)

#Hill boxplot on raw data:
hill_boxplot <- ggplot(alpha_3, aes(x = diagnose, y = hill)) +
  geom_boxplot(aes(fill = diagnose), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, size = 1.5, alpha = 0.7, show.legend = FALSE) +
  scale_color_manual(name = "Group:",
                     labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_shape_manual(name = "Group:",
                     labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
                     values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Group:",
                    labels = c("HC", "IBS-C", "IBS-M", "IBS-D"),
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  ylab("Hill numbers") +
  xlab("") +
  #remove the legend
  guides(fill = FALSE, color = FALSE) +
  theme_bw()
hill_boxplot

#combine the plots you want: 
alpha_plots <- estimated_species_richness_boxplot + evenness_boxplot + simpson_boxplot + plot_layout(nrow = 1) #ncol = 1
alpha_plots
