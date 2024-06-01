#Principle coordinate analysis (POcA):

#from: https://www.rpubs.com/roalle/mres_2019

set.seed(19981103)

#Bray-Curtis Dissimilarity: abundance-based dissimilarity metric
bray <- vegdist(filtered_raw_baseline_data[,-1:-6], method = "bray")

#SIMPER: similarity percentage. Identifying which species contribute most to beta diversity (Bray-Curtis dissimilarity)
simper <- simper(filtered_raw_baseline_data[,-1:-6], filtered_raw_baseline_data$diagnose, permutations = 99) 
summary(simper, ordered = TRUE, digits = max(3,getOption("digits") - 3)) 

#Ordinations: 
#PCoA: unconstrained ordination. 

#1: 
?cmdscale #Classical (Metric) Multidimensional Scaling (MDS) 
pcoa_bray <- cmdscale(bray, k = 2, eig = T) 

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa_bray_plotting <- as.data.frame(pcoa_bray$points)
colnames(pcoa_bray_plotting) <- c("axis_1", "axis_2")
pcoa_bray_plotting$diagnose <- raw_baseline_data$diagnose 

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa_bray$eig[1]/(sum(pcoa_bray$eig))
pcoa_bray$eig[2]/(sum(pcoa_bray$eig))

# create a PCoA plot
pcoa_bray_plot <- ggplot(pcoa_bray_plotting, aes(x = axis_1, y = axis_2, color = diagnose, shape = diagnose, fill = diagnose)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.75, geom = "polygon", alpha = 0.1, show.legend = FALSE) + #not exact fitting of the data. ggforce. 
  scale_color_manual(name = "Group:", 
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_shape_manual(name = "Group:",
                     values = c(21, 22, 23, 24)) +
  scale_fill_manual(name = "Group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  theme(legend.key.height = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black"), 
        legend.position = c(0.89, 0.87)) + 
  coord_fixed(xlim = c(-0.50, 0.50), ylim = c(-0.50, 0.50)) + 
  xlab("PCoA 1 (20.4%)") +
  ylab("PCoA 2 (8.1%)")
pcoa_bray_plot 

ggsave("pcoa_bray.png", plot = pcoa_bray_plot, 
        width = 18, height = 12, units = "cm", dpi = 300) 

#test significance off the PCoA plot: 
#PERMAVOVA: 
adonis2(filtered_raw_baseline_data[ ,-c(1:6)] ~ filtered_raw_baseline_data$diagnose, 
        permutations = 999, method = "bray")
#filtered_raw_baseline_data$diagnose: Pr(>F): 0.487. 

