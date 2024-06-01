#Species Accumulation Curves 

set.seed(19981103)

#The function "specaccum" finds species accumulation curves or the number of species 
#for a certain number of sampled sites or individuals.

#make one plot for each diagnose, and then combine the four plots: 

#the curves needs to be constructed on minimally filtered rarefied data: 
min_filt_rarefied_data <- rrarefy(filtered_comm, sample = 2836) |>
  as.data.frame() #minimal filtering, rarefied data.

min_filt_rarefied_data <- min_filt_rarefied_data |> 
  as_tibble(rownames = "p_ID") 
#p_ID in meta_data need to be character:
meta_data <- meta_data |> 
  mutate(p_ID = as.character(p_ID))

#join the minimally filtered raraefied data with the meta data: 
data <- min_filt_rarefied_data |> 
  left_join(meta_data, by = c("p_ID" = "p_ID")) |> 
  relocate(p_ID, diagnose) |> 
  select(-ONT_code, -group, -age, -sex, -BMI)

#HC:  
HC <- data |> 
  filter(diagnose == "HC") |> 
  select(-c(1:2)) |> 
  as.data.frame() 

HC_exa <- specaccum(HC, "exact", permutations = 999) #finds the expected (mean) species richness
HC_data <- data.frame(sites=HC_exa$sites, richness=HC_exa$richness, SD=HC_exa$sd)

#IBS-D: 
IBS_D <- data |> 
  filter(diagnose == "IBS-D") |> 
  select(-c(1:2)) |> 
  as.data.frame()

IBS_D_exact <- specaccum(IBS_D, "exact", permutations = 999) #finds the expected (mean) species richness
IBS_D_data <- data.frame(sites=IBS_D_exact$sites, richness=IBS_D_exact$richness, SD=IBS_D_exact$sd)

#IBS-C:
IBS_C <- data |> 
  filter(diagnose == "IBS-C") |> 
  select(-c(1:2)) |> 
  as.data.frame() 

IBS_C_exact <- specaccum(IBS_C, "exact", permutations = 999) #finds the expected (mean) species richness
IBS_C_data <- data.frame(sites=IBS_C_exact$sites, richness=IBS_C_exact$richness, SD=IBS_C_exact$sd)

#IBS-M:
IBS_M <- data |> 
  filter(diagnose == "IBS-M") |> 
  select(-c(1:2)) |> 
  as.data.frame()

IBS_M_exact <- specaccum(IBS_M, "exact") #finds the expected (mean) species richness
IBS_M_data <- data.frame(sites=IBS_M_exact$sites, richness=IBS_M_exact$richness, SD=IBS_M_exact$sd)

#combine the data:
all_data <- bind_rows(
  HC = HC_data,
  `IBS-C` = IBS_C_data,
  `IBS-D` = IBS_D_data,
  `IBS-M` = IBS_M_data,
  .id = "diagnose")

#visualize the data:
sac_all_plot <- ggplot(all_data, aes(x=sites, y=richness, fill = diagnose) ) +
  geom_line(aes(color = diagnose)) +
  geom_ribbon(aes(ymin=(richness-2*SD),ymax=(richness+2*SD)), alpha=0.2) +
  scale_color_manual(name = "Group:", 
                     values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_fill_manual(name = "Group:",
                    values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x ="Number of participants", y= "Number of OTUs") +
  xlim(0, 40) +
  theme(legend.position = c(0.90, 0.20)) +
  theme(legend.key.height = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black"))
sac_all_plot

#save the plot:
ggsave("rare_sac_plot.png", plot = sac_all_plot, 
       width = 18, height = 10, units = "cm", dpi = 300)

