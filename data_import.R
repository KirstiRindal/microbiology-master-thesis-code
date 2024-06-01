
#load packages:  
library(readxl)
library(tidyverse)
library(patchwork)
library(vegan)
library(ggvegan)
library(ggsignif)
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
conflict_prefer("map", winner = "purrr")

#set default theme:
theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#import data: 

#the identification key: 
key <- read_excel("data") 

#meta data: 
meta_data <- read_excel("data") |> 
  select(p_id, diagnose, age, sex, BMI) |> 
  left_join(key, by = c("p_id" = "p_id")) |> 
  #column_to_rownames("p_id") |>
  relocate(p_id, ONT_code, group) |> 
  drop_na()

#demographic_data$p_id <- as.character(demographic_data$p_id)  #change vector to character. 
#meta_data$age <- round(meta_data$age, 1) #round Age value to zero digits. 
#meta_data$BMI <- round(meta_data$BMI, 1) #round BMI value to zero digits.

#raw OTU count data: 
raw_baseline_data <- read.table("data") 
raw_baseline_data <- as.data.frame(t(raw_baseline_data))  
raw_baseline_data <- rownames_to_column(raw_baseline_data, "ONT_code")   

raw_baseline_data <- raw_baseline_data |> 
  left_join(meta_data, by = c("ONT_code" = "ONT_code")) |> 
  select(-ONT_code) |>
  relocate(p_id, group, diagnose, sex, age, BMI) |> 
  select(where(~ !all(. == 0))) |> #remove the columns with only zeros.
  drop_na() #92 patients have complete meta data and microbiota data. 

comm <- raw_baseline_data |> 
  column_to_rownames("p_id") |> 
  select(-(group:BMI))

#filtering: 
#before filtering: 555 OTUs. 

#pivot_longer the data and remove rare taxa:  
filtered_long_raw_baseline_data <- pivot_longer(raw_baseline_data, 
                                                cols = -c(p_id, group, diagnose, sex, age, BMI), 
                                                names_to = "bacteria_OTU",
                                                values_to = "abundance") |> 
  #remove bacteria OTUs that have a maximum abundance across all 92 samples >=3
  filter(max(abundance) >= 3, .by = "bacteria_OTU") 
#after filtering: 447 OTUs, which means that 108 OTUs had a maximum abundance >=3
  
#back to wide format: 
filtered_raw_baseline_data <- filtered_long_raw_baseline_data |> 
  pivot_wider(names_from = bacteria_OTU, values_from = abundance, values_fill = 0) |> 
  as.data.frame()

filtered_comm <- filtered_raw_baseline_data |> 
  column_to_rownames("p_id") |> 
  select(-(group:BMI))
  
# #pivot longer the data and 10% prevalence filtering:  
# prev_filtered_long_raw_baseline_data <- pivot_longer(raw_baseline_data, 
#                                                 cols = -c(p_id, group, diagnose, sex, age, BMI), 
#                                                 names_to = "bacteria_OTU",
#                                                 values_to = "abundance") |> 
#   #remember to filter subgroup before prevalence filtering! 
#   mutate(n_patients = n_distinct(p_id)) |>
#   group_by(bacteria_OTU) |> 
#   filter(sum(abundance > 0) > 0.1 * n_patients) |> #236 
#   select(-n_patients) 

# #back to wide format: 
# prev_filtered_raw_baseline_data <- prev_filtered_long_raw_baseline_data |> 
#   pivot_wider(names_from = bacteria_OTU, values_from = abundance, values_fill = 0) |> 
#   as.data.frame()  

#import IBS symthome severity score (IBS-SSS) data: 
IBS_SSS_data <- read_excel("data") |> 
  mutate(IBS_SSS = score) |> 
  select(p_id, IBS_SSS) |> 
  drop_na() |> #94 patients have IBS-SSS data. 
  #categorize the IBS-SSS score into high and low groups:
  mutate(IBS_SSS_cat = case_when(
    IBS_SSS >= 300 ~ "high",
    IBS_SSS < 300 & IBS_SSS >= 175 ~ "low",
    IBS_SSS < 175 ~ "HC")) |> 
  #join with raw_baseline_data
  left_join(raw_baseline_data, by = c("p_id" = "p_id")) |> 
  column_to_rownames("p_id") |>
  drop_na() #84 

#count the number of the different groups: 
IBS_SSS_data |> 
  count(diagnose) 

#count the number in IBS_SSS_cat high: 
IBS_SSS_data |> 
  count(IBS_SSS_cat) #high: 16, low: 68

#which diagnosis make out the high group? 
IBS_SSS_data |> 
  filter(IBS_SSS_cat == "high") |> 
  count(diagnose) 

#import psychological data: 

#HADS (depression, anxiety):
HADS <- read_excel("data") |> 
  drop_na() |> #94 patients have HADS data. 
  #categorize the HADS total score into high and low: 
  mutate(HADS_dep_cat = case_when(
    HADS_dep >= 8 ~ "case",
    HADS_dep < 8 ~ "non-case")) |> 
  mutate(HADS_anx_cat = case_when(
    HADS_anx >= 8 ~ "case",
    HADS_anx < 8 ~ "non-case")) |>  #94 
  #join with raw_baseline_data
  left_join(raw_baseline_data, by = c("p_id" = "p_id")) |>
  drop_na() #85 

#count the number of the different groups:
HADS |> 
  count(diagnose)

#count the number in HADS_dep_cat high:
HADS |> 
  count(HADS_dep_cat)

HADS |> 
  count(HADS_anx_cat) 

#CFQ-11/FSS (fatigue) 
FSS <- read_excel("data") |> 
  drop_na() |> 
  mutate(FSS_cat = case_when(
    FSS_score_BL >= 4 ~ "case",
    FSS_score_BL < 4 ~ "non_case")) |> 
  #join with raw_baseline_data
  left_join(raw_baseline_data, by = c("p_id" = "p_id")) |>
  drop_na() 

#count the number in each diagnose group: 
FSS |>
  count(diagnose)

#count the number in FSS_cat high:
FSS |> 
  count(FSS_cat) #high: 45, low: 43 

#max score: 
max(FSS$FSS_score_BL) #11

# just_psych_data <- HADS |> 
#   inner_join(VSI, by = c("p_id" = "p_id")) |> 
#   inner_join(FSS, by = c("p_id" = "p_id")) |> 
#   #inner_join(BIS, by = c("p_id" = "p_id")) |> 
#   drop_na()
# 
# #meta data with all the psychological variables:
# all_psych_data <- meta_data |> 
#   inner_join(HADS, by = c("p_id" = "p_id")) |> 
#   inner_join(VSI, by = c("p_id" = "p_id")) |> 
#   inner_join(FSS, by = c("p_id" = "p_id")) |>
#   #inner_join(BIS, by = c("p_id" = "p_id")) |> 
#   drop_na()
# 
# all_data <- raw_baseline_data |> 
#   inner_join(IBS_SSS_data, by = c("p_id" = "p_id")) |>
#   inner_join(HADS, by = c("p_id" = "p_id")) |> 
#   inner_join(VSI, by = c("p_id" = "p_id")) |> 
#   inner_join(FSS, by = c("p_id" = "p_id")) |>
#   #inner_join(BIS, by = c("p_id" = "p_id")) |> 
#   #relocate(p_id, group, diagnose, sex, age, BMI, IBS_SSS, IBS_SSS_cat) |> 
#   relocate(p_id, group, diagnose, sex, age, BMI, IBS_SSS, IBS_SSS_cat, 
#            HADS_anx, HADS_dep, HADS_total, HADS_dep_cat, HADS_anx_cat, HADS_total_cat, 
#           VSI_score_BL, VSI_cat,  FSS_score_BL, FSS_cat) |> 
#   drop_na() #77
