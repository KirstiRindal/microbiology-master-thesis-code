#ANOVA-Like Differential Expression (ALDEx2)  

#vignette: https://www.bioconductor.org/packages/devel/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html

library(ALDEx2)

#example from vignette:
#subset only the last 400 features for efficiency
data(selex)
selex.sub <- selex[1200:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
x.aldex <- aldex(selex.sub, conds, mc.samples=128, test="t", effect=TRUE, 
                 include.sample.summary=TRUE, denom="all", verbose=TRUE, 
                 paired.test=FALSE, gamma=NULL)

# note default is FDR of 0.05
par(mfrow=c(1,3))
aldex.plot(x.aldex, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(x.aldex, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(x.aldex, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot') 

#identify the species that are differential abundant between the two conditions. 

#-------------------------------------------------------------------------------#

#1: the total IBS group vs HC:
#construct otu count table and condition information for the analysis:
OTU_counts <- raw_baseline_data|>
  column_to_rownames("p_ID") |>
  select(-(group:BMI)) |>
  t()

group <- meta_data$group 
diagnose <- meta_data$diagnose

# #compute an aldex object:
aldex_obj <- aldex(OTU_counts, group, mc.samples=128, test="t", effect=TRUE, 
                 include.sample.summary=TRUE, denom="all", verbose=TRUE, 
                 paired.test=FALSE, gamma=NULL)

#plot:
?aldex.plot
#Bland-Altman (MA plot): shows the relationship between (relative) abundance and difference.
#Effect plot: shows the relationship between difference and dispersion: the lines are equal difference and dispersion.
#Volcano plot: the lines represent a posterior predictive p-value of 0.001 and 1.5 fold difference.
#In all plots features that are not significant are in grey or black.
#Features that are statistically significant are in red.
#The log-ratio abundance axis is the clr value for the feature.

par(mfrow=c(1,3))
MA_plot <- aldex.plot(aldex_obj, type="MA", test="welch", main='Bland-Altman plot')
MW_plot <- aldex.plot(aldex_obj, type="MW", test="welch", main='Effect plot plot')
aldex.plot(aldex_obj, type="volcano", test="welch", main='volcano plot')
#no sig diff OTUs.

#2:IBS-D vs HC:
OTUs_IBS_D <- raw_baseline_data |> 
  filter(diagnose %in% c("IBS-D", "HC")) |> 
  tibble::column_to_rownames("p_ID") |>
  select(-(group:BMI)) |> 
  t()

diagnose_D <- meta_data |>  
  filter(diagnose %in% c("IBS-D", "HC")) |> 
  select(diagnose) 

#compute an aldex object:
alde_ojc_D <- aldex(OTUs_IBS_D, diagnose_D$diagnose, mc.samples=128, test="t", 
                    effect = TRUE, include.sample.summary = TRUE, denome = "all", 
                    paired.test = FALSE, verbose = TRUE) 
par(mfrow=c(1,3))
aldex.plot(alde_ojc_D, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(alde_ojc_D, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(alde_ojc_D, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot') 

#3: IBS-C vs HC:
diagnose_C <- meta_data |>  
  filter(diagnose %in% c("IBS-C", "HC")) |> 
  select(diagnose)

OTUs_IBS_C <- raw_baseline_data |> 
  filter(diagnose %in% c("IBS-C", "HC")) |> 
  select(-c(1:7)) |> 
  t()

#compute an aldex object:
alde_ojc_C <- aldex(OTUs_IBS_C, diagnose_C$diagnose, mc.samples=128, test="t", effect = TRUE, 
                    include.sample.summary = TRUE, paired.test = FALSE) 

par(mfrow=c(1,3))
aldex.plot(alde_ojc_C, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(OTUs_IBS_C, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(OTUs_IBS_C, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot')
#no sig diff OTUs.

#4: IBS-M
#IBS-M vs HC:
diagnose_M <- meta_data |>  
  filter(diagnose %in% c("IBS-M", "HC")) |> 
  select(diagnose)

OTUs_IBS_M <- raw_baseline_data |> 
  filter(diagnose %in% c("IBS-M", "HC")) |> 
  select(-c(1:7)) |> 
  t()

#compute an aldex object:
alde_ojc_M <- aldex(OTUs_IBS_M, diagnose_M$diagnose, mc.samples=2, ntest="t", effect = TRUE, 
                    include.sample.summary = TRUE, paired.test = FALSE) 

par(mfrow=c(1,3))
aldex.plot(alde_all_M, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(alde_all_M, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(alde_all_M, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot')
#no sig diff OTUs.

#high IBS-SSS (>300) vs HC 
IBS_cat <- IBS_SSS_data |>  #threshold = 300
  select(IBS_SSS_cat) |> 
  filter(IBS_SSS_cat %in% c("HC", "high")) 

OTUs_SSS <- IBS_SSS_data |> 
  filter(IBS_SSS_cat %in% c("HC", "high")) |>
  select(-c(IBS_SSS:BMI)) |> 
  t()

#compute an aldex object:
alde_ojc_SSS <- aldex(OTUs_SSS, IBS_cat$IBS_SSS_cat, mc.samples=2, ntest="t", effect = TRUE, 
                    include.sample.summary = TRUE, paired.test = FALSE) 

par(mfrow=c(1,3))
aldex.plot(alde_all_SSS, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(alde_all_SSS, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(alde_all_SSS, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot')
#no sig diff OTUs.

#low vs high HADS depression score 
HADS_dep_cat <- HADS |>  
  select(HADS_dep_cat)

OTUs_HADS_dep <- HADS |> 
  select(-c(p_ID:BMI)) |> 
  t()

#compute an aldex object:
alde_ojc_HADS_dep_cat <- aldex(OTUs_HADS_dep, HADS_dep_cat$HADS_dep_cat, mc.samples=2, ntest="t", effect = TRUE, 
                      include.sample.summary = TRUE, paired.test = FALSE) 

par(mfrow=c(1,3))
aldex.plot(alde_all_HADS_dep, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(alde_all_HADS_dep, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(alde_all_HADS_dep, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot')

#identify the one significant OTU: 
sig_OTU <-  alde_all_HADS_dep |> 
  #filter((alde_all_HADS_dep$we.ep< 0.05)) |> #12. Expected p-value of Welch's t-test, a posterior predictive p-value
  #filter((alde_all_HADS_dep$we.eBH < 0.05)) |> #5: Expected Benjamini-Hochberg corrected p-value of Welch's t test
  #filter((alde_all_HADS_dep$wi.ep < 0.05)) |> #4: Expected p-value of Wilcoxon rank test
  filter((alde_all_HADS_dep$wi.eBH < 0.05)) #0: Expected Benjamini-Hochberg corrected p-value of Wilcoxon test

#high vs low HADS anxiety score
HADS_anx_cat <- HADS |>  
  select(HADS_anx_cat)

OTUs_HADS_anx <- HADS |> 
  select(-c(p_ID:BMI)) |> 
  t()

#compute an aldex object:
alde_ojc_HADS_anx <- aldex(OTUs_HADS_anx, HADS_anx_cat$HADS_anx_cat, mc.samples=2, ntest="t", effect = TRUE, 
                      include.sample.summary = TRUE, paired.test = FALSE) 

par(mfrow=c(1,3))
aldex.plot(alde_all_HADS_anx, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(alde_all_HADS_anx, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(alde_all_HADS_anx, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot')
#no sig diff OTUs.

#high vs low fatigue score 
FSS_cat <- FSS |>  #threshold = 4
  select(FSS_cat)

OTUs_FSS <- FSS |> 
  select(-c(p_ID:BMI)) |> 
  t()

#compute an aldex object:
alde_ojc_FSS <- aldex(OTUs_FSS, FSS_cat$FSS_cat, mc.samples=2, ntest="t", effect = TRUE, 
                      include.sample.summary = TRUE, paired.test = FALSE) 

par(mfrow=c(1,3))
aldex.plot(alde_all_FSS, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(alde_all_FSS, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(alde_all_FSS, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot')
#no sig diff OTUs.


