#Differential abundance analysis (DAA) using ZicoSeq

#from: https://cran.r-project.org/web/packages/GUniFrac/vignettes/ZicoSeq.html 
#method article: https://link.springer.com/article/10.1186/s40168-022-01320-0
#study article: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015216 
#example: https://search.r-project.org/CRAN/refmans/GUniFrac/html/SimulateMSeq.html 

library(GUniFrac)

set.seed(19980311)

#1: IBS group vs HC:

#raw OTU count data:
ZicoSeq_OTU_data <- raw_baseline_data |>
  select(-c(p_ID:BMI)) |>
  select(where(~ !all(. == 0))) |>
  drop_na()

#meta data:
ZicoSeq_meta_data_g <- meta_data |>
  select(p_ID, group, diagnose, sex, age, BMI) |>
  mutate(group = factor(group)) |>
  as.data.frame()

ZicoSeq.obj_g <- ZicoSeq(meta.dat = ZicoSeq_meta_data_g, feature.dat = t(ZicoSeq_OTU_data),
                         grp.name = 'group', adj.name = c("sex", "age", "BMI"), feature.dat.type = "count",
                         # Filter to remove rare taxa
                         prev.filter = 0.1, mean.abund.filter = 0,
                         max.abund.filter = 0, min.prop = 0,
                         # Winsorization to replace outliers
                         is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                         # Posterior sampling
                         is.post.sample = TRUE, post.sample.no = 25,
                         # Use the square-root transformation
                         link.func = list(function (x) x^0.5), stats.combine.func = max,
                         # Permutation-based multiple testing correction
                         perm.no = 999,  strata = NULL,
                         # Reference-based multiple stage normalization
                         ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                         # Family-wise error rate control
                         is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_g$p.adj.fdr <= 0.05)

#visualization:
ZicoSeq_g_plot <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_g, pvalue.type = 'p.adj.fdr',
                               cutoff = 0.05, text.size = 10, out.dir = NULL) +
  scale_y_continuous(limits = c(0, 1.35))
  #labs(tag = "IBS vs HC")
  #ggtitle("")
ZicoSeq_g_plot

#2: IBS subgroups vs HC:
#meta data:
ZicoSeq_meta_data <- meta_data |>
  select(p_ID, group, diagnose, sex, age, BMI) |>
  mutate(diagnose = factor(diagnose)) |>
  as.data.frame()

ZicoSeq.obj_d <- ZicoSeq(meta.dat = ZicoSeq_meta_data, feature.dat = t(ZicoSeq_OTU_data),
                         grp.name = 'diagnose', adj.name = c("sex", "age", "BMI"),
                         feature.dat.type = "count",
                         # Filter to remove rare taxa
                         prev.filter = 0.1, mean.abund.filter = 0,
                         max.abund.filter = 0, min.prop = 0,
                         # Winsorization to replace outliers
                         is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                         # Posterior sampling
                         is.post.sample = TRUE, post.sample.no = 25,
                         # Use the square-root transformation
                         link.func = list(function (x) x^0.5), stats.combine.func = max,
                         # Permutation-based multiple testing correction
                         perm.no = 999,  strata = NULL,
                         # Reference-based multiple stage normalization
                         ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                         # Family-wise error rate control
                         is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_d$p.adj.fdr <= 0.05)

#visualization:
ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_d, pvalue.type = 'p.adj.fdr', cutoff = 0.05,
             text.size = 10, out.dir = NULL)

#3: IBS-D vs HC: 
#raw OTU count data: 
ZicoSeq_data_D <- raw_baseline_data |>  
  filter(diagnose %in% c("IBS-D", "HC")) |> 
  select(-c(p_ID:BMI)) |> 
  select(where(~ !all(. == 0))) |> 
  drop_na() 

#meta data: 
ZicoSeq_meta_data_D <- meta_data |> 
  select(p_ID, group, diagnose, sex, age, BMI) |> 
  filter(diagnose %in%  c("IBS-D", "HC")) |>
  #mutate(diagnose = factor(diagnose)) |> 
  as.data.frame() 

ZicoSeq.obj_D <- ZicoSeq(meta.dat = ZicoSeq_meta_data_D, feature.dat = t(ZicoSeq_data_D), 
                       grp.name = 'diagnose', adj.name = c("sex", "age", "BMI"), feature.dat.type = "count", 
                       # Filter to remove rare taxa
                       prev.filter = 0.1, mean.abund.filter = 0,  
                       max.abund.filter = 0, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_D$p.adj.fdr <= 0.05)

#visualization: 
ZicoSeq.obj_D_plot <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_D, pvalue.type = 'p.adj.fdr', 
                                   cutoff = 0.05, text.size = 10, out.dir = NULL) +
  ggtitle("HC vs IBS-D") +
  scale_y_continuous(limits = c(0, 1.35)) +
  scale_x_continuous(limits = c(-0.05, 0.1)) +
  theme(legend.position = "none") 
ZicoSeq.obj_D_plot
  
#4: IBS-C vs HC: 
#raw OTU count data: 
ZicoSeq_data_C <- raw_baseline_data |>  
  filter(diagnose %in% c("IBS-C", "HC")) |> 
  select(-c(p_ID:BMI)) |> 
  select(where(~ !all(. == 0))) |> 
  drop_na() 

#meta data: 
ZicoSeq_meta_data_C <- meta_data |> 
  select(p_ID, group, diagnose, sex, age, BMI) |> 
  filter(diagnose %in%  c("IBS-C", "HC")) |>
  #mutate(diagnose = factor(diagnose)) |> 
  as.data.frame() 

ZicoSeq.obj_C <- ZicoSeq(meta.dat = ZicoSeq_meta_data_C, feature.dat = t(ZicoSeq_data_C), 
                       grp.name = 'diagnose', adj.name = c("sex", "age", "BMI"), feature.dat.type = "count", 
                       # Filter to remove rare taxa
                       prev.filter = 0.1, mean.abund.filter = 0,  
                       max.abund.filter = 0, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_C$p.adj.fdr <= 0.05)

#visualization: 
ZicoSeq.obj_C_plot <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_C, pvalue.type = 'p.adj.fdr', 
                                   cutoff = 0.05, text.size = 10, out.dir = NULL) +
  ggtitle("HC vs IBS-C") +
  scale_y_continuous(limits = c(0, 1.35)) +
  scale_x_continuous(limits = c(-0.05, 0.1)) +
  theme(legend.position = "none") 
ZicoSeq.obj_C_plot

#5: IBS-M vs HC: 
#raw OTU count data: 
ZicoSeq_data_M <- raw_baseline_data |>  
  filter(diagnose %in% c("IBS-M", "HC")) |> 
  select(-c(p_ID:BMI)) |> 
  select(where(~ !all(. == 0))) |> 
  drop_na() 

#meta data: 
ZicoSeq_meta_data_M <- meta_data |> 
  select(p_ID, group, diagnose, sex, age, BMI) |> 
  filter(diagnose %in%  c("IBS-M", "HC")) |>
  #mutate(diagnose = factor(diagnose)) |> 
  as.data.frame() 

ZicoSeq.obj_M <- ZicoSeq(meta.dat = ZicoSeq_meta_data_M, feature.dat = t(ZicoSeq_data_M), 
                         grp.name = 'diagnose', adj.name = c("sex", "age", "BMI"), feature.dat.type = "count", 
                         # Filter to remove rare taxa
                         prev.filter = 0.1, mean.abund.filter = 0,  
                         max.abund.filter = 0, min.prop = 0, 
                         # Winsorization to replace outliers
                         is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                         # Posterior sampling 
                         is.post.sample = TRUE, post.sample.no = 25, 
                         # Use the square-root transformation
                         link.func = list(function (x) x^0.5), stats.combine.func = max,
                         # Permutation-based multiple testing correction
                         perm.no = 999,  strata = NULL, 
                         # Reference-based multiple stage normalization
                         ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                         # Family-wise error rate control
                         is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_M$p.adj.fdr <= 0.05)

#visualization: 
ZicoSeq.obj_M_plot <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_M, pvalue.type = 'p.adj.fdr', 
                                   cutoff = 0.05, text.size = 10, out.dir = NULL) +
  ggtitle("HC vs IBS-M") +
  scale_y_continuous(limits = c(0, 1.35)) +
  scale_x_continuous(limits = c(-0.05, 0.1)) + 
  theme(legend.position = "none") 
ZicoSeq.obj_M_plot

combined_diagnose_plot <- ZicoSeq.obj_D_plot | ZicoSeq.obj_M_plot | ZicoSeq.obj_C_plot
combined_diagnose_plot

#-------------------------------------------------------------------------------#

#IBS-SSS >= 300 vs HC: 

#6: IBS-SSS on a continuous scale: 
comm_IBS_SSS <- IBS_SSS_data |> 
  column_to_rownames(var = "p_ID") |> 
  #only include the "high" and "HC" category: 
  filter(IBS_SSS_cat %in% c("high", "HC")) |>
  select(-c(IBS_SSS:BMI)) |> 
  select(where(~ !all(. == 0))) |> 
  drop_na() 

#meta data: 
ZicoSeq_meta_data_SSS <- IBS_SSS_data |> 
  column_to_rownames(var = "p_ID") |>
  filter(IBS_SSS_cat %in% c("high", "HC")) |>
  select(diagnose, sex, age, BMI, IBS_SSS, IBS_SSS_cat) |> 
  as.data.frame() 

#IBS-SSS on a continuous scale: 
ZicoSeq.obj_SSS_con <- ZicoSeq(meta.dat = ZicoSeq_meta_data_SSS, 
                               feature.dat = t(comm_IBS_SSS), 
                           grp.name = 'IBS_SSS', 
                           adj.name = c("diagnose", "sex", "age", "BMI"), 
                           feature.dat.type = "count", 
                           # Filter to remove rare taxa
                           prev.filter = 0.1, mean.abund.filter = 0,  
                           max.abund.filter = 0, min.prop = 0, 
                           # Winsorization to replace outliers
                           is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                           # Posterior sampling 
                           is.post.sample = TRUE, post.sample.no = 25, 
                           # Use the square-root transformation
                           link.func = list(function (x) x^0.5), stats.combine.func = max,
                           # Permutation-based multiple testing correction
                           perm.no = 999,  strata = NULL, 
                           # Reference-based multiple stage normalization
                           ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                           # Family-wise error rate control
                           is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_SSS_con$p.adj.fdr <= 0.05)

#visualization: 
ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_SSS_con, pvalue.type = 'p.adj.fdr', 
             cutoff = 0.05, text.size = 10, out.dir = NULL) +
  scale_y_continuous(limits = c(0, 1.35)) +
  scale_x_continuous(limits = c(-0.05, 0.1)) +
  #remove title 
  theme(title = element_blank()) 

# #IBS-SSS in categories:
# ZicoSeq.obj_SSS_cat <- ZicoSeq(meta.dat = ZicoSeq_meta_data_SSS, feature.dat = t(ZicoSeq_data_SSS),
#                            grp.name = 'IBS_SSS_cat', adj.name = c("sex", "age", "BMI"), feature.dat.type = "count",
#                            # Filter to remove rare taxa
#                            prev.filter = 0.1, mean.abund.filter = 0,
#                            max.abund.filter = 0.002, min.prop = 0,
#                            # Winsorization to replace outliers
#                            is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
#                            # Posterior sampling
#                            is.post.sample = TRUE, post.sample.no = 25,
#                            # Use the square-root transformation
#                            link.func = list(function (x) x^0.5), stats.combine.func = max,
#                            # Permutation-based multiple testing correction
#                            perm.no = 99,  strata = NULL,
#                            # Reference-based multiple stage normalization
#                            ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
#                            # Family-wise error rate control
#                            is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# 
# #detected differential OTUs
# which(ZicoSeq.obj_SSS_cat$p.adj.fdr <= 0.05)
# 
# #visualization:
# ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_SSS_cat, pvalue.type = 'p.adj.fdr', cutoff = 0.05,
#              text.size = 10, out.dir = NULL)

#-------------------------------------------------------------------------------#

#HADS: 
#raw OTU count data: 
comm_HADS <- HADS |>  
  column_to_rownames(var = "p_ID") |> 
  select(-c(HADS_anx:BMI)) |> 
  select(where(~ !all(. == 0))) |> 
  drop_na() 

#meta data: 
ZicoSeq_meta_data_HADS <- HADS |> 
  select(p_ID, diagnose, sex, age, BMI, HADS_dep, HADS_dep_cat, HADS_anx, HADS_anx_cat) |> 
  as.data.frame() 

#8: HADS depression on a continuous scale:
ZicoSeq.obj_HADS_d_con <- ZicoSeq(meta.dat = ZicoSeq_meta_data_HADS, feature.dat = t(comm_HADS),
                           grp.name = 'HADS_dep', adj.name = c("diagnose", "sex", "age", "BMI"), 
                           feature.dat.type = "count",
                           # Filter to remove rare taxa
                           prev.filter = 0.1, mean.abund.filter = 0,
                           max.abund.filter = 0, min.prop = 0,
                           # Winsorization to replace outliers
                           is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                           # Posterior sampling
                           is.post.sample = TRUE, post.sample.no = 25,
                           # Use the square-root transformation
                           link.func = list(function (x) x^0.5), stats.combine.func = max,
                           # Permutation-based multiple testing correction
                           perm.no = 999,  strata = NULL,
                           # Reference-based multiple stage normalization
                           ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                           # Family-wise error rate control
                           is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_HADS_d_con$p.adj.fdr <= 0.05)

#visualization:
HADS_dep_plot <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_HADS_d_con, pvalue.type = 'p.adj.fdr', cutoff = 0.05,
             text.size = 10, out.dir = NULL) +
  scale_y_continuous(limits = c(0, 1.35)) +
  scale_x_continuous(limits = c(-0.05, 0.1)) +
  ggtitle("HADS depression") +
  theme(legend.position  = "none")
HADS_dep_plot

# #9: HADS depression in categories:
# ZicoSeq.obj_HADS_d_cat <- ZicoSeq(meta.dat = ZicoSeq_meta_data_HADS, feature.dat = t(ZicoSeq_data_HADS),
#                          grp.name = 'HADS_dep_cat', adj.name = c("diagnose", "sex", "age", "BMI"), feature.dat.type = "count",
#                          # Filter to remove rare taxa
#                          prev.filter = 0.1, mean.abund.filter = 0,
#                          max.abund.filter = 0, min.prop = 0,
#                          # Winsorization to replace outliers
#                          is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
#                          # Posterior sampling
#                          is.post.sample = TRUE, post.sample.no = 25,
#                          # Use the square-root transformation
#                          link.func = list(function (x) x^0.5), stats.combine.func = max,
#                          # Permutation-based multiple testing correction
#                          perm.no = 999,  strata = NULL,
#                          # Reference-based multiple stage normalization
#                          ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
#                          # Family-wise error rate control
#                          is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# 
# #detected differential OTUs
# which(ZicoSeq.obj_HADS_d_cat$p.adj.fdr <= 0.05)
# 
# #visualization:
# HADS_d_plot_cat <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_HADS_d_cat, pvalue.type = 'p.adj.fdr', cutoff = 0.05,
#              text.size = 10, out.dir = NULL) +
#   ggtitle("HADS depression") +
#   theme(legend.position  = "none") +
#   scale_y_continuous(limits = c(0, 1.35)) + 
#   scale_x_continuous(limits = c(-0.05, 0.1))
# HADS_d_plot_cat

#10: HADS anxiety on a continuous scale:
ZicoSeq.obj_HADS_a_con <- ZicoSeq(meta.dat = ZicoSeq_meta_data_HADS, feature.dat = t(comm_HADS),
                                  grp.name = 'HADS_anx', adj.name = c("diagnose", "sex", "age", "BMI"), 
                                  feature.dat.type = "count",
                                  # Filter to remove rare taxa
                                  prev.filter = 0.1, mean.abund.filter = 0,
                                  max.abund.filter = 0.002, min.prop = 0,
                                  # Winsorization to replace outliers
                                  is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                                  # Posterior sampling
                                  is.post.sample = TRUE, post.sample.no = 25,
                                  # Use the square-root transformation
                                  link.func = list(function (x) x^0.5), stats.combine.func = max,
                                  # Permutation-based multiple testing correction
                                  perm.no = 999,  strata = NULL,
                                  # Reference-based multiple stage normalization
                                  ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                                  # Family-wise error rate control
                                  is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_HADS_a_con$p.adj.fdr <= 0.05)

#visualization:
HADS_anx_plot <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_HADS_a_con, pvalue.type = 'p.adj.fdr', cutoff = 0.05,
             text.size = 10, out.dir = NULL) +
  scale_y_continuous(limits = c(0, 1.35)) +
  scale_x_continuous(limits = c(-0.05, 0.1)) +
  ggtitle("HADS anxiety") +
  #remove legend
  theme(legend.position = "none")
HADS_anx_plot

# #11: HADS anxiety in categories:
# ZicoSeq.obj_HADS_a_cat <- ZicoSeq(meta.dat = ZicoSeq_meta_data_HADS, feature.dat = t(ZicoSeq_data_HADS),
#                                   grp.name = 'HADS_anx_cat', adj.name = c("diagnose", "sex", "age", "BMI"), feature.dat.type = "count",
#                                   # Filter to remove rare taxa
#                                   prev.filter = 0.1, mean.abund.filter = 0,
#                                   max.abund.filter = 0, min.prop = 0,
#                                   # Winsorization to replace outliers
#                                   is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
#                                   # Posterior sampling
#                                   is.post.sample = TRUE, post.sample.no = 25,
#                                   # Use the square-root transformation
#                                   link.func = list(function (x) x^0.5), stats.combine.func = max,
#                                   # Permutation-based multiple testing correction
#                                   perm.no = 999,  strata = NULL,
#                                   # Reference-based multiple stage normalization
#                                   ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
#                                   # Family-wise error rate control
#                                   is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# 
# #detected differential OTUs
# which(ZicoSeq.obj_HADS_a_cat$p.adj.fdr <= 0.05)
# 
# #visualization:
# HADS_a_plot_cat <- ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_HADS_a_cat, pvalue.type = 'p.adj.fdr', cutoff = 0.05,
#              text.size = 10, out.dir = NULL) +
#   ggtitle("HADS anxiety") +
#   theme(legend.position  = "none") 
#   # scale_y_continuous(limits = c(0, 1.35)) +
#   # scale_x_continuous(limits = c(-0.05, 0.1))
# HADS_a_plot_cat 


#Fatigue: 
#raw OTU count data: 
comm_FSS <- FSS |>  
  column_to_rownames(var = "p_ID") |>
  select(-c(FSS_score_BL:BMI)) |> 
  select(where(~ !all(. == 0))) |> 
  drop_na() 

#meta data: 
ZicoSeq_meta_data_FSS <- FSS |> 
  select(p_ID, diagnose, sex, age, BMI, FSS_score_BL, FSS_cat) |> 
  as.data.frame() 

#14: fatugue on a continuous scale:
ZicoSeq.obj_FSS_con <- ZicoSeq(meta.dat = ZicoSeq_meta_data_FSS, feature.dat = t(comm_FSS),
                               grp.name = 'FSS_score_BL', adj.name = c("diagnose", "sex", "age", "BMI"), feature.dat.type = "count",
                               # Filter to remove rare taxa
                               prev.filter = 0.1, mean.abund.filter = 0,
                               max.abund.filter = 0, min.prop = 0,
                               # Winsorization to replace outliers
                               is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                               # Posterior sampling
                               is.post.sample = TRUE, post.sample.no = 25,
                               # Use the square-root transformation
                               link.func = list(function (x) x^0.5), stats.combine.func = max,
                               # Permutation-based multiple testing correction
                               perm.no = 999,  strata = NULL,
                               # Reference-based multiple stage normalization
                               ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                               # Family-wise error rate control
                               is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)

#detected differential OTUs
which(ZicoSeq.obj_FSS_con$p.adj.fdr <= 0.05)

#visualization:
ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_FSS_con, pvalue.type = 'p.adj.fdr', cutoff = 0.05,
             text.size = 10, out.dir = NULL) +
  scale_y_continuous(limits = c(0, 1.35))
  

# #15: FSS in categories, with threshold at 4:
# ZicoSeq.obj_FSS_cat <- ZicoSeq(meta.dat = ZicoSeq_meta_data_FSS, feature.dat = t(ZicoSeq_data),
#                                grp.name = 'FSS_cat', adj.name = c("diagnose", "sex", "age", "BMI"), feature.dat.type = "count",
#                                # Filter to remove rare taxa
#                                prev.filter = 0.1, mean.abund.filter = 0,
#                                max.abund.filter = 0.002, min.prop = 0,
#                                # Winsorization to replace outliers
#                                is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
#                                # Posterior sampling
#                                is.post.sample = TRUE, post.sample.no = 25,
#                                # Use the square-root transformation
#                                link.func = list(function (x) x^0.5), stats.combine.func = max,
#                                # Permutation-based multiple testing correction
#                                perm.no = 99,  strata = NULL,
#                                # Reference-based multiple stage normalization
#                                ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
#                                # Family-wise error rate control
#                                is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# 
# #detected differential OTUs
# which(ZicoSeq.obj_FSS_cat$p.adj.fdr <= 0.05)
# 
# #visualization:
# ZicoSeq.plot(ZicoSeq.obj = ZicoSeq.obj_FSS_cat, pvalue.type = 'p.adj.fdr', 
#              cutoff = 0.05, text.size = 10, out.dir = NULL)
