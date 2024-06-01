#Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOMBC2)

#https://www.nature.com/articles/s41467-022-28034-z#citeas
#https://www.nature.com/articles/s41467-020-17041-7
#tutorial: https://microbiome.github.io/course_2022_radboud/differential-abundance-analysis-demo.html#ancom-bc
#tutorial: https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html 
#example from R: https://assets.researchsquare.com/files/rs-4144233/v1/dfa6a3a11cfd0027ea62538b.pdf 

library(ANCOMBC)
library(phyloseq)
library(ggtext)
library(ggrepel)
library(patchwork)

#the data needs to be phyloseq object; consisting of abundance table and sample table.
#https://joey711.github.io/phyloseq/import-data
#https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

#abundance table: OTUs as rows and samples as columns 
otu_abund <- comm |> 
  t() 

#sample info table: samples as rows. 
sample_info <- meta_data |> 
  column_to_rownames("p_ID") |> 
  as.data.frame() #is class "data.frame" so good. 

#taxonomy table: optional to include? 
#two identical rows "OTUs" and "OTUs"
#needs to have rownames! 

OTU = otu_table(otu_abund, taxa_are_rows = TRUE)
samples = sample_data(sample_info)
#TAX = tax_table(otu_tax)

phyloseq_obj <- phyloseq(OTU, samples)

set.seed(19981103)

ANCOMBC_analysis = ancombc2(data = phyloseq_obj, 
                            assay_name = "counts",
                            tax_level = NULL, 
                            fix_formula = "diagnose + age + sex + BMI", 
                            rand_formula = NULL, #I don't have random effects. 
                            p_adj_method = "fdr", #FDR correction.
                            pseudo_sens = TRUE,
                            prv_cut = 0.10, #10% prevalence cut-off. 
                            lib_cut = 2500,  #library size cut-off. The smallest library size is 2386? 
                            s0_perc = 0.05, #5% sparsity cut-off.
                            group = "diagnose", 
                            struc_zero = TRUE, #expected, so TRUE. 
                            neg_lb = FALSE, #negative lower bound. 
                            alpha = 0.05, #significance level.
                            n_cl = 1, #number of clusters: 1 = no parallel computing.
                            verbose = TRUE, 
                            global = FALSE, #global test.
                            pairwise = FALSE, #multiple pairwise comparisons. 
                            dunnet = TRUE, #Dunnett's type of test: Multiple pairwise comparisons against a pre-specified group. 
                            #trend = TRUE, #trend test???, not working.  
                            iter_control = list(tol = 1e-02, max_iter = 20, verbose = FALSE), #number of interations. Note that setting max_iter = 1 and B = 1 is only for the sake of speed. Use default or larger values for max_iter and B for better performance
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(), 
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                            trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)) 


#The results encompass: 
#1) log fold changes, 
#2) standard errors, 
#3) test statistics, 
#4) p-values, 
#5) adjusted p-values, 
#6) indicators denoting whether the taxon is differentially abundant (TRUE) or not (FALSE)
#7) indicators denoting whether the taxon passed the sensitivity analysis (TRUE) or not (FALSE).

#results:
#the primary results: identifies taxa with differential abundance based on the chosen covariate.
res_prim = ANCOMBC_analysis$res 

# #1: global test:
# res_global = ANCOMBC_analysis$res_global
# 
# sig_OTUs_global <- res_global |>
#   filter(`diff_abn` == TRUE) |> #95
#   filter(`passed_ss` == TRUE)
# #11 sign DAA OTUs.

# #2: multiple pairwise comparisons
# res_pair = ANCOMBC_analysis$res_pair

#3: Dunnett's type of test
res_dunn <- ANCOMBC_analysis$res_dunn 

sig_OTUs_IBS_C_dunn <- res_dunn |> 
  filter(`diff_diagnoseIBS-C` == TRUE) |> #2
  filter(`passed_ss_diagnoseIBS-C` == TRUE) #0

sig_OTUs_IBS_D_dunn <- res_dunn |> 
  filter(`diff_diagnoseIBS-D` == TRUE) |> #9
  filter(`passed_ss_diagnoseIBS-D` == TRUE) #1
#Enterobacter_himalayensis

sig_OTUs_IBS_M_dunn <- res_dunn |> 
  filter(`diff_diagnoseIBS-M` == TRUE) |> #16
  filter(`passed_ss_diagnoseIBS-M` == TRUE) #2
#Anaerotignum_lactatifermentans
#Pseudomonas_sp002966775

# #4: trend test
# res_trend = ANCOMBC_analysis$res_trend #pattern analysis

#whit sensitivity score is important for keeping the power but keeping a low FDR. 
#filter the results based on the passed sensitivity score: 
res_OTU_sens_IBSC = res_dunn |> 
  filter(`passed_ss_diagnoseIBS-C` == TRUE)

res_OTU_sens_IBSM = res_dunn |>
  filter(`passed_ss_diagnoseIBS-M` == TRUE)

res_OTU_sens_IBSD = res_dunn |>
  filter(`passed_ss_diagnoseIBS-D` == TRUE)
  
#visualize: volcano plots
ancom_vplot_C <- ggplot(data=res_OTU_sens_IBSC, aes(x=`lfc_diagnoseIBS-C`, y=-log10(`p_diagnoseIBS-C`), col=`diff_diagnoseIBS-C`)) + 
  geom_point() + 
  labs(x=NULL, y=NULL) +
  labs(tag = "C", title = "IBS-C vs HC") + 
  theme(plot.title=element_text(hjust=0.5)) +
  scale_color_manual(name = "Differential \nabundance:", values = c("grey", "red")) +
  #add horizontal dashed line at p = 0.05:
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "black") +
  theme(legend.position = "none") 
ancom_vplot_C

ancom_vplot_M <- ggplot(data=res_OTU_sens_IBSM, aes(x=`lfc_diagnoseIBS-M`, y=-log10(`p_diagnoseIBS-M`), col=`diff_diagnoseIBS-M`)) + 
  geom_point() + 
  labs(x="Log change in abundance", y=NULL) +
  labs(tag = "B", title = "IBS-M vs HC") + 
  theme(plot.title=element_text(hjust=0.5)) +
  #add species names:
  geom_text_repel(data = subset(res_OTU_sens_IBSM, `p_diagnoseIBS-M` < 0.003), aes(label = taxon)) +
  #geom_text(data = subset(res_OTU_sens_IBSM, `p_diagnoseIBS-M` < 0.003), aes(label = taxon)) +
  scale_color_manual(name = "Differential \nabundance:", values = c("grey", "red")) +
  #add horizontal dashed line at p = 0.05:
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "black") +
  theme(legend.position = "none") 
ancom_vplot_M

ancom_vplot_D <- ggplot(data=res_OTU_sens_IBSD, aes(x=`lfc_diagnoseIBS-D`, y=-log10(`p_diagnoseIBS-D`), col=`diff_diagnoseIBS-D`)) + 
  geom_point() + 
  labs(x=NULL, y="-log10 p-value") +
  labs(tag = "A", title = "IBS-D vs HC") + 
  theme(plot.title=element_text(hjust=0.5)) +
  #add species names:
  geom_text_repel(data = subset(res_OTU_sens_IBSD, `p_diagnoseIBS-D` < 0.003), aes(label = taxon)) +
  #geom_text(data = subset(res_OTU_sens_IBSD, `p_diagnoseIBS-D` < 0.005), aes(label = taxon)) +
  scale_color_manual(name = "Differential \nabundance:", values = c("grey", "red")) +
  #add horizontal dashed line at p = 0.05:
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "black") +
  #remove legend: 
  theme(legend.position = "none")
ancom_vplot_D

#combine the three plots, so they have one legend and share y and x axis: 
ancom_vplot_D + ancom_vplot_M + ancom_vplot_C + plot_layout(nrow = 1)

ggsave("ancom_vplot_IBS.png", width = 18, height = 12, units = "cm", dpi = 300)

# #investigate the effects of the co-variates: age, sex, BMI:
# 
# #age:
# sig_OTUs_age <- res_prim |> 
#   filter(`diff_age` == TRUE) |> #3
#   filter(`passed_ss_age` == TRUE) #1 
# 
# #sex male: what about female? 
# sig_OTUs_sexMale <- res_prim |> 
#   filter(`diff_sexMale` == TRUE) |> #116
#   filter(`passed_ss_sexMale` == TRUE) #27 
# 
# #BMI: 
# sig_OTUs_BMI <- res_prim |>
#   filter(`diff_BMI` == TRUE) |> #0
#   filter(`passed_ss_BMI` == TRUE) #0

# #are there any structural zeros? 
# struct_zero <- ANCOMBC_analysis$zero_ind

# #bias corrected log ratio table with NAs as zero. 
# bias_correct_log_table <- ANCOMBC_analysis$bias_correct_log_table
# bias_correct_log_table[is.na(bias_correct_log_table)] = 0

#-------------------------------------------------------------------------------#

#high IBS-SSS (>300) vs HC: 

#abundance table: OTUs as rows and samples as columns. 
SSS_otu_abund <- IBS_SSS_data |> 
  filter(IBS_SSS_cat %in% c("HC", "high")) |> 
  as.data.frame() |> 
  select(-c(IBS_SSS:BMI)) |> 
  t() 

#sample info table: samples as rows. 
sample_info <- IBS_SSS_data |> 
  filter(IBS_SSS_cat %in% c("HC", "high")) |> 
  select(diagnose, age, sex, BMI, IBS_SSS, IBS_SSS_cat) |>
  as.data.frame()

#taxonomy table: optional to include? 
#two identical rows "OTUs" and "OTUs"
#needs to have rownames! 

OTU = otu_table(SSS_otu_abund, taxa_are_rows = TRUE)
samples = sample_data(sample_info)
#TAX = tax_table(otu_tax)

phyloseq_obj <- phyloseq(OTU, samples)

set.seed(19981103)
# Note that setting max_iter = 1 and B = 1 is only for the sake of speed
# Use default or larger values for max_iter and B for better performance

ANCOMBC_analysis_SSS = ancombc2(data = phyloseq_obj, 
                            assay_name = "counts",
                            tax_level = NULL, 
                            fix_formula = "IBS_SSS_cat + age + sex + BMI", 
                            rand_formula = NULL, #I don't have random effects. 
                            p_adj_method = "fdr", #fdr is good.  
                            pseudo_sens = TRUE,
                            prv_cut = 0.10, #10% prevalence cut-off. 
                            lib_cut = 2500,  #library size cut-off. The smallest library size is 2836.  
                            s0_perc = 0.05, #5% sparsity cut-off.
                            group = "IBS_SSS_cat", 
                            struc_zero = TRUE, #expected. 
                            neg_lb = FALSE, #I don't understand this...
                            alpha = 0.05, #significance level.
                            n_cl = 1, #number of clusters: 1 = no parallel computing.
                            verbose = TRUE, 
                            #which of these tests are relevant for me? 
                            global = TRUE, #global test: 
                            pairwise = TRUE, #multiple pairwise comparisons: 
                            dunnet = TRUE, #!!Dunnett's type of test: Multiple pairwise comparisons against a pre-specified group. 
                            #trend = TRUE, #trend test: not working.  
                            iter_control = list(tol = 1e-02, max_iter = 20, verbose = FALSE), #this is default 
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(), 
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                            #pattern analysis: not relevant for me.
                            trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)) #I don't understand this...

res_prim_SSS = ANCOMBC_analysis_SSS$res 

# res_dunn_SSS <- ANCOMBC_analysis_SSS$res_dunn 
# res_OTU_sens_SSS = res_dunn_SSS |>
#   filter(`passed_ss_IBS_SSS_cathigh` == TRUE)

sig_OTUs_SSS_high <- res_prim_SSS |>
  filter(`diff_IBS_SSS_cathigh` == TRUE) |> #53 
  filter(`passed_ss_IBS_SSS_cathigh` == TRUE) #6

# #the intercept is the IBS-SSS high category? 
# sig_OTUs_SSS_HC <- res_prim_SSS |>
#   filter(`diff_(Intercept)` == TRUE) |> #0
#   filter(`passed_ss_(Intercept)` == TRUE) #0

res_OTU_sens_SSS = res_prim_SSS |> 
  filter(`passed_ss_IBS_SSS_cathigh` == TRUE)

ancom_vplot_SSS <- ggplot(data=res_OTU_sens_SSS, aes(x=`lfc_IBS_SSS_cathigh`, y=-log10(`p_IBS_SSS_cathigh`), col=`diff_IBS_SSS_cathigh`)) + 
  geom_point() + 
  labs(x="Log fold change in abundance", y="-log10 p-value") +
  #add species names:
  #geom_text_repel(data = subset(res_OTU_sens_SSS, `p_IBS_SSS_cathigh` < 0.003), aes(label = taxon)) +
  geom_text(data = subset(res_OTU_sens_SSS, `p_IBS_SSS_cathigh` < 0.01), aes(label = taxon), show.legend = FALSE) +
  scale_color_manual(name = "Differential \nabundance:", values = c("grey", "red")) +
  #add horizontal dashed line at p = 0.05:
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "black") +
  theme(legend.position = c(0.91, 0.12)) +
  theme(legend.key.height = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black"))  
  #remove the whole legend
  #theme(legend.position = "none") 
ancom_vplot_SSS

#save plot: 
ggsave("ancom_vplot_SSS.png", plot = ancom_vplot_SSS, device = "png", width = 18, height = 12, units = "cm", dpi = 300)

#------------------------------------------------------------------------------#

#HADS 

#anxiety 

#abundance table: OTUs as rows and samples as columns. 
HADS_otu_abund <- HADS |> 
  as.data.frame() |> 
  column_to_rownames("p_ID") |> 
  select(-c(HADS_anx:BMI)) |>
  t() 

#sample info table: samples as rows. 
sample_info <- HADS |> 
  select(BG_ID, diagnose, age, sex, BMI, HADS_anx, HADS_dep, HADS_dep_cat, HADS_anx_cat) |>
  column_to_rownames("p_ID") |> 
  as.data.frame()

#taxonomy table: optional to include? 
#two identical rows "OTUs" and "OTUs"
#needs to have rownames! 

OTU = otu_table(HADS_otu_abund, taxa_are_rows = TRUE)
samples = sample_data(sample_info)
#TAX = tax_table(otu_tax)

phyloseq_obj <- phyloseq(OTU, samples)

set.seed(19981103)
# Note that setting max_iter = 1 and B = 1 is only for the sake of speed
# Use default or larger values for max_iter and B for better performance

ANCOMBC_analysis_HADS_anx = ancombc2(data = phyloseq_obj, 
                                assay_name = "counts",
                                tax_level = NULL, 
                                fix_formula = "HADS_anx_cat + diagnose + age + sex + BMI", 
                                rand_formula = NULL, #I don't have random effects. 
                                p_adj_method = "fdr", #fdr is good.  
                                pseudo_sens = TRUE,
                                prv_cut = 0.10, #10% prevalence cut-off. 
                                lib_cut = 2500,  #library size cut-off. The smallest library size is 2836.  
                                s0_perc = 0.05, #5% sparsity cut-off.
                                group = "HADS_anx_cat", 
                                struc_zero = TRUE, #expected. 
                                neg_lb = FALSE, #I don't understand this...
                                alpha = 0.05, #significance level.
                                n_cl = 1, #number of clusters: 1 = no parallel computing.
                                verbose = TRUE, 
                                #which of these tests are relevant for me? 
                                global = TRUE, #global test: 
                                pairwise = TRUE, #multiple pairwise comparisons: 
                                dunnet = TRUE, #!!Dunnett's type of test: Multiple pairwise comparisons against a pre-specified group. 
                                trend = TRUE, #trend test: not working.  
                                iter_control = list(tol = 1e-02, max_iter = 20, verbose = FALSE), #this is default 
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(), 
                                mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                #pattern analysis: not relevant for me.
                                trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)) #I don't understand this...

res_prim_HADS_anx = ANCOMBC_analysis_HADS_anx$res  

sig_OTUs_HADS_low <- res_prim_HADS_anx |>
  filter(`diff_HADS_anx_catnon-case` == TRUE) |> #0
  filter(`passed_ss_HADS_anx_catnon-case` == TRUE) #0

#the intercept is the IBS-SSS high category? 
sig_OTUs_HADS_high <- res_prim_HADS_anx |>
  filter(`diff_(Intercept)` == TRUE) |> #0
  filter(`passed_ss_(Intercept)` == TRUE) #0

res_OTU_sens_HADS_anx = res_prim_HADS_anx |>
  filter(`passed_ss_HADS_anx_catnon-case` == TRUE)

ancom_vplot_HADS_anx <- ggplot(data=res_OTU_sens_HADS_anx, aes(x=`lfc_HADS_anx_catnon-case`, y=-log10(`p_HADS_anx_catnon-case`), col=`diff_HADS_anx_catnon-case`)) + 
  geom_point() + 
  labs(x="Log fold change in abundance", y="-log10 p-value") +
  #add species names:
  #geom_text(data = subset(res_OTU_sens_SSS, `p_IBS_SSS_cathigh` < 0.01), aes(label = taxon)) +
  scale_color_manual(name = "Differential \nabundance:", values = c("grey", "red")) +
  #add horizontal line at p = 0.05:
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  theme(legend.position = c(0.92, 0.1)) +
  theme(legend.key.height = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) 
ancom_vplot_HADS_anx

#save the plot:
ggsave("ancom_vplot_HADS_anx.png", plot = ancom_vplot_HADS_anx, device = "png", width = 18, height = 12, units = "cm", dpi = 300)

#depression
ANCOMBC_analysis_HADS_dep = ancombc2(data = phyloseq_obj, 
                                     assay_name = "counts",
                                     tax_level = NULL, 
                                     fix_formula = "HADS_dep_cat + diagnose + age + sex + BMI", 
                                     rand_formula = NULL, #I don't have random effects. 
                                     p_adj_method = "fdr", #fdr is good.  
                                     pseudo_sens = TRUE,
                                     prv_cut = 0.10, #10% prevalence cut-off. 
                                     lib_cut = 2500,  #library size cut-off. The smallest library size is 2836.  
                                     s0_perc = 0.05, #5% sparsity cut-off.
                                     group = "HADS_dep_cat", 
                                     struc_zero = TRUE, #expected. 
                                     neg_lb = FALSE, #I don't understand this...
                                     alpha = 0.05, #significance level.
                                     n_cl = 1, #number of clusters: 1 = no parallel computing.
                                     verbose = TRUE, 
                                     #which of these tests are relevant for me? 
                                     global = TRUE, #global test: 
                                     pairwise = TRUE, #multiple pairwise comparisons: 
                                     dunnet = TRUE, #!!Dunnett's type of test: Multiple pairwise comparisons against a pre-specified group. 
                                     trend = TRUE, #trend test: not working.  
                                     iter_control = list(tol = 1e-02, max_iter = 20, verbose = FALSE), #this is default 
                                     em_control = list(tol = 1e-5, max_iter = 100),
                                     lme_control = lme4::lmerControl(), 
                                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                     #pattern analysis: not relevant for me.
                                     trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)) #I don't understand this...

res_prim_HADS_dep = ANCOMBC_analysis_HADS_dep$res  

sig_OTUs_HADS_dep_low <- res_prim_HADS_dep |>
  filter(`diff_HADS_dep_catnon-case` == TRUE) |> #54
  filter(`passed_ss_HADS_dep_catnon-case` == TRUE) #0

#the intercept is the IBS-SSS high category? 
sig_OTUs_HADS_dep_high <- res_prim_HADS_dep |>
  filter(`diff_(Intercept)` == TRUE) |> #0
  filter(`passed_ss_(Intercept)` == TRUE) #0

res_OTU_sens_HADS_dep = res_prim_HADS_dep |>
  filter(`passed_ss_HADS_dep_catnon-case` == TRUE)

ancom_vplot_HADS_dep <- ggplot(data=res_OTU_sens_HADS_dep, aes(x=`lfc_HADS_dep_catnon-case`, y=-log10(`p_HADS_dep_catnon-case`), col=`diff_HADS_dep_catnon-case`)) + 
  geom_point() + 
  labs(x="Log fold change in abundance", y="-log10 p-value") +
  #add species names:
  #geom_text(data = subset(res_OTU_sens_SSS, `p_IBS_SSS_cathigh` < 0.01), aes(label = taxon)) +
  scale_color_manual(name = "Differential \nabundance:", values = c("grey", "red")) + 
  #add horizontal line at p = 0.05:
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  theme(legend.position = c(0.92, 0.1)) +
  theme(legend.key.height = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) 
ancom_vplot_HADS_dep

#save plot: 
ggsave("ancom_vplot_HADS_dep.png", plot = ancom_vplot_HADS_dep, device = "png", width = 18, height = 12, units = "cm", dpi = 300)

#CFQ-11/FSS

#abundance table: OTUs as rows and samples as columns. 
FSS_otu_abund <- FSS |> 
  as.data.frame() |> 
  column_to_rownames("p_ID") |> 
  select(-c(FSS_score_BL:BMI)) |>
  t() 

#sample info table: samples as rows. 
sample_info <- FSS |> 
  select(BG_ID, diagnose, age, sex, BMI, FSS_score_BL, FSS_cat) |>
  column_to_rownames("p_ID") |> 
  as.data.frame()

#taxonomy table: optional to include? 
#two identical rows "OTUs" and "OTUs"
#needs to have rownames! 

OTU = otu_table(FSS_otu_abund, taxa_are_rows = TRUE)
samples = sample_data(sample_info)
#TAX = tax_table(otu_tax)

phyloseq_obj <- phyloseq(OTU, samples)

set.seed(19981103)
# Note that setting max_iter = 1 and B = 1 is only for the sake of speed
# Use default or larger values for max_iter and B for better performance

ANCOMBC_analysis_FSS = ancombc2(data = phyloseq_obj, 
                                assay_name = "counts",
                                tax_level = NULL, 
                                fix_formula = "FSS_cat + diagnose + age + sex + BMI", 
                                rand_formula = NULL, #I don't have random effects. 
                                p_adj_method = "fdr", #fdr is good.  
                                pseudo_sens = TRUE,
                                prv_cut = 0.10, #10% prevalence cut-off. 
                                lib_cut = 2500,  #library size cut-off. The smallest library size is 2836.  
                                s0_perc = 0.05, #5% sparsity cut-off.
                                group = "FSS_cat", 
                                struc_zero = TRUE, #expected. 
                                neg_lb = FALSE, #I don't understand this...
                                alpha = 0.05, #significance level.
                                n_cl = 1, #number of clusters: 1 = no parallel computing.
                                verbose = TRUE, 
                                #which of these tests are relevant for me? 
                                global = TRUE, #global test: 
                                pairwise = TRUE, #multiple pairwise comparisons: 
                                dunnet = TRUE, #!!Dunnett's type of test: Multiple pairwise comparisons against a pre-specified group. 
                                trend = TRUE, #trend test: not working.  
                                iter_control = list(tol = 1e-02, max_iter = 20, verbose = FALSE), #this is default 
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(), 
                                mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                #pattern analysis: not relevant for me.
                                trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)) #I don't understand this...

res_prim_FSS = ANCOMBC_analysis_FSS$res  

sig_OTUs_FSS_low <- res_prim_FSS |>
  filter(`diff_FSS_catnon_case` == TRUE) |> #0
  filter(`passed_ss_FSS_catnon_case` == TRUE) #0

#the intercept is the IBS-SSS high category? 
sig_OTUs_FSS_high <- res_prim_FSS |>
  filter(`diff_(Intercept)` == TRUE) |> #0
  filter(`passed_ss_(Intercept)` == TRUE) #0

res_OTU_sens_FSS = res_prim_FSS |>
  filter(`passed_ss_FSS_catnon_case` == TRUE)

ancom_vplot_FSS <- ggplot(data=res_OTU_sens_FSS, aes(x=`lfc_FSS_catnon_case`, y=-log10(`p_FSS_catnon_case`), col=`diff_FSS_catnon_case`)) + 
  geom_point() + 
  labs(x="Log fold change in abundance", y="-log10 p-value") +
  #add species names:
  #geom_text(data = subset(res_OTU_sens_SSS, `p_IBS_SSS_cathigh` < 0.01), aes(label = taxon)) +
  scale_color_manual(name = "Differential \nabundance", values = c("grey", "red")) +
  #add horizontal line at p = 0.05:
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  theme(legend.position = c(0.92, 0.1)) + 
  theme(legend.key.height = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", linewidth = 0.2, linetype = "solid", color = "black")) 
ancom_vplot_FSS

#save plot:
ggsave("ancom_vplot_FSS.png", plot = ancom_vplot_FSS, device = "png", width = 18, height = 12, units = "cm", dpi = 300)
