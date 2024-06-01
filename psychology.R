#Psychological data: HADS, fatigue. 

# #HADS by group:
# 
# #anxiety: 
# HADS_anx_g_bplot <- HADS |> 
#   ggplot(aes(x = group, y = HADS_anx, fill = group)) +
#   geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) +
#   geom_jitter(aes(color = group), show.legend = FALSE, width = 0.3) +
#   geom_hline(aes(yintercept = 11), linetype="dashed", linewidth = 0.5) +
#   geom_hline(aes(yintercept = 8), linetype="dashed", linewidth = 0.5) +
#   scale_fill_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   scale_color_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   labs(x = "", y = "HADS anxiety score") +
#   ylim(-1, 20) + 
#   theme_classic()
# HADS_anx_g_bplot 
# 
# #depression: 
# HADS_dep_g_bplot <- HADS |> 
#   ggplot(aes(x = group, y = HADS_dep, fill = group)) +
#   geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) +
#   geom_jitter(aes(color = group), width = 0.3, show.legend = FALSE) +
#   geom_hline(aes(yintercept = 11), linetype="dashed", linewidth = 0.5) +
#   geom_hline(aes(yintercept = 8), linetype="dashed", linewidth = 0.5) +
#   scale_fill_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   scale_color_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   labs(x = "", y = "HADS depression score") +
#   ylim(-1, 20) +
#   theme_classic()
# HADS_dep_g_bplot
# 
# combine_HADS_g_bp <- HADS_dep_g_bplot | HADS_anx_g_bplot
# combine_HADS_g_bp

#HADSby diagnose:

#anxiety: 
HADS_anx_d_bplot <- HADS |> 
  ggplot(aes(x = diagnose, y = HADS_anx, fill = diagnose)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, show.legend = FALSE) +
  geom_hline(aes(yintercept = 11), linetype="dashed", linewidth = 0.7) +
  geom_hline(aes(yintercept = 8), linetype="dashed", linewidth = 0.7) +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x = "", y = "HADS anxiety score") +
  ylim(-1, 21) 
HADS_anx_d_bplot

#significance test:
#kruksal-wallis test:
kruskal.test(HADS_anx ~ diagnose, data = HADS)
#pairwise wilcoxon test:
pairwise.wilcox.test(HADS$HADS_anx, HADS$diagnose, p.adjust.method = "bonferroni")

#depression:
HADS_dep_d_bplot <- HADS |> 
  ggplot(aes(x = diagnose, y = HADS_dep, fill = diagnose)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(aes(color = diagnose), width = 0.3, show.legend = FALSE) +
  geom_hline(aes(yintercept = 11), linetype="dashed", linewidth = 0.7) +
  geom_hline(aes(yintercept = 8), linetype="dashed", linewidth = 0.7) +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x = "", y = "HADS depression score") +
  ylim(-1, 21)
HADS_dep_d_bplot

#significance test:
#kruksal-wallis test:
kruskal.test(HADS_dep ~ diagnose, data = HADS)
#pairwise wilcoxon test:
pairwise.wilcox.test(HADS$HADS_dep, HADS$diagnose, p.adjust.method = "bonferroni")

combine_HADS_d_bp <- HADS_dep_d_bplot | HADS_anx_d_bplot
combine_HADS_d_bp

ggsave("combine_HADS_d_bp.png", plot = combine_HADS_d_bp, width = 18, height = 12, units = "cm", dpi = 300)

# #HADS total by group:
# HADS_tot_plot_g <- HADS |> 
#   ggplot(aes(x = group, y = HADS_total, fill = group)) +
#   geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) +
#   geom_jitter(aes(color = group), width = 0.2, show.legend = FALSE) +
#   geom_hline(aes(yintercept = 11), linetype="dashed", linewidth = 0.5) +
#   geom_hline(aes(yintercept = 8), linetype="dashed", linewidth = 0.5) +
#   scale_fill_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   scale_color_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   labs(x = "", y = "HADS total score") +
#   ylim(-1, 30) +
#   theme_classic()
# HADS_tot_plot_g
# 
# #HADS total by diagnose:
# HADS_tot_plot_d <- HADS |> 
#   ggplot(aes(x = diagnose, y = HADS_total, fill = diagnose)) +
#   geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape =NA) +
#   geom_jitter(aes(color = diagnose), width = 0.2, show.legend = FALSE) +
#   geom_hline(aes(yintercept = 11), linetype="dashed", linewidth = 0.5) +
#   geom_hline(aes(yintercept = 8), linetype="dashed", linewidth = 0.5) +
#   scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
#   labs(x = "", y = "HADS total score") +
#   ylim(-1, 30) +
#   theme_classic()
# HADS_tot_plot_d

# #CFQ-11/FSS (fatigue) by group:
# FSS_plot_g <- FSS |> 
#   ggplot(aes(x = group, y = FSS_score_BL, fill = group)) +
#   geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) + 
#   geom_jitter(aes(color = group), width = 0.3, show.legend = FALSE) +
#   #geom_hline(aes(yintercept = 11), linetype="dashed", linewidth = 0.7) +
#   scale_fill_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   scale_color_manual(values = c("IBS" = "#0080FF", "HC" = "#FF9900")) +
#   labs(x = "", y = "BIS score") +
#   ylim(-1, 12) 
# FSS_plot_g

#CFQ-11/FSS (fatigue)  by diagnose:
FSS_plot_d <- FSS |> 
  ggplot(aes(x = diagnose, y = FSS_score_BL, fill = diagnose)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) + 
  geom_jitter(aes(color = diagnose), width = 0.3, alpha = 0.7, show.legend = FALSE) +
  geom_hline(aes(yintercept = 4), linetype="dashed", linewidth = 0.7) +
  scale_fill_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  scale_color_manual(values = c("HC" = "#FF9900", "IBS-C" = "#39cce3", "IBS-M" = "#3974e3", "IBS-D" = "#131863")) +
  labs(x = "", y = "Chalder fatigue scale (CFQ-11) score") +
  scale_y_continuous(limits = c(-1, 12))
FSS_plot_d

ggsave("FSS_plot_d.png", plot = FSS_plot_d, width = 18, height = 12, units = "cm", dpi = 300)
