# group-wise effects gfap on glu:
em_results_glu <- lmer_GFAP_glu_groupwise %>% emtrends(~ diagnose_group, var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% confint()
em_results_p_glu <-lmer_GFAP_glu_groupwise %>% emtrends(~ diagnose_group, var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% test()
coef_GFAP_glu_group_1 <- paste0(round((em_results_glu$GFAP.trend[1:4])*1000,2), " [", round((em_results_glu$lower.CL[1:4])*1000,2), "; ", sprintf("%.2f", (em_results_glu$upper.CL[1:4])*1000), "]) \u00d7 10\u00b3")
coef_GFAP_glu_group_2 <- paste0("p = ", round(em_results_p_glu$p.value[1:4], 3))

# group-wise effects gfap on GABA:
em_results_gaba <- lmer_GFAP_gaba_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% confint()
em_results_p_gaba <-lmer_GFAP_gaba_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% test()
coef_GFAP_gaba_group_1 <- paste0(round((em_results_gaba$GFAP.trend[1:4])*1000,2), " [", round((em_results_gaba$lower.CL[1:4])*1000,2), "; ", round((em_results_gaba$upper.CL[1:4])*1000,2), "]) \u00d7 10\u00b3")
coef_GFAP_gaba_group_2 <- paste0("p = ",sprintf("%.3f",(em_results_p_gaba$p.value[1:4])))

# group-wise effects gfap on glu/GABA:
em_results_glu_gaba <- lmer_GFAP_glu_gaba_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% confint()
em_results_p_glu_gaba <-lmer_GFAP_glu_gaba_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% test()
coef_GFAP_glu_gaba_group_1 <- paste0(sprintf("%.2f",(em_results_glu_gaba$GFAP.trend[1:4])*1000), " [", round((em_results_glu_gaba$lower.CL[1:4])*1000,2), "; ", round((em_results_glu_gaba$upper.CL[1:4])*1000,2), "]) \u00d7 10\u00b3")
coef_GFAP_glu_gaba_group_2 <- paste0("p = ",sprintf("%.3f",em_results_p_glu_gaba$p.value[1:4], 3))

# group-wise effects gfap on gln:
em_results_gln <- lmer_GFAP_gln_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% confint()
em_results_p_gln <-lmer_GFAP_gln_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% test()
coef_GFAP_gln_group_1 <- paste0(round((em_results_gln$GFAP.trend[1:4])*1000,2), " [", round((em_results_gln$lower.CL[1:4])*1000,2), "; ", round((em_results_gln$upper.CL[1:4])*1000,2), "]) \u00d7 10\u00b3")
coef_GFAP_gln_group_2 <- paste0("p = ",round(em_results_p_gln$p.value[1:4], 3))

# group-wise effects gfap on glu/gln:
em_results_glu_gln <- lmer_GFAP_glu_gln_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% confint()
em_results_p_glu_gln <-lmer_GFAP_glu_gln_groupwise %>% emtrends(~ diagnose_group,var = "GFAP", data = NeuroMET %>% filter(!is.na(GFAP))) %>% test()
coef_GFAP_glu_gln_group_1 <- paste0(sprintf("%.2f",((em_results_glu_gln$GFAP.trend[1:4])*1000)), " [", sprintf("%.2f",(em_results_glu_gln$lower.CL[1:4])*1000,2), "; ", sprintf("%.2f",(em_results_glu_gln$upper.CL[1:4])*1000,2), "]) \u00d7 10\u00b3")
coef_GFAP_glu_gln_group_2 <- paste0("p = ",sprintf("%.3f",(em_results_p_glu_gln$p.value[1:4])))

# Model Glutamate
# Define new data grid
new_data_gln_year <- expand.grid(
  year = seq(min(NeuroMETs$year, na.rm = TRUE), max(NeuroMETs$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMETs$diagnose_group)),
  age_bl = mean(NeuroMETs$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1)   # or change to "male" or both if needed
)

# Predict values
preds <- predict(lmer_long_gln_groupwise_s, newdata = new_data_gln_year, re.form = NA, se.fit = TRUE)

label_positions <- new_data_gln_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 7, .groups = "drop")
coef_gln <- data.frame(emtrends(lmer_long_gln_groupwise_s, ~ diagnose_group, var = "year")) %>%
  mutate(pval = (emtrends(lmer_long_gln_groupwise_s, ~ diagnose_group, var = "year") %>% test())$p.value,
         pval_txt = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3))),
         label = paste0("β = ", format(year.trend, nsmall = 1, digits = 1), 
                        " [", round(lower.CL, 2), "; ", round(upper.CL, 2), "]\n", 
                        pval_txt)) %>%
  select(diagnose_group, label) %>%
  left_join(label_positions, by = "diagnose_group")

# Store predictions and confidence intervals
new_data_gln_year$pred <- preds$fit
new_data_gln_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_gln_year$upper <- preds$fit + 1.96 * preds$se.fit

pA1 <- ggplot(new_data_gln_year, aes(x = year, y = pred, color = diagnose_group)) +
  
  geom_line(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(gln)), 
            aes(x = year, y = gln, group=record_id, color=diagnose_group), alpha=0.5) +
  # geom_point(data=NeuroMET %>% filter(!is.na(diagnose_group) & !is.na(gln)), 
  #            aes(x = year, y = gln, group=record_id, color=diagnose_group, size = gln_CRLB), alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(pTau_group)), aes(x = year, y = gln, color = pTau_group, size = gln_CRLB),  alpha = 0.8)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  theme_minimal() +
  labs(x = "Year", y = "Glutamine [mmol/L]",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis at visit 1") +
  facet_wrap(~diagnose_group, nrow = 1) + 
  geom_text(data = coef_gln, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
  scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  scale_size_continuous(range = c(0.8, 2.5), name = "CRLB-based size \nfor Glutamine") +  # optional legend
  guides(fill= "none", color = "none")+
  theme(legend.position = 'top',
        axis.title.x=element_blank(),
        strip.background = element_rect(fill = "grey85",colour = NA),
        axis.line = element_line(colour = "black"),
        axis.ticks.x= element_line(colour = "black"))

# Model Glutamate/Glutamine
# Define new data grid
new_data_glu_gln_year <- expand.grid(
  year = seq(min(NeuroMETs$year, na.rm = TRUE), max(NeuroMETs$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMETs$diagnose_group)),
  age_bl = mean(NeuroMETs$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1)   # or change to "male" or both if needed
)

# Predict values
preds <- predict(lmer_long_glu_gln_groupwise_s, newdata = new_data_glu_gln_year, re.form = NA, se.fit = TRUE)

label_positions <- new_data_glu_gln_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 1.7, .groups = "drop")
coef_glu_gln <- data.frame(emtrends(lmer_long_glu_gln_groupwise_s, ~ diagnose_group, var = "year")) %>%
  mutate(pval = (emtrends(lmer_long_glu_gln_groupwise_s, ~ diagnose_group, var = "year") %>% test())$p.value,
         pval_txt = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3))),
         label = paste0("β = ", format(year.trend, nsmall = 1, digits = 1), 
                        " [", round(lower.CL, 2), "; ", round(upper.CL, 2), "]\n", 
                        pval_txt)) %>%
  select(diagnose_group, label) %>%
  left_join(label_positions, by = "diagnose_group")

# Store predictions and confidence intervals
new_data_glu_gln_year$pred <- preds$fit
new_data_glu_gln_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_glu_gln_year$upper <- preds$fit + 1.96 * preds$se.fit

pA2 <- ggplot(new_data_glu_gln_year, aes(x = year, y = pred, color = diagnose_group)) +
  
  geom_line(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(glu_gln)), 
            aes(x = year, y = glu_gln, group=record_id, color=diagnose_group), alpha=0.5) +
  # geom_point(data=NeuroMET %>% filter(!is.na(diagnose_group) & !is.na(glu_gln)), 
  #            aes(x = year, y = glu_gln, group=record_id, color=diagnose_group, size = glu_gln_CRLB), alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(pTau_group)), aes(x = year, y = glu_gln, color = pTau_group),  alpha = 0.8)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  theme_minimal() +
  labs(x = "Year", y = "Log(Glutamate/Glutamine)",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis at visit 1") +
  facet_wrap(~diagnose_group, nrow = 1) + 
  geom_text(data = coef_glu_gln, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
  scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  guides(fill= "none", color="none", size = "none")+
  theme(legend.position = 'top',
        strip.background = element_rect(fill = "grey85",colour = NA),
        axis.line = element_line(colour = "black"),
        axis.ticks.x= element_line(colour = "black"))

### group-wise coefficients
#Glutamate
em_results_glu$diagnose_group <- fct_relevel(
    em_results_glu$diagnose_group, 
    "AD Aβ+", "MCI Aβ+", "CN Aβ+", "CN Aβ-"
)
p_glu = data.frame (coef_GFAP_glu_group_2,  y = as.factor(c("CN Aβ-","CN Aβ+", "MCI Aβ+","AD Aβ+")), x = rep(0.005, 4))
errbar_ggplot_glu <- ggplot(em_results_glu, aes(diagnose_group, (GFAP.trend*1000), ymin = (lower.CL*1000), ymax= (upper.CL*1000), color = diagnose_group ))+
  geom_point(size = 1.8)+
  geom_linerange(size = 1) +
  coord_flip() +
  labs(x = "", y = (expression( "β"["Glutamate~GFAP"]~"\u{00d7}"~"10"~"\u{00b3}")), color = "Clinical Diagnosis at Visit 1")+
  geom_hline(yintercept = 0)+
  scale_color_manual(values= c("CN Aβ-" =diagnose_colors[1], "CN Aβ+" = diagnose_colors[2], "MCI Aβ+" = diagnose_colors[3], "AD Aβ+" = diagnose_colors[4]), guide = guide_legend(reverse = TRUE)) +
  theme_minimal()+
  theme(legend.position = "none")
#add p-values
F <- errbar_ggplot_glu+ geom_text(data = p_glu, aes(x = em_results_glu$diagnose_group, y = em_results_glu$GFAP.trend, label = coef_GFAP_glu_group_2), inherit.aes = FALSE,  vjust = -1.5, hjust=1.2, size = 3.5)

#GABA
em_results_gaba$diagnose_group <- fct_relevel(
    em_results_gaba$diagnose_group, 
    "AD Aβ+", "MCI Aβ+", "CN Aβ+", "CN Aβ-"
)
p_gaba = data.frame (coef_GFAP_gaba_group_2,  y = as.factor(c("CN Aβ-","CN Aβ+", "MCI Aβ+","AD Aβ+")), x = rep(0.005, 4))
errbar_ggplot_gaba <- ggplot(em_results_gaba, aes(diagnose_group, (GFAP.trend*1000), ymin = (lower.CL*1000), ymax= (upper.CL*1000), color = diagnose_group ))+
  geom_point(size = 1.8)+
  geom_linerange(size = 1) +
  coord_flip() +
  labs(x = "", y = expression("β"["GABA~GFAP"]~"\u{00d7}"~"10"~"\u{00b3}"), color = "Clinical Diagnosis at Visit 1")+
  geom_hline(yintercept = 0)+
  scale_color_manual(values= c("CN Aβ-" =diagnose_colors[1], "CN Aβ+" = diagnose_colors[2], "MCI Aβ+" = diagnose_colors[3], "AD Aβ+" = diagnose_colors[4]), guide = guide_legend(reverse = TRUE)) +
  theme_minimal()
#add p-values
G <- errbar_ggplot_gaba+ geom_text(data = p_gaba, aes(x = em_results_gaba$diagnose_group, y = em_results_gaba$GFAP.trend, label = coef_GFAP_gaba_group_2), inherit.aes = FALSE,  vjust = -1.5, hjust=1.2, size = 3.5)

#Glutamin
em_results_gln$diagnose_group <- fct_relevel(
    em_results_gln$diagnose_group, 
    "AD Aβ+", "MCI Aβ+", "CN Aβ+", "CN Aβ-"
)
p_gln = data.frame (coef_GFAP_gln_group_2,  y = as.factor(c("CN Aβ-","CN Aβ+", "MCI Aβ+","AD Aβ+")), x = rep(0.005, 4))
errbar_ggplot_gln <- ggplot(em_results_gln, aes(diagnose_group, (GFAP.trend*1000), ymin = (lower.CL*1000), ymax= (upper.CL*1000), color = diagnose_group ))+
  geom_point(size = 1.8)+
  geom_linerange(size = 1) +
  coord_flip() +
  labs(x = "", y = expression( "β"["Glutamine~GFAP"]~"\u{00d7}"~"10"~"\u{00b3}"), color = "Clinical Diagnosis at Visit 1")+
  geom_hline(yintercept = 0)+
  scale_color_manual(values= c("CN Aβ-" =diagnose_colors[1], "CN Aβ+" = diagnose_colors[2], "MCI Aβ+" = diagnose_colors[3], "AD Aβ+" = diagnose_colors[4]), guide = guide_legend(reverse = TRUE)) +
  theme_minimal()+
  theme(legend.position = "none")
#add p-values
H <- errbar_ggplot_gln + geom_text(data = p_gln, aes(x = em_results_gln$diagnose_group, y = em_results_gln$GFAP.trend, label = coef_GFAP_gln_group_2), inherit.aes = FALSE,  vjust = -1.5, hjust=1.2, size = 3.5)

#Glutamate/Glutamine
em_results_glu_gln$diagnose_group <- fct_relevel(
    em_results_glu_gln$diagnose_group, 
    "AD Aβ+", "MCI Aβ+", "CN Aβ+", "CN Aβ-"
)
p_glu_gln = data.frame (coef_GFAP_glu_gln_group_2,  y = as.factor(c("CN Aβ-","CN Aβ+", "MCI Aβ+","AD Aβ+")), x = rep(0.005, 4))
errbar_ggplot_glu_gln <- ggplot(em_results_glu_gln, aes(diagnose_group, (GFAP.trend*1000), ymin = (lower.CL*1000), ymax= (upper.CL*1000), color = diagnose_group ))+
  geom_point(size = 1.8)+
  geom_linerange(size = 1) +
  coord_flip() +
  labs(x = "", y = expression( "β"["log(Glutamate/Glutamine)~GFAP"]~"\u{00d7}"~"10"~"\u{00b3}"), color = "Clinical Diagnosis at Visit 1")+
  geom_hline(yintercept = 0)+
  scale_color_manual(values= c("CN Aβ-" =diagnose_colors[1], "CN Aβ+" = diagnose_colors[2], "MCI Aβ+" = diagnose_colors[3], "AD Aβ+" = diagnose_colors[4]), guide = guide_legend(reverse = TRUE)) +
  theme_minimal()+
  theme(legend.position = "none")
#add p-values
I <- errbar_ggplot_glu_gln+ geom_text(data = p_glu_gln, aes(x = em_results_glu_gln$diagnose_group, y = em_results_glu_gln$GFAP.trend, label = coef_GFAP_glu_gln_group_2), inherit.aes = FALSE,  vjust = -1.5, hjust=1.0, size = 3.5)

Plot4A1 <- pA1/pA2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")
Plot4A <- wrap_elements(Plot4A1 + plot_annotation(title = "A Group-wise Longitudinal Trajectory of 7T MRS Glutamine and Log(Glutamate/Glutamine) ratio"))

Plot4C <- (F+G+H+I)+
  plot_layout(nrow = 1, guides = 'collect') &
  theme(legend.position = "top")
Plot4C <- wrap_elements(Plot4C + plot_annotation(title = "B Group-wise Associations between Plasma GFAP and Neurotransmitters adjusted for PCC+Precuneus Volume"))

Figure4<- Plot4A / Plot4C + plot_layout(heights = c(1.8, 1)) 
Figure4
ggsave(
  filename = "Figure4.jpeg",
  plot = Figure4,
  width = 10,  
  height = 10,
  dpi = 300
)
