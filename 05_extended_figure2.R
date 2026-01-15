# Model pTau
# Define new data grid
new_data_pTau_year <- expand.grid(
  year = seq(min(NeuroMETs$year, na.rm = TRUE), max(NeuroMETs$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMETs$diagnose_group)),
  age_bl = mean(NeuroMETs$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1) 
)

# Predict values
preds <- predict(lmer_long_pTau_groupwise_s, newdata = new_data_pTau_year, re.form = NA, se.fit = TRUE)

# Store predictions and confidence intervals
new_data_pTau_year$pred <- preds$fit
new_data_pTau_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_pTau_year$upper <- preds$fit + 1.96 * preds$se.fit

label_positions <- new_data_pTau_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 9.5, .groups = "drop")
coef_pTau <- data.frame(emtrends(lmer_long_pTau_groupwise_s, ~ diagnose_group, var = "year")) %>%
  mutate(pval = (emtrends(lmer_long_pTau_groupwise_s, ~ diagnose_group, var = "year") %>% test())$p.value,
         pval_txt = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3))),
         label = paste0("β = ", format(year.trend, nsmall = 1, digits = 1), " [", format(lower.CL, nsmall = 1, digits = 1), "; ", format(upper.CL, nsmall = 2, digits = 1), "]\n", pval_txt)) %>%
  select(diagnose_group, label) %>%
  left_join(label_positions, by = "diagnose_group")

p4 <- ggplot(new_data_pTau_year, aes(x = year, y = pred, color = diagnose_group)) +
  
  geom_line(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(pTau)), 
            aes(x = year, y = pTau, group=record_id, color=diagnose_group), alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(pTau)), 
             aes(x = year, y = pTau, group=record_id, color=diagnose_group), size=1, alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(pTau_group)), aes(x = year, y = pTau, color = pTau_group), size = 1.6, alpha = 0.9)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  geom_hline(yintercept = 2.08, linetype = "longdash", color = diagnose_colors[4], linewidth = 1) +
  theme_minimal() +
  labs(x = "Year", y = "Plasma p-Tau 181 [pg/mL]",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis at visit 1") +
  facet_wrap(~diagnose_group, nrow = 1) +  
  geom_text(data = coef_pTau, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
  scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  guides(fill= "none") +
  theme(legend.position = 'top', 
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x= element_line(colour = "black"),
        axis.line = element_line(colour = "black"))


# Model NMM
# Define new data grid
new_data_nmm_year <- expand.grid(
  year = seq(min(NeuroMETs$year, na.rm = TRUE), max(NeuroMETs$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMETs$diagnose_group)),
  age_bl = mean(NeuroMETs$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1) 
)

# Predict values
preds <- predict(lmer_long_nmm_groupwise_s, newdata = new_data_nmm_year, re.form = NA, se.fit = TRUE)

# Store predictions and confidence intervals
new_data_nmm_year$pred <- preds$fit
new_data_nmm_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_nmm_year$upper <- preds$fit + 1.96 * preds$se.fit

label_positions <- new_data_nmm_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 3.5, .groups = "drop")
coef_nmm <- data.frame(emtrends(lmer_long_nmm_groupwise_s, ~ diagnose_group, var = "year")) %>%
  mutate(pval = (emtrends(lmer_long_nmm_groupwise_s, ~ diagnose_group, var = "year") %>% test())$p.value,
         pval_txt = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3))),
         label = paste0("β = ", format(year.trend, nsmall = 1, digits = 1), " [", format(lower.CL, nsmall = 1, digits = 1), "; ", format(upper.CL, nsmall = 1, digits = 1), "]\n", pval_txt)) %>%
  select(diagnose_group, label) %>%
  left_join(label_positions, by = "diagnose_group")

p5 <- ggplot(new_data_nmm_year, aes(x = year, y = pred, color = diagnose_group)) +
  
  geom_line(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(nmm)), 
            aes(x = year, y = nmm, group=record_id, color=diagnose_group), alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(nmm)), 
             aes(x = year, y = nmm, group=record_id, color=diagnose_group), size=1, alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(pTau_group)), aes(x = year, y = nmm, color = pTau_group), size = 1.6, alpha = 0.9)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  theme_minimal() +
  labs(x = "Year", y = "Memory Ability (NMM)",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis atvisit 1") +
  facet_wrap(~diagnose_group, nrow = 1) +    
  geom_text(data = coef_nmm, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
  scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  guides(fill= "none") +
  theme(legend.position = 'none', 
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        strip.text.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.x= element_line(colour = "black"))

# Define new data grid
new_data_gfap_year <- expand.grid(
  year = seq(min(NeuroMET$year, na.rm = TRUE), max(NeuroMET$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMET$diagnose_group)),
  age_bl = mean(NeuroMET$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1)   # or change to "male" or both if needed
)

# Predict values
preds <- predict(lmer_long_gfap_groupwise, newdata = new_data_gfap_year, re.form = NA, se.fit = TRUE)

label_positions <- new_data_gfap_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 525, .groups = "drop")
coef_gfap <- data.frame(emtrends(lmer_long_gfap_groupwise, ~ diagnose_group, var = "year")) %>%
  mutate(pval = (emtrends(lmer_long_gfap_groupwise, ~ diagnose_group, var = "year") %>% test())$p.value,
         pval_txt = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3))),
         label = paste0("β = ", format(year.trend, nsmall = 2, digits = 2), 
                        " [", round(lower.CL, 2), "; ", sprintf("%.2f", upper.CL), "]\n", 
                        pval_txt)) %>%
  select(diagnose_group, label) %>%
  left_join(label_positions, by = "diagnose_group")

# Store predictions and confidence intervals
new_data_gfap_year$pred <- preds$fit
new_data_gfap_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_gfap_year$upper <- preds$fit + 1.96 * preds$se.fit

A <- ggplot(new_data_gfap_year, aes(x = year, y = pred, color = diagnose_group)) +
  geom_line(data=NeuroMET %>% filter(!is.na(diagnose_group) & !is.na(GFAP)), 
            aes(x = year, y = GFAP, group=record_id, color=diagnose_group), alpha=0.5) +
  # geom_point(data=NeuroMET %>% filter(!is.na(diagnose_group) & !is.na(GFAP)), 
  #            aes(x = year, y = GFAP, group=record_id, color=diagnose_group), alpha=0.5) +
  geom_point(data=NeuroMET %>% filter(!is.na(pTau_group)), aes(x = year, y = GFAP, color = pTau_group),  alpha = 0.8)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  theme_minimal() +
  labs(x = "Year", y = "Plasma GFAP [pg/mL]",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis at visit 1") +
  facet_wrap(~diagnose_group, nrow = 1) + 
  geom_text(data = coef_gfap, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
  scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  guides(fill= "none", color = guide_legend(nrow=2))+
  theme(legend.position = 'none',
        strip.background = element_rect(fill = "grey85",colour = NA),
        axis.line = element_line(colour = "black"),
        axis.ticks.x= element_line(colour = "black"),
        strip.text.x = element_blank())

figureE1 <- p4/p5/A

ggsave(
  filename = "extended.figure1.jpeg",
  plot = figureE1,
  width = 10,  
  height = 8,
  dpi = 300
)
