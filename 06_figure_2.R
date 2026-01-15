Boxplot2 <- NeuroMET %>%
  filter(visit == "t1" )%>%
  select(record_id, diagnose, glu, gaba, glu_gaba, pTau_v1_group, diagnose_group, glu_CRLB, gaba_CRLB)%>%
  pivot_longer(cols = 3:5, names_to = "biomarker", values_drop_na = TRUE)%>%
  mutate(biomarker = factor(biomarker, levels = c("glu", "gaba", "glu_gaba"))) %>% # Set order: Glutamate first
  group_by(diagnose_group)%>% 
  mutate(CRLB = ifelse(biomarker == "glu", glu_CRLB , ifelse(biomarker == "gaba", gaba_CRLB,  NA)))%>%
  ggplot()+ 
  aes(x= diagnose_group, y= value, fill = diagnose_group)+
  geom_boxplot(outlier.shape = NA, lwd = 0.3, alpha = 0.9)+
  geom_jitter(aes(size = ifelse(biomarker == "glu", glu_CRLB, ifelse(biomarker == "gaba", gaba_CRLB, 0.8))), 
              position=position_jitter(0.2), alpha = 0.25)+
  xlab(NULL)+
  ylab(NULL) +
  labs(fill = "Diagnosis at Visit 1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ ####remove background grid
  scale_fill_manual(values=met.brewer ("Hiroshige", direction = -1, n = 4))+
  #scale_color_manual(values=met.brewer ("Hiroshige", direction = -1, n = 4))+
  scale_size_continuous(range = c(0.8, 2.5), name = "CRLB-based size for Glutamate and GABA") +  # optional legend
  facet_wrap(~biomarker, scales = "free", labeller = labeller(biomarker = c(glu = "Glutamate [mmol/L]", gaba = "GABA [mmol/L]", glu_gaba = "log(Glutamate/GABA)")))+
  theme(legend.position = "top")+
  guides(fill = "none")


Figure2A <- wrap_elements(Boxplot2 + plot_annotation(title = "A Group Differences for Levels of Glutamate, GABA and Log(Glutamate/GABA) at Visit 1")) 

# Model Glutamate
# Define new data grid
new_data_glu_year <- expand.grid(
  year = seq(min(NeuroMETs$year, na.rm = TRUE), max(NeuroMETs$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMETs$diagnose_group)),
  age_bl = mean(NeuroMETs$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1)   # or change to "male" or both if needed
)

# Predict values
preds <- predict(lmer_long_glu_groupwise_s, newdata = new_data_glu_year, re.form = NA, se.fit = TRUE)

label_positions <- new_data_glu_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 17.75, .groups = "drop")
coef_glu <- data.frame(emtrends(lmer_long_glu_groupwise_s, ~ diagnose_group, var = "year")) %>%
  mutate(pval = (emtrends(lmer_long_glu_groupwise_s, ~ diagnose_group, var = "year") %>% test())$p.value,
         pval_txt = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3))),
         label = paste0("β = ", format(year.trend, nsmall = 1, digits = 1), 
                        " [", round(lower.CL, 2), "; ", round(upper.CL, 2), "]\n", 
                        pval_txt)) %>%
  select(diagnose_group, label) %>%
  left_join(label_positions, by = "diagnose_group")

# Store predictions and confidence intervals
new_data_glu_year$pred <- preds$fit
new_data_glu_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_glu_year$upper <- preds$fit + 1.96 * preds$se.fit

p1 <- ggplot(new_data_glu_year, aes(x = year, y = pred, color = diagnose_group)) +
  
  geom_line(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(glu)), 
            aes(x = year, y = glu, group=record_id, color=diagnose_group), alpha=0.5) +
  # geom_point(data=NeuroMET %>% filter(!is.na(diagnose_group) & !is.na(glu)), 
  #            aes(x = year, y = glu, group=record_id, color=diagnose_group, size = glu_CRLB), alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(pTau_group)), aes(x = year, y = glu, color = pTau_group, size = glu_CRLB),  alpha = 0.8)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  theme_minimal() +
  labs(x = "Year", y = "Glutamate [mmol/L]",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis at visit 1") +
  facet_wrap(~diagnose_group, nrow = 1) + 
  geom_text(data = coef_glu, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
  scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  scale_size_continuous(range = c(0.8, 2.5), name = "CRLB-based size \nfor Glutamate and GABA") +  # optional legend
  guides(fill= "none")+
  theme(legend.position = 'top',
        axis.title.x=element_blank(),
        strip.background = element_rect(fill = "grey85",colour = NA),
        axis.line = element_line(colour = "black"),
        axis.ticks.x= element_line(colour = "black"))
# Model GABA

# Define new data grid
new_data_gaba_year <- expand.grid(
  year = seq(min(NeuroMETs$year, na.rm = TRUE), max(NeuroMETs$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMETs$diagnose_group)),
  age_bl = mean(NeuroMETs$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1) 
)

# Predict values
preds <- predict(lmer_long_gaba_groupwise_s, newdata = new_data_gaba_year, re.form = NA, se.fit = TRUE)

# Store predictions and confidence intervals
new_data_gaba_year$pred <- preds$fit
new_data_gaba_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_gaba_year$upper <- preds$fit + 1.96 * preds$se.fit

label_positions <- new_data_gaba_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 5.2, .groups = "drop")
coef_gaba <- data.frame(emtrends(lmer_long_gaba_groupwise_s, ~ diagnose_group, var = "year")) %>% mutate(label = paste0("β = ", round(year.trend, 2), " [",  round(lower.CL, 2), "; ", round(upper.CL, 2), "]\n", "p = ", round((emtrends(lmer_long_gaba_groupwise_s, ~ diagnose_group, var = "year")%>% test())$p.value, 3))) %>% select(diagnose_group, label)%>%
  left_join(label_positions, by = "diagnose_group")

p2 <- ggplot(new_data_gaba_year, aes(x = year, y = pred, color = diagnose_group)) +
  
  geom_line(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(gaba)), 
            aes(x = year, y = gaba, group=record_id, color=diagnose_group), alpha=0.5) +
  # geom_point(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(gaba)), 
  #            aes(x = year, y = gaba, group=record_id, color=diagnose_group), size=1, alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(pTau_group)), aes(x = year, y = gaba, color = pTau_group, size = gaba_CRLB), alpha = 0.8)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  theme_minimal() +
  labs(x = "Year", y = "GABA [mmol/L]",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis at visit 1") +
  facet_wrap(~diagnose_group, nrow = 1) +  
  geom_text(data = coef_gaba, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
  scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) + 
  scale_size_continuous(range = c(0.8, 2.5), name = "CRLB-based size \nfor Glutamate and GABA") +  # optional legend
  guides(fill= "none", size = "none") +
  theme(legend.position = 'none', 
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x= element_line(colour = "black"),
            axis.line = element_line(colour = "black"))

# Glutamate / GABA
# Define new data grid
new_data_glu_gaba_year <- expand.grid(
  year = seq(min(NeuroMETs$year, na.rm = TRUE), max(NeuroMETs$year, na.rm = TRUE), length.out = 100),
  diagnose_group = levels(factor(NeuroMETs$diagnose_group)),
  age_bl = mean(NeuroMETs$age_bl, na.rm = TRUE), # Fixed age at mean
  sex = as.factor(1) 
)

# Predict values
preds <- predict(lmer_long_glu_gaba_groupwise_s  , newdata = new_data_glu_gaba_year, re.form = NA, se.fit = TRUE)

# Store predictions and confidence intervals
new_data_glu_gaba_year$pred <- preds$fit
new_data_glu_gaba_year$lower <- preds$fit - 1.96 * preds$se.fit
new_data_glu_gaba_year$upper <- preds$fit + 1.96 * preds$se.fit

label_positions <- new_data_glu_gaba_year %>%
  group_by(diagnose_group) %>%
  summarise(x = max(year), y = 1.9, .groups = "drop")
coef_glu_gaba <- data.frame(emtrends(lmer_long_glu_gaba_groupwise_s, ~ diagnose_group, var = "year")) %>% mutate(label = paste0("β = ", round(year.trend, 2), " [",  round(lower.CL, 2), "; ", sprintf("%.2f", upper.CL), "]\n", "p = ", round((emtrends(lmer_long_glu_gaba_groupwise_s, ~ diagnose_group, var = "year")%>% test())$p.value, 3))) %>% select(diagnose_group, label)%>%
  left_join(label_positions, by = "diagnose_group")

p3 <- ggplot(new_data_glu_gaba_year, aes(x = year, y = pred, color = diagnose_group)) +
 
  geom_line(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(glu_gaba)), 
            aes(x = year, y = glu_gaba, group=record_id, color=diagnose_group), alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(diagnose_group) & !is.na(glu_gaba)), 
             aes(x = year, y = glu_gaba, group=record_id, color=diagnose_group), size=1, alpha=0.5) +
  geom_point(data=NeuroMETs %>% filter(!is.na(pTau_group)), aes(x = year, y = glu_gaba, color = pTau_group), size = 1.6, alpha = 0.9)+
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = diagnose_group), linetype = 0, alpha = 0.2) +
  theme_minimal() +
  labs(x = "Year", y = "log(Glutamate/GABA)",
       color = "Clinical Diagnosis at visit 1",
       fill = "Clinical Diagnosis at visit 1") +
  facet_wrap(~diagnose_group, nrow = 1) +  
  geom_text(data = coef_glu_gaba, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3)+
scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  guides(fill= "none") +
  theme(legend.position = 'none', 
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.ticks.x= element_line(colour = "black"),
            axis.line = element_line(colour = "black"))


Figure2B <-wrap_elements(( p1/p2/p3 ) + plot_annotation(title = "B Group-wise Longitudinal Trajectories of Glutamate, GABA and Log(Glutamate/GABA)"))

Figure2 <- Figure2A /Figure2B+ 
  plot_layout(heights = c(2.3, 5))
Figure2

ggsave(
  filename = "Figure2.jpeg",
  plot = Figure2,
  width = 10,  
  height = 13,
  dpi = 300
)
