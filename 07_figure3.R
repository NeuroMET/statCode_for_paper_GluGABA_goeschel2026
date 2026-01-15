ylim_min_glu <- 9.9
ylim_max_glu <- 16.85
ylim_min_gaba <- 1.75 
ylim_max_gaba <- 4.80
ylim_min_glu_gaba <- 0.95
ylim_max_glu_gaba <- 1.90


coef_pTau_glu <- data.frame(label = coef_glu$label[1], x = 10, y= 17)
coef_pTau_gaba <- data.frame(label = coef_gaba$label[1], x = 10, y= 4.9)
coef_pTau_glu_gaba <- data.frame(label = coef_glu_gaba$label[1], x =10, y=1.9)

plot_from_model <- function(model, data, ylim_min, ylim_max, outcome_var, y_lab, size_var = NULL, coef_pTau) {
  # Create new data for predictions
  new_data <- expand.grid(
    pTau = seq(min(data$pTau, na.rm = TRUE), max(data$pTau, na.rm = TRUE), length.out = 100),
    year = 0,
    age_bl = mean(data$age_bl, na.rm = TRUE),
    sex = as.factor(1)   # or change to "male" or both if needed
  )

  # Predict values and confidence intervals
  preds <- predict(model, newdata = new_data, re.form = NA, se.fit = TRUE)
  new_data$pred <- preds$fit
  new_data$lower <- preds$fit - 1.96 * preds$se.fit
  new_data$upper <- preds$fit + 1.96 * preds$se.fit

  # Create ggplot
  plot <- ggplot(new_data, aes(x = pTau, y = pred)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2A9D8F", alpha = 0.2) +
    geom_line(size = 1.2, color = "#2A9D8F")+
    geom_vline(xintercept = 2.08, linetype = "longdash", color = diagnose_colors[4], linewidth = 1) +
    labs(x = "Plasma pTau181 [pg/mL]", y = y_lab)+
    geom_text(data = coef_pTau, aes(x = x, y = y, label = label), color = "black", hjust = 1, vjust = 1, size = 3.5)+
    theme_minimal() 
  
  if (!is.null(size_var)) {
    plot <- plot +
      geom_point(
        data = data %>% filter(!is.na(get(outcome_var)) & visit == "t1"),
        aes(x = pTau, y = get(outcome_var), group = record_id, size = !!sym(size_var)),
        color = "#2A9D8F", alpha = 0.5
      ) +
      scale_size_continuous(range = c(0.8, 2.5), name = "CRLB-based size for \nGlutamate and GABA")
  } else {
    plot <- plot +
      geom_point(
        data = data %>% filter(!is.na(get(outcome_var)) & visit == "t1"),
        aes(x = pTau, y = get(outcome_var), group = record_id),
        size = 1.5, color = "#2A9D8F", alpha = 0.5
      )
  }
    
  plot+
    guides(fill = "none", size = "none") +
    coord_cartesian(ylim = c(ylim_min, ylim_max))

  return(list(ABC = plot))
}
A_data <- plot_from_model(model_glu, NeuroMET, outcome_var= "glu", y_lab = "Glutamate [mmol/L]",  ylim_min = ylim_min_glu, ylim_max= ylim_max_glu, size_var = "glu_CRLB", coef_pTau = coef_pTau_glu)
B_data <- plot_from_model(model_gaba, NeuroMET, "gaba", "GABA [mmol/L]", ylim_min = ylim_min_gaba , ylim_max= ylim_max_gaba, size_var = "gaba_CRLB", coef_pTau = coef_pTau_gaba)
C_data <- plot_from_model(model_glu_gaba, NeuroMET, "glu_gaba", "log(Glutamate/GABA)", ylim_min = ylim_min_glu_gaba , ylim_max= ylim_max_glu_gaba, coef_pTau =coef_pTau_glu_gaba)

A <- A_data$ABC+ guides (size = "none")
B <- B_data$ABC 
C <- C_data$ABC 


##### plot 3D

new_data <- expand.grid(
    diagnose_group = levels(NeuroMET$diagnose_group)[1:4],
    glu = seq(min(NeuroMET$glu, na.rm = TRUE), max(NeuroMET$glu, na.rm = TRUE), length.out = 100),
    year = 0,
    age_bl = mean(NeuroMET$age_bl, na.rm = TRUE),
    edu = mean(NeuroMET$edu, na.rm = TRUE),
    sex = as.factor(1)  # or change to "male" or both if needed
  )

  # Predict values and confidence intervals
  preds <- predict(lmer_nmm_glu_groupwise, newdata = new_data, re.form = NA, se.fit = TRUE)
  new_data$pred <- preds$fit
  new_data$lower <- preds$fit - 1.96 * preds$se.fit
  new_data$upper <- preds$fit + 1.96 * preds$se.fit

  # Create ggplot
  D <- ggplot(new_data, aes(x = glu, y = pred, color =  diagnose_group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill =  diagnose_group), alpha = 0.2, colour = NA) +
    geom_line(size = 1.2) +
    geom_point(data = NeuroMET %>% filter(!is.na(glu) & visit == "t1"),
               aes(x = glu, y = nmm, group = diagnose_group, size = glu_CRLB), 
               alpha = 0.5) +
    labs(x = "Glutamate [mmol/L]", y = "Memory Ability (NMM)", color = "", fill = "") +
      scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  scale_size_continuous(range = c(0.8, 2.5), name = "CRLB-based size for \nGlutamate and GABA")+
    theme_minimal()+ 
  guides(color = "none", fill = "none", size = "none") 
  # theme(legend.position = c(0.4, 0.95), legend.direction = "horizontal",legend.box = "horizontal")
  
##### plot 3E

new_data <- expand.grid(
    diagnose_group = levels(NeuroMET$diagnose_group)[1:4],
    gaba = seq(min(NeuroMET$gaba, na.rm = TRUE), max(NeuroMET$gaba, na.rm = TRUE), length.out = 100),
    year = 0,
    age_bl = mean(NeuroMET$age_bl, na.rm = TRUE),
    edu = mean(NeuroMET$edu, na.rm = TRUE),
    sex = as.factor(1)  # or change to "male" or both if needed
  )

  # Predict values and confidence intervals
  preds <- predict(lmer_nmm_gaba_groupwise, newdata = new_data, re.form = NA, se.fit = TRUE)
  new_data$pred <- preds$fit
  new_data$lower <- preds$fit - 1.96 * preds$se.fit
  new_data$upper <- preds$fit + 1.96 * preds$se.fit

  # Create ggplot
  E <- ggplot(new_data, aes(x = gaba, y = pred, color =  diagnose_group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill =  diagnose_group), alpha = 0.2, colour = NA) +
    geom_line(size = 1.2) +
    geom_point(data = NeuroMET %>% filter(!is.na(gaba) & visit == "t1"),
               aes(x = gaba, y = nmm, group = diagnose_group, size = gaba_CRLB), 
                alpha = 0.5) +
    labs(x = "GABA [mmol/L]", y = "Memory Ability (NMM)", color = "", fill = "") +
      scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  scale_size_continuous(range = c(0.8, 2.5), name = "CRLB-based size for \nGlutamate and GABA")+
  theme_minimal()
  # guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE), size = "none") +
  # theme(legend.position = c(0.4, 0.95), legend.direction = "horizontal",legend.box = "horizontal")
  
##### plot 3F

new_data <- expand.grid(
    diagnose_group = levels(NeuroMET$diagnose_group)[1:4],
    glu_gaba = seq(min(NeuroMET$glu_gaba, na.rm = TRUE), max(NeuroMET$glu_gaba, na.rm = TRUE), length.out = 100),
    year = 0,
    age_bl = mean(NeuroMET$age_bl, na.rm = TRUE),
    edu = mean(NeuroMET$edu, na.rm = TRUE),
    sex = as.factor(1)  # or change to "male" or both if needed
  )

  # Predict values and confidence intervals
  preds <- predict(lmer_nmm_glu_gaba_groupwise, newdata = new_data, re.form = NA, se.fit = TRUE)
  new_data$pred <- preds$fit
  new_data$lower <- preds$fit - 1.96 * preds$se.fit
  new_data$upper <- preds$fit + 1.96 * preds$se.fit

  # Create ggplot
  F <- ggplot(new_data, aes(x = glu_gaba, y = pred, color =  diagnose_group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill =  diagnose_group), alpha = 0.2, colour = NA) +
    geom_line(size = 1.2) +
    geom_point(data = NeuroMET %>% filter(!is.na(glu_gaba) & visit == "t1"),
               aes(x = glu_gaba, y = nmm, group = diagnose_group), 
               size = 1.5, alpha = 0.5) +
    labs(x = "log(Glutamate/GABA)", y = "Memory Ability (NMM)", color = "", fill = "") +
      scale_color_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +
  scale_fill_manual(values=c(diagnose_colors, diagnose_colors[4], diagnose_colors[1])) +  
  theme_minimal()+
  guides(color = "none", fill = "none", size = "none") 
  # theme(legend.position = c(0.4, 0.95), legend.direction = "horizontal",legend.box = "horizontal")
  
Plot3A <- (A+B+C) &
  theme(legend.position = "top")
Plot3A <- wrap_elements(Plot3A + plot_annotation(title = "A Associations with Continuous Plasma p-Tau181 at Visit 1"))
  
Plot3B <- (D+E+F) &
  theme(legend.position = "top") 
Plot3B <- wrap_elements(Plot3B + plot_annotation(title = "B Associations with Memory Ability at Visit 1"))
  
Plot3 <-   Plot3A / Plot3B+
  plot_layout(heights = c(1, 1.15))+ 
  theme(plot.tag.position = "topleft")

Plot3 

ggsave(
  filename = "Figure3.jpeg",
  plot = Plot5,
  width = 9,  
  height = 8,
  dpi = 300
)
