# 1. Sample Size per Group

n_counts <- NeuroMET %>%
  filter (visit == "t1")%>%
  count(diagnose_group) %>%
  mutate(label = paste0("n = ", n)) %>%
  select(-n) %>%
  pivot_wider(names_from = diagnose_group, values_from = label) %>%
  mutate(variable = "Sample size")


# 2. Continuous variable summaries with 95% CIs

vars_cont <- c("age", "edu", "mmse", "nmm", "glu", "gaba", "glu_gaba", "pTau", "GFAP")

cont_summary <- vars_cont %>%
  map_dfr(function(var) {
    NeuroMET %>%
      filter (visit == "t1")%>%
      group_by(diagnose_group) %>%
      summarise(
        mean = mean(.data[[var]], na.rm = TRUE),
        se = sd(.data[[var]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[var]]))),
        n = sum(!is.na(.data[[var]])),
        .groups = "drop"
      ) %>%
            mutate(
        ci_low = mean - qt(0.975, df = n - 1) * se,
        ci_high = mean + qt(0.975, df = n - 1) * se,
        mean_ci = if (var %in% c("age", "edu", "mmse", "glu")) {
          sprintf("%.1f [%.1f;%.1f]", mean, ci_low, ci_high)
        } else {
          sprintf("%.2f [%.2f;%.2f]", mean, ci_low, ci_high)
        },
        variable = var
      ) %>%
      select(variable, diagnose_group, mean_ci)
  })

cont_table_wide <- cont_summary %>%
  pivot_wider(names_from = diagnose_group, values_from = mean_ci)

# 3. ANCOVA models and adjusted p-values
models <- list(
  mmse     = lm(mmse ~ diagnose_group + age + sex + edu, data = NeuroMET%>%
      filter (visit == "t1")),
  nmm      = lm(nmm ~ diagnose_group + age + sex + edu, data = NeuroMET%>%
      filter (visit == "t1")),
  glu      = lm(glu ~ diagnose_group + age + sex, data = NeuroMET%>%
      filter (visit == "t1"), weights = glu_CRLB),
  gaba     = lm(gaba ~ diagnose_group + age + sex, data = NeuroMET%>%
      filter (visit == "t1"), weights = gaba_CRLB),
  glu_gaba = lm(glu_gaba ~ diagnose_group + age + sex, data = NeuroMET%>%
      filter (visit == "t1")),
  pTau     = lm(pTau ~ diagnose_group + age + sex, data = NeuroMET%>%
      filter (visit == "t1")),
  GFAP     = lm(GFAP ~ diagnose_group + age + sex, data = NeuroMET%>%
      filter (visit == "t1")),
  age      = lm(age ~ diagnose_group, data = NeuroMET%>%
      filter (visit == "t1")),
  edu      = lm(edu ~ diagnose_group, data = NeuroMET%>%
      filter (visit == "t1"))
)

pvals_cont <- map_dfr(names(models), function(var) {
  model <- models[[var]]
  tibble(
    variable = var,
    `p-value` = anova(model)["diagnose_group", "Pr(>F)"] %>% sprintf("%.3f", .)
  )
})

cont_table_wide <- cont_table_wide %>%
    left_join(pvals_cont, by = "variable")


# 4. Categorical variable summaries and p-values

categorical_vars <- c("sex", "APOEe4")

cat_summary <- categorical_vars %>%
  map_dfr(function(var) {
    NeuroMET%>%
      filter (visit == "t1") %>%
      group_by(diagnose_group, value = .data[[var]]) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(diagnose_group) %>%
      mutate(pct = n / sum(n) * 100) %>%
      mutate(label = paste0(n, " (", round(pct, 0), "%)")) %>%
      mutate(variable = paste0(var, " = ", value)) %>%
      select(variable, diagnose_group, label)
  })

cat_table_wide <- cat_summary %>%
  pivot_wider(names_from = diagnose_group, values_from = label)


# p-values (Fisher test if needed)
pvals_cat <- categorical_vars %>%
  map_dfr(function(var) {
    tbl <- table(NeuroMET[[var]], NeuroMET$diagnose_group)
    pval <- tryCatch({
      if (any(tbl < 5)) {
        fisher.test(tbl)$p.value
      } else {
        chisq.test(tbl)$p.value
      }
    }, error = function(e) NA)

    tibble(base_var = var, `p-value` = ifelse(is.na(pval), NA, sprintf("%.3f", pval)))
  })
pvals_cat_2 <- tibble(
  variable = c("sex = 0", "sex = 1", "APOEe4 = non-carrier", "APOEe4 = carrier"),
  `p-value` = c(
    pvals_cat$`p-value`[1],  # sex p-value repeated
    pvals_cat$`p-value`[1],
    pvals_cat$`p-value`[2],  # APOEe4 p-value repeated
    pvals_cat$`p-value`[2]
  )
)


# Expand p-values for each categorical variable level
cat_table_wide <- cat_table_wide %>%
    left_join(pvals_cat_2, by = "variable")


# 5. Combine everything


# Combine sample size, continuous, categorical
final_table <- bind_rows(
  n_counts,
  cont_table_wide,
  cat_table_wide
) %>%
  select(variable, "CN Aβ-", "CN Aβ+", "MCI Aβ+", "AD Aβ+", `p-value`) %>%
  mutate(variable = case_when(
    variable == "n" ~ "Sample size",
    variable == "sex = 1" ~ "Sex = Female",
    variable == "sex = 0" ~ "Sex = Male",
    variable == "APOEe4 = 0" ~ "APOEe4 = 0",
    variable == "APOEe4 = 1" ~ "APOEe4 = 1",
    TRUE ~ variable
  )) %>%
  filter (!(variable == "Sex = Male" | variable == "APOEe4 = non-carrier"))

# Rename variables to display-friendly format
final_table_cleaned <- final_table %>%
  mutate(variable = case_when(
    variable == "age" ~ "Age [years]",
    variable == "Sex = Female" ~ "Sex (Female)",
    variable == "edu" ~ "Education [years]",
    variable == "APOEe4 = carrier" ~ "APOEε4 (carrier)",
    variable == "mmse" ~ "MMSE",
    variable == "nmm" ~ "Memory (NMM)",
    variable == "glu" ~ "MRS Glutamate [mmol/L]",
    variable == "gaba" ~ "MRS GABA [mmol/L]",
    variable == "glu_gaba" ~ "log(Glutamate/GABA)",
    variable == "pTau" ~ "Plasma pTau-181 [pg/mL]",
    variable == "GFAP" ~ "GFAP [pg/mL]",
    TRUE ~ variable
  ))

desired_order <- c(
  "Sample size",
  "Age [years]",
  "Sex (Female)",
  "Education [years]",
  "APOEε4 (carrier)",
  "MMSE",
  "Memory (NMM)",
  "MRS Glutamate [mmol/L]",
  "MRS GABA [mmol/L]",
  "log(Glutamate/GABA)",
  "Plasma pTau-181 [pg/mL]",
  "GFAP [pg/mL]"
)

final_table_ordered <- final_table_cleaned %>%
  mutate(variable = factor(variable, levels = desired_order)) %>%
  arrange(variable) %>%
  mutate(`p-value` = ifelse(as.numeric(`p-value`) <= 0.001, "<0.001", `p-value`))

final_table_ordered %>%
  kable(
    format = "html",  # use "latex" for PDF or "markdown" for markdown docs
    digits = 2,
    col.names = c("Variable", "CN Aβ-", "CN Aβ+", "MCI Aβ+", "AD Aβ+", "p-value"),
    caption = "<b>Table 1. Descriptive statistics and group comparisons of the characteristics of the study population at visit 1.</b> For continuous variables (e.g., age, MMSE, metabolite concentrations), the mean and 95% confidence interval are shown. P-values are derived from linear models: models for MRS data and plasma pTau181 and GFAP are adjusted for age and sex; models for MMSE and memory (NMM) are additionally adjusted for education. For categorical variables (sex and APOEε4 status), absolute and relative frequencies are presented, and p-values are based on Chi-square tests. <i>Abbreviations: AD = Alzheimer’s Disease, APOEε4 = Apolipoprotein E epsilon 4 allele, CN = Cognitively Normal, GABA = Gamma-Aminobutyric Acid, Glu = Glutamate, MCI = Mild Cognitive Impairment, MMSE = Mini-Mental State Examination, MRS = Magnetic Resonance Spectroscopy, pTau181 = Plasma phosphorylated tau at threonine 181</i>", align = "c"
  ) %>%
  kable_paper()
