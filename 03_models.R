# Models which excluded one individual with disproportional influence on the model. Model fit was tested for several spline dfs using anova() and chosing accordin best AIC/BIC

NeuroMETs <- NeuroMET %>% filter (!(record_id == "NeuroMetXXX" & visit == "t4"))

# Model 1, 1 individual excluded

lmer_long_glu_groupwise_s       <- lmer(glu ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 2) + sex + (1|record_id),
       data = NeuroMETs,REML = TRUE, na.action = na.omit, weight = glu_CRLB)

lmer_long_gaba_groupwise_s <- lmer(gaba ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex+ (1|record_id),
       data = NeuroMETs,REML = TRUE, na.action = na.omit, weight = gaba_CRLB)

lmer_long_glu_gaba_groupwise_s   <- lmer(glu_gaba ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 2) + sex+ (1|record_id),
       data = NeuroMETs,REML = TRUE, na.action = na.omit)

lmer_long_pTau_groupwise_s <- lmer(pTau ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex +(1|record_id),
       data = NeuroMETs,REML = TRUE, na.action = na.omit)

lmer_long_nmm_groupwise_s <- lmer(nmm ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex+ (1|record_id),
       data = NeuroMETs,REML = TRUE, na.action = na.omit)

lmer_long_gfap_groupwise_s <- lmer(GFAP ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex + (1|record_id), data = NeuroMETs,REML = TRUE, na.action = na.omit)

lmer_long_gln_groupwise_s       <- lmer(gln ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 2) + sex + (1|record_id), data = NeuroMETs,REML = TRUE, na.action = na.omit, weight = gln_CRLB)

lmer_long_glu_gln_groupwise_s       <- lmer(glu_gln ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex + (1|record_id), data = NeuroMETs,REML = TRUE, na.action = na.omit)


# Model 1, all data

lmer_long_glu_groupwise       <- lmer(glu ~ diagnose_group * ns(year, df = 2) + ns(age_bl, df = 2) + sex + (1|record_id), data = NeuroMET,REML = TRUE, na.action = na.omit, weight = glu_CRLB)

lmer_long_gaba_groupwise <- lmer(gaba ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex+ (1|record_id),
       data = NeuroMET,REML = TRUE, na.action = na.omit, weight = gaba_CRLB)

lmer_long_glu_gaba_groupwise   <- lmer(glu_gaba ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 2) + sex+ (1|record_id),
       data = NeuroMET,REML = TRUE, na.action = na.omit)

lmer_long_pTau_groupwise <- lmer(pTau ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex +(1|record_id),
       data = NeuroMET,REML = TRUE, na.action = na.omit)

lmer_long_nmm_groupwise <- lmer(nmm ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex+ (1|record_id),
       data = NeuroMET,REML = TRUE, na.action = na.omit)

lmer_long_gfap_groupwise <- lmer(GFAP ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex + (1|record_id), data = NeuroMET,REML = TRUE, na.action = na.omit)

lmer_long_gln_groupwise       <- lmer(gln ~ diagnose_group * ns(year, df = 2) + ns(age_bl, df = 2) + sex + (1|record_id), data = NeuroMET,REML = TRUE, na.action = na.omit, weight = gln_CRLB)

lmer_long_glu_gln_groupwise       <- lmer(glu_gln ~ diagnose_group * ns(year, df = 1) + ns(age_bl, df = 1) + sex + (1|record_id), data = NeuroMET,REML = TRUE, na.action = na.omit)


# Model 2

model_glu <- lmer(glu ~ ns(pTau, df = 1) * ns(year, df = 1) + ns(age_bl, df = 2) + sex + (1 | record_id),
                  data = NeuroMET %>% filter(!is.na(pTau)),
                  REML = TRUE, na.action = na.omit, weight = glu_CRLB)

model_gaba <- lmer(gaba ~ ns(pTau, df = 1) * ns(year, df = 1) + ns(age_bl, df = 1) + sex + (1 | record_id),
                   data = NeuroMET %>% filter(!is.na(pTau)),
                   REML = TRUE, na.action = na.omit, weight = gaba_CRLB)

model_glu_gaba <- lmer(glu_gaba ~ ns(pTau, df = 1) * ns(year, df = 1) + ns(age_bl, df = 2) + sex + (1 | record_id),
                       data = NeuroMET %>% filter(!is.na(pTau)),
                       REML = TRUE, na.action = na.omit)

# Model 3

lmer_nmm_glu_groupwise <- lmer(nmm ~ diagnose_group * ns(glu, df = 1) + ns(year, df = 2) + ns(age_bl, df = 1) + sex + edu + (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit, weight = glu_CRLB)

lmer_nmm_gaba_groupwise <- lmer(nmm ~ diagnose_group * ns(gaba, df = 1) + ns(year, df = 2) + ns(age_bl, df = 1) + sex + edu + (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit, weight = gaba_CRLB)

lmer_nmm_glu_gaba_groupwise <- lmer(nmm ~ diagnose_group * ns(glu_gaba, df = 1) + ns(year, df = 2) + ns(age_bl, df = 1) + sex + edu + (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit)

# Model 4

NeuroMET <- NeuroMET %>% mutate(pcc_precuneus_vol= Whole_pcc_cs_long + Whole_precuneus_cs_long)

lmer_GFAP_glu_groupwise <- lmer(glu ~ diagnose_group * ns(GFAP, df = 1) + ns(pcc_precuneus_vol ,df = 1)+ ns(year, df = 1) + ns(age_bl, df = 2) + sex+ (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit, weight = glu_CRLB)

lmer_GFAP_gaba_groupwise <- lmer(gaba ~ diagnose_group * ns(GFAP, df = 1)+ ns(pcc_precuneus_vol,df = 1) + ns(year, df = 1) + ns(age_bl, df = 1) + sex+ (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit, weight = gaba_CRLB)

lmer_GFAP_glu_gaba_groupwise <- lmer(glu_gaba ~ diagnose_group * ns(GFAP, df = 1)+ ns(pcc_precuneus_vol,df = 1)  + ns(year, df = 1) + ns(age_bl, df = 2) + sex+ (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit)

lmer_GFAP_gln_groupwise <- lmer(gln ~ diagnose_group * ns(GFAP, df = 1)+ ns(pcc_precuneus_vol,df = 1)  + ns(year, df = 1) + ns(age_bl, df = 2) + sex+ (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit, weight = gln_CRLB)

lmer_GFAP_glu_gln_groupwise <- lmer(glu_gln ~ diagnose_group * ns(GFAP, df = 1) + ns(pcc_precuneus_vol,df = 1) + ns(year, df = 1) + ns(age_bl, df = 1) + sex+ (1|record_id),
                  data = NeuroMET, REML = TRUE, na.action = na.omit)

# Model Age Groupwise

model_age_glu_groupwise       <- lmer(glu ~ diagnose_group* ns(year, df = 2) * ns(age_bl, df = 1) + sex + (1 | record_id), data = NeuroMET,REML = TRUE, na.action = na.omit, weight = glu_CRLB)

model_age_gaba_groupwise       <- lmer(gaba ~ diagnose_group * ns(year, df = 1) * ns(age_bl, df = 1)  + sex + (1 | record_id), data = NeuroMET,REML = TRUE, na.action = na.omit, weight = gaba_CRLB)

model_age_glu_gaba_groupwise       <- lmer(glu_gaba ~ diagnose_group * ns(year, df = 1) * ns(age_bl, df = 2) + sex + (1 | record_id), data = NeuroMET,REML = TRUE, na.action = na.omit)

em_sex_glu <- emmeans(lmer_long_glu_groupwise, pairwise ~ sex, cov.reduce = list(year = function(x) 0))
em_sex_gaba <- emmeans(lmer_long_gaba_groupwise, pairwise ~ sex, cov.reduce = list(year = function(x) 0))
em_sex_glu_gaba <- emmeans(lmer_long_glu_gaba_groupwise, pairwise ~ sex, cov.reduce = list(year = function(x) 0))
