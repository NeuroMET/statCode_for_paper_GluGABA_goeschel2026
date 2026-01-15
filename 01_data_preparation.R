library(readxl)
library(tidyverse)
library(kableExtra)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(knitr)
library(ggeffects)
library(splines)
library(cowplot)
library(broom)
library(glue)
library(pbkrtest)
library(MetBrewer)

diagnose_colors <- c("#376795","#72bcd5","#ffd06f","#ef8a47")

NeuroMET <- read_excel("OpenData_NeuroMET.xlsx")

NeuroMET$diagnose <- factor(NeuroMET$diagnose, labels=c('HC', 'SCD', 'MCI','AD'))
NeuroMET <- NeuroMET %>% filter(!is.na(diagnose))
NeuroMET <- NeuroMET %>% filter(!(visit == 'smartage'))
NeuroMET$visit <- factor(NeuroMET$visit, labels=c('t1', 't2','t3', 't4', 't5'))
#create a variable diagnose_bl
Diagnose_1 <- NeuroMET %>% filter(visit == 't1') %>% mutate(diagnose_bl= diagnose) %>% select(record_id, diagnose_bl)
NeuroMET <- NeuroMET %>%
  full_join(Diagnose_1, by = "record_id")

NeuroMET$APOE <- factor(NeuroMET$genotype_APOE, labels=c('22','23','24','33','34','44'))
NeuroMET$APOEe4 <- ifelse(NeuroMET$APOE=="22" | NeuroMET$APOE=="23" | NeuroMET$APOE=="33", "0", "1")
NeuroMET$APOEe4 <- factor(NeuroMET$APOEe4, labels=c('non-carrier','carrier'))

NeuroMET$edu <- as.numeric(NeuroMET$education_years_cer)
NeuroMET$sex <- factor(NeuroMET$sex)

#global cognition
NeuroMET$mmse <- as.numeric(NeuroMET$cer_mmse_sum)

#memory 
NeuroMET$nmm <- as.numeric(NeuroMET$NMM_value)

#biomarkers
NeuroMET$Ab40 <- as.numeric(NeuroMET$plasma_ab_40)
NeuroMET$Ab42 <- as.numeric(NeuroMET$plasma_ab_42)
NeuroMET$Ab42_40 <- NeuroMET$Ab42/NeuroMET$Ab40
NeuroMET$pTau <- as.numeric(NeuroMET$plasma_ptau181)
NeuroMET$GFAP <- as.numeric(NeuroMET$plasma_gfap)
NeuroMET$nfl <- as.numeric(NeuroMET$plasma_nfl)

#Calculate Slopes
ptau_slopes <- NeuroMET %>%
  group_by(record_id) %>%
  do(tidy(lm(pTau ~ year, data = .))) %>%
  filter(term == "year") %>%
  select(record_id, ptau_slope = estimate)
NeuroMET <- NeuroMET %>%
  left_join(ptau_slopes, by = "record_id")%>%
  rename(ptau_slopes = ptau_slope)

gfap_slopes <- NeuroMET %>%
  group_by(record_id) %>%
  do(tidy(lm(GFAP ~ year, data = .))) %>%
  filter(term == "year") %>%
  select(record_id, gfap_slope = estimate)
NeuroMET <- NeuroMET %>%
  left_join(gfap_slopes, by = "record_id")%>%
  rename(gfap_slopes = gfap_slope)

nfl_slopes <- NeuroMET %>%
  group_by(record_id) %>%
  do(tidy(lm(nfl ~ year, data = .))) %>%
  filter(term == "year") %>%
  select(record_id, nfl_slope = estimate)
NeuroMET <- NeuroMET %>%
  left_join(nfl_slopes, by = "record_id")%>%
  rename(nfl_slopes = nfl_slope)


pTau_sd <- sd(NeuroMET$pTau, na.rm = TRUE)
#create the two variables pTau group and at visit one
NeuroMET <- NeuroMET %>%
  group_by(record_id) %>%  # Group by participant
  mutate(pTau_v1_group = as.factor(ifelse(pTau[visit == "t1"] >= 2.08, "Aβ+", "Aβ-")),
         pTau_group= as.factor(ifelse(pTau >= 2.08, "Aβ+", "Aβ-")),
         pTau_group2 = as.factor(ifelse (pTau > 2.08 + 0.5 * pTau_sd, "high", ifelse(pTau < 2.08 - 0.3 * pTau_sd, "low", "intermediate"))),
         age_bl = age[visit == "t1"]) %>%
  ungroup()

#test
#NeuroMET %>% filter(pTau_group == "Aβ-" & visit == "t1") %>% select (record_id, visit, diagnose, diagnose_bl, pTau, pTau_group) %>% print (n = 61)

#Find cases which converted back from Aβ+ to Aβ-
cases_t1_Aß_pos <- NeuroMET %>%
  filter(pTau_v1_group == "Aβ+" & pTau <= 2.08) %>%
  select(record_id, visit, diagnose_bl, pTau_v1_group, pTau)

# Exclude one case HC pTau at visit1 of 3.42 with follow-up pTau of 1.55 and 1.64
NeuroMET <- NeuroMET %>% filter(!(record_id == "NeuroMetXXX"))

# Exclude one case SCD pTau at visit1 of 2.76 with follow-up pTau of 1.66
NeuroMET <- NeuroMET %>% filter(!(record_id == "NeuroMetXXX"))

# Exclude one case HC pTau at visit1 of 2.14 with follow-up pTau of 1.84 (no Glutamate at follow-ups available)
NeuroMET <- NeuroMET %>% filter(!(record_id == "NeuroMetXXX"))


#pTau as factor
NeuroMET$pTau_group <- as.factor(NeuroMET$pTau_group)
NeuroMET$pTau_v1_group <- as.factor(NeuroMET$pTau_v1_group)

#alternative grouping
NeuroMET <- NeuroMET %>%
  mutate(diagnose_group = case_when(
    pTau_v1_group == "Aβ-" & (diagnose_bl == "HC" | diagnose_bl == "SCD")  ~ "CN Aβ-",
    pTau_v1_group == "Aβ+" & (diagnose_bl == "HC" | diagnose_bl == "SCD") ~ "CN Aβ+", 
    pTau_v1_group == "Aβ+" & diagnose_bl == "MCI" ~ "MCI Aβ+",
    record_id == "NeuroMetXXX" ~ "MCI Aβ+",   ## define as Aβ+ because pTau at visit 1 was > 2.00 and for follow ups Aβ+
    record_id == "NeuroMetXXX" ~ "AD Aβ+",    ## define as Aβ+ because pTau at visit 1 was > 2.00 and for follow ups Aβ+
    pTau_v1_group == "Aβ+" & diagnose_bl == "AD" ~ "AD Aβ+",
    TRUE ~ "other"
  ))%>% 
  mutate(diagnose_group = factor(diagnose_group, levels = c("CN Aβ-", "CN Aβ+", "MCI Aβ+","AD Aβ+", "other" )))
  
### EXCLUDE MCIs and AD who are Aβ-, except for the MCIs NeuroMETXXX and XXX who became Aβ+ in the follow-up: "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX"
NeuroMET <- NeuroMET %>% filter(!(pTau_v1_group == "Aβ-" & (diagnose == "MCI" | diagnose == "AD") & !(record_id == "NeuroMetXXX")& !(record_id == "NeuroMetXXX")))


#########

NeuroMET$gaba <- as.numeric(NeuroMET$GABA_conc)
NeuroMET$gaba[NeuroMET$gaba == "NaN"] <- NA
NeuroMET$glu <- as.numeric(NeuroMET$Glu_conc)
NeuroMET$glu[NeuroMET$glu == "NaN"] <- NA
NeuroMET$gln <- as.numeric(NeuroMET$Gln_conc)
NeuroMET$gln[NeuroMET$gln == "NaN"] <- NA

#####EXCLUDE 8 participants without MRS data at visit 1. but keep XXX because has follow-up data! "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX" "NeuroMetXXX"
invalid_ids <- NeuroMET %>%
  filter(visit == "t1", is.na(glu), !(record_id == "NeuroMetXXX")) %>%
  pull(record_id) %>%
  unique() 
NeuroMET <- NeuroMET %>% filter (!(record_id %in% invalid_ids))


###z-scores to calc glu/gaba
mean_z <- NeuroMET %>% filter(diagnose == 'HC' & visit == 't1') %>%  summarise(across(c(glu, gaba), base::mean, na.rm = TRUE))

sd_z <- NeuroMET %>% filter(diagnose == 'HC' & visit == 't1') %>%  summarise(across(c(glu, gaba), stats::sd, na.rm = TRUE))

cols <- list("glu", "gaba")

NeuroMET[, paste0("z_score_", cols)] <- lapply(cols, function(col){ 
  x <- NeuroMET[col]
  y <- (x - as.vector(mean_z[col])) / as.vector(sd_z[col])
  y[[1]]
  })
#NeuroMET$glu_gaba <- as.numeric(NeuroMET$z_score_glu / NeuroMET$z_score_gaba)

NeuroMET <- NeuroMET %>%
  mutate(glu_gaba = log(glu / gaba))

NeuroMET <- NeuroMET %>%
  mutate(glu_gln = log(glu / gln))

#calculate rate of change (slope) for those participants with 2 or more visits
glu_slopes <- NeuroMET %>%
  filter(!is.na(glu)) %>%
  group_by(record_id) %>%
  filter(n() >= 2) %>%  # Only if they have 2+ visits
  nest() %>%
  mutate(model = map(data, ~ lm(glu ~ year, data = .x)),
         slope = map_dbl(model, ~ coef(.x)["year"])) %>%
  select(record_id, slope)
NeuroMET <- NeuroMET %>%
  left_join(glu_slopes, by = "record_id")%>%
  rename(glu_slopes = slope)

gaba_slopes <- NeuroMET %>%
  filter(!is.na(gaba)) %>%
  group_by(record_id) %>%
  filter(n() >= 2) %>%  # Only if they have 2+ visits
  nest() %>%
  mutate(model = map(data, ~ lm(gaba ~ year, data = .x)),
         slope = map_dbl(model, ~ coef(.x)["year"])) %>%
  select(record_id, slope)
NeuroMET <- NeuroMET %>%
  left_join(gaba_slopes, by = "record_id")%>%
  rename(gaba_slopes = slope)

glu_gaba_slopes <- NeuroMET %>%
  filter(!is.na(glu_gaba)) %>%
  group_by(record_id) %>%
  filter(n() >= 2) %>%  # Only if they have 2+ visits
  nest() %>%
  mutate(model = map(data, ~ lm(glu_gaba ~ year, data = .x)),
         slope = map_dbl(model, ~ coef(.x)["year"])) %>%
  select(record_id, slope)
NeuroMET <- NeuroMET %>%
  left_join(glu_gaba_slopes, by = "record_id")%>%
  rename(glu_gaba_slopes = slope)

# For weighted analyses of MRS data, CRLB are mutated to weights (mean = 1)
NeuroMET$gaba_aCRLB <- as.numeric(NeuroMET$GABA_CRLB_abs) # absolute CRLB
NeuroMET$gaba_inv_aCRLB <- 1/((NeuroMET$gaba_aCRLB)^2) # reciprocal and squared aCRLB
NeuroMET$gaba_CRLB <- ((NeuroMET$gaba_inv_aCRLB - mean(NeuroMET$gaba_inv_aCRLB, na.rm =TRUE))/ (max(NeuroMET$gaba_inv_aCRLB, na.rm =TRUE) - min(NeuroMET$gaba_inv_aCRLB, na.rm =TRUE)))+1
mean(NeuroMET$gaba_CRLB, na.rm = TRUE)   # should be 1 
sum(NeuroMET$gaba_CRLB, na.rm=TRUE)      # should be number of participants

NeuroMET$glu_aCRLB <- as.numeric(NeuroMET$Glu_CRLB_abs) # absolute CRLB
NeuroMET$glu_inv_aCRLB <- 1/((NeuroMET$glu_aCRLB)^2) # reciprocal and squared aCRLB
NeuroMET$glu_CRLB <- ((NeuroMET$glu_inv_aCRLB - mean(NeuroMET$glu_inv_aCRLB, na.rm =TRUE))/ (max(NeuroMET$glu_inv_aCRLB, na.rm =TRUE) - min(NeuroMET$glu_inv_aCRLB, na.rm =TRUE)))+1
mean(NeuroMET$glu_CRLB, na.rm = TRUE)   # should be 1 
sum(NeuroMET$glu_CRLB, na.rm=TRUE)      # should be number of participants

NeuroMET$gln_aCRLB <- as.numeric(NeuroMET$Gln_CRLB_abs) # absolute CRLB
NeuroMET$gln_inv_aCRLB <- 1/((NeuroMET$gln_aCRLB)^2) # reciprocal and squared aCRLB
NeuroMET$gln_CRLB <- ((NeuroMET$gln_inv_aCRLB - mean(NeuroMET$gln_inv_aCRLB, na.rm =TRUE))/ (max(NeuroMET$gln_inv_aCRLB, na.rm =TRUE) - min(NeuroMET$gln_inv_aCRLB, na.rm =TRUE)))+1
mean(NeuroMET$gln_CRLB, na.rm = TRUE)   # should be 1 
sum(NeuroMET$gln_CRLB, na.rm=TRUE)      # should be number of participants

#Add drop-out reasons to the dataset

dropout_reasons <- read_excel("drop-out_reasons_new.xlsx") %>% mutate(record_id = Pseudonym)

# Join to NeuroMET by record_id
NeuroMET <- NeuroMET %>%
  left_join(dropout_reasons, by = "record_id")

#Summarize drop-out reasons
#NeuroMET %>% filter(visit == "t1") %>% group_by(dropout_reasons)%>%summarize(n= n()) 

#NeuroMET %>% filter(visit == "t1") %>% group_by(diagnose_group)%>%summarize(n= n()) 

#exclude MRI volumes with insufficient quality
NeuroMET <- NeuroMET %>% mutate(pcc_lh_cs = lh_posteriorcingulate_cs*qc_lh_posteriorcingulate_cs,
                                pcc_rh_cs = rh_posteriorcingulate_cs*qc_rh_posteriorcingulate_cs,
                                pcc_lh_long = lh_posteriorcingulate_volume_long* qc_lh_posteriorcingulate_long,
                                pcc_rh_long = rh_posteriorcingulate_volume_long* qc_rh_posteriorcingulate_long,
                                precuneus_lh_cs = lh_precuneus_cs*qc_lh_posteriorcingulate_cs,
                                precuneus_rh_cs = rh_precuneus_cs*qc_rh_posteriorcingulate_cs,
                                precuneus_lh_long = lh_precuneus_volume_long*qc_lh_precuneus_long,
                                precuneus_rh_long = rh_precuneus_volume_long*qc_rh_precuneus_long,
                                etiv_cs= EstimatedTotalIntraCranialVol_cs,
                                etiv_long= EstimatedTotalIntraCranialVol_long) %>% 
                         mutate(pcc_lh_cs = replace(pcc_lh_cs, pcc_lh_cs==0, NA),
                                pcc_rh_cs = replace(pcc_rh_cs, pcc_rh_cs==0, NA),
                                pcc_lh_long = replace(pcc_lh_long, pcc_lh_long==0, NA),
                                pcc_rh_long = replace(pcc_rh_long, pcc_rh_long==0, NA),
                                precuneus_lh_cs = replace(precuneus_lh_cs, precuneus_lh_cs==0, NA),
                                precuneus_rh_cs= replace(precuneus_rh_cs, precuneus_rh_cs==0, NA),
                                precuneus_lh_long = replace(precuneus_lh_long, precuneus_lh_long==0, NA),
                                precuneus_rh_long = replace(precuneus_rh_long, precuneus_rh_long==0, NA),
                                etiv_cs = replace(etiv_cs, etiv_cs==0, NA)
                                )

# Calculate mean or copy unilateral volume data if other side did not pass quality control
mean_lh_rh <- function (lh, rh){
  if (! is.na(lh)){
    if (! is.na(rh)) {
      (lh + rh ) /2
    }
    else if(is.na(rh)){
      (lh + lh ) /2
    }
  }
  else if (is.na(lh)){
    if (! is.na(rh)){
      (rh + rh ) /2
    }
  } else {NA}
}

NeuroMET <- NeuroMET %>%
  mutate(Whole_pcc_cs = mapply(mean_lh_rh, pcc_lh_cs, pcc_rh_cs),
         Whole_precuneus_cs = mapply(mean_lh_rh, precuneus_lh_cs, precuneus_rh_cs),
         Whole_pcc_long = mapply(mean_lh_rh, pcc_lh_long, pcc_rh_long),
         Whole_precuneus_long = mapply(mean_lh_rh, precuneus_lh_long, precuneus_rh_long))

# mapply gives out lists with NULLs instead of NA. Replace NULL by NA and unlist volumetric data
NeuroMET$Whole_pcc_cs[sapply(NeuroMET$Whole_pcc_cs, is.null)] <- NA
NeuroMET$Whole_precuneus_cs[sapply(NeuroMET$Whole_precuneus_cs, is.null)] <- NA
NeuroMET$Whole_pcc_long[sapply(NeuroMET$Whole_pcc_long, is.null)] <- NA
NeuroMET$Whole_precuneus_long[sapply(NeuroMET$Whole_precuneus_long, is.null)] <- NA


NeuroMET <- NeuroMET %>%
  mutate(Whole_pcc_cs = unlist(Whole_pcc_cs),
         Whole_precuneus_cs = unlist(Whole_precuneus_cs),
         Whole_pcc_long = unlist(Whole_pcc_long),
         Whole_precuneus_long = unlist(Whole_precuneus_long)
         )

#adjust to eTIV cross-sectional (adjusted to HC in T1 only)
ROI_list <- list("Whole_pcc_cs", "Whole_precuneus_cs")
NeuroMET_HC_t1 <- NeuroMET %>%
  filter (diagnose_group=="CN Aβ-" & year ==0)
meanTIV_HC = mean(NeuroMET_HC_t1$EstimatedTotalIntraCranialVol_cs, na.rm = TRUE)
adj_ROIs <- sapply(ROI_list, function (ROI) {
  NeuroMET_HC_t1$ROI_HC <- unlist(NeuroMET_HC_t1[, ROI])
  regr <- lm(ROI_HC ~ EstimatedTotalIntraCranialVol_cs, NeuroMET_HC_t1)
  slope <- summary(regr)$coefficients[2,1]
  adjVol = unlist(NeuroMET[, ROI]) - slope * (NeuroMET$EstimatedTotalIntraCranialVol_cs - meanTIV_HC)
})
adj_ROIs_cs <- data.frame(record_id = NeuroMET$record_id, year = NeuroMET$year, as.data.frame(adj_ROIs))%>%
  mutate(Whole_pcc_cs_adj = V1, Whole_precuneus_cs_adj = V2)%>%
  select(-V1, -V2)

#adjust to eTIV longitudinal (adjusted to HC in T1 only)
ROI_list <- list("Whole_pcc_long", "Whole_precuneus_long")
NeuroMET_HC_t1 <- NeuroMET %>%
  filter (diagnose_group=="CN Aβ-" & year ==0)
meanTIV_HC = mean(NeuroMET_HC_t1$EstimatedTotalIntraCranialVol_long, na.rm = TRUE)
adj_ROIs <- sapply(ROI_list, function (ROI) {
  NeuroMET_HC_t1$ROI_HC <- unlist(NeuroMET_HC_t1[, ROI])
  regr <- lm(ROI_HC ~ EstimatedTotalIntraCranialVol_long, NeuroMET_HC_t1)
  slope <- summary(regr)$coefficients[2,1]
  adjVol = unlist(NeuroMET[, ROI]) - slope * (NeuroMET$EstimatedTotalIntraCranialVol_long - meanTIV_HC)
})
adj_ROIs_long <- data.frame(record_id = NeuroMET$record_id, year = NeuroMET$year, as.data.frame(adj_ROIs)) %>%
  mutate(Whole_pcc_long_adj = V1, Whole_precuneus_long_adj = V2)%>%
  select(-V1, -V2)

NeuroMET <- NeuroMET%>%
  full_join(adj_ROIs_cs, by = c("record_id", "year"))%>%
  full_join(adj_ROIs_long, by = c("record_id", "year"))%>%
  # merge cross sectional and longitudinal data. When there is no Whole_hip_long_adj, then select cross-sectional data 
  mutate(Whole_pcc_cs_long = ifelse(is.na(Whole_pcc_cs_adj) != is.na(Whole_pcc_long_adj), Whole_pcc_cs_adj, Whole_pcc_long_adj),
         Whole_precuneus_cs_long = ifelse(is.na(Whole_precuneus_cs_adj) != is.na(Whole_precuneus_long_adj), Whole_precuneus_cs_adj, Whole_precuneus_long_adj)) 
