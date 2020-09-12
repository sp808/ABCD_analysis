## Install packages
library(dplyr)
library(ggplot2)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)

# Import and subset data
setwd("~/Projects/DEAPext")
data <- readRDS("nda2.0.2.Rds")
data_subset <- select(data, src_subject_id,eventname,ksads_ptsd_raw_754_p:ksads_ptsd_raw_770_p,
                      nihtbx_cryst_agecorrected,nihtbx_fluidcomp_agecorrected,nihtbx_totalcomp_agecorrected,mrif_score,fsqc_qc)

# Filters
data_subset <- filter(data_subset, eventname == "baseline_year_1_arm_1")
data_subset <- filter(data_subset, mrif_score == "No abnormal findings" | mrif_score == "Normal anatomical variant of no clinical significance")
data_subset <- filter(data_subset, fsqc_qc == "accept")
# Convert Yes/No to binary (hidden
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_754_p = ifelse(ksads_ptsd_raw_754_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_755_p = ifelse(ksads_ptsd_raw_755_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_756_p = ifelse(ksads_ptsd_raw_756_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_757_p = ifelse(ksads_ptsd_raw_757_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_758_p = ifelse(ksads_ptsd_raw_758_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_759_p = ifelse(ksads_ptsd_raw_759_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_760_p = ifelse(ksads_ptsd_raw_760_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_761_p = ifelse(ksads_ptsd_raw_761_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_762_p = ifelse(ksads_ptsd_raw_762_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_763_p = ifelse(ksads_ptsd_raw_763_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_764_p = ifelse(ksads_ptsd_raw_764_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_765_p = ifelse(ksads_ptsd_raw_765_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_766_p = ifelse(ksads_ptsd_raw_766_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_767_p = ifelse(ksads_ptsd_raw_767_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_768_p = ifelse(ksads_ptsd_raw_768_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_769_p = ifelse(ksads_ptsd_raw_769_p == "Yes",1,0))
data_subset <- data_subset %>%
  mutate(ksads_ptsd_raw_770_p = ifelse(ksads_ptsd_raw_770_p == "Yes",1,0))
df_names = colnames(data_subset) # grab column names
ptsd_summated <- rowSums(data_subset[,c(df_names[3:18])], na.rm=TRUE) # row-wise sums of columns 3-18
data_subset$summated <- ptsd_summated

# Histogram
plot <- data_subset %>%
  filter( summated>2) %>%
  ggplot( aes(x=summated)) +
  geom_histogram( binwidth=0.5, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Bin size = 1") +
  theme(
    plot.title = element_text(size=15)
  )
plot
