## Housekeeping #######
library(dplyr)
library(ggplot2)
library(ggExtra)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)
library(RColorBrewer)
library(reshape2)
library(Hmisc)
library(readr)
library(tidyr)
library(viridis)
library(grid)
library(gridExtra)
library(hexbin)
library(paletteer)
library(cowplot)

# Import dataframe #######
data <- readRDS("nda2.0.2.Rds")

# Subset and filter data #########
data_subset <- select(data, src_subject_id,eventname,sports_activity_ss_read_hours_p, sports_activity_read_years_p, screentime_wkdy_1,
                      nihtbx_cryst_agecorrected,nihtbx_fluidcomp_agecorrected,nihtbx_totalcomp_agecorrected,mrif_score,fsqc_qc)
# convert to daily reading to match television 
data_subset$sports_activity_ss_read_hours_p <- data_subset$sports_activity_ss_read_hours_p/7

# Filters
data_subset <- filter(data_subset, eventname == "baseline_year_1_arm_1")
data_subset <- filter(data_subset, mrif_score == "No abnormal findings" | mrif_score == "Normal anatomical variant of no clinical significance")
data_subset <- filter(data_subset, fsqc_qc == "accept")
# filter top 15% outliers, max() won't work unless you add na.rm=T to not count NAs
data_subset <- filter(data_subset, sports_activity_ss_read_hours_p < (max(data_subset$sports_activity_ss_read_hours_p, na.rm=T) * 0.5)) 
data_subset <- data_subset[!is.na(data_subset$screentime_wkdy_1), ] # filter NA


# Cognition scatter plots ######
# Reorder screentime factor correctly and colorcode scatter by screentime category with scale_color_viridis
a <- ggplot(data_subset %>% mutate(screentime_wkdy_1 = factor(screentime_wkdy_1, levels=c("None", "< 30 minutes", ".5", "1 hour", "2 hours", "3 hours", "4+ hours"))),
  aes(x=sports_activity_ss_read_hours_p, y=nihtbx_cryst_agecorrected, color = screentime_wkdy_1)) + 
  geom_point(shape=16, size = 2) + scale_color_viridis(discrete = TRUE) + labs(y =NULL, x = NULL) + 
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, lwd=1, color="orchid", linetype="dashed") + theme(legend.position = "none")

b <- ggplot(data_subset %>% mutate(screentime_wkdy_1 = factor(screentime_wkdy_1, levels=c("None", "< 30 minutes", ".5", "1 hour", "2 hours", "3 hours", "4+ hours"))),
  aes(x=sports_activity_ss_read_hours_p, y=nihtbx_fluidcomp_agecorrected, color = screentime_wkdy_1)) + 
  geom_point(shape=16, size = 2) + scale_color_viridis(discrete = TRUE) + labs(y =NULL, x = NULL) + theme(legend.position = "none") +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, color="orchid", linetype="dashed")

c <- ggplot(data_subset %>% mutate(screentime_wkdy_1 = factor(screentime_wkdy_1, levels=c("None", "< 30 minutes", ".5", "1 hour", "2 hours", "3 hours", "4+ hours"))),
  aes(x=sports_activity_ss_read_hours_p, y=nihtbx_totalcomp_agecorrected, color = screentime_wkdy_1)) + 
  geom_point(shape=16, size = 2) + scale_color_viridis(discrete = TRUE) + labs(y =NULL, x = NULL) + 
  theme(legend.justification = c("right"), legend.margin = margin(0, 0, 0, 0), legend.position ="none", legend.text=element_text(size=8)) + 
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, color="orchid", linetype="dashed")

grid.arrange(a,b,c, nrow = 1)
v <- arrangeGrob(a,b,c, nrow=1) # ggsave won't save a grid.arrange
ggsave("reading_cognitive.pdf", plot=v, path="./plots", dpi=150)

# Television violin plots ######
violin_theme <- theme(axis.text.x=element_blank(),
                      axis.ticks.x = element_blank(),)

v1 <- ggplot(data_subset %>% mutate(screentime_wkdy_1 = factor(screentime_wkdy_1, levels=c("None", "< 30 minutes", ".5", "1 hour", "2 hours", "3 hours", "4+ hours"))),
  aes(x=screentime_wkdy_1, y=nihtbx_cryst_agecorrected, fill=screentime_wkdy_1)) + # color= sets line, fill=sets fill
  geom_violin(trim=TRUE, alpha=0.6, lwd = 0.3) + scale_fill_viridis(discrete=TRUE) + scale_color_viridis(discrete=TRUE) + labs(y =NULL, x = NULL) +
  theme(legend.position ="none") + violin_theme + stat_summary(fun.data=mean_sdl, mult=1, geom="point")

v2 <- ggplot(data_subset %>% mutate(screentime_wkdy_1 = factor(screentime_wkdy_1, levels=c("None", "< 30 minutes", ".5", "1 hour", "2 hours", "3 hours", "4+ hours"))),
             aes(x=screentime_wkdy_1, y=nihtbx_fluidcomp_agecorrected, fill=screentime_wkdy_1)) + # color= sets line, fill=sets fill
  geom_violin(trim=TRUE, alpha=0.6, lwd = 0.3) + scale_fill_viridis(discrete=TRUE) + scale_color_viridis(discrete=TRUE) + labs(y =NULL, x = NULL) +
  theme(legend.position ="none") + violin_theme + stat_summary(fun.data=mean_sdl, mult=1, geom="point")

v3 <- ggplot(data_subset %>% mutate(screentime_wkdy_1 = factor(screentime_wkdy_1, levels=c("None", "< 30 minutes", ".5", "1 hour", "2 hours", "3 hours", "4+ hours"))),
             aes(x=screentime_wkdy_1, y=nihtbx_totalcomp_agecorrected, fill=screentime_wkdy_1)) + # color= sets line, fill=sets fill
  geom_violin(trim=TRUE, alpha=0.6, lwd = 0.3) + scale_fill_viridis(discrete=TRUE) + scale_color_viridis(discrete=TRUE) + labs(y =NULL, x = NULL) +
  theme(legend.position ="none") + violin_theme + stat_summary(fun.data=mean_sdl, mult=1, geom="point")

grid.arrange(v1,v2,v3, nrow = 1)
violins <- arrangeGrob(v1,v2,v3, nrow=1) # ggsave won't save a grid.arrange
ggsave("television_cognitive.pdf", plot=violins, path="./plots", dpi=150)

# Television vs. Reading #####

violin_theme <- theme(axis.text.x=element_blank(),
                      axis.ticks.x = element_blank(),)

k1 <- ggplot(data_subset %>% mutate(screentime_wkdy_1 = factor(screentime_wkdy_1, levels=c("None", "< 30 minutes", ".5", "1 hour", "2 hours", "3 hours", "4+ hours"))),
             aes(x=screentime_wkdy_1, y=sports_activity_ss_read_hours_p, fill=screentime_wkdy_1)) + # color= sets line, fill=sets fill
  geom_violin(trim=TRUE, alpha=0.6, lwd = 0.3) + scale_fill_viridis(discrete=TRUE, name="Daily Television") + scale_color_viridis(discrete=TRUE) + labs(y =NULL, x = NULL) +
  theme(legend.position ="bottom") + violin_theme + stat_summary(fun.data=mean_sdl, mult=1, geom="point")
k1
ggsave("television_vs_reading.pdf", plot=k1, path="./plots", dpi=150)
