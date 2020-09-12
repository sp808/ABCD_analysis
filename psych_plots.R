## Housekeeping #########
library(dplyr)
library(ggplot2)
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

# Import and subset data ##########----------
setwd("~/Projects/DEAPext")
data <- readRDS("nda2.0.2.Rds")
data_subset <- select(data, src_subject_id,eventname,prodrom_psych_ss_number,prodrom_psych_ss_severity_score,
                      nihtbx_cryst_agecorrected,nihtbx_fluidcomp_agecorrected,nihtbx_totalcomp_agecorrected,mrif_score,fsqc_qc)

# Filters
data_subset <- filter(data_subset, eventname == "baseline_year_1_arm_1")
data_subset <- filter(data_subset, mrif_score == "No abnormal findings" | mrif_score == "Normal anatomical variant of no clinical significance")
data_subset <- filter(data_subset, fsqc_qc == "accept")

### Double raincloud ############################# 
# define theme
raincloud_theme = theme(
  text = element_text(size = 10, family = "Helvetica"),
  axis.title.y=element_blank(),
  axis.title.x = element_text(size = 10),
  axis.text = element_text(size = 10),
  axis.text.y=element_blank(),
  axis.ticks.y = element_blank(),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

# Calculate summary statistics for box plots 
x=data_subset$prodrom_psych_ss_number
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

# Plot
source( './plots/geom_flat_violin.R' )
g <- ggplot(data = data_subset, aes(y = prodrom_psych_ss_number, x ="", fill ="")) +
  geom_flat_violin(position = position_nudge(x = .55, y = 0), alpha = 0.5, width = 3) +
  geom_point(aes(y = prodrom_psych_ss_number, x = "", color=""), position = position_jitter(width = .5), size = .005, alpha = 0.3) +
  geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.5) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  labs(y =NULL, x = NULL) +
  scale_fill_manual(values = c("snow4")) +
  scale_colour_manual(values = c("snow4")) +
  coord_flip() +
  theme_bw() +
  raincloud_theme

####### Raincloud severity score and subplot
x=data_subset$prodrom_psych_ss_severity_score
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

# Plot
h <- ggplot(data = data_subset, aes(y = prodrom_psych_ss_severity_score, x ="", fill ="")) +
  geom_flat_violin(position = position_nudge(x = .55, y = 0), alpha = 0.6, width = 3, lwd=0.5) +
  geom_point(aes(y = prodrom_psych_ss_severity_score, x = "", color=""), position = position_jitter(width = .5), size = .005, alpha = 0.4) +
  geom_boxplot(width = .15, outlier.shape = NA, alpha = 0.5, lwd=0.5) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  labs(y =NULL, x = NULL) +
  scale_fill_manual(values = c("coral4")) +
  scale_colour_manual(values = c("coral4")) +
  coord_flip() +
  theme_bw() +
  raincloud_theme

grid.arrange(g,h, nrow = 2)
j <- arrangeGrob(g,h, nrow=2) # ggsave won't save a grid.arrange
ggsave("psych_rainclouds.pdf", plot=j, path="./plots", dpi=150)

########### Cognition scatter plots ######
scatter_theme = theme(
  text = element_text(size = 12, family = "Helvetica"),
  axis.title.y=element_blank(),
  axis.title.x = element_text(size = 12),
  axis.text = element_text(size = 12),
  axis.text.y=element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

a <- ggplot(data_subset, aes(x=prodrom_psych_ss_number, y=nihtbx_cryst_agecorrected, color=prodrom_psych_ss_severity_score)) + 
  geom_point(shape=16, size = 3) + scale_color_continuous(type='viridis') + labs(y =NULL, x = NULL) + 
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, lwd=1, color="orchid", linetype="dashed") + theme(legend.position = "none")

b <- ggplot(data_subset, aes(x=prodrom_psych_ss_number, y=nihtbx_fluidcomp_agecorrected, color=prodrom_psych_ss_severity_score)) + 
  geom_point(shape=16, size = 3) + scale_color_continuous(type='viridis') + labs(y =NULL, x = NULL) + theme(legend.position = "none") +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, color="orchid", linetype="dashed")

c <- ggplot(data_subset, aes(x=prodrom_psych_ss_number, y=nihtbx_totalcomp_agecorrected, color=prodrom_psych_ss_severity_score)) + 
  geom_point(shape=16, size = 3) + scale_color_continuous(type='viridis', name = "Severity\nScore") + labs(y =NULL, x = NULL) + 
  theme(legend.justification = c("top"), legend.margin = margin(0, 0, 0, 0), legend.position = c(0.85, 0.85)) + 
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, color="orchid", linetype="dashed")

grid.arrange(a,b,c, nrow = 1)
v <- arrangeGrob(a,b,c, nrow=1) # ggsave won't save a grid.arrange
ggsave("psych_cognitive.pdf", plot=v, path="./plots", dpi=150)

