#plots results from norm.wig files made with runCN
#V1 is chromosome, V2 position, V3 CN

library(dplyr)
library(ggplot2)

#hetData <- read.delim("~/Dropbox/Dunham Lab/yEvo/Caffeine/220526_seq/caffeine_CNV/RB_caffeine_cnv_all.txt")
#hetData <- subset(hetData, clone != "3gray1")
#hetData <- subset(hetData, clone != "6yellow1")
hetDataA <- read.delim("~/Documents/Research/dunham/freeze_thaw_evolution/march2023_FTevo/CNV/A/A24_S14_1000bp_norm.txt")

plotA <- ggplot() +
  geom_jitter(data=hetDataA,
              aes(x=pos,#or V2
                  y=cnv, #or V3
                  shape='1',
                  size='0.25',
                  color = ifelse(chr %% 2 == 0, "A", "B"), #or chr is V1
                  fill = ifelse(chr %% 2 == 0, "A", "B")),
              alpha = 0.5, show.legend=F) +
  scale_color_manual(
    name = NULL,
    values=c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  scale_fill_manual(
    name = NULL,
    values = c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  xlab("genome position (1000bp windows)") +
  ylab("copy number") +
  ylim(0,2) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18)) +
  coord_fixed(2000000)

hetDataE <- read.delim("~/Documents/Research/dunham/freeze_thaw_evolution/march2023_FTevo/CNV/E/E15_S22_1000bp_norm.txt")

plotE <- ggplot() +
  geom_jitter(data=hetDataE,
              aes(x=pos,#or V2
                  y=cnv, #or V3
                  shape='1',
                  size='0.25',
                  color = ifelse(chr %% 2 == 0, "A", "B"), #or chr is V1
                  fill = ifelse(chr %% 2 == 0, "A", "B")),
              alpha = 0.5, show.legend=F) +
  scale_color_manual(
    name = NULL,
    values=c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  scale_fill_manual(
    name = NULL,
    values = c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  xlab("genome position (1000bp windows)") +
  ylab("copy number") +
  ylim(0,2) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18)) +
  coord_fixed(2000000)

hetDataG <- read.delim("~/Documents/Research/dunham/freeze_thaw_evolution/march2023_FTevo/CNV/G/G15_S24_1000bp_norm.txt")

plotG <- ggplot() +
  geom_jitter(data=hetDataG,
              aes(x=pos,#or V2
                  y=cnv, #or V3
                  shape='1',
                  size='0.25',
                  color = ifelse(chr %% 2 == 0, "A", "B"), #or chr is V1
                  fill = ifelse(chr %% 2 == 0, "A", "B")),
              alpha = 0.5, show.legend=F) +
  scale_color_manual(
    name = NULL,
    values=c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  scale_fill_manual(
    name = NULL,
    values = c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  xlab("genome position (1000bp windows)") +
  ylab("copy number") +
  ylim(0,2) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18)) +
  coord_fixed(2000000)

CNVdataA <- read.delim("~/Documents/Research/dunham/freeze_thaw_evolution/march2023_FTevo/CNV/A/A24_S14_1000bp_norm.txt")
CNVdataE <- read.delim("~/Documents/Research/dunham/freeze_thaw_evolution/march2023_FTevo/CNV/E/E15_S22_1000bp_norm.txt")
CNVdataG <- read.delim("~/Documents/Research/dunham/freeze_thaw_evolution/march2023_FTevo/CNV/G/G15_S24_1000bp_norm.txt")

CNVdataA$dataset <- "A"
CNVdataE$dataset <- "E"
CNVdataG$dataset <- "G"

combined_data <- bind_rows(CNVdataA, CNVdataE, CNVdataG)

ggplot() +
  facet_wrap(~dataset, dir="v") +
  geom_jitter(data=combined_data,
              aes(x=pos,#or V2
                  y=cnv, #or V3
                  shape='1',
                  size='0.25',
                  color = ifelse(chr %% 2 == 0, "A", "B"), #or chr is V1
                  fill = ifelse(chr %% 2 == 0, "A", "B")),
              alpha = 0.5, show.legend=F) +
  scale_color_manual(
    name = NULL,
    values=c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  scale_fill_manual(
    name = NULL,
    values = c("A"='#33b4ff',"B"='#0072b2',"C"='#e69f00',"D"='#ffbe33')) +
  xlab("genome position (1000bp windows)") +
  ylab("copy number") +
  ylim(0,2) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18), plot.margin = margin(2, 10, 4, 4, "mm")) +
  coord_fixed(2000000)

