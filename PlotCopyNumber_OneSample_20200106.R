# Adapted by Zilong Zeng

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Args
# Args[1] = WorkDirectory
# Args[2] = Sample_name
# Args[3] = Size

library(ggplot2)
library(colorspace)

pdf(paste(args[1],'/',args[2],'/',args[2],'_CopyNumber_', args[3],'bp.pdf',sep=""),height=5, width=18)

hetData <- read.table(paste(args[1],'/',args[2],'/CNV_new_', args[3], 'bp/',args[2],'_', args[3],'bp_norm.wig',sep=""), sep='\t', head=FALSE)
#cutoff <- data.frame( x = c(-Inf, Inf), y = 4, cutoff = factor(4) )

hetPlot <- ggplot() + 
  geom_jitter(data=hetData,
              aes(x=V2,
                  y=V3,
                  shape='1',
                  size=1,                  
                  color = ifelse(V1 %% 2 == 0, "A", "B"),
                  fill = ifelse(V1 %% 2 == 0, "A", "B")),
                  alpha = 0.5) +
  scale_color_manual(
    name = NULL,
    values=darken(c("A"='#33b4ff',
                    "B"='#0072b2',
                    "C"='#e69f00',
                    "D"='#ffbe33')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("A"='#33b4ff',
               "B"='#0072b2',
               "C"='#e69f00',
               "D"='#ffbe33')) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 230218, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 1043402, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 1360022, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 2891955, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 3468829, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 3738990, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 4829930, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 5392573, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 5832461, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 6578212, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 7245028, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 8323205, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 9247636, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 10031969, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 11123260, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 12071326, linetype = "longdash", color = "grey80") +
  scale_shape(solid = TRUE)+
  scale_y_continuous(name = "Copy Number",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
                     limits = c(0, 7.5)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0, 0),
                     limits = c(-100000, 12171326),
                     breaks = c(-100000,
                                115109,
                                636810,
                                1201712,
                                2125988.5,
                                3180392,
                                3603909.5,
                                4284460,
                                5111251.5,
                                5612517,
                                6205336.5,
                                6911620,
                                7784116.5,
                                8785420.5,
                                9639802.5,
                                10577614.5,
                                11597293),
                     labels = c("Chr.",
                                "I",
                                "II",
                                "III",
                                "IV",
                                "V",
                                "VI",
                                "VII",
                                "VIII",
                                "IX",
                                "X",
                                "XI",
                                "XII", 
                                "XIII",
                                "XIV", 
                                "XV", 
                                "XVI")) +
  theme_bw(base_line_size = 1) +
  theme(line = element_line(color="grey80", linewidth=1),
        axis.line = element_line(color="grey80", linewidth=1),
        axis.line.x.bottom = element_line(color="grey80", linewidth=1),
        axis.line.y.left = element_line(color="grey80", linewidth=1),
        axis.line.y.right = element_line(color="grey80", linewidth=1),
        axis.line.x.top = element_line(color="grey80", linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", linewidth = 1),
        axis.ticks.x = element_blank(),
        text = element_text(color = "black", size = 20),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

print(hetPlot)

dev.off()


