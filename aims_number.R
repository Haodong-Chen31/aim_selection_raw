#!/usr/bin/R

library(ggplot2)
library(ggthemes)

arg=commandArgs(T)
aim_number <- arg[1] #top的aim文件
show_num <- as.integer(arg[2])

#统计AIMs数
aims_num <- read.table(aim_number,col.names = c('top','na','num'))
aims_num_p <- head(aims_num,show_num)
aims_num_p$top_num <- as.numeric(gsub("top", "", aims_num_p$top))
aims_num_p <- aims_num_p[order(aims_num_p$top_num), ]
aims_num_p$top <- factor(aims_num_p$top, levels = aims_num_p$top)

# 定义要显示的标签
# 创建序列并添加"top"前缀
label_seq <- paste0("top", seq(0, show_num, 50))
# 将top0替换为top1
label_seq[label_seq == "top0"] <- "top1"

p2 <- ggplot(aims_num_p) +
  geom_line(aes(x = top, y = num, group = 1), color = "#577590", linewidth = 0.6) +
  geom_point(aes(x = top, y = num), color = "#577590", size = 0.7) +
  geom_hline(yintercept = 0) +
  # 设置Y轴范围为[0,30000]，并调整次要坐标轴
  scale_y_continuous(
    limits = c(0, 9000),
    breaks = seq(0, 9000, by = 1000),
  ) +
  scale_x_discrete(
    breaks = label_seq)+
  labs(x = "Top SNPs", y='the number of AIMs') +
  theme_hc() +
  theme(
    axis.text.x = element_text(size = 12, angle = 300,vjust = 1,hjust = 0),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10))  # top/right/bottom/left，单位 pt
ggsave('aim_number.pdf',p2,width = 19,height = 10)
