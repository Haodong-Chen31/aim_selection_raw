library(ggthemes)
library(ggplot2)
library(dplyr)

aac_mcc <- read.csv("aac_mcc_top1to1000.csv")

# 添加 top_num 和 factor
aac_mcc_p <- aac_mcc %>%
  mutate(
    top_num = seq(1, n()),
    top = factor(paste0("top", top_num), levels = paste0("top", top_num))
  )

# 定义要显示的X轴标签
label_seq <- paste0("top", seq(0, 1000, 100))
label_seq[label_seq == "top0"] <- "top1"

# reshape 数据，让 AAC 和 MCC 在同一列，用一个变量区分
aac_mcc_long <- aac_mcc_p %>%
  select(top, aac, mcc) %>%
  tidyr::pivot_longer(cols = c("aac", "mcc"), names_to = "metric", values_to = "value")

# 添加线型信息
aac_mcc_long <- aac_mcc_long %>%
  mutate(
    line_type = ifelse(metric == "aac", "solid", "dashed"),
    point_shape = ifelse(metric == "aac", 16, 17)  # 实心圆 / 三角
  )

# 为x添加索引位置
x_pos_targ <- which(levels(aac_mcc_long$top) == "top539")

# 绘图
p2 <- ggplot(aac_mcc_long, aes(x = top, y = value, color = metric, linetype = metric, shape = metric, group = metric)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.4) +
  geom_hline(yintercept = 0, color = "black") +
  #geom_hline(yintercept = 0.98, linetype = "dashed", color = "grey") + # y上加一条阈值虚线
  # 在y左侧为阈值添加标签
  #annotate(
  #  "text",
  #  x = -Inf, y = 0.98,
  #  label = "0.98",
  #  hjust = -0.2,
  #  vjust = -0.5,
  #  color = "grey40",
  #  size = 3
  #) +
  geom_vline(xintercept = x_pos_targ,linetype = "dashed",color = "grey") + # x上加一条阈值虚线
  annotate(
    "text",
    x = x_pos_targ,
    y = -Inf,
    label = "top539",
    angle = 90,
    vjust = 1.2,
    hjust = 1.1,
    color = "grey40",
    size = 3
  ) + # x上阈值标签
  scale_color_manual(
    name = "Metric",
    values = c("aac" = "#577590", "mcc" = "#F94144"),
    labels = c("AAC", "MCC")
  ) +
  scale_linetype_manual(
    name = "Metric",
    values = c("aac" = "solid", "mcc" = "dashed"),
    labels = c("AAC", "MCC")
  ) +
  scale_shape_manual(
    name = "Metric",
    values = c("aac" = 16, "mcc" = 17),
    labels = c("AAC", "MCC")
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_discrete(breaks = label_seq) +
  labs(x = "Top SNPs", y = "", title = "AAC and MCC values") +
  theme_hc() +
  theme(
    axis.text.x = element_text(size = 10, angle = 300, vjust = 1, hjust = 0),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)  # top/right/bottom/left，单位 pt
  )+
  coord_cartesian(clip = "off")

ggsave("aac_mcc_top1to1000_more0.98.pdf", p2, bg = "white", width = 10, height = 7)

