library(ggthemes)
library(ggplot2)
library(dplyr)

# 1. 读取数据
aac_mcc <- read.csv("aac_mcc_top1to1000.csv")

# 2. 处理数据：只保留 aac 相关信息
aac_mcc_p <- aac_mcc %>%
  mutate(
    top_num = seq(1, n()),
    top = factor(paste0("top", top_num), levels = paste0("top", top_num))
  )

# 3. 定义要显示的 X 轴标签
label_seq <- paste0("top", seq(0, 1000, 100))
label_seq[label_seq == "top0"] <- "top1"

# 4. 绘图：只绘制 aac 曲线
p2 <- ggplot(aac_mcc_p, aes(x = top, y = aac, group = 1)) +
  # 绘制折线
  geom_line(color = "#577590", linewidth = 0.6) +
  # 绘制散点
  geom_point(color = "#577590", size = 0.8, shape = 16) +
  
  # 设置坐标轴
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    expand = c(0, 0.02)
  ) +
  scale_x_discrete(breaks = label_seq) +
  
  # 添加参考线（例如 0.95 或 0.98）
  #geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey60") +
  #annotate(
  #  "text", 
  #  x = 1, y = 0.96, 
  #  label = "Threshold: 0.95", 
  #  hjust = 0, color = "grey40", size = 3
  #) +
  
  # 主题设置
  theme_hc() +
  labs(
    title = "Average Accuracy for Top SNPs",
    x = "Number of SNPs (Top N)",
    y = "Average Accuracy (AAC)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90")
  )

# 5. 保存图片
ggsave("nb_aac_line_top1to1000.pdf", p2, width = 10, height = 6)

# 打印结果查看
print(p2)
