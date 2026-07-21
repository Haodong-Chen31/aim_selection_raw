library(sysfonts)
library(showtextdb)
library(showtext)
library(tidyverse)
library(ggthemr)
library(patchwork)
library(ggthemes)
library(ggrepel)
library(ggpubr)

setwd("E:/Rworking/popgenetics/aim_1205_more_extract/")

showtext::showtext_auto()
sf::sf_use_s2(FALSE)

input_file <- "tagSNP_union.top61_number_snp.csv"
snp <- read.csv(input_file)
color_set <- c("#801e91","#344fa8","#f7cb34","#a8a8aa","#ffe4c7")
snp$chr <- factor(snp$chrom, levels = c(1:22,'Y','MT'))

# 从输入文件名提取前缀
# 移除路径（如果有）
file_base <- basename(input_file)
# 移除扩展名（.csv或其他）
file_prefix <- sub("\\.[^.]*$", "", file_base)

#bar
p_bar <- ggplot(snp, aes(chr, total, fill = chr)) +
  geom_col(width = 0.84) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = rep(color_set, length.out = length(unique(snp$chr)))) +
  geom_text(aes(x=chr,label=total),vjust = -0.5,size=4)+
  scale_x_discrete(limits = as.character(c(1:22,'Y','MT'))) +
  # Hide color legend
  guides(fill = "none", color = "none")+
  labs(x = "Chromosome",y='count')+
  theme_pubclean()+
  theme(
    #panel.grid = element_blank(),
    plot.title = element_text(size = 32),
    axis.text.y = element_text(size = 14, hjust = 1),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 15, face = 'plain'),  # x轴标题大小
    axis.title.y = element_text(size = 15, face = 'plain')
  )
  
ggsave(paste0(file_prefix, ".pdf"),p_bar,width = 14,height = 10)
