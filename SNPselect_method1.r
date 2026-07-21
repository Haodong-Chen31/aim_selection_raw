# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 设置默认值
input_file <- "combined_fst_test.txt"
mean_threshold <- 0.1
cv_threshold <- 0.5

# 如果提供了参数，则使用参数值
if (length(args) >= 1) {
  input_file <- args[1]
}
if (length(args) >= 2) {
  mean_threshold <- as.numeric(args[2])
}
if (length(args) >= 3) {
  cv_threshold <- as.numeric(args[3])
}

# 检查输入文件是否存在
if (!file.exists(input_file)) {
  stop("错误: 输入文件不存在 - ", input_file)
}

cat("使用的参数:\n")
cat("输入文件 =", input_file, "\n")
cat("mean_threshold =", mean_threshold, "\n")
cat("cv_threshold =", cv_threshold, "\n")

library(qqman)

# 读取数据
data <- read.table(input_file, header=T)

# 生成输出文件的基础名称（去掉路径和扩展名）
base_name <- tools::file_path_sans_ext(basename(input_file))

# 将FST列中小于0的值替换为0
fst_cols <- 4:ncol(data)
data[, fst_cols] <- apply(data[, fst_cols], 2, function(x) {
  x[x < 0] <- 0
  return(x)
})

#过滤掉全NA行
fst_cols <- 4:ncol(data)
na_count <- rowSums(is.na(data[, fst_cols]))
data_clean <- data[na_count < length(fst_cols), ]

# 计算均值和变异系数
mean_fst <- rowMeans(data_clean[, 4:ncol(data_clean)], na.rm = TRUE)
sd_fst <- apply(data_clean[, 4:ncol(data_clean)], 1, sd, na.rm = TRUE)
cv_fst <- sd_fst / mean_fst
data_clean$Mean_FST <- mean_fst
data_clean$SD_FST <- sd_fst
data_clean$CV_FST <- cv_fst

# 输出Mean_FST的范围
mean_range <- range(data_clean$Mean_FST, na.rm = TRUE)
cat(paste("Mean_FST范围: [", round(mean_range[1], 4), ", ", round(mean_range[2], 4), "]\n", sep = ""))

# 统计超过阈值的SNP数量
snps_above_threshold <- sum(data_clean$Mean_FST > mean_threshold, na.rm = TRUE)
total_snps <- nrow(data_clean)
percentage_above <- round(snps_above_threshold / total_snps * 100, 2)

cat(paste("总SNP数量: ", total_snps, "\n", sep = ""))
cat(paste("Mean_FST > ", mean_threshold, " 的SNP数量: ", snps_above_threshold, 
          " (占 ", percentage_above, "%)\n", sep = ""))

# 筛选候选位点
candidate_snps <- data_clean[data_clean$Mean_FST > mean_threshold & data_clean$CV_FST > cv_threshold, ]

#manhattan_plot
color_set <- c("#801e91","#344fa8","#f7cb34","#a8a8aa","#ffe4c7")
pdf(paste0(base_name, '_manhattan_mean', mean_threshold, '_cv', cv_threshold, '.pdf'), width =12,height =8)
manhattan(data_clean, col=color_set, logp = F, ylab='Fst', abline(lty = 2 , h = mean_threshold),
          chr="CHROM", bp="POS", snp="SNPID", p="Mean_FST", ylim=c(0,1), annotateTop = FALSE)
dev.off()

# 输出结果
write.table(candidate_snps$SNPID, 
            paste0(base_name, '_candidate_snps_mean', mean_threshold, '_cv', cv_threshold, '.snpid'), 
            row.names=F, col.names = FALSE, quote=F, sep="\t")

# 创建颜色列，根据阈值条件标记
data_clean$color <- ifelse(data_clean$Mean_FST > mean_threshold & data_clean$CV_FST > cv_threshold, 
                           "blue", "gray")

# 创建散点图，查看阈值范围内的位点
pdf(paste0(base_name, '_meanfst_cv_dot_mean', mean_threshold, '_cv', cv_threshold, '.pdf'), width =12,height =8)
plot(data_clean$Mean_FST, data_clean$CV_FST,
     col = data_clean$color,
     pch = 16,  # 实心圆点
     cex = 1, # 点的大小
     xlab = "Mean FST",
     ylab = "Coefficient of Variation (CV) of FST",
     main = paste0("FST Mean vs CV (Mean > ", mean_threshold, ", CV > ", cv_threshold, ")"))

# 添加阈值参考线
abline(v = mean_threshold, lty = 2, col = "black")  # 垂直虚线
abline(h = cv_threshold, lty = 2, col = "black")    # 水平虚线

# 添加图例
legend("topright", 
       legend = c(paste("High Mean FST & CV SNPs (n =", nrow(candidate_snps), ")"), "Other SNPs"),
       col = c("blue", "gray"), 
       pch = 16,
       bty = "n")
dev.off()

cat("分析完成!\n")
cat("找到", nrow(candidate_snps), "个候选SNP\n")
cat("输出文件:\n")
cat("- Manhattan图:", paste0(base_name, '_manhattan_mean', mean_threshold, '_cv', cv_threshold, '.pdf'), "\n")
cat("- 散点图:", paste0(base_name, '_meanfst_cv_dot_mean', mean_threshold, '_cv', cv_threshold, '.pdf'), "\n")
cat("- 候选SNP列表:", paste0(base_name, '_candidate_snps_mean', mean_threshold, '_cv', cv_threshold, '.snpid'), "\n")
