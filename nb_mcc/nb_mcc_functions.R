#!/usr/bin/env Rscript

##########################
## 0. 服务器安全设置
##########################
Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)

##########################
## 1. 加载包
##########################
suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(doRNG)

  library(caret)
  library(naivebayes)
  library(ggplot2)
  library(reshape2)
  library(genio)
})
source('applied_maM.R')

##########################
## 0. 处理命令行参数
##########################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  # 如果从命令行接收到参数
  n_from <- as.numeric(args[1])
  n_to   <- as.numeric(args[2])
} else {
  # 默认值（如果没有传参）
  n_from <- 1
  n_to   <- 1000
  message("当前将使用默认设置：n_from = 1, n_to = 1000\n")
}

# 动态生成文件名
output_filename <- sprintf("acc_mcc_top%dto%d.csv", n_from, n_to)

# 封装函数：运行完整的朴素贝叶斯分析流程
run_naive_bayes_analysis <- function(genotype_data_path, 
                                     n_folds = 5, 
                                     n_repeats = 10,
                                     conf_level = 0.95) {
  
  # 加载必要的包
  suppressPackageStartupMessages({
    library(caret)
    library(naivebayes)
    library(ggplot2)
    library(reshape2)
    library(genio)
  })
  
  # 1. 读取PLINK数据
  cat("步骤1: 读取PLINK数据...\n")
  plink_data <- read_plink(genotype_data_path)
  fn <- basename(genotype_data_path)
  
  # 从plink_data中提取基因型和人群标签
  genotypes_from_plink <- t(as.data.frame(plink_data$X))
  population_labels_from_plink <- plink_data$fam$fam
  
  # 组合成数据框
  df_actual_plink <- data.frame(
    population = population_labels_from_plink, 
    genotypes_from_plink,
    check.names = FALSE
  )
  
  # 2. 数据预处理
  cat("步骤2: 数据预处理...\n")
  
  # 保留至少能满足交叉验证的人群样本数
  pop_counts <- table(df_actual_plink$population)
  keep_pops <- names(pop_counts[pop_counts >= n_folds])
  df_filtered <- df_actual_plink[df_actual_plink$population %in% keep_pops, ]
  
  # 去除掉有NA的列(即位点)
  col_with_na <- which(apply(df_filtered, 2, anyNA))
  
  # 只有当有包含NA的列时才进行过滤
  if (length(col_with_na) > 0) {
    df_filtered <- df_filtered[, -col_with_na]
    cat(sprintf("已移除 %d 个包含NA值的列(位点)\n", length(col_with_na)))
  } else {
    cat("未发现包含NA值的列(位点)\n")
  }
  
  # 转换数据类型
  df_filtered_factor <- df_filtered
  df_filtered_factor$population <- as.factor(df_filtered_factor$population)
  
  # 将基因型列转换为因子
  # 没有ID列，从第2列开始
  df_filtered_factor[, -1] <- lapply(df_filtered_factor[, -1], factor)
  
  # 移除零方差/近零方差特征
  cat("步骤3: 移除零方差/近零方差特征...\n")
  
  # 确定population列的位置
  population_col <- which(names(df_filtered_factor) == "population")
  
  # 使用nearZeroVar
  if (population_col == 1) {
    zv_nzv_indices <- nearZeroVar(df_filtered_factor[, -1], saveMetrics = FALSE)
    if (length(zv_nzv_indices) > 0) {
      df_filtered_factor_cleaned_by_nzv <- df_filtered_factor[, -(zv_nzv_indices + 1)]
      cat(sprintf("已移除 %d 个零方差/近零方差特征\n", length(zv_nzv_indices)))
    } else {
      df_filtered_factor_cleaned_by_nzv <- df_filtered_factor
      cat("未发现零方差/近零方差特征\n")
    }
  } else {
    # 如果population不是第一列，需要特殊处理
    zv_nzv_indices <- nearZeroVar(df_filtered_factor[, -population_col], saveMetrics = FALSE)
    if (length(zv_nzv_indices) > 0) {
      # 调整索引，跳过population列
      df_filtered_factor_cleaned_by_nzv <- df_filtered_factor[, -(zv_nzv_indices + population_col)]
      cat(sprintf("已移除 %d 个零方差/近零方差特征\n", length(zv_nzv_indices)))
    } else {
      df_filtered_factor_cleaned_by_nzv <- df_filtered_factor
      cat("未发现零方差/近零方差特征\n")
    }
  }
  
  # 3. 训练朴素贝叶斯模型
  # 朴素贝叶斯 + 分层 K 折交叉验证 + AC 评估实现
  cat("步骤4: 训练朴素贝叶斯模型...\n")
  
  set.seed(43)
  
  # 定义交叉验证参数
  # classProbs = TRUE 确保模型输出类别概率，虽然这里不直接用，但有时对调试有用
  train_control_nb <- trainControl(
    method = "repeatedcv",
    number = n_folds,
    repeats = n_repeats,
    classProbs = TRUE,
    savePredictions = "final",
    allowParallel = FALSE
  )
  
  # 训练模型
  # m折交叉验证 × 重复n次 = m × n次独立的训练/测试循环
  x_data <- df_filtered_factor_cleaned_by_nzv[, -population_col] # 提取特征
  y_data <- df_filtered_factor_cleaned_by_nzv[[population_col]] # 提取标签
  model_nb_cv <- train(
    x = x_data,
    y = y_data,
    #population ~ .,
    #data = df_filtered_factor_cleaned_by_nzv,
    method = "naive_bayes",
    trControl = train_control_nb,
    metric = "Accuracy"
  )
  
  # 4. 计算性能指标
  cat("步骤5: 计算性能指标...\n")
  print("朴素贝叶斯模型训练结果:")
  print(model_nb_cv)
  
  # 准确率相关统计
  accuracy_values <- model_nb_cv$resample$Accuracy
  mean_accuracy_values <- mean(accuracy_values)
  
  # 准确率的95%置信区间 增加对恒定数据的处理
  if (length(unique(accuracy_values)) > 1) {
    # 如果数据有波动，正常运行 t.test
    accuracy_values_ci <- t.test(accuracy_values)$conf.int
    accuracy_ci_lower <- accuracy_values_ci[1]
    accuracy_ci_upper <- accuracy_values_ci[2]
  } else {
    # 如果所有值都相同（例如全是 1），CI 就等于该值本身
    cat("警告：所有重采样准确率完全一致，手动设置置信区间。\n")
    accuracy_ci_lower <- mean_accuracy_values
    accuracy_ci_upper <- mean_accuracy_values
  }

  # 变异系数
  cv_percentage <- sd(accuracy_values) / mean(accuracy_values) * 100

  # 5. 计算宏观平均MCC
  cat("步骤6: 计算宏观平均MCC...\n")
  
  # 创建混淆矩阵的函数
  create_confusion_matrix <- function(actual, predicted) {
    all_classes <- sort(union(unique(actual), unique(predicted)))
    n <- length(all_classes)
    cm <- matrix(0, nrow = n, ncol = n,
                 dimnames = list(Predicted = all_classes, 
                                 Actual = all_classes))
    
    for(i in 1:length(actual)) {
      pred_idx <- which(all_classes == predicted[i])
      act_idx <- which(all_classes == actual[i])
      cm[pred_idx, act_idx] <- cm[pred_idx, act_idx] + 1
    }
    
    return(cm)
  }
  
  # 获取预测结果
  all_predictions_nb_mcc <- model_nb_cv$pred
  
  # 构建概率矩阵
  confusion_matrix <- create_confusion_matrix(
    all_predictions_nb_mcc$obs, 
    all_predictions_nb_mcc$pred
  )
  prob_matrix <- confusion_matrix / sum(confusion_matrix)
  
  # 计算宏观MCC和置信区间
  #source('applied_maM.R')
  macro_mcc_result <- compute_macroMcc_CI(prob_matrix, nrow(df_filtered_factor_cleaned_by_nzv), conf_level)
  
  # 6. 可视化概率矩阵
  cat("步骤7: 生成可视化...\n")
  
  # 可视化函数
  visualize_prob_matrix <- function(prob_matrix, title = "Matrix of probabilities") {
    # 转换为长格式
    prob_melted <- melt(prob_matrix)
    colnames(prob_melted) <- c("Predicted", "Actual", "Probability")
    
    # 创建热图
    p <- ggplot(prob_melted, aes(x = Actual, y = Predicted, fill = Probability)) +
      geom_tile(color = "white", size = 1) +
      geom_text(aes(label = round(Probability, 3)),
                color = "black", size = 1.4) +
      scale_fill_gradient2(low = "white", high = "steelblue",
                           midpoint = max(prob_melted$Probability)/2) +
      labs(title = title,
           x = "Actual",
           y = "Predicted") +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
	    axis.text.y = element_text(size = 5),
            plot.title = element_text(hjust = 0.5, size = 10),
      	    # 2. 设置 X 和 Y 轴标题大小 (默认 ~11)
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            
            # 3. 设置图例标题大小 (默认 ~11)
            legend.title = element_text(size = 8),
            # 如果图例里的数字(刻度)也想变小，可以加下面这行
            legend.text = element_text(size = 6))
    
    # 对角线分析
    diag_values <- diag(prob_matrix)
    cat("\n对角线值（正确预测概率）:\n")
    diag_summary <- data.frame(
      类别 = names(diag_values),
      正确率 = round(diag_values, 4)
    )
    print(diag_summary)
    cat(sprintf("平均正确率: %.4f\n", mean(diag_values)))
    
    return(p)
  }
  
  # 生成可视化
  prob_matrix_plot <- visualize_prob_matrix(prob_matrix, "Matrix of probabilities")
  
  # 7. 汇总结果
  cat("步骤8: 汇总分析结果...\n")
  
  # 删除函数内不再需要的超大型中间变量
  rm(plink_data, df_actual_plink, genotypes_from_plink, x_data)

  # 强制执行垃圾回收
  gc()

  # 数据概览
  data_summary <- list(
    sample_num = nrow(df_filtered_factor_cleaned_by_nzv),
    pop_num = length(unique(df_filtered_factor_cleaned_by_nzv$population)),
    aim_num = ncol(df_filtered_factor_cleaned_by_nzv) - 1,
    cv_set = sprintf("%d CV with %d repeats", n_folds, n_repeats)
  )
  
  # 模型性能
  model_performance <- list(
    mean_accuracy = mean_accuracy_values,
    accuracy_ci_lower = accuracy_ci_lower,
    accuracy_ci_upper = accuracy_ci_upper,
    accuracy_sd = sd(accuracy_values),
    cv_percentage = cv_percentage,
    macro_mcc = macro_mcc_result["maM"],
    mcc_ci_lower = macro_mcc_result["confidence_interval1"],
    mcc_ci_upper = macro_mcc_result["confidence_interval2"]
  )
  
  # 模型详细信息
  model_details <- list(
    optimal_parameters = model_nb_cv$bestTune,
    resample_times = nrow(model_nb_cv$resample)
  )
  
  # 8. 创建最终输出列表
  results <- list(
    # 数据信息
    data_summary = data_summary,
    
    # 模型性能
    model_performance = model_performance,
    
    # 详细结果
    model_details = model_details,
    
    # 原始对象
    #trained_model = model_nb_cv,
    probability_matrix = prob_matrix,
    
    # 可视化
    plot_prob_matrix = prob_matrix_plot
    
    # 原始数据
    #processed_data = df_filtered_factor_cleaned_by_nzv,
    
    # 完整预测结果
    #all_predictions = all_predictions_nb_mcc
  )
  
  # 9. 打印总结报告
  cat("\n=================== 分析完成 ===================\n")
  cat("数据概览:\n")
  cat(sprintf("  总样本数: %d\n", data_summary$sample_num))
  cat(sprintf("  群体数: %d\n", data_summary$pop_num))
  cat(sprintf("  SNP特征数: %d\n", data_summary$aim_num))
  cat(sprintf("  交叉验证: %s\n", data_summary$cv_set))
  
  cat("\n模型性能:\n")
  cat(sprintf("  平均准确率: %.4f\n", model_performance$mean_accuracy))
  cat(sprintf("  准确率95%% CI: [%.4f, %.4f]\n", 
              model_performance$accuracy_ci_lower, model_performance$accuracy_ci_upper))
  cat(sprintf("  准确率标准差: %.4f\n", model_performance$accuracy_sd))
  cat(sprintf("  变异系数: %.2f%%\n", model_performance$cv_percentage))
  cat(sprintf("  宏观MCC: %.4f\n", model_performance$macro_mcc))
  cat(sprintf("  MCC 95%% CI: [%.4f, %.4f]\n", 
              model_performance$mcc_ci_lower, model_performance$mcc_ci_upper))
  
  cat(sprintf("\n模型稳定性: %s\n", 
              ifelse(cv_percentage < 10, "稳定", "需关注")))
  cat("================================================\n")
  
  return(results)
}

# 使用示例函数
print_results_summary <- function(results, fn, n_folds, n_repeats) {
  cat("\n============== 结果摘要：", fn, n_folds, "折", n_repeats,"次重复 ==============\n")
  
  # 从results中提取关键信息
  perf <- results$model_performance
  
  cat(sprintf("重采样结果的平均准确率: %.4f (95%% CI: [%.4f, %.4f])\n",
              perf$mean_accuracy, perf$accuracy_ci_lower, perf$accuracy_ci_upper))
  cat(sprintf("宏观MCC: %.4f (95%% CI: [%.4f, %.4f])\n",
              perf$macro_mcc, perf$mcc_ci_lower, perf$mcc_ci_upper))
  cat(sprintf("变异系数: %.2f%%\n", perf$cv_percentage))
  
  cat("\n")
  cat("朴素贝叶斯分类器在测试数据上获得了")
  cat(sprintf(" %.1f%% (95%% CI: %.1f%%-%.1f%%)", 
              perf$mean_accuracy*100, perf$accuracy_ci_lower*100, perf$accuracy_ci_upper*100))
  cat("的准确率，宏观MCC为")
  cat(sprintf(" %.3f (95%% CI: %.3f-%.3f)。\n",
              perf$macro_mcc, perf$mcc_ci_lower, perf$mcc_ci_upper))
  
  cat("=====================================================\n")
}

# 快速使用函数（一键运行所有分析）
quick_analysis_nb_mcc <- function(genotype_data_path, 
                                  n_folds=10, 
                                  n_repeats=10, 
                                  conf_level =0.95, 
                                  output_plot = TRUE) {
  cat("开始朴素贝叶斯分类分析...\n")
  fn <- basename(genotype_data_path)
  # 运行完整分析
  results <- run_naive_bayes_analysis(genotype_data_path, n_folds, n_repeats, conf_level = conf_level)
  
  # 打印结果摘要
  print_results_summary(results, fn, n_folds, n_repeats)
  
  # 显示可视化图形
  if (output_plot) {
    print(results$plot)
  }
  
  # 返回完整结果
  return(results)
}

# 使用示例：
# results <- quick_analysis("test_random100snp")
# 
# 访问特定结果：
# results$model_performance$平均准确率
# results$model_performance$宏观MCC
# results$plot


##########################
## 3. 并行设置
##########################
n_cores <- min(100, parallel::detectCores() - 1)
cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

cat("========================================\n")
cat("Linux parallel NB-MCC analysis started\n")
cat("Cores used:", n_cores, "\n")
cat("========================================\n")

##########################
## 4. 并行运行分析
##########################
results_list <- foreach(
  n = seq(n_from, n_to, 1),
  .packages = c("caret","naivebayes","ggplot2","reshape2","genio"),
  .export = c(
    "quick_analysis_nb_mcc",
    "run_naive_bayes_analysis",
    "print_results_summary"
  ),
  #.export = NULL,
  .options.RNG = 20240126
) %dorng% {
  options(expressions = 100000)  # 关键：在子进程中再次提高限制（双重保险）

  sink(sprintf("log_top%d.txt", n))   # 重定向输出
  res <- quick_analysis_nb_mcc(
    paste0("tagSNP_union/tagSNP_union.top", n),
    n_folds = 10,
    n_repeats = 10,
    conf_level = 0.95,
    output_plot = FALSE
  )
  sink()  # 结束重定向

  list(
    n   = n,
    aac = res$model_performance$mean_accuracy,
    mcc = res$model_performance$macro_mcc,
    res = res
  )
}
names(results_list) <- paste0("top", sapply(results_list, `[[`, "n"))

##########################
## 5. 汇总结果
##########################
all_nb_mcc_results <- vector("list", max(1:100))
aac_mcc <- data.frame(
  top = integer(),
  aac = numeric(),
  mcc = numeric()
)

for (x in results_list) {
  all_nb_mcc_results[[x$n]] <- x$res
  aac_mcc <- rbind(
    aac_mcc,
    data.frame(top = x$n, aac = x$aac, mcc = x$mcc)
  )
}

aac_mcc <- aac_mcc[order(aac_mcc$top), ]

write.csv(
  aac_mcc,
  file = output_filename,
  row.names = FALSE
)

##########################
## 6. 清理
##########################
stopCluster(cl)
registerDoSEQ()

cat("========================================\n")
cat("Analysis finished successfully\n")
cat(sprintf("Output: %s\n", output_filename))
cat("========================================\n")






#批量运行
#all_nb_mcc_results <- list()
#aac_mcc <- data.frame(aac = numeric(),mcc = numeric())
#for (n in seq(101,200,1)) {
#  all_nb_mcc_results[[n]] <-
#    quick_analysis_nb_mcc(paste0('tagSNP_union/tagSNP_union.top',n),
#                          n_folds = 10,n_repeats = 10,conf_level = 0.95,output_plot = F)
#  # 统计值数据
#  new_row <- data.frame(
#    aac = all_nb_mcc_results[[n]]$model_performance$mean_accuracy,
#    mcc = all_nb_mcc_results[[n]]$model_performance$macro_mcc
#  )
#  aac_mcc <- rbind(aac_mcc,new_row)
#}
#结果保存csv文件
#write.csv(aac_mcc,'acc_mcc_top1to100.csv')

####################################### 概率矩阵图 ############################################
library(patchwork)

plot_prob_matrix_list <- lapply(
  results_list,
  function(x) x$res$plot_prob_matrix
)

# 假设 plot_prob_matrix_list 已经存在，有 10,000 个 ggplot 元素
# 每 100 个一组
plots_per_page <- 100
num_pages <- ceiling(length(plot_prob_matrix_list) / plots_per_page)

# 将 list 分组
plot_chunks <- split(plot_prob_matrix_list, ceiling(seq_along(plot_prob_matrix_list)/plots_per_page))

# 循环保存每一页
for (i in seq_along(plot_chunks)) {
  chunk <- plot_chunks[[i]]
  
  # 创建大图 十列
  p <- wrap_plots(chunk, ncol = 10) &
    theme(
      plot.margin = margin(2, 2, 2, 2, "pt")
    )
  
  # 文件名按顺序编号
  filename <- sprintf("prob_matrix_top%dto%d_page%03d.pdf", n_from, n_to, i)
  ggsave(filename, p, width = 39, height = 30)
  
  cat(sprintf("Saved page %d: %s\n", i, filename))
}

#plot_prob_matrix_list <- lapply(
#  results_list,
#  function(x) x$res$plot_prob_matrix
#)

#plot_prob_matrix_list <- list()
#for (n in seq(1,100,1)) {
#  plot_prob_matrix_list[[n]] <- all_nb_mcc_results[[n]]$plot_prob_matrix
#}

# 画图时一行几个图
#p1 <- wrap_plots(plot_prob_matrix_list, ncol = 10) &
#  theme(
#    plot.margin = margin(2, 2, 2, 2, "pt")  # 设置每个子图的边距为1
#  )
#ggsave('prob_matrix_top500.pdf',p1,width = 39,height = 30)
####################################### 折现图 ############################################
# library(ggthemes)
# 
# aac_mcc_p <- aac_mcc
# aac_mcc_p$top_num <- as.numeric(seq(1,nrow(aac_mcc),1))
# aac_mcc_p <- aac_mcc_p[order(aac_mcc_p$top_num), ]
# aac_mcc_p$top <- paste0('top',aac_mcc_p$top_num)
# aac_mcc_p$top <- factor(aac_mcc_p$top, levels = aac_mcc_p$top)
# 
# # 定义要显示的标签
# # 创建序列并添加"top"前缀
# label_seq <- paste0("top", seq(0, 100, 5))
# # 将top0替换为top1
# label_seq[label_seq == "top0"] <- "top1"
# p2 <- ggplot(aac_mcc_p) +
#   geom_line(aes(x = top, y = aac, group = 1), color = "#577590", linewidth = 0.6) +
#   geom_point(aes(x = top, y = aac), color = "#577590", size = 1) +
#   geom_hline(yintercept = 0) +
#   geom_line(aes(x = top, y = mcc, group = 1),
#             color = "#F94144", linewidth = 0.6, linetype = "dashed") +
#   geom_point(aes(x = top, y = mcc),
#              color = "#F94144", size = 1, shape = 17) +
#   # 设置Y轴范围为[0,1]，并调整次要坐标轴
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, by = 0.1),
#   ) +
#   scale_x_discrete(
#     breaks = label_seq)+
#   # 加0.95的线
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey")+
#   # 在左侧添加标签
#   annotate(
#     "text",
#     x = -Inf, y = 0.95,  # 放在图形最左侧
#     label = "0.95",
#     hjust = -0.2,  # 稍微向左偏移
#     vjust = -0.5,  # 稍微向上
#     color = "grey40",
#     size = 3
#   )+
#   labs(x = "Top SNPs", y='', title = "AAC and MCC values") +
#   theme_hc() +
#   theme(
#     axis.text.x = element_text(size = 10, angle = 300,vjust = 1,hjust = 0),
#     axis.title.y.left = element_text(color = "#577590"),
#     plot.title = element_text(hjust = 0.5),
#     axis.text.y.left = element_text(color = "#577590"),
#     axis.title.y.right = element_text(color = "#F94144"),
#     axis.text.y.right = element_text(color = "#F94144")
#   )
# ggsave(paste0('aac_mcc.pdf'),p2,bg = 'white', width = 10, height = 7)
