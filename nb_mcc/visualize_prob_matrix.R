library(ggplot2)

# 可视化函数
  visualize_prob_matrix <- function(prob_matrix, title = "Matrix of probabilities") {
    # 转换为长格式
    prob_melted <- melt(prob_matrix)
    colnames(prob_melted) <- c("Predicted", "Actual", "Probability")

    # 创建热图
    p <- ggplot(prob_melted, aes(x = Actual, y = Predicted, fill = Probability)) +
      geom_tile(color = "white", size = 1) +
      geom_text(aes(label = round(Probability, 3)),
                color = "black", size = 2) +
      scale_fill_gradient2(low = "white", high = "steelblue",
                           midpoint = max(prob_melted$Probability)/2) +
      labs(title = title,
           x = "Actual",
           y = "Predicted") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5))

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

