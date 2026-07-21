#!/usr/bin/R

library(vegan)
library(ggplot2)
library(patchwork)
library(ggthemes)

arg=commandArgs(T)
origin_evec <- arg[1]  # 原始PCA evec结果的绝对路径

# 原始PCA evec结果的绝对路径
#origin_evec <- "/mnt/hdd1/HaodongChen/1kg/smartpca_EAS/exclude_outlier_true/merged_1kg_EAS_autosomes_qc_pruned.evec"
# 普鲁克分析结果存放目录
output_dir  <- "procrustes_results"

source('procrustes_analysis.R')

results   <- vector("list", 2000)
plot_list <- vector("list", 2000)
m2_cor    <- data.frame()

failed_n <- c()  # 记录失败的 n

for (n in 1:2000) {

  message("Running top", n)

  res <- tryCatch(
    run_procrustes_analysis(
      origin_evec = origin_evec,
      sub_evec = sprintf(
        "pca_results/smartpca_aim%d.evec",
        n
      ),
      n_pcs_input = 5,
      output_dir = output_dir
    ),
    error = function(e) {
      message("❌ Failed at top", n, ": ", e$message)
      return(NULL)
    }
  )

  # ---------- 如果失败，跳过 ----------
  if (is.null(res)) {
    failed_n <- c(failed_n, n)
    next
  }

  # ---------- 成功才继续 ----------
  results[[n]] <- res
  plot_list[[n]] <- res$plot

  m2_cor <- rbind(
    m2_cor,
    data.frame(
      top = res$top_value,
      m2 = res$m2,
      correlation = res$correlation,
      p_value = res$p_value
    )
  )
}

write.csv(m2_cor, "procrustes_m2_cor.csv", row.names = FALSE)

write.table(failed_n,
  file = "failed_procrustes_top.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

