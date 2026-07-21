run_procrustes_analysis <- function(origin_evec,
                                    sub_evec,
                                    output_dir,
                                    n_pcs_input = 10,
                                    save_plot = TRUE) {

  stopifnot(n_pcs_input >= 2)

  # ---------- 路径 ----------
  origin_evec <- normalizePath(origin_evec, mustWork = TRUE)
  sub_evec    <- normalizePath(sub_evec, mustWork = TRUE)
  output_dir  <- normalizePath(output_dir, mustWork = FALSE)

  # ---------- 提取 top ----------
  top_value <- sub(".*(top[0-9]+).*", "\\1", basename(sub_evec))

  # ---------- 构造列名 ----------
  coln <- c("Sample", paste0("PC", 1:n_pcs_input), "Pop")

  # ---------- 读取 PCA ----------
  evec_all <- read.table(origin_evec, col.names = coln, stringsAsFactors = FALSE)
  evec_sub <- read.table(sub_evec,    col.names = coln, stringsAsFactors = FALSE)

  rownames(evec_all) <- evec_all$Sample
  rownames(evec_sub) <- evec_sub$Sample

  # ---------- 样本对齐 ----------
  common_samples <- intersect(rownames(evec_all), rownames(evec_sub))
  evec_all <- evec_all[common_samples, ]
  evec_sub <- evec_sub[common_samples, ]

  # ---------- 明确只使用 PC1 & PC2 ----------
  X <- as.matrix(evec_all[, c("PC1", "PC2")])
  Y <- as.matrix(evec_sub[, c("PC1", "PC2")])

  # ---------- Procrustes ----------
  pro_result <- procrustes(X, Y, symmetric = TRUE)

  set.seed(1949)
  pro_test <- protest(X, Y, permutations = 10000)

  # ---------- 作图数据 ----------
  plot_df <- data.frame(
    full_PC1 = pro_result$X[, 1],
    full_PC2 = pro_result$X[, 2],
    sub_PC1  = pro_result$Yrot[, 1],
    sub_PC2  = pro_result$Yrot[, 2]
  )

  # ---------- 作图 ----------
  p <- ggplot(plot_df) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_segment(
      aes(x = sub_PC1, y = sub_PC2, xend = full_PC1, yend = full_PC2),
      arrow = arrow(length = unit(0.15, "cm")),
      color = "grey70"
    ) +
    geom_point(aes(x = sub_PC1, y = sub_PC2),
               color = "#46F0F0", size = 1.8) +
    geom_point(aes(x = full_PC1, y = full_PC2),
               color = "#577590", size = 1.8, shape = 17) +
    theme_bw() +
    labs(
      subtitle = paste0("M² = ", round(pro_test$ss, 4),
                        "  p-value = ", format(pro_test$signif, scientific = TRUE, digits = 3)),
      x = "PC1",
      y = "PC2"
    ) +
    theme(
      panel.grid = element_blank(),
      plot.subtitle = element_text(hjust = 0.5)
    )+
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste("correlation =", round(pro_test$scale, 3)),
      hjust = 1.1,
      vjust = 1.5
    )

  # ---------- 保存 ----------
  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    ggsave(
      file.path(
        output_dir,
        paste0("procrustes_", top_value, ".pdf")
      ),
      p,
      width = 8,
      height = 6
    )
  }

  # ---------- 返回 ----------
  return(list(
    top_value = top_value,
    input_pcs = n_pcs_input,
    m2 = pro_test$ss,
    correlation = pro_test$scale,
    p_value = pro_test$signif,
    plot = p
  ))
}

