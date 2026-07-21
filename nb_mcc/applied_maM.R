# 用于计算多分类问题的宏观平均MCC（Matthews Correlation Coefficient）及其渐近置信区间
# https://github.com/tamuraj/Multiclass-MCC/blob/main/applied_maM.R

###############################################################
# Point estimates and asymptotic confidence intervals for maM #
###############################################################

# The code derives point estimates and asymptotic confidence intervals for maM. The function arguments are shown below.

# prob_matri: Matrix of probabilities
# n         : samplesize
# conf.level: confidence interval level


compute_macroMcc_CI <- function(prob_matrix, n, conf.level) {
  # Function to calculate TP, FP, FN, and TN
  calculate_measures <- function(prob_matrix) {
    r <- nrow(prob_matrix)
    
    TP <- numeric(r)
    FP <- numeric(r)
    FN <- numeric(r)
    TN <- numeric(r)
    
    for (i in 1:r) {
      # True Positive (TP)
      TP[i] <- prob_matrix[i, i]
      
      # False Positive (FP)
      FP[i] <- sum(prob_matrix[, i]) - prob_matrix[i, i]
      
      # False Negative (FN)
      FN[i] <- sum(prob_matrix[i, ]) - prob_matrix[i, i]
      
      # True Negative (TN)
      TN[i] <- sum(prob_matrix) - (TP[i] + FP[i] + FN[i])
    }
    
    return(list(TP = TP, FP = FP, FN = FN, TN = TN))
  }
  
  
  calculate_mcc <- function(TP, FP, FN, TN) {
    numerator <- TP * TN - FP * FN
    denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    
    
    mcc <- ifelse(denominator == 0, 0, numerator / denominator)
    
    return(mcc)
  }
  
  
  macro_averaged_mcc <- function(prob_matrix) {
    measures <- calculate_measures(prob_matrix)
    TP <- measures$TP
    FP <- measures$FP
    FN <- measures$FN
    TN <- measures$TN
    
    r <- length(TP)
    mcc_values <- numeric(r)
    
    for (i in 1:r) {
      mcc_values[i] <- calculate_mcc(TP[i], FP[i], FN[i], TN[i])
    }
    
    macro_mcc <- mean(mcc_values)
    
    return(macro_mcc)
  }
  
  # Function to compute A_{ii}^{maM}
  compute_Aii_maM <- function(p_i_dot, p_dot_i, p_ii, r) {
    term1 <- (1 - p_i_dot - p_dot_i) / sqrt(p_i_dot * p_dot_i * (1 - p_i_dot) * (1 - p_dot_i))
    
    numerator_inner <- (p_dot_i + p_i_dot) * (1 - p_dot_i) - p_i_dot * p_dot_i
    inner_expression <- numerator_inner * (1 - p_i_dot) - p_i_dot * p_dot_i * (1 - p_dot_i)
    numerator <- (p_ii - p_dot_i * p_i_dot) * inner_expression
    denominator <- 2 * (p_dot_i * p_i_dot * (1 - p_dot_i) * (1 - p_i_dot))^(3/2)
    term2 <- numerator / denominator
    
    Aii_maM <- (1 / r) * (term1 - term2)
    
    return(Aii_maM)
  }
  
  # Function to compute A_{ij}^{maM}
  compute_Aij_maM <- function(p_j_dot, p_dot_j, p_i_dot, p_dot_i, p_jj, p_ii, r) {
    term1_part1 <- -p_j_dot / sqrt(p_j_dot * p_dot_j * (1 - p_j_dot) * (1 - p_dot_j))
    
    numerator1_inner <- p_j_dot * (1 - p_j_dot) * (1 - p_dot_j) - p_j_dot * p_dot_j * (1 - p_j_dot)
    numerator1 <- (p_jj - p_j_dot * p_dot_j) * numerator1_inner
    denominator1 <- 2 * (p_j_dot * p_dot_j * (1 - p_j_dot) * (1 - p_dot_j))^(3/2)
    term1_part2 <- numerator1 / denominator1
    
    term1 <- term1_part1 - term1_part2
    
    term2_part1 <- -p_dot_i / sqrt(p_i_dot * p_dot_i * (1 - p_i_dot) * (1 - p_dot_i))
    
    numerator2_inner <- p_dot_i * (1 - p_i_dot) * (1 - p_dot_i) - p_i_dot * p_dot_i * (1 - p_dot_i)
    numerator2 <- (p_ii - p_i_dot * p_dot_i) * numerator2_inner
    denominator2 <- 2 * (p_i_dot * p_dot_i * (1 - p_dot_i) * (1 - p_i_dot))^(3/2)
    term2_part2 <- numerator2 / denominator2
    
    term2 <- term2_part1 - term2_part2
    
    Aij_maM <- (1 / r) * (term1 + term2)
    
    return(Aij_maM)
  }
  
  
  compute_gradient <- function(prob_matrix) {
    r <- nrow(prob_matrix)
    
    gradient <- numeric(r * r)
    
    for (i in 1:r) {
      for (j in 1:r) {
        if (i == j) {
          gradient[(i - 1) * r + j] <- compute_Aii_maM(
            sum(prob_matrix[i, ]), sum(prob_matrix[, j]), prob_matrix[i, j], r
          )
        } else {
          gradient[(i - 1) * r + j] <- compute_Aij_maM(
            sum(prob_matrix[j, ]), sum(prob_matrix[, j]), sum(prob_matrix[i, ]), sum(prob_matrix[, i]), prob_matrix[j, j], prob_matrix[i, i], r
          )
        }
      }
    }
    
    return(gradient)
  }
  
  # Calculate maM
  mcc <- macro_averaged_mcc(prob_matrix)
  
  
  d <- compute_gradient(prob_matrix)
  pa <- c(t(prob_matrix))
  d1 <- c(d)
  p <- t(t(pa))
  dp <- diag(pa)
  ma <- rbind(d1)
  
  
  s1 <- sum(pa * (d1^2)) - (sum(pa * d1))^2
  
  # Calculate confidence interval
  z <- qnorm(1 - (1 - conf.level) / 2)
  v <- s1 / n 
  s <- sqrt(v)
  ci <- mcc + c(-1, 1) * z * s
  
  output <- c("maM" = mcc, "confidence_interval" = ci)
  return(output)
}


# Example 
# prob_matrix <- matrix(c(
#   0.1, 0.05, 0.05,
#   0.1, 0.3, 0.1,
#   0.05, 0.1, 0.15
# ), nrow = 3, byrow = TRUE)
# compute_macroMcc_CI(prob_matrix,1000,0.95)