function (count_data, metadata, numerator, denominator, response_variable, 
          Factor, log_transform, p_adjust, test_type) 
{
count_data = DF
metadata = c57_nos2KO_mouse_metadata
numerator = "Strep"
denominator = "Mock"
response_variable = "Metabolite"
Factor = "Treatment"
log_transform = TRUE
p_adjust = "BH"
test_type = "welch"
  rownames(count_data) <- count_data[, response_variable]
  count_data[, response_variable] <- NULL
  data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]
  data_Transpose <- as.data.frame(t(data_Int))
  data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), 
                                        data_Transpose))
  Factor = metadata[, Factor]
  data_Transpose$Factor = Factor[match(metadata$Sample, data_Transpose$Sample)]
  data_Subset <- filter(data_Transpose, Factor == numerator | 
                          Factor == denominator)
  rownames(data_Subset) <- data_Subset[, "Sample"]
  data_Subset[, "Sample"] <- NULL
  data_Numeric <- data_Subset[sapply(data_Subset, function(x) is.numeric(x))]
  data_Numeric <- data.frame(lapply(data_Numeric, function(x) as.numeric(as.character(x))), 
                             check.names = F, row.names = rownames(data_Numeric))
# transpose and rownames are now plot names
  
  data_Log <- as.data.frame(log(data_Numeric))
  data_Log$Factor = data_Subset$Factor
  cols_to_test <- data_Log[sapply(data_Log, function(x) is.numeric(x))]
  Vect = colnames(cols_to_test)
  
# Vect is plot names as an array
  
  data_Numeric$Factor <- data_Subset$Factor
  if (log_transform == FALSE) {
    data_mod = data_Numeric
  }
  else if (log_transform == TRUE) {
    data_mod = data_Log
  }
  model = data_mod[, "Factor"]

# model = Factor array (ours would be control and heated)
s <- (data_mod[Vect])

    Run_Tests <- function(data_Log, Vect, model) {
      results <- ldply(Vect, function(Metabolite) {
        t_val = t.test(data_mod[[Metabolite]] ~ model, 
                       var.equal = FALSE)$statistic
        p_val = t.test(data_mod[[Metabolite]] ~ model, 
                       var.equal = FALSE)$p.value
        return(data.frame(Metabolite = Metabolite, t_value = t_val, 
                          pval = p_val))
      })
    }
    
library(plyr)
  results <- Run_Tests(data_Log = data_Log, Vect = Vect, model = model)
  data_Log$Factor <- factor(data_Log$Factor, levels = c(numerator, 
                                                        denominator))
  data_Numeric$Factor = data_Log$Factor
  Means <- data_Numeric %>% group_by(Factor) %>% summarise_all(funs(mean))
  Means_T = as.data.frame(t(Means))
  colnames(Means_T) <- as.character(unlist(Means_T[1, ]))
  Means_T = Means_T[-1, ]
  colnames(Means_T)[1] <- numerator
  colnames(Means_T)[1] <- paste(colnames(Means_T)[1], "mean", 
                                sep = ".")
  colnames(Means_T)[2] <- denominator
  colnames(Means_T)[2] <- paste(colnames(Means_T)[2], "mean", 
                                sep = ".")
  Means_T[, 3] <- as.numeric(as.character(Means_T[, 1]))/as.numeric(as.character(Means_T[, 
                                                                                         2]))
  colnames(Means_T)[3] <- "Fold_Change"
  Means_T[, 4] = log2(as.numeric(as.character(Means_T[, 3])))
  colnames(Means_T)[4] = "log2FoldChange"
  Means_T <- cbind(rownames(Means_T), data.frame(Means_T, 
                                                 row.names = NULL))
  colnames(Means_T)[1] <- response_variable
  st.err <- function(x) sd(x)/sqrt(length(x))
  stdev <- data_Log %>% group_by(Factor) %>% summarise_all(funs(sd))
  stdev <- as.data.frame(stdev)
  rownames(stdev) <- stdev[, "Factor"]
  stdev[, "Factor"] <- NULL
  stdev_t = as.data.frame(t(stdev))
  stdev_t <- cbind(rownames(stdev_t), data.frame(stdev_t, 
                                                 row.names = NULL))
  colnames(stdev_t)[1] <- response_variable
  colnames(stdev_t)[2] <- numerator
  colnames(stdev_t)[2] <- paste(colnames(stdev_t)[2], "stdev", 
                                sep = ".")
  colnames(stdev_t)[3] <- denominator
  colnames(stdev_t)[3] <- paste(colnames(stdev_t)[3], "stdev", 
                                sep = ".")
  SE <- data_Log %>% group_by(Factor) %>% summarise_all(funs(st.err))
  SE <- as.data.frame(SE)
  rownames(SE) <- SE[, "Factor"]
  SE[, "Factor"] <- NULL
  SE_t = as.data.frame(t(SE))
  SE_t <- cbind(rownames(SE_t), data.frame(SE_t, row.names = NULL))
  colnames(SE_t)[1] <- response_variable
  colnames(SE_t)[2] <- numerator
  colnames(SE_t)[2] <- paste(colnames(SE_t)[2], "std.err", 
                             sep = ".")
  colnames(SE_t)[3] <- denominator
  colnames(SE_t)[3] <- paste(colnames(SE_t)[3], "std.err", 
                             sep = ".")
  colnames(results)[1] <- response_variable
  results$padj = p.adjust(results$pval, method = p_adjust)
  count_data = cbind(rownames(count_data), data.frame(count_data, 
                                                      row.names = NULL))
  colnames(count_data)[1] <- "Metabolite"
  results = merge(results, count_data, by = response_variable, 
                  all = TRUE)
  results = merge(Means_T, results, by = response_variable, 
                  all.y = TRUE)
  results = merge(SE_t, results, by = response_variable, all.y = TRUE)
  results = merge(stdev_t, results, by = response_variable, 
                  all.y = TRUE)
  class(results) = append(class(results), "cpd")
  return(results)
}
