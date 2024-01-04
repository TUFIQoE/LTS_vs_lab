library(tidyverse)
library(ordinal)
library(optimx)

rm(list = ls())

data_all <- readRDS("all_data_1204.RDS")
#data_all$rejected <- ifelse(data_all$external_id %in% c("61", "65", "42", "32", "41", "25", "20"), 'yes', 'no')
data_all <- data_all[!(data_all$external_id %in% c("61", "65", "42", "32", "41", "25", "20")), ]
f_mean <- function(x) {
  tmp <- x[,1]
  tmp[is.na(x[,1])] <- 0
  norm <- rep(1, nrow(x))
  norm[is.na(x[,1])] <- 0
  for (i in c(2:6)){
    tmp_norm <- rep(1, nrow(x))
    tmp_norm[is.na(x[,i])] <- 0
    norm <- norm + tmp_norm
    tmp_na <- x[,i]
    tmp_na[is.na(tmp_na)] <- 0
    tmp <- tmp + tmp_na
  }
  return(tmp / norm)
}

f_min <- function(x) {
  tmp <- x[,1]
  tmp[is.na(x[,1])] <- Inf  # Initialize with positive infinity
  for (i in c(2:6)){
    tmp_na <- x[,i]
    tmp_na[is.na(tmp_na)] <- Inf  # Treat missing values as positive infinity
    tmp <- pmin(tmp, tmp_na)  # Calculate the element-wise minimum
  }
  tmp[tmp == Inf] <- NA  # Replace positive infinity with NA
  return(tmp)
}

f_max <- function(x) {
  tmp <- x[,1]
  tmp[is.na(x[,1])] <- Inf  # Initialize with positive infinity
  for (i in c(2:6)){
    tmp_na <- x[,i]
    tmp_na[is.na(tmp_na)] <- Inf  # Treat missing values as positive infinity
    tmp <- pmax(tmp, tmp_na)  # Calculate the element-wise minimum
  }
  tmp[tmp == Inf] <- NA  # Replace positive infinity with NA
  return(tmp)
}

predict_and_update <- function(model, data, weights, model_column) {
  predicted_probs <- predict(model, newdata = data, type = "prob")
  predicted_probs_matrix <- do.call(rbind, predicted_probs)
  weights <- matrix(rep(c(1,2,3,4,5), 13), nrow = 13, ncol = 5, byrow = TRUE)
  data[[model_column]] <- rowSums(predicted_probs_matrix * weights)
  return(data)
}

normalize_beta <- function(model) {
  beta_values <- model$beta
  sum_of_coefficients <- sum(beta_values)
  normalized_beta <- beta_values / sum_of_coefficients
  return(normalized_beta)
}

data_all$f_mean = f_mean(as.matrix(data_all[,c("mo", "tu", "we", "th", "fr", "sa")]))
data_all$f_max = f_max(as.matrix(data_all[,c("mo", "tu", "we", "th", "fr", "sa")]))
data_all$f_min = f_min(as.matrix(data_all[,c("mo", "tu", "we", "th", "fr", "sa")]))
data_all$os <- factor(data_all$q2, ordered = TRUE, levels = c(1, 2, 3, 4, 5))

steps <- c(0, 1/6, 2/6, 3/6, 4/6, 5/6, 6/6)

kumaraswamy_weights <- function(steps, a, b){
  kum_weights <- rep(0, length(steps) - 1)
  for (i in 2:length(steps)){
    kum_weights[i-1] = (1 - steps[i-1]^a)^b - (1 - steps[i]^a)^b
  }
  return(kum_weights)
}

f_kum <- function(x, steps, a, b) {
  norm_vec <- t(replicate(nrow(x), kumaraswamy_weights(steps, a, b)))
  norm_vec[is.na(x)] <- 0
  rowSums(x * norm_vec, na.rm = TRUE) / rowSums(norm_vec)
}

objective_fun <- function(params, df, steps) {
  a <- params[1]
  b <- params[2]
  
  df$tmp <- f_kum(as.matrix(df[,c("mo", "tu", "we", "th", "fr", "sa")]), 
                  steps, a, b)
  
  model <- clm(os ~ tmp, data = df)  
  
  return(AIC(model))  
}



kum_par <- c(a = 1, b = 1)
lower_bounds <- c(0.00001, 0.00001)
opt_result_lts <- optimx(kum_par, objective_fun, method="L-BFGS-B", 
                     lower=lower_bounds, df = data_all[data_all$group == "lts",], steps = steps)
print(opt_result_lts)
opt_result_lab_pl <- optimx(kum_par, objective_fun, method="L-BFGS-B", 
                         lower=lower_bounds, df = data_all[data_all$group == "labPl",], steps = steps)
print(opt_result_lab_pl)
opt_result_lab_nor <- optimx(kum_par, objective_fun, method="L-BFGS-B", 
                         lower=lower_bounds, df = data_all[data_all$group == "labNor",], steps = steps)
print(opt_result_lab_nor)
opt_result_lab <- optimx(kum_par, objective_fun, method="L-BFGS-B", 
                             lower=lower_bounds, df = data_all[data_all$group != "lts",], steps = steps)
print(opt_result_lab)

ggplot() + 
  geom_point(aes(steps[2:7], kumaraswamy_weights(steps, opt_result_lts$a, opt_result_lts$b)), color = "green") +
  geom_point(aes(steps[2:7], kumaraswamy_weights(steps, opt_result_lab$a, opt_result_lab$b)), color = "cyan") +
  geom_point(aes(steps[2:7], kumaraswamy_weights(steps, opt_result_lab_pl$a, opt_result_lab_pl$b)), color = "red") +
  geom_point(aes(steps[2:7], kumaraswamy_weights(steps, opt_result_lab_nor$a, opt_result_lab_nor$b)), color = "blue")


# ============ Testing GLZ ==================

data_all <- data_all  %>% 
  mutate(motu = (mo + tu)/2,
         weth = (we + th)/2,
         frsa = (fr + sa)/2
         )
data_lts <- data_all %>% filter(data_all$group == "lts")
data_lts <- data_lts %>%  mutate(f_best_fun = f_kum(as.matrix(
    data_lts[,c("mo", "tu", "we", "th", "fr", "sa")]), 
    steps, opt_result_lts$a, opt_result_lts$b))

data_nlts <- data_all %>% filter(data_all$group != "lts")
data_nlts <- data_nlts %>%  mutate(f_best_fun = f_kum(as.matrix(
  data_nlts[,c("mo", "tu", "we", "th", "fr", "sa")]), 
  steps, opt_result_lab$a, opt_result_lab$b))

data_labPl <- data_all %>% filter(data_all$group == "labPl")
data_labPl <- data_labPl %>%  mutate(f_best_fun = f_kum(as.matrix(
  data_labPl[,c("mo", "tu", "we", "th", "fr", "sa")]), 
  steps, opt_result_lab_pl$a, opt_result_lab_pl$b))

data_labNor <- data_all %>% filter(data_all$group == "labNor")
data_labNor <- data_labNor %>%  mutate(f_best_fun = f_kum(as.matrix(
  data_labNor[,c("mo", "tu", "we", "th", "fr", "sa")]), 
  steps, opt_result_lab_nor$a, opt_result_lab_nor$b))

d1 <- rbind(data_lts, data_nlts)
d2 <- rbind(data_labPl, data_labNor)


best_lts <- clm(os ~ mo + we + fr, data = data_lts)  #best
AIC(best_lts)

best_labPl <- clm(os ~ mo + th + sa, data = data_labPl)  #best
AIC(best_labPl)

best_labNor <- clm(os ~  mo + tu + we + th + fr + sa , data = data_labNor)  #best
AIC(best_labNor) 

