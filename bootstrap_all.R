library(tidyverse)
library(ggpubr)
library(ordinal)
library(optimx)
library(stats)
library(data.table)

rm(list = ls())

# functions
# Different classes of weighting functions ============================
kumaraswamy_weights <- function(steps, a, b){
  kum_weights <- rep(0, length(steps) - 1)
  for (i in 2:length(steps)){
    kum_weights[i-1] = (1 - steps[i-1]^a)^b - (1 - steps[i]^a)^b
  }
  return(kum_weights)
}

kum_par <- c(1,1, 0.6,1, 0.9,1, 1,0.5, 5,0.9, 0.5,0.5, 0.9,0.9, 5,5, 3,6, 2,6)
dim(kum_par) <- c(2, 10)
kum_par <- t(kum_par)
colnames(kum_par) <- c("a", "b")
kum_par <- kum_par[1, ]

steps <- c(0, 1/6, 2/6, 3/6, 4/6, 5/6, 6/6)

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

N <- 10000 


# specify conditions =======================================================

all_data <- read_csv("lab_and_lts.csv")
clean_data <- all_data[!(all_data$external_id %in% c("61", "65", "42", "32", "41", "25", "20")), ]

conditions <- list(
  full_lts = list( all_data$group == "lts"),
  full_labAll = list(all_data$group != "lts"),
  full_labPl = list(all_data$group == "labPl"),
  full_labNor = list(all_data$group == "labNor")
)

# Create DataFrames using a loop
data_frames <- list()
for (condition_name in names(conditions)) {
  condition <- conditions[[condition_name]]
  data_frames[[condition_name]] <- all_data %>% filter(!!!condition)
}


conditions_clean <- list(
  clean_lts = list( clean_data$group == "lts"),
  clean_labAll = list(clean_data$group != "lts"),
  clean_labPl = list(clean_data$group == "labPl"),
  clean_labNor = list(clean_data$group == "labNor")
)

clean_data_frames <- list()
for (condition_name in names(conditions_clean)) {
  condition <- conditions_clean[[condition_name]]
  clean_data_frames[[condition_name]] <- clean_data %>% filter(!!!condition)
}

all_data_frames <- c(data_frames, clean_data_frames)
df_names <- names(all_data_frames)

# Specify the file path where you want to save the output
output_file <- "optimal_models_parameters.txt"

# Create an empty file or open an existing file in append mode
file_con <- file(output_file, "a")

for (i in seq_along(all_data_frames)){
  data_by_weeks <- all_data_frames[[i]]
  data_by_weeks$os <- factor(data_by_weeks$q2, ordered = TRUE, 
                             levels = c(1, 2, 3, 4, 5))
  file_name <- paste0(df_names[i], ".rds")
  
  init_params <- kum_par
  lower_bounds <- c(0.00001, 0.00001)
  opt_data <- data_by_weeks 
  opt_result <- optimx(init_params, objective_fun, method="L-BFGS-B", 
                       lower=lower_bounds, df = opt_data, steps = steps)
  print(opt_result)
  output <- capture.output(print(opt_result))
  writeLines(output, con = file_con)
  
  data_by_weeks$f_best_fun <- f_kum(as.matrix(
    data_by_weeks[,c("mo", "tu", "we", "th", "fr", "sa")]), 
    steps, opt_result$a, opt_result$b)
  opt_model <- clm(os ~ f_best_fun, data = data_by_weeks)
  summary(opt_model)

  ggplot(NULL, aes(steps[2:7], kumaraswamy_weights(steps, opt_result$a, opt_result$b))) + 
    geom_point()

# 1. Predict probabilities
  new_data <- data.frame(f_best_fun = data_by_weeks$f_best_fun)
  predicted_probs <- predict(opt_model, newdata = new_data, type = "prob")

# Convert predicted_probs from list to matrix
  predicted_probs_matrix <- do.call(rbind, predicted_probs)

# Pre-allocate a matrix to store the samples
  random_samples_matrix <- matrix(0, nrow = nrow(predicted_probs_matrix), ncol = N)

# Draw random levels based on predicted probabilities for each observation
  for(i in 1:nrow(predicted_probs_matrix)) {
    random_samples_matrix[i, ] <- sample(1:ncol(predicted_probs_matrix), 
                                        size = N, 
                                        replace = TRUE, 
                                        prob = predicted_probs_matrix[i, ])
  }

  tmp_data <- data_by_weeks[,c("mo", "tu", "we", "th", "fr", "sa", "q2")]
  n_steps <- length(steps) - 1
  points_bootstrap <- matrix(0, nrow = n_steps*ncol(random_samples_matrix), ncol = 3)
  for(i in 1:ncol(random_samples_matrix)) {
    tmp_data$q2 <- random_samples_matrix[ ,i]
    tmp_data$os <- factor(tmp_data$q2, ordered = TRUE, levels = c(1, 2, 3, 4, 5))
    init_params <- kum_par # uniform distribution
    lower_bounds <- c(0.00001, 0.00001)
    opt_result <- optimx(init_params, objective_fun, method="L-BFGS-B", 
                        lower=lower_bounds, df = tmp_data, steps = steps)

    cat(i, kumaraswamy_weights(steps, opt_result$a, opt_result$b), "\n")
  # done protect against estimation error
    if (is.na(opt_result$kkt1) == FALSE & is.na(opt_result$kkt2) == FALSE & opt_result$kkt1 == TRUE & opt_result$kkt2 == TRUE) {
 
      points_bootstrap[((i - 1)*n_steps + 1):(i*n_steps), 1] = i
      points_bootstrap[((i - 1)*n_steps + 1):(i*n_steps), 2] = steps[2:(n_steps + 1)]
      points_bootstrap[((i - 1)*n_steps + 1):(i*n_steps), 3] = kumaraswamy_weights(steps, opt_result$a, opt_result$b)
    }
    else {
      points_bootstrap[((i - 1)*n_steps + 1):(i*n_steps), 1] = i
      points_bootstrap[((i - 1)*n_steps + 1):(i*n_steps), 2] = steps[2:(n_steps + 1)]
      points_bootstrap[((i - 1)*n_steps + 1):(i*n_steps), 3] = NA 
      }
  }
  points_bootstrap <- na.omit(points_bootstrap)
  bootstrap_data <- tibble(rep = points_bootstrap[ ,1], step = points_bootstrap[ ,2], weight = points_bootstrap[ ,3])
  saveRDS(bootstrap_data, file_name)
}

close(file_con)