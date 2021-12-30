if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(Peptides)) install.packages("Peptides", repos = "http://cran.us.r-project.org")
if(!require(ggpubr)) install.packages("ggpubr", repos = "http://cran.us.r-project.org")
if(!require(broom)) install.packages("broom", repos = "http://cran.us.r-project.org")
if(!require(glue)) install.packages("glue", repos = "http://cran.us.r-project.org")

library(tidyverse)
library(caret)
library(data.table)
library(Peptides)
library(ggpubr)
library(broom)
library(glue)


# Calculate the charge of amino acids at pH 2.7 (default pH for peptide chromatography)
sort(c(sapply(aaList(), function(x) charge(x, pH=2.7))))

# Calculate hydrophobicity for amino acids
sort(c(sapply(aaList(), hydrophobicity)))

# Calculate molecular weight
sort(c(sapply(aaList(), mw)))

# Composition of amino acids in data set
aaComp(paste(peptide_rt$sequence, collapse = ""))

# Calculating and predicting peptide features
peptide_rt <- peptide_rt %>%
  mutate(length = str_length(sequence), # adding lenght of peptide
         charge = charge(sequence, pH=2.7), # adding charge
         hydrophobicity = hydrophobicity(sequence), # adding hydrophobicity
         molecular_weight = mw(sequence), # adding molecular weight
         aliphatic_index = aIndex(sequence), #adding aliphatic index
  )
head(peptide_rt)

median(peptide_rt$length)

pep_rt <- peptide_rt %>% 
  ggplot(aes(x = RT/60)) +
  geom_histogram(bins = 30) +
  ggtitle("Retention time (min)")

pep_length <- peptide_rt %>%
  ggplot(aes(x = length)) +
  geom_histogram(bins = 30) +
  ggtitle("Length")

pep_charge <- peptide_rt %>%
  ggplot(aes(x = charge)) +
  geom_histogram(bins = 12) +
  ggtitle("Charge")

pep_hydrophobicity <- peptide_rt %>%
  ggplot(aes(x = hydrophobicity)) +
  geom_histogram(bins = 30) +
  ggtitle("Hydrophobicity")

pep_mw <- peptide_rt %>%
  ggplot(aes(x = molecular_weight)) +
  geom_histogram(bins = 30) +
  ggtitle("Molecular weight")

pep_aliphatic <- peptide_rt %>%
  ggplot(aes(x = aliphatic_index)) +
  geom_histogram(bins = 30) +
  ggtitle("Aliphatic index")

ggarrange(pep_rt, pep_length, pep_charge, 
          pep_hydrophobicity, pep_mw, pep_aliphatic,
          ncol = 3, nrow = 2)

pep_rt_length <- peptide_rt %>% 
  ggplot(aes(x = RT/60, y = length)) +
  geom_point(alpha = 0.05) +
  ggtitle("RT vs. length") +
  geom_smooth(method = "loess")

pep_rt_charge <- peptide_rt %>% 
  ggplot(aes(x = RT/60, y = charge)) +
  geom_point(alpha = 0.05) +
  ggtitle("RT vs. charge") +
  geom_smooth(method = "loess")

pep_rt_hydrophobicity <- peptide_rt %>% 
  ggplot(aes(x = RT/60, y = hydrophobicity)) +
  geom_point(alpha = 0.05) +
  ggtitle("RT vs. hydrophobicity") +
  geom_smooth(method = "loess")

pep_rt_mw <- peptide_rt %>% 
  ggplot(aes(x = RT/60, y = molecular_weight)) +
  geom_point(alpha = 0.05) +
  ggtitle("RT vs. MW") +
  geom_smooth(method = "loess")

pep_rt_aliphatic <- peptide_rt %>% 
  ggplot(aes(x = RT/60, y = aliphatic_index)) +
  geom_point(alpha = 0.05) +
  ggtitle("RT vs. aliphatic index") +
  geom_smooth(method = "loess")

ggarrange(pep_rt_length, pep_rt_charge, pep_rt_hydrophobicity, 
          pep_rt_mw, pep_rt_aliphatic,
          ncol = 3, nrow = 2)

# Split data in train and test set. Test set is 10% of peptide_rt data.
suppressWarnings(set.seed(1, sample.kind="Rounding"))
test_index <- createDataPartition(y = peptide_rt$RT, times = 1, p = 0.1, list = FALSE)
train_data <- peptide_rt[-test_index,]
dim(train_data)

test_data <- peptide_rt[test_index,]
dim(test_data)

# Calculate the RMSE for predicting RT with peptide length
fit_l <- lm(RT ~ length, data = train_data)
RMSE(predict(fit_l, test_data), test_data$RT)/60

# Calculate the RMSE for predicting RT with peptide hydrophobicity
fit_h <- lm(RT ~ hydrophobicity, data = train_data)
RMSE(predict(fit_h, test_data), test_data$RT)/60

# Calculate the RMSE for predicting RT with peptide molecular weight
fit_m <- lm(RT ~ molecular_weight, data = train_data)
RMSE(predict(fit_m, test_data), test_data$RT)/60

# Calculate the RMSE for predicting RT with peptide aliphatic index
fit_a <- lm(RT ~ aliphatic_index, data = train_data)
RMSE(predict(fit_a, test_data), test_data$RT)/60

# Calculate the RMSE for predicting RT with peptide charge
fit_c <- lm(RT ~ charge, data = train_data)
RMSE(predict(fit_c, test_data), test_data$RT)/60

# Calculate the RMSE for predicting RT with all available peptide features
fit_lhmac <- lm(RT ~ length + hydrophobicity + molecular_weight + aliphatic_index + charge, data = train_data)
RMSE(predict(fit_lhmac, test_data), test_data$RT)/60

# Define models to be tested and create empty vectors to record data
models <- c("glm", "knn", "svmLinear", "rf")
model_rmse <- numeric()
sec_per_model <- numeric()

# Train all models, take the time for each training and calculate RMSE on
# predictions for the test set
for (model in models){
  tic <- Sys.time()
  suppressWarnings(fit <- train(RT ~ length + hydrophobicity + 
                                  molecular_weight + aliphatic_index + charge,
                                method = model,
                                data = train_data))
  toc <- Sys.time()
  duration <- round(as.numeric(difftime(toc, tic, units="secs")),1)
  sec_per_model[model] <- duration
  
  rmse <- RMSE(predict(fit, test_data), test_data$RT)/60
  model_rmse[model] <- rmse
}

# Overview of time needed to train each model in seconds
sec_per_model

# Overview RMSE in min for RT predictions for the test set for each model
model_rmse

