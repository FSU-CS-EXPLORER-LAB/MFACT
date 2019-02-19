# ============================= Recommendation Model =======================================
# Purpose: Train and test performance of classification model and whether or
#          not to recommend for further study
# Date: 06/2017 

# ==========================================================================================


### Setting up environment =====================================
setwd("~/Google Drive/MFACT/Classification")

# --- Install and load relevant packages
packages <- c("tibble", "dplyr",  "ggplot2", "glmnet")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

lapply(packages, require, character.only = TRUE) # load packages

# --- Read data from your source
df <- read.csv(file = "data.csv", skip = 1, header = TRUE)


# --- Take a look at the data, examine variable types and formats
glimpse(df)

df <- df %>% select(-DIFF_comm)


# --- Correct the columns read with a percentage sign
perc_ind <- which(apply(df, 2, function(x) length(grep("%", x))>1))
df[, perc_ind] <- lapply(df[,perc_ind], function(x) as.numeric(sub("%", "", x))) 

# --- Correct columns read with a comma in the middle
comma_ind <- which(apply(df, 2, function(x) length(grep(",", x))>1))
df[, comma_ind] <- lapply(df[,comma_ind], function(x) as.numeric(gsub(",", "", x))) 

# --- a quick overall check of levels and statistics
summary(df)




### Building a Classification Model ==============================================


# --- create c sets of Monte Carlo Cross Validation indices by bootstrap sampling
n <- nrow(df)
c <- 100 # number of cross validation samples
cv_ind <- vector("list", c)
for (i in 1:c){
  set.seed(3*i+1)
  cv_ind[[i]]<- sample(1:n, round(.8*n), replace = F)
} # 100 sets of monte carlo training indices


# --- create 0/1 binary variable to serve as response (optional)
df$Response1 <- model.matrix(~ df$Response - 1)[,2]


# --- Define metric log loss for evaluation
logloss <- function(actual, predicted, eps = 1e-15) {
  predicted = pmin(pmax(predicted, eps), 1-eps)
  - (sum(actual * log(predicted) + (1 - actual) * log(1 - predicted))) / length(actual)
}

# --- Define lower (trivial model) and upper (full model) formula for forward selection
f_lo <- as.formula("Response1 ~ 1")
f_hi <- as.formula(paste0("Response1 ~ ", paste(colnames(df%>% select(-c(DIFF, Response1, Response, App, CLTag))), collapse = "+")))


# --- initialize evaluation metrics for the cross validated process
lloss_step <- vector("numeric", c) # log loss
mclass_step <- vector("numeric", c) # mis-classification
FPR_step <-  vector("numeric", c) # false positive rate
FNR_step <-  vector("numeric", c) # false negative rate
beta_step <- vector("list", c) # beta - coefficients
yhat <- vector("list", c) # predicted y (on CV test set)
yy <- vector("list", c) # true y (on CV test set)

mclass_step_t <- vector("numeric", c) # mis-classification after post-modeling filter
FPR_step_t <-  vector("numeric", c) # false positive rate after post-modeling filter
FNR_step_t <-  vector("numeric", c) # false negative rate after post-modeling filter

for (i in 1:c){ # looping through the cross validated set
  m_step <- step(glm(f_lo, data = df[cv_ind[[i]],], family = "binomial"), 
                 scope =  f_hi,
                 direction = "forward", 
                 trace = 0, 
                 steps = 5)  
  
  beta_step[[i]] <- m_step$coefficients
  
  yhat[[i]] <- predict(m_step, df[-cv_ind[[i]],], type = "response")
  yy[[i]] <- df$Response1[-cv_ind[[i]]]
  # add a rule to further screen 
  yhat2 <- yhat[[i]]
  yhat2[which(df$CL[-cv_ind[[i]]] != "bw")] = 0
  
  
  lloss_step[i] <- logloss(yy[[i]], yhat[[i]])
  mclass_step[i] <- mean((yhat[[i]]> .5) !=  yy[[i]])
  FPR_step[i] <- sum((yhat[[i]] > .5)*1 -  yy[[i]]==1)/sum(yy[[i]]== 0)
  FNR_step[i] <- sum((yhat[[i]] > .5)*1 -  yy[[i]]== -1)/sum(yy[[i]]== 1)
  
  
  # another set of metric computed with tweeked yhat
  mclass_step_t[i] <- mean((yhat2 > .5) !=  yy[[i]])
  FPR_step_t[i] <- sum((yhat2 > .5)*1 -  yy[[i]]==1)/sum(yy[[i]]== 0)
  FNR_step_t[i] <- sum((yhat2 > .5)*1 -  yy[[i]]== -1)/sum(yy[[i]]== 1)
}



### Model evaluation =========================================================================

# --- aggregate coefficients on all models
for (i in 1:c){
  if ( i == 1){
    res <- as.data.frame(t(beta_step[[i]]))
  } else{
    res <- bind_rows(res, as.data.frame(t(beta_step[[i]])))
  }
}

write.table(t(res),
            row.names = rownames(t(res)),
            file = "~/Google Drive/Com_Mat_Plot/Classification/StepSelection/Coef100.csv",
            sep = ",")



# --- evaluation metrics

metric_step <- cbind(lloss_step, mclass_step, FPR_step, FNR_step, mclass_step_t, FPR_step_t, FNR_step_t)
apply(metric_step, 2, function(x) mean(x, trim=.02, na.rm = TRUE))
apply(metric_step, 2, summary)
# write.table(tmp, sep = "\t")

write.table(metric_step,
            file = "~/Google Drive/Com_Mat_Plot/Classification/StepSelection/metric.csv",
            sep = ",")

# False negative rate by number of true positives in the c number of CV sets
tmp <- data.frame(nPos = unlist(lapply(yy, function(x) sum(x==1))))
tmp$b <- cut(tmp$nPos, breaks = 4)
tmp$FNR <- FNR_step

tmp %>% group_by(b) %>%
  summarise(n_obs = n(),
            m_FNR = mean(FNR), 
            std = sd(FNR),
            min_FNR = min(FNR),
            max_FNR = max(FNR))


# --- which variables get selected the most?
beta_mtrx_step <- t(res)
beta_rank <- as.data.frame(apply(beta_mtrx_step, 1, function(x) sum(!is.na(x))/c) %>% sort(decreasing = TRUE))
beta_rank <- rownames_to_column(beta_rank, "var")
colnames(beta_rank) <- c("var", "rank")

beta_mean <- as.data.frame(apply(t(res), 1, function(x) mean(x, na.rm=T)))
beta_mean <- rownames_to_column(beta_mean, "var")
colnames(beta_mean) <- c("var", "beta")
result_step <- (beta_rank) %>% left_join(as.data.frame(beta_mean), by = "var")

write.table(result_step,
            file = "~/Google Drive/Com_Mat_Plot/Classification/StepSelection/Coef.csv",
            sep = ",")




### Fitting a final model with the selected model for future recommendation ===================
# f_final <- as.formula(paste0("Response1 ~ ", paste(result_step$var[2:(2+5)], collapse = "+")))
# 
# # --- use the entire data set
# m <- glm(Response1 ~ CL + PoSYN + R + Tasyn + CRComm, data = df, family = "binomial")
