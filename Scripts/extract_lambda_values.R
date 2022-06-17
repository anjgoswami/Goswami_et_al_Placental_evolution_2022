library(BTprocessR)
library(tidyverse)

logfilesA <- list.files(pattern = "*a_run.txt.Log.txt$", recursive = TRUE)

logfilesA <- logfilesA[grepl("var_rate", logfilesA)]
logfilesA <- logfilesA[grepl("lambda", logfilesA)]

my_results <- tibble(filename=logfilesA)

treelist <- c("70_75","75_80","80_85","85_90","90_95","95_100")
my_results<- my_results %>%
  mutate(., Age = str_extract(filename, paste0(treelist,collapse = "|")))


  extractor_func<-function(x){
  post<-BTprocessR::loadPosterior(x)
  mean_lambda<-mean(post$Lambda)
  sd_lambda<-sd(post$Lambda)
  return(tibble(mean_lambda=mean_lambda, sd_lambda=sd_lambda))
}

raw_lambda_results <- lapply(logfilesA, function(x) extractor_func(x))
lambda_results <- do.call(rbind, raw_lambda_results) %>% as_tibble()
my_results <- bind_cols(my_results,lambda_results)

my_results <- my_results %>% mutate(Topology = "Tree 3")

write.csv(my_results, file = "/mnt/shared/projects/nhm/goswamilab/BT_placental/lambda_summary_tree3.csv")
