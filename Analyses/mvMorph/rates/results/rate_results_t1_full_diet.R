library(tidyverse)

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_70_75/full/results/diet"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=Diet)) +
  geom_boxplot()


resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_75_80/full/results/diet"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=Diet)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_80_85/full/results/diet"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=Diet)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_85_90/full/results/diet"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=Diet)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_90_95/full/results/diet"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=Diet)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_95_100/full/results/diet"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=Diet)) +
  geom_boxplot()

ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
  geom_boxplot()

  