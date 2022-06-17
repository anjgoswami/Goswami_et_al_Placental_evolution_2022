library(dplyr)
library(tidyverse)



load("./list_of_results_diet.Rdata")


#simple plot, single tree
data<-bind_rows(lapply(list_of_results, '[[', "param"))
data <- data %>% mutate(bone = names(list_of_results))
data <- data %>%
        pivot_longer(!bone, names_to = "diet", values_to = "rate")
data<- data %>% mutate(log_rate = log(data$rate))
data$bone <- factor(data$bone, levels = names(list_of_results))
        
g <- ggplot(data, aes(x=bone, y=log_rate, color=diet, shape = diet))+geom_point()+
scale_shape_manual(values=c(21, 21,22,22,23,24,25,25))
g + theme_bw() + theme(axis.text.x = element_text(angle = 90))

#for multiple trees, julien's code

names_bones = c("Nasal", "Premaxilla (d)", "Maxilla (d)", "Jugal", "Frontal", "Parietal", "Squamosal (vault)", "Squamosal (zygomatic)", "Glenoid", "Supraoccipital", "Occipial Condyle", "Basioccipital", "Basisphenoid", "Pterygoid", "Palatine", "Maxilla (v)", "Premaxilla (v)")


res_concat <- list()
for(i in 1:length(names_bones)){
        # load results for bone "i"
        names_i <- paste(names_bones[i],".Rdata", sep = "")
        load(file=names_i)
        
        # compute statistics
        param <- sapply(sampletree, function(x) x$fit$param)
        
        res_concat[[i]] <- t(param)
        cat("Done",i,"\n")
}



# load the results (list with the 15 bones and results for each 100 trees) - relative (to aquatic) rates
load(file="list_of_results_diet.Rdata")

att = seq(0, 5*14, by=5)
spaces <- as.numeric(sapply(att, function(x) 1:5 + x))
table <- do.call("cbind",res_concat)
col <- c("red", "darkgreen","chocolate4","deepskyblue","lightgreen")
boxplot(table,boxwex=0.8, outline=FALSE, col=col, las=2, ylab="Rates (relative to carnivore)",
        xaxt = "n", at=spaces)
abline(h=1)
abline(v=seq(5,5*16,by=5), lwd=0.3)
legend("topleft", inset=.02, title="Diet",
       c("Herbivore","Omnivore","Invertivore", "Sanguinivore"), fill=col, horiz=TRUE, cex=0.8)

text(x = seq(1, 5*15,by=5)-5, y = par("usr")[3] - 1, srt = 45, adj = 0,
     labels = names_bones, xpd = TRUE)


# load the results (list with the 15 bones and results for each 100 trees) - absolute rates for the seven categories

load(file="results_abs_rates.Rdata")


att = seq(0, 9*14, by=9)
spaces <- as.numeric(sapply(att, function(x) 1:7 + x))
table <- do.call("cbind",res_concat)
col <- c("blue","darkgreen","chocolate4","deepskyblue","lightgreen","chocolate","red")
boxplot(table,boxwex=0.8, outline=FALSE, col=col, las=2, ylab="Evolutionary rates",
        xaxt = "n", at=spaces)
abline(v=seq(8,9*16,by=9), lwd=0.3)
legend("topleft", inset=.02, title="Habitat",
       c("aquatic","arboreal","burrower","semi aquatic","semi arboreal","semi burrower","terrestrial"), fill=col, horiz=TRUE, cex=0.8)

text(x = seq(1, 9*15,by=9)-9, y = par("usr")[3] - 1, srt = 45, adj = 0,
     labels = names_bones, xpd = TRUE)






## just change the order of the categories to make it easier to follow
simpleOrder <- c("aquatic","semi aquatic","arboreal","semi arboreal","burrower","semi burrower","terrestrial")
simpleOrder_nb <- c(1,4,2,5,3,6,7)
tot_grp <- sapply(seq(0,14*7,by=7), function(x) simpleOrder_nb + x)
att = seq(0, 9*14, by=9)
spaces <- as.numeric(sapply(att, function(x) 1:7 + x))
table <- do.call("cbind",res_concat)
col <- c("blue","deepskyblue","darkgreen","lightgreen","chocolate4","chocolate","red")
boxplot(table[,as.numeric(tot_grp)],boxwex=0.8, outline=FALSE, col=col, las=2, ylab="Evolutionary rates",
        xaxt = "n", at=spaces)

abline(v=seq(8,9*16,by=9), lwd=0.3)
legend("topleft", inset=.02, title="Habitat",
       simpleOrder, fill=col, horiz=TRUE, cex=0.8)

text(x = seq(1, 9*15,by=9)-9, y = par("usr")[3] - 1, srt = 45, adj = 0,
     labels = names_bones, xpd = TRUE)