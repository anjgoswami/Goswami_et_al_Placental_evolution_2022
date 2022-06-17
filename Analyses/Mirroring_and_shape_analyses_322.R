library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
library(phytools)
library(doParallel)
library(RRphylo)
library(patchwork)
library(gridExtra)
library(extrafont)
#devtools::install_github("willgearty/deeptime")
library(deeptime)
#BiocManager::install("ggtree")
#library(ggtree)
#devtools::install_github("rnfelice/hot.dots")
library(hot.dots)
#devtools::install_github("rnfelice/SURGE")
library(SURGE)
library(fs)
library(here)
library(data.table)
library(strap)
#devtools::install_github("wabarr/ggphylomorpho")
#library(ggphylomorpho)
library(plotly)
library(ggedit)
library(geoscale)
library(ggrepel)
library(gginnards)
library(ggtree)
library(treeio)
source("./Scripts/utility_functions.R")
library(remotes)
library(paleotree)
source("./dtt_mod.R")
library(ggnewscale)
library(magrittr)
library(furniture)
#install for BayesTraits analyses source(here("Scripts","utility_functions.R"))

species.data<-read.csv("./Raw_Data/full_species_data.csv")
load("./Data/shape.data.322.R")
load("./Data/size.data.322.R")

shape.temp<-geomorph::two.d.array(shape.data)
write.csv(shape.temp, file = "./Data/shape.data.322.csv", quote = FALSE)

extant.data<-shape.data[,,which(species.data$Extant_Extinct=="Extant")]
extant.temp<-geomorph::two.d.array(extant.data)
write.csv(extant.temp, file = "./Data/extant.data.csv", quote = FALSE)

extant.size<-cs[which(species.data$Extant_Extinct=="Extant")]
write.csv(extant.size, file = "./Data/extant.size.csv", quote = FALSE)

name.check(sample_trees[[1]], two.d.array(extant.data))

###############

curve_table <- read_csv("./Raw_Data/placental_curves.csv")
my_curves <- create_curve_info(curve_table, n_fixed = 66)


outliers <- plotOutliers(shape.data)
#Desmocyon_matthewsi is most average skull, followed by Vulpes_pallida

###########Plot modules on skull to check ids and shape data

mod_colours<-read.csv("./Raw_Data/module.color.table.csv")
modulecolors<-as.factor(mod_colours$Color)
color.palette <- c( "#6a3d9a","dimgrey","#fb9a99",  "gold", "#009E73",  "#D55E00", "#CC79A7", "cyan2",  "#e31a1c", "#0072B2", "#b2df8a", "#E69F00",  "whitesmoke" ,  "deeppink",   "#a6cee3",   "#F0E442","blue","red","brown", "black")
levels(modulecolors)<-color.palette


open3d()
shade3d(mesh1,col="white")
spheres3d(Mirrored[,,1],col=modulecolors,radius=.001)

###########split out regions
lm_premax_d <- my_curves$Curves[which(curve_table$Module%in%c("premax_d"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_max_d <- my_curves$Curves[which(curve_table$Module%in%c("max_d"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_frontal <- my_curves$Curves[which(curve_table$Module%in%c("frontal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_nasal <- my_curves$Curves[which(curve_table$Module%in%c("nasal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_jugal <- my_curves$Curves[which(curve_table$Module%in%c("jugal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_parietal <- my_curves$Curves[which(curve_table$Module%in%c("parietal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_premax_v <- my_curves$Curves[which(curve_table$Module%in%c("premax_v"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_max_v <- my_curves$Curves[which(curve_table$Module%in%c("max_v"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_squam_v <- my_curves$Curves[which(curve_table$Module%in%c("squam_v"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_glenoid <- my_curves$Curves[which(curve_table$Module%in%c("glenoid"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_squam_z <- my_curves$Curves[which(curve_table$Module%in%c("squam_z"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_supraocc <- my_curves$Curves[which(curve_table$Module%in%c("supraocc"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_occ_cond <- my_curves$Curves[which(curve_table$Module%in%c("occ_cond", "supraocc_cond"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_basiocc <- my_curves$Curves[which(curve_table$Module%in%c("basiocc"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_basisph <- my_curves$Curves[which(curve_table$Module%in%c("basisph"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_pterygoid <- my_curves$Curves[which(curve_table$Module%in%c("pterygoid"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_palatine <- my_curves$Curves[which(curve_table$Module%in%c("palatine"))]%>%unlist(.)%>%unique(.)%>%sort(.)


######PCA

pca_results <- gm.prcomp(shape.data)
summary(pca_results)
plot(pca_results, main = "PCA")

#                          Comp1       Comp2       Comp3      Comp4
#Eigenvalues            0.01388667 0.006117949 0.003728552 0.002666065
#Proportion of Variance 0.34049634 0.150009987 0.091422817 0.065370981
#Cumulative Proportion  0.34049634 0.490506329 0.581929146 0.647300127

PC95<-pca_results$x[,1:42]
PC99<-pca_results$x[,1:112]
extremes<-pca_results$shapes
loadings<-pca_results$rotation

#PCA without pterygoid (makes almost no difference)
shape.nopter<-shape.data[-lm_pterygoid,,]
pca.nopter_results <- gm.prcomp(shape.nopter)
summary(pca.nopter_results)
plot(pca.nopter_results, main = "PCA")



# Test allometry --------------------------------------------

allometry1 <- geomorph::procD.lm(shape.data~log(cs))
summary(allometry1)
alloplot<-plotAllometry(allometry1, size=log(cs),method="CAC", pch=19)



alloshapes <- picknplot.shape(alloplot)

open3d();spheres3d(alloshapes$shapes[[1]],radius=0.002,col=modulecolours)
open3d();spheres3d(alloshapes$shapes[[1]],radius=0.002,col=mod_colours)

#                 Df    SS      MS     Rsq      F      Z     Pr(>F)
#log(Y.gpa$Csize)    1  1.9898 1.98982 0.15247 57.387 8.4825  0.001 **
#Residuals        319 11.0609 0.03467 0.84753
#Total            320 13.0508

res_shape.data<-allometry1$residuals

####PCA on residuals of shape w/o allometry

allom_pca_results <- gm.prcomp(res_shape.data)
summary(allom_pca_results)
plot(allom_pca_results, main = "PCA")

#                          Comp1        Comp2       Comp3       Comp4       Comp5
#Eigenvalues            0.009613214 0.005135412 0.00361578 0.002644728 0.001564765
#Proportion of Variance 0.278116467 0.148570764 0.10460685 0.076513683 0.045269649
#Cumulative Proportion  0.278116467 0.426687231 0.53129408 0.607807760 0.653077409

allomPC95<-allom_pca_results$x[,1:47]
allomPC99<-allom_pca_results$x[,1:120]
allom_loadings<-allom_pca_results$rotation




####import phylogeny and conduct phylogenetic comparative analyses

tree1 <- read.nexus("./Analyses/Bayes_Traits/final_BT_analyses/tree2_whole_skull_all_trees_322/tree80_85_85.nex")


nc1<-name.check(tree1,two.d.array(shape.data))

#geoscale trees
#tree1 <- drop.tip(tree1,nc1$tree_not_data)
#tree1$root.time <- max(nodeHeights(tree1))
#class(tree1)
#geoscalePhylo(tree=tree1,cex.tip=0.3)

#if more than one tree
# treelist <- lapply(all_trees,drop.tip,nc1$tree_not_data)
# for (i in 1:length(treelist)){
#   treelist[[i]]$root.time <- max(nodeHeights(treelist[[i]]))
# }
# geoscalePhylo(tree=treelist[[2]],cex.tip=0.7)
#

#tree1_root<-root(tree1,1,r=TRUE) #if tree unrooted for some reason




# Test Phylogenetic Signal --------------------------------------------
#more than one tree
physiglist <- lapply(1:3, function(k) geomorph::physignal(shape.data, treelist[[k]]))

#Observed Phylogenetic Signal (K) tree 1: 0.1802
#P-value: 0.001
#Effect Size: 6.1465

#Observed Phylogenetic Signal (K) tree 2: 0.1794
#P-value: 0.001
#Effect Size: 6.3232

#Observed Phylogenetic Signal (K) tree 3: 0.1724
#P-value: 0.001
#Effect Size: 5.9188
#Based on 1000 random permutations

# Test phylogenetic allometry --------------------------------------------
allometry.phylo <- lapply(1:6, function (k) geomorph::procD.pgls(shape.data~log(cs), phy=tree1[[k]]))
summary(allometry.phylo)
alloplot<-plotAllometry(allometry.phylo[[1]], size=log(cs),method="CAC",pch=19)

alloshapes <- picknplot.shape(alloplot)

open3d();spheres3d(alloshapes$shapes[[1]],radius=0.001,col=mod_colours)

# #if more than one tree
# > summary(allometry.phylo[[1]])
#
# Analysis of Variance, using Residual Randomization
# Permutation procedure: Randomization of null model residuals
# Number of permutations: 1000
# Estimation method: Generalized Least-Squares (via OLS projection)
# Sums of Squares and Cross-products: Type I
# Effect sizes (Z) based on F distributions
#
# Df      SS        MS     Rsq     F     Z Pr(>F)
# cs          1 0.00760 0.0076002 0.02024 6.611 2.398  0.006 **
#   Residuals 320 0.36788 0.0011496 0.97976
# Total     321 0.37549
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Call: procD.lm(f1 = shape.data ~ cs, iter = iter, seed = seed, RRPP = TRUE,
#                SS.type = SS.type, effect.type = effect.type, int.first = int.first,      Cov = Cov, data = data, print.progress = print.progress)
# > summary(allometry.phylo[[2]])
#
# Analysis of Variance, using Residual Randomization
# Permutation procedure: Randomization of null model residuals
# Number of permutations: 1000
# Estimation method: Generalized Least-Squares (via OLS projection)
# Sums of Squares and Cross-products: Type I
# Effect sizes (Z) based on F distributions
#
# Df      SS        MS     Rsq    F      Z Pr(>F)
# cs          1 0.00826 0.0082582 0.02377 7.79 2.6778  0.004 **
#   Residuals 320 0.33923 0.0010601 0.97623
# Total     321 0.34749
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Call: procD.lm(f1 = shape.data ~ cs, iter = iter, seed = seed, RRPP = TRUE,
#                SS.type = SS.type, effect.type = effect.type, int.first = int.first,      Cov = Cov, data = data, print.progress = print.progress)
# > summary(allometry.phylo[[3]])
#
# Analysis of Variance, using Residual Randomization
# Permutation procedure: Randomization of null model residuals
# Number of permutations: 1000
# Estimation method: Generalized Least-Squares (via OLS projection)
# Sums of Squares and Cross-products: Type I
# Effect sizes (Z) based on F distributions
#
# Df      SS        MS     Rsq      F      Z Pr(>F)
# cs          1 0.00712 0.0071187 0.01861 6.0677 2.3894  0.007 **
#   Residuals 320 0.37543 0.0011732 0.98139
# Total     321 0.38255
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Call: procD.lm(f1 = shape.data ~ cs, iter = iter, seed = seed, RRPP = TRUE,
#                SS.type = SS.type, effect.type = effect.type, int.first = int.first,      Cov = Cov, data = data, print.progress = print.progress)
# > summary(allometry.phylo[[4]])
#
# Analysis of Variance, using Residual Randomization
# Permutation procedure: Randomization of null model residuals
# Number of permutations: 1000
# Estimation method: Generalized Least-Squares (via OLS projection)
# Sums of Squares and Cross-products: Type I
# Effect sizes (Z) based on F distributions
#
# Df      SS        MS     Rsq      F      Z Pr(>F)
# cs          1 0.00683 0.0068316 0.02083 6.8064 2.4984  0.004 **
#   Residuals 320 0.32118 0.0010037 0.97917
# Total     321 0.32801
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Call: procD.lm(f1 = shape.data ~ cs, iter = iter, seed = seed, RRPP = TRUE,
#                SS.type = SS.type, effect.type = effect.type, int.first = int.first,      Cov = Cov, data = data, print.progress = print.progress)
# > summary(allometry.phylo[[5]])
#
# Analysis of Variance, using Residual Randomization
# Permutation procedure: Randomization of null model residuals
# Number of permutations: 1000
# Estimation method: Generalized Least-Squares (via OLS projection)
# Sums of Squares and Cross-products: Type I
# Effect sizes (Z) based on F distributions
#
# Df      SS        MS     Rsq      F      Z Pr(>F)
# cs          1 0.00815 0.0081477 0.02087 6.8214 2.3446  0.012 *
#   Residuals 320 0.38222 0.0011944 0.97913
# Total     321 0.39036
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Call: procD.lm(f1 = shape.data ~ cs, iter = iter, seed = seed, RRPP = TRUE,
#                SS.type = SS.type, effect.type = effect.type, int.first = int.first,      Cov = Cov, data = data, print.progress = print.progress)
# > summary(allometry.phylo[[6]])
#
# Analysis of Variance, using Residual Randomization
# Permutation procedure: Randomization of null model residuals
# Number of permutations: 1000
# Estimation method: Generalized Least-Squares (via OLS projection)
# Sums of Squares and Cross-products: Type I
# Effect sizes (Z) based on F distributions
#
# Df      SS        MS     Rsq      F      Z Pr(>F)
# cs          1 0.00723 0.0072279 0.02241 7.3353 2.5024  0.007 **
#   Residuals 320 0.31532 0.0009854 0.97759
# Total     321 0.32255
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Call: procD.lm(f1 = shape.data ~ cs, iter = iter, seed = seed, RRPP = TRUE,
#                SS.type = SS.type, effect.type = effect.type, int.first = int.first,      Cov = Cov, data = data, print.progress = print.progress)

pres_shape.data_1<-allometry.phylo[[1]]$pgls.residuals


#####pretty PC plots

#symbol by superorder
# Afrotheria = 8, Euarchontoglires = 5, Laurasiatheria = 13, Stem = 3, Xenarthra = 2
#fill by order
#outline black/grey by extinct

PCA1<-gm.prcomp(shape.data)
write.table(PCAsummary$PC.summary, file = "./Analyses/pcasummary.txt")

#following uses 42 componenents, 95% total variance
pcscores.95<-PCA1$x[,1:42]
pcscores<-as_tibble(pcscores.95)%>%mutate(.,Tip_Label=rownames(pcscores.95))
#pcscores<-pcscores %>% select(., Comp1, Comp2, Comp3, Comp4, Tip_Label)
pcscores<-left_join(pcscores,species.data) #join my pca results with my clade and ecological data

colortable<- read_csv("./Raw_Data/new_order_colors.csv")
names(colortable$Order_colour)<-colortable$Order

phyloseq <- c("Zalambdalestidae", "Cimolesta", "Leptictida", "Cingulata", "Pilosa", "Afrosoricida", "Macroscelidea","Tubulidentata", "Hyracoidea", "Proboscidea", "Desmostylia", "Sirenia", "Embrithopoda", "Dermoptera", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", "Eulipotyphla", "Acreodi", "Amblypoda", "Artiodactyla", "Cetacea", "Litopterna", "Notoungulata", "Astrapotheria","Perissodactyla", "Creodonta", "Carnivora", "Pholidota")
neworder <- colortable$Order_colour
neworder[order(factor(names(neworder), levels = phyloseq))]
pcscores$Order<-factor(pcscores$Order, levels = phyloseq)
pcscores<-pcscores%>%mutate(order2=Order)
#pcscores$order2[which(pcscores$Extant_Extinct=="Extinct")]<-NA


PCAsummary<-summary(PCA1)

#First plot with the legend
g <- ggplot(pcscores, aes(x=Comp1, y=Comp2, shape=Superorder,fill=Order,color=Extant_Extinct)) +
  geom_point(size=1.5,stroke=.4)+
  scale_shape_manual(values=c(21, 22,23,25, 24))+
  scale_color_manual(values = c("black","grey60"))+
  scale_fill_manual(values = neworder)+
  labs(
    x = paste0("Principal Component 1 (", signif((PCAsummary$PC.summary[2,1]*100),3), "%)"),
    y = paste0("Principal Component 2 (", signif((PCAsummary$PC.summary[2,2]*100),3),"%)")
  )+
  scale_x_continuous(breaks = seq(-0.1,0.5,by=.1))+
  scale_y_continuous(breaks = seq(-0.2,0.3,by=.1))+
  theme_bw()+
  coord_fixed(ratio=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text=element_text(size=6),
        axis.title=element_text(size=7));g
#g+guides(color = guide_legend(override.aes = list(size = 5)))
hull_clade <- pcscores%>%
  group_by(Superorder) %>%
  slice(chull(Comp1,Comp2))%>%
  rename(x=Comp1)%>%
  rename(y=Comp2)

g+new_scale_color()+
geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Superorder), alpha = .1, inherit.aes = FALSE)+
scale_color_manual(values = neworder)#+
guides(color = guide_legend(override.aes = list(size = 5)))

#g+ geom_text_repel(label = sub("_", Tip_Label))
g

figfolder <- "./Analyses/Figures/Morphospace/Whole_Skull/"
library(ggrepel)
ggsave(filename = paste0(figfolder,"morphospace no legend_temp.pdf"),
       plot = g + theme(legend.position = "none"),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace no legend plus chulls light.pdf"),
       plot = g + theme(legend.position = "none") + geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order), alpha = .1, inherit.aes = FALSE),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace no legend plus labels_temp.pdf"),
       plot = g + theme(legend.position = "none") + geom_text_repel(label = sub("_", Tip_Label)),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")

#extract the legend by itself
library(ggpubr)
leg <- as_ggplot(get_legend(g))
ggsave(filename = paste0(figfolder,"morphospace just legend 34 big.pdf"),
       plot = leg,
       device = "pdf",
       width = 20,
       height=4,
       units = "in")

######################
#Repeat for PC3 vs PC4
######################

g <- ggplot(pcscores, aes(x=Comp3, y=Comp4, shape=Superorder,fill=Order,color=Extant_Extinct)) +
  geom_point(size=1.5,stroke=.3)+
  scale_shape_manual(values=c(21, 22,23,25, 24))+
  scale_color_manual(values = c("black","grey60"))+
  scale_fill_manual(values = neworder)+
  labs(
    x = paste0("Principal Component 3 (", signif((PCAsummary$PC.summary[2,3]*100),3), "%)"),
    y = paste0("Principal Component 4 (", signif((PCAsummary$PC.summary[2,4]*100),3),"%)")
  )+
  scale_x_continuous(breaks = seq(-0.1,.5,by=.1))+
  scale_y_continuous(breaks = seq(-0.2,.3,by=.1))+
  theme_bw()+
  coord_fixed(ratio=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text=element_text(size=6),
        axis.title=element_text(size=7));g
hull_clade <- pcscores%>%
  group_by(Order) %>%
  slice(chull(Comp3,Comp4))%>%
  rename(x=Comp3)%>%
  rename(y=Comp4)
g+geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order), alpha = .2, inherit.aes = FALSE)
#g+ geom_text_repel(label = sub("_", Tip_Label))
g

ggsave(filename = paste0(figfolder,"morphospace 34 no legend_temp.pdf"),
       plot = g + theme(legend.position = "none"),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace 34 no legend plus chulls_temp.pdf"),
       plot = g + theme(legend.position = "none") + geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order), alpha = .1,inherit.aes = FALSE),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace 34 no legend plus labels_temp.pdf"),
       plot = g + theme(legend.position = "none") + geom_text_repel(label = sub("_", Tip_Label)),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")

###### REPEAT WITHOUT CETACEA::


#First plot with the legend
g <- ggplot(pcscores %>% subset(Order != "Cetacea"), aes(x=Comp1, y=Comp2, shape=Superorder,fill=Order,color=Extant_Extinct)) +
  geom_point(size=1.5,stroke=.3)+
  scale_shape_manual(values=c(21, 22,23,25, 24))+
  scale_color_manual(values = c("black","grey60"))+
  scale_fill_manual(values = neworder)+
  labs(
    x = paste0("Principal Component 1 (", signif((PCAsummary$PC.summary[2,1]*100),3), "%)"),
    y = paste0("Principal Component 2 (", signif((PCAsummary$PC.summary[2,2]*100),3),"%)")
  )+
  scale_x_continuous(breaks = seq(-0.1,.5,by=.1))+
  scale_y_continuous(breaks = seq(-0.2,.3,by=.1))+
  theme_bw()+
  coord_fixed(ratio=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text=element_text(size=6),
        axis.title=element_text(size=7));g
hull_clade <- pcscores %>% subset(Order != "Cetacea") %>%
  group_by(Order) %>%
  slice(chull(Comp1,Comp2))%>%
  rename(x=Comp1)%>%
  rename(y=Comp2)

g+new_scale_color()+
  geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order,color=Order), alpha = .1, inherit.aes = FALSE)+
  scale_color_manual(values = neworder)

g+ geom_text_repel(label = sub("_", Tip_Label))
g

figfolder <- "/Volumes/MYPASSPORT/Eutheria/Analyses/Figures/Morphospace/Whole_Skull/"
library(ggrepel)
ggsave(filename = paste0(figfolder,"morphospace no legend_temp-no whale.pdf"),
       plot = g + theme(legend.position = "none"),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace no legend plus chulls_temp-no whale.pdf"),
       plot = g + theme(legend.position = "none") + geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order), alpha = .2, inherit.aes = FALSE),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace no legend plus labels_temp-no whale.pdf"),
       plot = g + theme(legend.position = "none") + geom_text_repel(label = sub("_", Tip_Label)),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")

#extrat the legend by itself
library(ggpubr)
leg <- as_ggplot(get_legend(g))
ggsave(filename = paste0(figfolder,"morphospace just legend-no whale.pdf"),
       plot = leg,
       device = "pdf",
       width = 20,
       height=4,
       units = "in")

######################
#Repeat for PC3 vs PC4
######################

g <- ggplot(pcscores, aes(x=Comp3, y=Comp4, shape=Superorder,fill=Order,color=Extant_Extinct, lable = Tip_Label)) +
  geom_point(size=1.5,stroke=.3)+
  scale_shape_manual(values=c(21, 22,23,25, 24))+
  scale_color_manual(values = c("black","grey60"))+
  scale_fill_manual(values = neworder)+
  labs(
    x = paste0("Principal Component 3 (", signif((PCAsummary$PC.summary[2,3]*100),3), "%)"),
    y = paste0("Principal Component 4 (", signif((PCAsummary$PC.summary[2,4]*100),3),"%)")
  )+
  scale_x_continuous(breaks = seq(-0.1,.5,by=.1))+
  scale_y_continuous(breaks = seq(-0.2,.3,by=.1))+
  theme_bw()+
  coord_fixed(ratio=1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        axis.text=element_text(size=6),
        axis.title=element_text(size=7));g
hull_clade <- pcscores%>%
  group_by(Order) %>%
  slice(chull(Comp3,Comp4))%>%
  rename(x=Comp3)%>%
  rename(y=Comp4)
g+geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order), alpha = .2, inherit.aes = FALSE)
g+ geom_text_repel(label = sub("_", Tip_Label))
g

ggsave(filename = paste0(figfolder,"morphospace 34 no legend_temp-no whale.pdf"),
       plot = g + theme(legend.position = "none"),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace 34 no legend plus chulls_temp-no whale.pdf"),
       plot = g + theme(legend.position = "none") + geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order), alpha = .1,inherit.aes = FALSE),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")
ggsave(filename = paste0(figfolder,"morphospace 34 no legend plus labels_temp-no whale.pdf"),
       plot = g + theme(legend.position = "none") + geom_text_repel(label = sub("_", Tip_Label)),
       device = "pdf",
       width = 12,
       height=8,
       units = "cm")

ggplotly(g)
l <- plotly::ggplotly(g)
####extract and plot extreme shapes for figure

snapshotfolder<-"./Analyses/Figures/Morphospace/Whole_Skull/"
#snapshotfolder<-"C:/Morphospace/"

open3d()
spheres3d(PCA1$shapes$shapes.comp1$min,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC1_min_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC1_min_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC1_min_ventral.png"))

open3d()
spheres3d(PCA1$shapes$shapes.comp1$max,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC1_max_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC1_max_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC1_max_ventral.png"))

open3d()
spheres3d(PCA1$shapes$shapes.comp2$min,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC2_min_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC2_min_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC2_min_ventral.png"))

open3d()
spheres3d(PCA1$shapes$shapes.comp2$max,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC2_max_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC2_max_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC2_max_ventral.png"))

#############

open3d()
spheres3d(PCA1$shapes$shapes.comp3$min,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC3_min_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC3_min_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC3_min_ventral.png"))

open3d()
spheres3d(PCA1$shapes$shapes.comp3$max,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC3_max_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC3_max_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC3_max_ventral.png"))

open3d()
spheres3d(PCA1$shapes$shapes.comp4$min,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC4_min_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC4_min_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC4_min_ventral.png"))

open3d()
spheres3d(PCA1$shapes$shapes.comp4$max,col=modulecolors,radius=0.001)
rgl.viewpoint(fov=0)
rgl.snapshot(paste0(snapshotfolder,"PC4_max_dorsal.png"))
rgl.snapshot(paste0(snapshotfolder,"PC4_max_lateral.png"))
rgl.snapshot(paste0(snapshotfolder,"PC4_max_ventral.png"))



# adaptive landscape lookin one
ggplot(pcscores, aes(x=Comp1, y=Comp2))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  scale_fill_continuous(type = "viridis")

#for the hexagon one
ggplot(pcscores, aes(x=Comp1, y=Comp2))+ geom_hex(bins=30) +    scale_fill_continuous(type = "viridis")

#3D PCA
x<-plot_ly(x=pcscores$Comp1, y=pcscores$Comp2, z=pcscores$Comp3, type="scatter3d", mode = "markers",color=pcscores$Superorder)
htmlwidgets::saveWidget(p, "./Analyses/Figures/3dpca.html")

supersyms<-c("Afrotheria" = "circle", "Xenarthra" = "triangle-up", "Stem" = "triangle-down", "Laurasiatheria" = "diamond", "Euarchontaglires" = "square")

p <- plot_ly(x=pcscores$Comp1, y=pcscores$Comp2, z=pcscores$Comp3,
             type = 'scatter3d',
             color = pcscores$Order,  colors = neworder,
             symbol = pcscores$Superorder, symbols = supersyms,
             stroke = pcscores$Extant_Extinct, strokes = c("black","grey60"),
             mode = 'markers',
             size = c(2,5))
p

##########allometry-removed PC plots
#symbol by superorder
# Afrotheria = 8, Euarchontoglires = 5, Laurasiatheria = 13, Stem = 3, Xenarthra = 2
#colour by order
#open/closed by extinct

allomPCA<-gm.prcomp(res_shape.data)
#following uses 41 componenents, 95% total variance
pcscores.95<-allomPCA$x[,1:47]
allom_pcscores<-as_tibble(pcscores.95)%>%mutate(.,Tip_Label=rownames(pcscores.95))
#pcscores<-pcscores %>% select(., Comp1, Comp2, Comp3, Comp4, Tip_Label)
allom_pcscores<-left_join(allom_pcscores,species.data) #join my pca results with my clade and ecological data

colortable<- read_csv("./Raw_Data/new_order_colors.csv")
names(colortable$Order_colour)<-colortable$Order

phyloseq <- c("Zalambdalestidae", "Cimolesta", "Leptictida", "Cingulata", "Pilosa", "Afrosoricida", "Macroscelidea","Tubulidentata", "Hyracoidea", "Proboscidea", "Desmostylia", "Sirenia", "Embrithopoda", "Dermoptera", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", "Eulipotyphla", "Acreodi", "Amblypoda", "Artiodactyla", "Cetacea", "Litopterna", "Notoungulata", "Astrapotheria","Perissodactyla", "Creodonta", "Carnivora", "Pholidota")
neworder <- colortable$Order_colour
neworder[order(factor(names(neworder), levels = phyloseq))]
allom_pcscores$Order<-factor(allom_pcscores$Order, levels = phyloseq)
allom_pcscores<-allom_pcscores%>%mutate(order2=Order)
allom_pcscores$order2[which(allom_pcscores$Extant_Extinct=="Extinct")]<-NA

allom_g <- ggplot(allom_pcscores, aes(x=Comp1, y=Comp2, shape=Superorder,color=Order,fill=order2,label = Tip_Label)) +
  geom_point(size=4,stroke=1.5)+
  scale_shape_manual(values=c(21, 22,23,3, 24))+
  scale_color_manual(values = neworder)+
  scale_fill_manual(values = neworder)+
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    title = "Eutheria Skull Morphospace"
  )+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())
hull_clade <- allom_pcscores%>%
  group_by(Order) %>%
  slice(chull(Comp1,Comp2))%>%
  rename(x=Comp1)%>%
  rename(y=Comp2)
allom_g+geom_polygon(data = hull_clade, aes(x=x,y=y,fill=Order), alpha = .1)
allom_g+ geom_text(hjust = 0, nudge_x = 0.008)

###PhyloPCAs
### Phylogenetic PCA - PCA based on GLS-centering and projection
# This is the same as the method described by Revell (2009)
phyloPCA <- gm.prcomp(shape.data, phy = tree, GLS=TRUE)
summary(phyloPCA)
pPCAscores95_1<-phyloPCA.1$x[,1:53]

#                          Comp1        Comp2       Comp3       Comp4        Comp5
#Eigenvalues            0.000926251 0.0006787885 0.000282448 0.000253681 0.0001992151
#Proportion of Variance 0.237927856 0.1743616781 0.072552944 0.065163518 0.0511727628
#Cumulative Proportion  0.237927856 0.4122895340 0.484842478 0.550005995 0.6011787580

plot(phyloPCA, phylo = TRUE, main = "phylo PCA")

phyloPCA.2 <- gm.prcomp(shape.data, phy = treelist[[2]], GLS = TRUE)
summary(phyloPCA.2)
pPCAscores95_2<-phyloPCA.2$x[,1:53]

phyloPCA.3 <- gm.prcomp(shape.data, phy = treelist[[3]], GLS = TRUE)
summary(phyloPCA.3)
pPCAscores95_3<-phyloPCA.3$x[,1:52]

write.table(pPCAscores95_1*1000, file = "./Data/pPCscores/pPCAscores95_1.txt", quote = FALSE,
            col.names = FALSE)
write.table(pPCAscores95_2*1000, file = "./Data/pPCscores/pPCAscores95_2.txt", quote = FALSE,
            col.names = FALSE)
write.table(pPCAscores95_3*1000, file = "./Data/pPCscores/pPCAscores95_3.txt", quote = FALSE,
            col.names = FALSE)



### 3D plot with a phylogeny and time on the z-axis
plot(phyloPCA, time.plot = TRUE)
plot(phyloPCA, time.plot = TRUE, bg = neworder, phylo.par = list(tip.labels = FALSE,
                                                                 tip.txt.cex = 2, edge.color = "blue", edge.width = 2))


# Repeated with new composite tree based on Upham in 75_80 bin
tree75_80_1<-read.nexus("./Analyses/mvMorph/tree75_80_1.nex")
phyloPCA.3 <- gm.prcomp(shape.data, phy = tree75_80_1, GLS=TRUE)
summary(phyloPCA.3)



# PhyloPCA on phylogenetic allometry-corrected data - looks very similar to uncorrected pPCA
# Slight differences if use phylogenetically corrected allometry residuals, pres_shape.data instead of res_shape.data
# with pres_shape.data, differences with res_shape.data primarily in Balaenoptera and Physeter

phyloPCA.allom.1 <- gm.prcomp(res_shape.data, phy = treelist[[1]], GLS=TRUE)
summary(phyloPCA.allom.1)
plot(phyloPCA.allom.1, phylo = TRUE, main = "phylo PCA")


####pPCAs per region

phyloPCA_nasal<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_nasal,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_nasal[[1]])


phyloPCA_premax_d<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_premax_d,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_premax_d[[1]])


phyloPCA_max_d<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_max_d,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_max_d[[1]])


phyloPCA_jugal<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_jugal,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_jugal[[1]])


phyloPCA_frontal<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_frontal,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_frontal[[1]])


phyloPCA_parietal<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_parietal,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_parietal[[1]])

phyloPCA_squam_v<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_squam_v,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_squam_v[[1]])


phyloPCA_squam_z<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_squam_z,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_squam_z[[1]])


phyloPCA_glenoid<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_glenoid,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_glenoid[[1]])


phyloPCA_supraocc<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_supraocc,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_supraocc[[1]])


phyloPCA_occ_cond<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_occ_cond,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_occ_cond[[1]])


phyloPCA_basiocc<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_basiocc,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_basiocc[[1]])


phyloPCA_basisph<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_basisph,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_basisph[[1]])


phyloPCA_palatine<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_palatine,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_palatine[[1]])


phyloPCA_max_v<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_max_v,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_max_v[[1]])


phyloPCA_premax_v<-lapply(1:3, function (k) geomorph::gm.prcomp(shape.data[lm_premax_v,,], phy=treelist[[k]], GLS = TRUE))
summary(phyloPCA_premax_v[[1]])



#######################################################################################################
############################################Visualisation############################################
warpmovie3d(x, y, n, col = "green", palindrome = FALSE,
            folder = NULL, movie = "warpmovie", add = FALSE, close = TRUE,
            countbegin = 0, ask = TRUE, radius = NULL, xland = NULL,
            yland = NULL, lmcol = "black", ...)

#not included in github repository due to file size
#mesh.vulpes<-ply2mesh(file="./Vulpes_pallida.ply")

max.mesh1 <- tps3d(mesh.vulpes,shape.data[,,318],PCA1$pc.shapes$PC1max,threads=1)
min.mesh1 <- tps3d(mesh.vulpes,shape.data[,,318],pc.all$pc.shapes$PC1min,threads=1)

warpmovie3d(min.mesh1,max.mesh1 ,col="white",n=50, movie = "warpmoviepc1")## create 50 images.

#make gif, check working directory

list.files(pattern = "*.png", full.names = T) %>%
  image_read() %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=10) %>% # animates, can opt for number of loops
  image_write("warped_PC1.gif") # write to current dir


#
######DISPARITY
disp_group<-morphol.disparity(shape.data~1,groups = species.data$Order)
sorted_disp<-disp_group$Procrustes.var[order(-disp_group$Procrustes.var)]
sorted_disp

disp_superorder<-morphol.disparity(shape.data~1,groups = species.data$Superorder)
sorted_disp_sup<-disp_superorder$Procrustes.var[order(-disp_superorder$Procrustes.var)]
sorted_disp_sup

tax_div<-read.csv("./Raw_Data/clade_diversity.csv")

colortable<- read_csv("./Raw_Data/new_order_colors.csv")
names(colortable$Order_colour)<-colortable$Order
phyloseq <- c("Cingulata", "Pilosa", "Afrosoricida", "Macroscelidea", "Proboscidea", "Sirenia", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", "Eulipotyphla", "Artiodactyla", "Cetacea", "Perissodactyla", "Carnivora")
neworder <- colortable$Order_colour[which[[]]]
neworder[order(factor(names(neworder), levels = phyloseq))]
#pcscores$order2[which(pcscores$Extant_Extinct=="Extinct")]<-NA
g<-ggplot(data=subset(tax_div,Order !="Cetacea"), aes(x=log(Extant_species), y=Disparity, shape=Superorder,fill=Colour, label = Order)) +
  geom_point(size=3,stroke=0)+
  scale_fill_manual(values = tax_div$Colour)+
  scale_shape_manual(values=c(21, 22,23,25, 24))+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor=element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      axis.text=element_text(size=8))
g+ geom_text_repel()


htmlwidgets::saveWidget(l, "./Analyses/Figures/anc.pca.html")

plotly(g)

dtt.results<-dtt.paleo(phy=multi2di(treelist[[3]]), data = two.d.array(shape.data), index = 'avg.sq',
              nsim = 2, plot = TRUE)

#######to put a time scale on it
# Use the geoscale package to plot the DTT results with timescale and reverse time axis so that it goes from past to present from left to right #
geoscalePlot(ages=((dtt.results$times)*branching.times(tree)[1]), data=dtt.results$dtt, units=c("Epoch","Period"), age.lim=(rev(range(dtt.results$times)*branching.times(tree)[1])), type="l", lwd=2, label="Disparity", data.lim=c(0,2), cex.age=2, cex.ts=1.5)
# The "poly" fxn below is from the "dtt_paleo.R" file # Alternatively, you can just run the fxn #
poly <- .dtt.polygon(dtt.results$sim, dtt.results$times*branching.times(tree)[1], alpha=0.05)  # "alpha" arg is confidence interval
polygon(poly[,"x"], poly[,"y"], col=.transparency("lightgray", 0.5), border=NA)


# Per-region disparity --------------------------------------------

module_defs <- read_csv("./Raw_Data/module_data.csv")

#variance
premax_v.disp<-morphol.disparity(two.d.array(shape.data[lm_premax_v,,])~1)
premax_d.disp<-morphol.disparity(two.d.array(shape.data[lm_premax_d,,])~1)
max_v.disp<-morphol.disparity(two.d.array(shape.data[lm_max_v,,])~1)
max_d.disp<-morphol.disparity(two.d.array(shape.data[lm_max_d,,])~1)
nasal.disp<-morphol.disparity(two.d.array(shape.data[lm_nasal,,])~1)
frontal.disp<-morphol.disparity(two.d.array(shape.data[lm_frontal,,])~1)
parietal.disp<-morphol.disparity(two.d.array(shape.data[lm_parietal,,])~1)
supraocc.disp<-morphol.disparity(two.d.array(shape.data[lm_supraocc,,])~1)
pterygoid.disp<-morphol.disparity(two.d.array(shape.data[lm_pterygoid,,])~1)
palatine.disp<-morphol.disparity(two.d.array(shape.data[lm_palatine,,])~1)
jugal.disp<-morphol.disparity(two.d.array(shape.data[lm_jugal,,])~1)
jaw_joint.disp<-morphol.disparity(two.d.array(shape.data[lm_glenoid,,])~1)
squam_v.disp<-morphol.disparity(two.d.array(shape.data[lm_squam_v,,])~1)
squam_z.disp<-morphol.disparity(two.d.array(shape.data[lm_squam_z,,])~1)
occ_cond.disp<-morphol.disparity(two.d.array(shape.data[lm_occ_cond,,])~1)
basiocc.disp<-morphol.disparity(two.d.array(shape.data[lm_basiocc,,])~1)
basisph.disp<-morphol.disparity(two.d.array(shape.data[lm_basisph,,])~1)


disparityvalues<-cbind(nasal.disp, premax_d.disp, premax_v.disp, max_d.disp, max_v.disp,
                            palatine.disp,jugal.disp, frontal.disp, parietal.disp,squam_v.disp,squam_z.disp,
                            jaw_joint.disp, supraocc.disp,occ_cond.disp,
                       basiocc.disp,basisph.disp,pterygoid.disp)
colnames(disparityvalues)<-c("nasal", "premax_d", "premax_v", "max_d", "max_v", "palatine", "jugal", "frontal",	"parietal", "squamosal_v",	"squamosal_z",	"jaw joint",	"supraoccipital", "occipital condyles",	"basiocc", "basisphenoid", "pterygoid")

disparity<-as.data.frame(disparityvalues,row.names = "Disparity")

disparity[order(-disparity)]

#             nasal       max_d     frontal    premax_d    parietal    premax_v      max_v     supraoccipital    jugal  squamosal_v
#Disparity 0.005459498 0.004808817 0.004693267 0.004540862 0.003416716 0.002939377 0.00274767     0.00267503 0.002450517  0.00148596
#            squamosal_z  palatine occipital condyles    jaw joint basisphenoid      basiocc    pterygoid
#Disparity 0.001223369 0.0011954      0.0007823639     0.0006836954 0.0006489057 0.0005754962 0.0004197404





# Compare Variation Between Groups ----------------------------------------
#eco <- read_csv(here("Data", "Ecology_data.csv"))

#eco <- eco %>% filter(., Tip_label %in% row.names(pca_results$x)) %>% arrange(., Tip_label)
#names(eco$extinct)<-eco$Tip_label
#brevirostres <- rep("other", 43)
#brevirostres[which(eco$Clade %in% c("Alligatoroidea" , "Crocodylidae"))] <- "Brevirostres"
#names(brevirostres) <- eco$Tip_label
#names(eco$Diet2)<-eco$Tip_label

diet <- species.data$Fine.Diet
loc <- species.data$Broad.Locomotion
#gdf <- geomorph.data.frame(coords = shape.data[,,names(eco$extinct)], is_extinct = eco$extinct, clade = brevirostres, diet = diet, size = size.data)

gdf <- geomorph.data.frame(coords = shape.data, diet = species.data$Fine.Diet, loc = species.data$Broad.Locomotion)

disparity_comp<-morphol.disparity(f1 = coords ~ diet, groups = diet, data = gdf)
disparity_comp2<-morphol.disparity(f1 = coords ~ loc, groups = loc, data = gdf)

# Per-landmark rate and variance --------------------------------------------

rateslist <- lapply(1:6, function(k) hot.dots::per_lm_rates(shape.data, phy[[k]]))
perlmvar <- hot.dots::per_lm_variance(shape.data)

library(Morpho)
mesh.vulpes<-ply2mesh(file="./Vulpes_pallida.ply")

load("./Data/Mirrored_322.R")
load("./Data/slid.lms.R")


#hotdots variance
shade3d(mesh.vulpes,col="#E4D1C0")
spheres3d(Mirrored[c(1:754),,"Vulpes_pallida"],radius = 0.5,col=perlmvar$Variance_Colors)
title3d(main="Variance")
rgl.snapshot(filename = "./Analyses/Figures/Vulpes_hotdot_variance_lateral.png")

#hotdots rates
shade3d(mesh.vulpes,col="#E4D1C0")
spheres3d(Mirrored[c(1:754),,"Vulpes_pallida"],radius = 0.5,col=rateslist[[3]]$Rate_Colors)
title3d(main="Rate")
rgl.snapshot(filename = "./Analyses/Figures/Vulpes_hotdot_rate_lateral.png")

plot(rateslist[[3]]$Per_Lm_Rates, perlmvar$Per_Lm_Variance)

#hotscats
#utility functions:
phylo.mat<-function(x,phy){
  C<-vcv.phylo(phy,anc.nodes=FALSE)
  C<-C[rownames(x),rownames(x)]
  invC <-fast.solve(C)
  eigC <- eigen(C)
  lambda <- zapsmall(eigC$values)
  if(any(lambda == 0)){
    warning("Singular phylogenetic covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect = eigC$vectors[,1:(length(lambda))]
  D.mat <- fast.solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect))
  rownames(D.mat) <- colnames(D.mat) <- colnames(C)
  list(invC = invC, D.mat = D.mat,C = C)
}
# fast.ginv
# same as ginv, but without traps (faster)
# used in any function requiring a generalized inverse
fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
  k <- ncol(X)
  Xsvd <- La.svd(X, k, k)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  rtu <-((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
  v <-t(Xsvd$vt)[, Positive, drop = FALSE]
  v%*%rtu
}
# fast.solve
# chooses between fast.ginv or qr.solve, when det might or might not be 0
# used in any function requiring a matrix inverse where the certainty of
# singular matrices is in doubt; mostly phylo. functions
fast.solve <- function(x) if(det(x) > 1e-8) qr.solve(x) else fast.ginv(x)
# sigma.d.multi
# multiple trait multivariate evolutionary rates
# used in: compare.multi.evol.rates
sigma.d.multi<-function(x,invC,D.mat,gps,Subset){
  sig.calc<-function(x.i,invC.i,D.mat.i,Subset){
    x.i<-as.matrix(x.i)
    N<-dim(x.i)[1];p<-dim(x.i)[2]
    ones<-matrix(1,N,N)
    x.c<- x.i - crossprod(ones,invC.i)%*%x.i/sum(invC.i)
    R<-crossprod(x.c, crossprod(invC.i,x.c))/N
    if(Subset==FALSE) sigma<-sigma<-sum((D.mat.i%*%x.c)^2)/N  else
      sigma<-sum((D.mat.i%*%x.c)^2)/N/p
    return(list(sigma=sigma,R=R))
  }
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  gps.combo <- combn(ngps, 2)
  global<-sig.calc(x,invC,D.mat,Subset)
  rate.global<-global$sigma; R<-global$R
  ngps<-nlevels(gps)
  rate.gps<-sapply(1:ngps, function(j){ sig.calc(x[,g==j],
                                                 invC,D.mat,Subset)$sigma  })
  sigma.d.ratio<-max(rate.gps)/min(rate.gps)
  sigma.d.rat <- sapply(1:ncol(gps.combo), function(j){
    rates<-c(rate.gps[levels(g)==gps.combo[1,j]],rate.gps[levels(g)==gps.combo[2,j]])
    max(rates)/min(rates)
  })
  if(length(sigma.d.rat) > 1) rate.mat <- dist(matrix(0, length(rate.gps),)) else
    rate.mat = 0
  for(i in 1:length(rate.mat)) rate.mat[[i]] <- sigma.d.rat[i]
  list(sigma.d.ratio = sigma.d.ratio, rate.global = rate.global,
       rate.gps = rate.gps, sigma.d.gp.ratio = rate.mat,R = R)
}
sig.calc<-function(x.i,invC.i,D.mat.i,Subset){
  x.i<-as.matrix(x.i)
  N<-dim(x.i)[1];p<-dim(x.i)[2]
  ones<-matrix(1,N,N)
  x.c<- x.i - crossprod(ones,invC.i)%*%x.i/sum(invC.i)
  R<-crossprod(x.c, crossprod(invC.i,x.c))/N
  if(Subset==FALSE) sigma<-sigma<-sum((D.mat.i%*%x.c)^2)/N  else
    sigma<-sum((D.mat.i%*%x.c)^2)/N/p
  return(list(sigma=sigma,R=R))
}
#HOTDOTS WITH THIS DATA:
##################
x.plac<-two.d.array(shape.data)
phy.parts <- phylo.mat(x.plac, phy[[3]])
invC <- phy.parts$invC
D.mat <- phy.parts$D.mat
C = phy.parts$C
global<-sig.calc(x.plac,invC,D.mat,Subset=TRUE)
################################
#simulate:
cov.mat<-matrix(0,nrow=length(diag(global$R)),ncol=length(diag(global$R)))
diag(cov.mat)<-diag(global$R)
#phy3<-phy[[3]]
iter=100
k=3
system.time(
  simdat <- sim.char(phy=phy[[3]],par=cov.mat, nsim=iter,model="BM")
)
expectedmat<-matrix(nrow=(dim(simdat)[2]/k),ncol=iter)
for (i in 1:iter){
  temp<-arrayspecs(simdat[,,i],p=754,k=3)
  expected<-rowSums(apply(temp,c(1,2),var))
  expectedmat[,i]<-expected
}
colnames(expectedmat)<-paste(c(1:iter))
simulatedvariance<-tibble::as_tibble(cbind(rates.vector,expectedmat))
simvarclean<-gather(simulatedvariance, "simnumber", "variance", 2:(iter+1))
# Create prediction interval
m1<-lm(variance~rates.vector, data=simvarclean)
newx <- seq(min(simvarclean$rates.vector), max(simvarclean$rates.vector), length.out = 100)
pred_interval <- predict(m1, newdata=data.frame(rates.vector=newx), interval="prediction", level = 0.95)
pred_interval <- as.data.frame(pred_interval)
pred_interval$rates = newx
#ggthemr_reset()

rates.and.vars


p <-ggplot(rates.and.vars, aes(x=rates.vector, y=variances, color=as.factor(modulecolors[1:754])))
p<- p + geom_point(size=1.8)+
                      scale_color_manual(values=color.palette,
                     labels = c("Nasal","Premaxilla (dorsal)","Maxilla (dorsal)","Jugal","Frontal","Parietal","Glenoid","Squamosal_Zyg","Squamosal_vault","Supraoccciput","Occipital_Condyles", "Basioccipital","Basisphenoid","Palatine","Maxilla (ventral)", "Premaxilla (ventral)"),
                     name = "Group")
p + theme_light()+ stat_smooth(method="lm",
                               se=TRUE,
                               fill="grey",
                               formula=y ~ x,
                               colour="red")+
  stat_smooth(mapping=aes(x=rates.vector,y=variance),
              linetype=2,
              se=TRUE,
              data=simvarclean,
              method=lm,fill="grey",
              formula=y ~ x, inherit.aes = FALSE)+
  geom_ribbon(mapping=aes(x=rates,ymin = lwr, ymax = upr),
              data=pred_interval,
              fill = "grey",
              alpha = 0.2,
              inherit.aes = FALSE)+
  theme(aspect.ratio = 1,legend.position="bottom")+
  labs(x="Evolutionary Rate (Sigma)",y="Procrustes Variance")
#
model1 <- lm(variances~rates.vector, data= rates.and.vars)
summary(model1)

# EMMLIi Test of Modularity --------------------------------------------------------
#import module definitions
module_defs <- read_csv("./Raw_Data/module_data.csv")
tree1 <- read.nexus("./Analyses/mvMorph/t1_list_of_trees.nex")

coords.2d <- two.d.array(shape.data)
treelist<-tree1

contrasts_tree1 <- apply(coords.2d, 2, FUN=pic, phy=multi2di(treelist[[1]])) %>% arrayspecs(., p=dim(shape.data)[1], k = 3)
contrasts_tree2 <- apply(coords.2d, 2, FUN=pic,  phy=multi2di(treelist[[2]])) %>% arrayspecs(., p=dim(shape.data)[1], k = 3)
contrasts_tree3 <- apply(coords.2d, 2, FUN=pic,  phy=multi2di(treelist[[3]])) %>% arrayspecs(., p=dim(shape.data)[1], k = 3)
contrasts_tree4 <- apply(coords.2d, 2, FUN=pic,  phy=multi2di(treelist[[4]])) %>% arrayspecs(., p=dim(shape.data)[1], k = 3)
contrasts_tree5 <- apply(coords.2d, 2, FUN=pic,  phy=multi2di(treelist[[5]])) %>% arrayspecs(., p=dim(shape.data)[1], k = 3)
contrasts_tree6 <- apply(coords.2d, 2, FUN=pic,  phy=multi2di(treelist[[6]])) %>% arrayspecs(., p=dim(shape.data)[1], k = 3)


corrmat_contrasts_t1<-dotcorr(contrasts_tree1)
corrmat_contrasts_t2<-dotcorr(contrasts_tree2)
corrmat_contrasts_t3<-dotcorr(contrasts_tree3)
corrmat_contrasts_t4<-dotcorr(contrasts_tree4)
corrmat_contrasts_t5<-dotcorr(contrasts_tree5)
corrmat_contrasts_t6<-dotcorr(contrasts_tree6)



EMMLI_PIC_tree1 <- EMMLi(corr = corrmat_contrasts_t1, N_sample = 321, mod = as.data.frame(module_defs[,-c(1:2,4)]))
EMMLI_PIC_tree2 <- EMMLi(corr = corrmat_contrasts_t2, N_sample = 321, mod = as.data.frame(module_defs[,-c(1:2,4)]))
EMMLI_PIC_tree3 <- EMMLi(corr = corrmat_contrasts_t3, N_sample = 321, mod = as.data.frame(module_defs[,-c(1:2,4)]))
EMMLI_PIC_tree4 <- EMMLi(corr = corrmat_contrasts_t4, N_sample = 321, mod = as.data.frame(module_defs[,-c(1:2,4)]))
EMMLI_PIC_tree5 <- EMMLi(corr = corrmat_contrasts_t5, N_sample = 321, mod = as.data.frame(module_defs[,-c(1:2,4)]))
EMMLI_PIC_tree6 <- EMMLi(corr = corrmat_contrasts_t6, N_sample = 321, mod = as.data.frame(module_defs[,-c(1:2,4)]))

source("./Scripts/plotnetwork2.R")
#mod_names<-
layout<-as.matrix(read.csv("./Raw_data/layout_modules.csv"))
plotNetwork2(EMMLI_PIC_tree1$rho$`Region.sep.Mod + sep.between`[2,], module_names = NULL, linecolour = "#56B4E9",title = NULL, layout = layout, minimum = 0)

nmodules=17
#rho list directs the code to the results. this will be the "rhos_best" thing for subsampled emmli
rholist<-t(EMMLI_PIC_tree1$rho$`Region.sep.Mod + sep.between`)
words<-strsplit(rownames(rholist), split = " ")
corrmat_new<-matrix(data = NA, nrow = nmodules, ncol = nmodules)
for (i in 1:(length(words)-1)){
  if (length(words[[i]]) == 2){
    corrmat_new[as.numeric(words[[i]][2]),as.numeric(words[[i]][2])] <- rholist[i,2]
  }
  if (length(words[[i]]) == 3){
    corrmat_new[as.numeric(words[[i]][1]),as.numeric(words[[i]][3])] <- rholist[i,2]
    corrmat_new[as.numeric(words[[i]][3]),as.numeric(words[[i]][1])] <- rholist[i,2]
  }
}
modnames<-c("Nasal","Premaxilla (d)","Maxilla (d)","Jugal","Frontal","Parietal","Squamosal (v)","Squamosal (z)","Glenoid","Supraoccipital","Occipital Condyle","Basioccipital","Basisphenoid","Palatine","Pterygoid","Maxilla (v)","Premaxilla (v)")
rownames(corrmat_new)<-colnames(corrmat_new)<-modnames
layout<-as.matrix(read.csv("./Raw_data/layout_modules.csv"))
within<-diag(corrmat_new)
between<-corrmat_new
qgraph(between, shape="circle", posCol="#056571",  labels = rownames(corrmat_new), label.cex=2, label.scale.equal=TRUE, vsize=(within)*10, diag = FALSE, title="phyloEMMLi_tree1_70_75",layout=layout, minimum=.4)


################################################
# Load and Plot Bayes Traits Results
#note that these BayesTraits output files are not incldued in the github repository dude to their large size and can be reproduced using the code in
#./Analyses/BayesTraits

color3<-colorRampPalette(c("#0c2c84","#225ea8","#31a354","#ffff00","#fe9929","#fc4e2a","red","darkred"))

library(BTRTools)

#Tree 1
##Check convergence of likelihood:
#test1 = tracePlots(file="C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_1_run_A.txt.VarRates.txt", plot = FALSE)
#test2 = tracePlots(file="C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_1_run_B.txt.VarRates.txt", plot = FALSE)
#my_list_of_chains = mcmc.list(list(test1[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")], test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))
#gelman.diag(my_list_of_chains)
#
#tree_1_BTraits<-BTRTools::rjpp(rjlog = "C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_1_run_A.txt.VarRates.txt",
#                                          rjtrees = "C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_1_run_A.txt.Output.trees",
#                                          tree = treelist[[1]]) #this is your time scaled tree that was used to input into bayestraits
#tree_1_BTraits.output <- tree_1_BTraits$data
#write.csv(tree_1_BTraits.output, "./Analyses/Bayes_Traits/tree1_BT.output.csv")
#
#
#tree_2_BTraits<-BTRTools::rjpp(rjlog = "C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_2_run_A.txt.VarRates.txt",
#                               rjtrees = "C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_2_run_A.txt.Output.trees",
#                               tree = treelist[[2]]) #this is your time scaled tree that was used to input into bayestraits
#tree_2_BTraits.output <- tree_2_BTraits$data
#write.csv(tree_2_BTraits.output, "./Analyses/Bayes_Traits/tree2_BT.output.csv")
#
#
#tree_3_BTraits<-BTRTools::rjpp(rjlog = "C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_3_run_A.txt.VarRates.txt",
#                               rjtrees = "C:/BayesTraits/btresults.tar/pPCscores/pPCAscores.wholeskull.95_3_run_A.txt.Output.trees",
#                               tree = treelist[[3]]) #this is your time scaled tree that was used to input into bayestraits
#tree_3_BTraits.output <- tree_3_BTraits$data
#write.csv(tree_3_BTraits.output, "./Analyses/Bayes_Traits/tree3_BT.output.csv")
#
#palatex<-c( "#D55E00",   "#009E73",  "#56B4E9", "brown", "black")
#tree_1_w_data <- add_rjpp_to_tree(tree_1_BTraits)
#tree_1_w_data <- full_join(tree_1_w_data, mutate(species.data, label = species.data$Tip_Label), by = "label")
#threshold <- .15 # the minimum posterior probability you want to plot a symbol for
#p<-ggtree(tree_1_w_data, aes(color = log(meanRate)), layout = "fan", size=1)+
#  scale_colour_gradientn(colours = color3(100))+
#  #geom_tippoint(aes(fill=Broad.Diet, x=x+4),color="black", shape=23)+
#  #scale_fill_manual(breaks = c("Carnivore","Herbivore","Invertivore", "Omnivore", "Sanguinivore"), values = palatex)+
#  #theme(legend.position="top")+
#  theme(legend.position=c(.15,.95),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
#  scale_size(range = c(1,2))+
#  labs(title="Tree 1",
#       color="log(Rate)")+
#  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
#  scale_size(range = c(1,2))+
#  geom_tiplab(label= sub("_", " ",tree_1_w_data@phylo$tip.label), size=3, offset=0, color = "black", fontface="italic")+
#  coord_cartesian(xlim = c(-90, 20), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
#                  ylim = c(-2, 325), #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
#                  expand = FALSE) +
#  scale_x_continuous(breaks=-periods$max_age[c(1:5)], labels=periods$max_age[c(1:5)]) +
#  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
#p <- revts(p);p
#ptree1 <-  gggeo_scale(p, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)
#ptree1

#p2<-p+geom_strip(taxa1 ="Crocodylus_mindorensis", taxa2 = "Prodiplocynodon_langi",  label="Crocodyloidea", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial")+
#  geom_strip(taxa1 ="Caiman_yacare", taxa2 =  "Leidyosuchus_canadensis",  label="Alligatoroidea", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial", color="grey")+
#  geom_strip(taxa1 ="Simosuchus_clarki" , taxa2 =  "Araripesuchus_wegeneri" ,  label="Notosuchia", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial", color = "grey")+
#  geom_strip(taxa1 ="Sarcosuchus_imperator", taxa2 =  "Pholidosaurus_sp",  label="Pholidosauridae", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial")+
#  geom_strip(taxa1 ="Pelagosaurus_typus" , taxa2 =  "Cricosaurus"  ,  label="Thalattosuchia", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial")


#ptree1b <-  gggeo_scale(p2, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)
#ptree1b

#ggsave(filename = "./Analyses/Figures/ws_rate_tree1.pdf",
#       plot = ptree1, device = cairo_pdf,
#       width = 11, height = 11, units = "cm")
#
#tree_2_w_data <- add_rjpp_to_tree(tree_2_BTraits)
#tree_2_w_data <- full_join(tree_2_w_data, mutate(species.data, label = species.data$Tip_Label), by = "label")
#threshold <- .15 # the minimum posterior probability you want to plot a symbol for
#p<-ggtree(tree_2_w_data, aes(color = log(meanRate)), size=1)+
#  scale_colour_gradientn(colours = color3(100))+
#  #geom_tippoint(aes(fill=Broad.Diet, x=x+4),color="black", shape=23)+
#  #scale_fill_manual(breaks = c("Carnivore","Herbivore","Invertivore", "Omnivore", "Sanguinivore"), values = palatex)+
#  #theme(legend.position="top")+
#  theme(legend.position=c(.15,.95),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
#  scale_size(range = c(1,2))+
#  labs(title="Tree 2",
#       color="log(Rate)")+
#  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
#  scale_size(range = c(1,2))+
#  geom_tiplab(label= sub("_", " ",tree_2_w_data@phylo$tip.label), size=3, offset=0, color = "black", fontface="italic")+
#  coord_cartesian(xlim = c(-90, 20), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
#                  ylim = c(-2, 325), #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
#                  expand = FALSE) +
#  scale_x_continuous(breaks=-periods$max_age[c(1:5)], labels=periods$max_age[c(1:5)]) +
#  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
#p <- revts(p);p
#ptree2 <-  gggeo_scale(p, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)
#ptree2
#
#tree_3_w_data <- add_rjpp_to_tree(tree_3_BTraits)
#tree_3_w_data <- full_join(tree_3_w_data, mutate(species.data, label = species.data$Tip_Label), by = "label")
#threshold <- .15 # the minimum posterior probability you want to plot a symbol for
#p<-ggtree(tree_3_w_data, aes(color = log(meanRate)), size=1)+
#  scale_colour_gradientn(colours = color3(100))+
#  #geom_tippoint(aes(fill=Broad.Diet, x=x+4),color="black", shape=23)+
#  #scale_fill_manual(breaks = c("Carnivore","Herbivore","Invertivore", "Omnivore", "Sanguinivore"), values = palatex)+
#  #theme(legend.position="top")+
#  theme(legend.position=c(.15,.95),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
#  scale_size(range = c(1,2))+
#  labs(title="Tree 3",
#       color="log(Rate)")+
#  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
#  scale_size(range = c(1,2))+
#  geom_tiplab(label= sub("_", " ",tree_3_w_data@phylo$tip.label), size=3, offset=0, color = "black", fontface="italic")+
#  coord_cartesian(xlim = c(-90, 20), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
#                  ylim = c(-2, 325), #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
#                  expand = FALSE) +
#  scale_x_continuous(breaks=-epochs$max_age[c(1:8)], labels=epochs$max_age[c(1:8)]) +
#  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
#p <- revts(p);p
#ptree3 <-  gggeo_scale(p, dat = "epochs", neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"),  size=3)
#ptree3



#Plot the diets on the tree
diet_pal<-c("royalblue2","red2","springgreen3","tan2","gold2","deeppink","tan4")
loc_pal<-c("royalblue2","springgreen4","tan4","grey","turquoise","springgreen2","tan2", "brown","skyblue")
hab_pal<-c("white","turquoise4", "darkgreen", "tan1", "springgreen2","gold2")



treedataX <- full_join(tree1, mutate(species.data, label = species.data$Tip_Label), by = "label")


p<-ggtree(treedataX, layout = "circular")
p + geom_tippoint(aes(fill=Fine_Diet, x=x+4),color="black", shape=23) +
  scale_fill_manual(breaks = c("Bulk invertivore", "Carnivore","Herbivore","Invertivore", "Omnivore", "Piscivore", "Social insectivore"), values = diet_pal)+
  new_scale("fill") +

  geom_tippoint(aes(fill=Broad_Locomotion, x=x+8),color="black", shape=24) +
  scale_fill_manual(breaks = c("Aquatic", "Arboreal", "Fossorial",  "Saxicolous", "Semi-aquatic",  "Semi-arboreal", "Semi-fossorial", "Terrestrial", "Volant"), values = loc_pal)+

  new_scale("fill") +
  geom_tippoint(aes(fill=Broad_Habitat, x=x+12),color="black", shape=21) +
  scale_fill_manual(breaks = c( "Extinct", "Aquatic", "Closed", "Desert", "Mixed_terrestrial", "Open"), values = hab_pal)


soc_pal <- c("springgreen2","red2")
devo_pal <- c("salmon2","turquoise2")
active_pal <- c("springgreen2","royalblue2","orange2", "blue4")

p<-ggtree(treedataX, layout = "circular")
p + geom_tippoint(aes(fill=Social, x=x+4),color="black", shape=23) +
  scale_fill_manual(breaks = c("Social","Solitary","","unknown"), values = c(soc_pal,"white","white"))+
  new_scale("fill")

  geom_tippoint(aes(fill=Broad_Locomotion, x=x+8),color="black", shape=24) +
  scale_fill_manual(breaks = c("Aquatic", "Arboreal", "Fossorial",  "Saxicolous", "Semi-aquatic",  "Semi-arboreal", "Semi-fossorial", "Terrestrial", "Volant"), values = loc_pal)+

  new_scale("fill") +
  geom_tippoint(aes(fill=Broad_Habitat, x=x+12),color="black", shape=21) +
  scale_fill_manual(breaks = c( "Extinct", "Aquatic", "Closed", "Desert", "Mixed_terrestrial", "Open"), values = hab_pal)

##############
#per-group skull means
func1<-function(group_name){
taxa<-species.data %>% filter(Fine_Diet==group_name) %>% pull(Tip_Label)
skulls<-shape.data[,,taxa]
ancs<-Rphylopars::anc.recon(two.d.array(skulls), keep.tip(tree1,taxa))
mean<-arrayspecs(ancs, k=3, p=dim(shape.data)[1])[,,1]
return(mean)
}
carnivoreanc<-func1("Carnivore")
clear3d();spheres3d(carnivoreanc,radius = 0.001,col=modulecolors)
rgl.snapshot(paste0("./plots/", "Carnivore mean.png"))
herbivoreanc<-func1("Herbivore")
clear3d();spheres3d(herbivoreanc,radius = 0.001,col=modulecolors)
rgl.snapshot(paste0("./plots/", "Herbivore mean.png"))
invertanc<-func1("Invertivore")
clear3d();spheres3d(invertanc,radius = 0.001,col=modulecolors)
rgl.snapshot(paste0("./plots/", "Invertivore mean.png"))
Sinsectanc<-func1("Social insectivore")
clear3d();spheres3d(Sinsectanc,radius = 0.001,col=modulecolors)
rgl.snapshot(paste0("./plots/", "social insectivore mean.png"))
Piscivoreanc<-func1("Piscivore")
clear3d();spheres3d(Piscivoreanc,radius = 0.001,col=modulecolors)
rgl.snapshot(paste0("./plots/", "Piscivore mean.png"))
Omnivoreanc<-func1("Omnivore")
clear3d();spheres3d(Omnivoreanc,radius = 0.001,col=modulecolors)
rgl.snapshot(paste0("./plots/", "Omnivore mean.png"))
BULKanc<-func1("Bulk invertivore")
clear3d();spheres3d(BULKanc,radius = 0.001,col=modulecolors)
rgl.snapshot(paste0("./plots/", "Bulk invertivore mean.png"))



########
