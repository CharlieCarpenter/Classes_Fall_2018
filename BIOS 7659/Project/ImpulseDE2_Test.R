#############################
## Playing with ImpulseDE2 ##
#############################

## Date of package
## 2018-10-30

#BiocManager::install("ImpulseDE2")
library(ImpulseDE2); library(ComplexHeatmap); library(RUVSeq); library(tidyverse)
library(RColorBrewer)

## Has 3 reps at each time point under each condition
## 
imp.dat <- read.table("~/Documents/Classes_Fall_2018/BIOS 7659/Project/Data/GSE69822_RNA-Seq_Raw_Counts.txt",
                      sep = ",", header = T)

## Pulling out negative controls
neg.control <- imp.dat[, grep("konost", names(imp.dat))]
rownames(neg.control) <- imp.dat[,"Gene_ID"]

## PTEN knock outs and wild type both pre-incubation with DMSO for 20 min 
## before EGF stimulation
imp.set <- imp.dat[,c(41:76)] %>% as.matrix 
rownames(imp.set) <- imp.dat[,"Gene_ID"]

#######################
## RUV Normalization ##
#######################

## Filtering out genes without read counts of at least 5 in 6 spots
## Hoping for one in each time point...
## 63677 genes to 18739
filtered <- imp.set[apply(imp.set, 1, function(x) length(x[x>5]) >= 6),]

colors <- brewer.pal(3, "Set2")

x <- as.factor(rep(c("case", "control"), each = 18))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))

set_bt <- betweenLaneNormalization(set, which = "full")

biasPlot(set, "gc", log = T)

set_bt[,grep("_0", colnames(set_bt))] %>%
  plotRLE(col = colors[as.factor(rep(c("case", "control"), each = 3))])

set_bt[,grep("_15", colnames(set_bt))] %>%
  plotRLE(col = colors[as.factor(rep(c("case", "control"), each = 3))])

set_bt[,grep("_40", colnames(set_bt))] %>%
  plotRLE(col = colors[as.factor(rep(c("case", "control"), each = 3))])

set_bt[,grep("_90", colnames(set_bt))] %>%
  plotRLE(col = colors[as.factor(rep(c("case", "control"), each = 3))])

set_bt[,grep("_300", colnames(set_bt))] %>%
  plotRLE(col = colors[as.factor(rep(c("case", "control"), each = 3))],
          main = "RLE Plot of Normalized Counts\nTime = 300 min")
#legend("bottom", c("PTEN KO", "WT"), col = colors[c(1,2)], pch = 15)

## Neg control genes that are in set
neg.filt <- rownames(neg.control)[rownames(neg.control) %in% rownames(set)]

## Only have neg controls at t = 300

set_ruv <- set_bt[,grep("_300", colnames(set_bt))] %>% RUVg(neg.filt, k = 1)

plotRLE(set_ruv)

###################
## IMPULSE MODEL ##
###################

imp.annotation <- data.frame(Sample = colnames(imp.set),
                             Condition = rep(c("case", "control"), each = 18),
                             Time = rep(c(0,15,40,90,180,300), each = 3, times = 2)
                             )

objectImpulseDE2 <- runImpulseDE2(
  matCountData = counts(set_bt),
  dfAnnotation = imp.annotation,
  boolCaseCtrl = TRUE, ## Case Control analysis
  vecConfounders = NULL,
  scaNProc = 1 ## only have one processor like a chump :(
)

#save(objectImpulseDE2,file = "~/Documents/Classes_Fall_2018/BIOS 7659/Project/Output/ImpulseModel.RData")

(objectImpulseDE2$dfImpulseDE2Results$padj < 0.0001) %>% sum

objectImpulseDE2$dfImpulseDE2Results %>% arrange(padj) %>% View

############
## DESeq2 ##
############



######################
## TRAJECTORY PLOTS ##
######################

plotImpulse <- plotGenes(
  vecGeneIDs       = NULL,
  scaNTopIDs       = 10,
  objectImpulseDE2 = objectImpulseDE2,
  boolCaseCtrl     = TRUE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL
)

gl <- scale_color_discrete(labels = c("PTEN", "combined", "Wild Type"))

print(plotImpulse[[1]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[2]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[3]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[4]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[5]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[6]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[7]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[8]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[9]] + theme_bw() + labs(x = "Time (min)") + gl)
print(plotImpulse[[10]] + theme_bw() + labs(x = "Time (min)") + gl)

##############
## HEATMAPS ##
##############

caseHeat <- plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           = "case",
  boolIdentifyTransients = FALSE,
  scaQThres              = 0.00001
)

draw(caseHeat$complexHeatmapFit)

controlHeat <- plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           = "control",
  boolIdentifyTransients = FALSE,
  scaQThres              = 0.00001
)

draw(controlHeat$complexHeatmapFit)
