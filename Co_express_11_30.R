rm(list = ls())

library(mouse430a2.db)
library(Biobase)

getwd()
# setwd("C:/Users/yrt05/Desktop/Systems biology project/GSE35338_RAW") # \\ abs path
# setwd("C:/Users/yrt05/Desktop/Systems Biology_Stroke-20210522T190533Z-001/Systems Biology_Stroke") # \\ abs path
setwd("H:\\Users\\yrt05\\Desktop\\Systems biology project\\GSE137482_RAW")
options(stringsAsFactors = F)

### load data
# load('GSE35338_eSet.Rdata')
# load('exprSet_mydata_GSE35338_rma.Rdata')
# # load('RMA_norm_data.Rdata')

load('data_sample_probe.Rdata')
exprSet = data_sample_probe
# boxplot(exprSet,las=2)

dat = exprSet
head(dat)

### Second: gene annotation
# ann=read.delim(file='GPL1261-56135.txt',comment.char = '#',stringsAsFactors = F)
# colnames(ann)
# 
# dat=as.data.frame(dat)
# 
# # dat=dat[,rownames(pdata)]
# dat$ID=rownames(dat)
# 
# 
# # annotation of probe
# id_symbol_expr<-na.omit(merge(x=dat,y=ann[c('ID','Gene.Symbol')],by='ID',all.x=T))
# symbol<-lapply(id_symbol_expr$Gene.Symbol,FUN = function(x){strsplit(x,'///')[[1]][1]})
# id_symbol_expr$Gene.Symbol<-as.character(symbol)
# # delete probes that not appear in GPL annotation file (few)
# 
# ids=id_symbol_expr[id_symbol_expr$Gene.Symbol != 'NA',]
# ids=ids[ids$ID %in% rownames(dat),]
# dat=dat[ids$ID,] 
# 
# dat=dat[,-22] # delete last column "ID"


#############################  
if(F){
  # ids = ids[,-33]
  ids$mean=apply(dat,1,mean) # mean probe value
  # ids$median=apply(dat,1,median) # mean probe value
  ids=ids[order(ids$Gene.Symbol,ids$mean,decreasing = T),] # rank with means
  # ids=ids[order(ids$Gene.Symbol,ids$median,decreasing = T),] # rank with means
  ids=ids[!duplicated(ids$Gene.Symbol),] #delete repeated probes
  dat=dat[ids$ID,]
  rownames(dat)=ids$Gene.Symbol
  # dat[1:4,1:4] 
  save(dat,ann,pdata, file = "all_21_days.Rdata")
}







### load data, neglect above
###
###

# setwd("C:/Users/yrt05/Desktop/Systems biology project/GSE35338_RAW") # \\ abs path
# load("all_21_days.Rdata")
dat_col = length(dat[1,])


setwd("H:\\Users\\yrt05\\Desktop\\Systems biology project\\GSE137482_RAW")
load('data_sample_probe.Rdata')
dat = data_sample_probe
dat_col = length(dat[1,])
# boxplot(exprSet,las=2)


for (i in 1:dat_col) {
  dat=dat[which(dat[,i]>0),]}
dim(dat) #10613  24


############### Add 'mean' column, get the first 5000 genes according to means

dat$mean=apply(dat,1,mean)
# dat$var=apply(dat,1,var)
# dat$median=apply(dat,1,median)
descend_order = order(dat$mean,decreasing = TRUE)
# descend_order = order(dat$median,decreasing = TRUE)
# descend_order = order(dat$var,decreasing = TRUE)

# load("all_21_days.Rdata")
# dat=dat[descend_order,]
dat=dat[descend_order,]
# dat=dat[order(dat$var,decreasing = TRUE),]
# dat=dat[order(dat$median,decreasing = TRUE),]


dim(dat)
# round(1/2*dim(dat)[1])

five = dat[1:round(1/2*dim(dat)[1]),]
# five = dat[1:5000,]
# save(five,file = 'five_11_30.Rdata')



### Dendrogram
five=t(five)
five=as.data.frame(five)
# datExpr=five[-23,]
datExpr=five[-nrow(five),]
# datExpr_1 = datExpr
# datExpr=five
sampleTree = hclust(dist(datExpr), method = "average");
# plot(sampleTree)
datExpr_tree <- hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub = "", xlab = "",
     cex.axis = 0.9, cex.main = 1.2, cex.lab = 1, cex = 0.7)




#### data trait
if(F){
  # BiocManager::install("GEOquery")
  library(GEOquery)
  a = getGEO('GSE137482')
  look_at_a = pData(a[[1]])
  group_info = look_at_a[,43]
  # look_at_a$treatment = apply(look_at_a,1,group_info)
  # test = sapply(metadata, function(x) { x$b <- rep(8,10);return(x)})
  
  # metadata=pData(a[[1]])[,c(2,8,43)]
  datTraits = data.frame(gsm=rownames(datExpr),
                         treatment=group_info)

  # save(datTraits, file = "datTraits_11_30.RData")
}

load("datTraits_11_30.RData")
datTraits = datTraits

sampleNames = rownames(datExpr);
traitRows = match(sampleNames, datTraits$gsm)  
rownames(datTraits) = datTraits[traitRows, 1] 
library(stringr)
datTraits$treatment = gsub("\\s+", "_", str_trim(datTraits$treatment))
datTraits$treatment <- gsub("-", "_", str_trim(datTraits$treatment))

table(datTraits$treatment)



# Load the WGCNA package
if(F){
  install.packages("BiocManager")
  BiocManager::install("WGCNA")
}  # install 'WGCNA'

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GO.db")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("impute")

library(WGCNA)
options(stringsAsFactors = FALSE)

dim(datExpr)

#### soft-thresholding
powers = c(seq(1,10, by=1), seq(12, 20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


sft[[2]]

# Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.8;

# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit (signed R^2)",type="n",main = paste("Scale independence"), ylim = c(0,0.9));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.8,col="red")
# 
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# # R_2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]





### Soft threshold Plot
x <- sft$fitIndices$Power
y1 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
y2 <- sft$fitIndices[,5]
sizeGrWindow(9, 5)
par(mfrow = c(1,1));

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit (signed R^2)",type="n",cex.lab=1.5, cex.axis=1.5, ylim = c(0,0.9));
lines(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type = "b", pch = 19, lwd=4)
text(sft$fitIndices[,1] + 0.2, -sign(sft$fitIndices[,3])*sft$fitIndices[,2] - 0.02,labels=powers, cex=1.5, col="red")
abline(h=0.8,col="red", lwd=3)

par(new = TRUE)   
plot(x, y2, axes = FALSE, col="blue", xlab = "", ylab = "")
lines(x, y2, type = "b", pch = 2, lwd=3, col="blue")
axis(side = 4, at = pretty(range(y2)))      
mtext("Mean Connectivity", side = 4, line = 3, cex=1.5)
legend(9.5 ,1000, legend=c("Signed R^2", "Mean Connectivity"), pch = c(19, 2),
       col=c("black", "blue"), lwd = 4, lty=1:2, cex=1, text.font=2, text.width = 4, box.lty=5, box.lwd=2, border = None,)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 2)      # Grid line width







###
### soft thresholding power Beta
power = sft$powerEstimate

mean_connectivity <- sft$fitIndices[,5]
mean_connectivity[power]

median_connectivity <- sft$fitIndices[,6]
median_connectivity[power]


###
cor <- WGCNA::cor
if(T){
  net = blockwiseModules(
    datExpr,
    power = power,
    TOMType = "unsigned", minModuleSize = 100,
    reassignThreshold = 0, mergeCutHeight = 0.2,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F,
    verbose = 3 )
  
  table(net$colors)
}


test <- net$colors


# plotDendroAndColors(
#   # dendro = net$dendrograms[[1]],
#   dendro = hierTOM,
#   colors = cbind(dynamicColors,labels2colors(net$colors)),
#   dendroLabels = FALSE, hang = 0.03,
#   addGuide = TRUE, guideHang = 0.05, main = "", cex.dendroLabels = 1, cex.colorLabels = 1.1, cex.lab = 1.5, cex.axis = 1.5
# )









Connectivity = softConnectivity(datExpr = datExpr, power = power) -1
par(mfrow = c(1,1))
scaleFreePlot(Connectivity, main = paste("soft threshold, power = ", power), truncated = F)
ADJ = adjacency(datExpr, power = power)
gc()

# TOM based on adjacency matrix
TOM = TOMsimilarity(ADJ)
dissTOM = 1 - TOM

### Hierarchical clustering with TOM matrix
hierTOM = hclust(as.dist(dissTOM), method = "average")
par(mfrow = c(1,1))
plot(hierTOM, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);






merged <- labels2colors(net$colors)
before_merge <- labels2colors(net$unmergedColors)


plotDendroAndColors(hierTOM, cbind(before_merge, merged),
                    c("Original Tree", "After merging"), 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "", cex.dendroLabels = 1, cex.colorLabels = 1.1, cex.lab = 1.5, cex.axis = 1.5)






#### Cut tree show

### dynamic tree cutting algorithm
dynamicMods = cutreeDynamic(dendro = hierTOM, distM = dissTOM, cutHeight = 0.8,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 100);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
dynamicColors = labels2colors(net$colors)

# table(dynamicColors)

# Calculate eigengenes
# MEList = moduleEigengenes(datExpr[, restConnectivity], colors = dynamicColors)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)

# MEs = MEList$eigengenes
MEs = MEList$averageExpr

# Calculate dissimilarity of module eigengenes
dim(cor(MEs))
MEDiss = 1-cor(MEs);  # 21 samples * 20 modules
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");


#### Plot the tree
par("mar")
par(mar=c(1,5,1,1))
plot(METree, 
     xlab = "", sub = "", cex.lab=1.5,cex=1.5, main="")

MEDissThres = 0.2   ###  fuse modules that has 80% correlation
abline(h=MEDissThres, col = "red")








merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
table(mergedColors)
mergedMEs = merge$newMEs;


####  4
if(F){
  mergedColors = labels2colors(net$colors)
  # table(mergedColors)
  moduleColors=mergedColors
  
  plotDendroAndColors(net$dendrogram[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}





if(F){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  par(mar = c(0,5,2,0))
  plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)

  sample_colors <- numbers2colors(as.numeric(factor(datTraits$treatment)),
                                  colors = c("red","green"),signed = FALSE)

  par(mar = c(1,4,3,1),cex=0.8)
  
  # png("sample-subtype-cluster.png",width = 800,height = 600)
  plotDendroAndColors(datExpr_tree, sample_colors,
                      groupLabels = colnames(sample),
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait heatmap")
  # dev.off()
}




par(cex=0.3, mar=c(5, 8, 4, 1))
# plotDendroAndColors(hierTOM, cbind(dynamicColors, mergedColors),

test <- labels2colors(net$colors)


# plotDendroAndColors(hierTOM, cbind(dynamicColors, test),
#                     c("Original Tree", "After merging"), 
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05, main = "", cex.dendroLabels = 1, cex.colorLabels = 1.1, cex.lab = 1.5, cex.axis = 1.5)








#### 5 module significance

table(datTraits$treatment)
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  design=model.matrix(~0+ datTraits$treatment)
  # colnames(design)=levels(datTraits$treatment)
  colnames(design)= c("MCAO_induced_stroke", "sham_surgery")
  
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); 
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  # png("Module-trait-relationships.png",width = 800,height = 1200,res = 120)
  par(mar = c(8, 7, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 1,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  # dev.off()

  ### histogram of gene significance
  Stroke = as.data.frame(design[,1]);
  names(Stroke) = "Stroke"
  y=Stroke
  # GS1=as.numeric(cor(y,datExpr[, restConnectivity], use="p"))
  # test = cor(y,datExpr, use="p")
  
  GS1=as.numeric(cor(y,datExpr, use="p"))
  
  GeneSignificance=abs(GS1)
  # Next module significance is defined as average gene significance.
  ModuleSignificance=tapply(GeneSignificance,
                            moduleColors, mean, na.rm=T)
  # sizeGrWindow(8,7)
  # par(mfrow = c(4,1))
  plotModuleSignificance(GeneSignificance,moduleColors)
  # plot_module_sig <- plotModuleSignificance(GeneSignificance, moduleColors, ps=1.5,cex.lab=1.5, cex.axis=1.5, main = NULL, sub = "", cex.sub=1, cex=1.5)
}



### "modu_sig_accord_gene_sig.png"
plotModuleSignificance(GeneSignificance, moduleColors, ps=1.5,cex.lab=1.5, cex.axis=1.5,  sub = "", cex.sub=1, cex=1.5)





### Cut Tree

before_merge <- labels2colors(net$unmergedColors)


labels2colors(net$colors)

table(labels2colors(net$colors))

# moduleColors <- labels2colors(net$colors)
MEList = moduleEigengenes(datExpr, colors = before_merge)
# MEs = MEList$eigengenes
MEs = MEList$averageExpr
# Calculate dissimilarity of module eigengenes
dim(cor(MEs))
MEDiss = 1-cor(MEs);  # 21 samples * 20 modules
# Cluster module eigengenes

rownames(MEDiss) <- gsub("AE","", str_trim(rownames(MEDiss)) )

METree = hclust(as.dist(MEDiss), method = "average");

#### Plot the tree
par("mar")
par(mar=c(1,5,1,1))
plot(METree, 
     xlab = "", sub = "", cex.lab=1.5,cex=1.5, main="")

MEDissThres = 0.2   ###  fuse modules that has 80% correlation
abline(h=MEDissThres, col = "red")




















### Plot "soft_threshold.png"
x <- sft$fitIndices$Power
y1 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
y2 <- sft$fitIndices[,5]
# sizeGrWindow(9, 5)
par(mfrow = c(1,1));
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit",type="n",cex.lab=1.5, cex.axis=1.5, ylim = c(0,0.9));
lines(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type = "b", pch = 19, lwd=4)
text(sft$fitIndices[,1] + 0.2, -sign(sft$fitIndices[,3])*sft$fitIndices[,2] - 0.02,labels=powers, cex=1.5, col="red")
abline(h=0.8,col="red", lwd=3)

par(new = TRUE)   
plot(x, y2, axes = FALSE, col="blue", xlab = "", ylab = "",cex.axis=1.5, cex.lab=1.5, cex = 1.5)
lines(x, y2, type = "b", pch = 2, lwd=3, col="blue")
axis(side = 4, at = pretty(range(y2)), cex = 1.5)      
mtext("Mean Connectivity", side = 4, line = 3, cex=1.5)
legend(8 ,1000, legend=c("Signed R^2", "Mean Connectivity"), pch = c(19, 2),
       col=c("black", "blue"), lwd = 4, lty=1:2, cex=1, text.font=2,  text.width = 4, box.lty=0, box.lwd=2, border = None,)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 2)      # Grid line width







# A = adjacency(t(datExprFemale), type = "distance")
kIM = intramodularConnectivity(ADJ, moduleColors, scaleByMax = TRUE) 


####  6 
Stroke = as.data.frame(design[,1]);
names(Stroke) = "Stroke"
module = "turquoise"
if(T){
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
 
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneModuleMembership[1:4,1:4]
  
  geneTraitSignificance = as.data.frame(cor(datExpr, Stroke, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(Stroke), sep="");
  names(GSPvalue) = paste("p.GS.", names(Stroke), sep="");
  
  # module = "black"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  # png("step6-Module_membership-gene_significance.png",width = 800,height = 600)
  #sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Stroke",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")
}


## step 7 
# intramodular connectivity
if(T){
  
  plotTOM = dissTOM^power;
  diag(plotTOM) = NA; 
  # TOMplot(plotTOM, hierTOM, moduleColors, main = "Network heatmap plot, all genes")
  nSelect = 400

  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # Theres no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  sizeGrWindow(9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^power;
  diag(plotDiss) = NA;
  
  # png("step7-Network-heatmap.png",width = 800,height = 600)
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, all 16000 genes")

  
  #### relations btw modules and stroke
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  Stroke = as.data.frame(design[,1]);
  names(Stroke) = "Stroke"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, Stroke))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(8,10);
  
  par(cex = 1)
  # png("step7-Eigengene-dendrogram.png",width = 800,height = 600)
  plotEigengeneNetworks(MET, "", marDendro = c(1,4,1,3), marHeatmap = c(5,4,1,3), cex.lab = 0.8, xLabelsAngle
                        = 90)
}




## step 8 
if(T){
  # Select module
  probes = colnames(datExpr)
  inModule = (moduleColors==module);
  modProbes = probes[inModule];  ########### gene in one module
  
  
  length(modProbes)
  
  
  ############### intra-modular connectivity
  intra_connect = intramodularConnectivity(ADJ, moduleColors, scaleByMax = FALSE)  # ADJ = 16000 x 16000

  name_record = list()
  module_intra_record <- matrix(1:length(modProbes), nrow = length(modProbes), ncol = )
  dim(module_intra_record)
  module_intra_record = as.data.frame(module_intra_record);
  for (i in 1:length(modProbes)){
    # print(i)
    probe_this = (rownames(intra_connect) == modProbes[i])
    ind = which(probe_this)
    module_intra_record[i,1] = intra_connect[ind,2]   # (intramodular) connectivity
    name_record = append(name_record,(rownames(intra_connect))[probe_this] ) # (intramodular) genes
  }
  rownames(module_intra_record) = name_record
  
  desc_ord = order(module_intra_record[,"V1"],decreasing = TRUE )
  
  library(wrapr)
  # desc_ord  =  module_intra_record[orderv(module_intra_record,decreasing = TRUE), , drop = FALSE] %.>%knitr::kable(.)
  
  
  module_intra_record_rearra = as.data.frame(module_intra_record[desc_ord,])
  
  rownames(module_intra_record_rearra) = name_record[desc_ord]
  
  first_ten_percent = round(0.1 * length(module_intra_record_rearra[,1]))
  first_ten_name = (rownames(module_intra_record_rearra))[1: first_ten_percent]  # first 20% gene names
  
  first_ten_record <- matrix(1:first_ten_percent, nrow = length(first_ten_percent), ncol = )
  first_ten_record = t(first_ten_record)
  first_ten_record = as.data.frame(first_ten_record);
  rownames(first_ten_record) = first_ten_name
  
  # first_ten_record_red = first_ten_record
  # save(first_ten_record_red,file = 'first_ten_record_brown.Rdata')
  
  first_ten_record_yellow = first_ten_record
  # save(first_ten_record_yellow,file = 'first_20_Coexp.Rdata')
  
  head(first_ten_record_yellow)
}




#### 9 
# export the selected module to Cytoscape
if(T){
  # Recalculate topological overlap
  TOM = TOMsimilarityFromExpr(datExpr, power = power); 
  # Select module
  # module = "brown";
  # module = "green";
  # Select module probes
  probes = colnames(datExpr) ## ??????????????probe???ǻ???
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  ## Ҳ????ȡָ??ģ???Ļ?????
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  ## ģ????Ӧ?Ļ?????ϵ??
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.5,   # adjacency threshold for including edges in the output.
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  );
}




# Construct numerical labels corresponding to the colors
colorOrder = c(module, standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
# save(MEs, moduleLabels, moduleColors, hierTOM, file = "NetworkConstruction.RData")


modNames = substring(names(mergedMEs), 3)
# module = "turquoise"
moduleGenes = moduleColors==module;



geneModule = as.data.frame(cor(datExpr[, ], mergedMEs, use = "p"))

# selected_genes = geneModule[moduleGenes, match(module, modNames)]

selected_genes_name = rownames(geneModule[moduleGenes, ])


# save(selected_genes_name,file = 'selected_genes_name_WGCNA_1130.Rdata')

# load('selected_genes_name_WGCNA_1130.Rdata')



















###############################

# # cor <- WGCNA::cor
# # # minModuleSize = 50,  mergeCutHeight = 0.25 correspond to correlation of 0.75
# # net = blockwiseModules(datExpr, power = 8,
# #                        TOMType = "unsigned", minModuleSize = 100,
# #                        reassignThreshold = 0, mergeCutHeight = 0.25,
# #                        numericLabels = TRUE, pamRespectsDendro = FALSE,
# #                        saveTOMs = TRUE,saveTOMFileBase = "MouseTOM", 
# #                        verbose = 3)
# # cor<-stats::cor
# # 
# # class(net)
# # names(net)
# # table(net$colors)
# # # 0    1    2    3    4    5 
# # # 88 1792 1169 1034  642  275 
# # 
# # # clustering result: 5 modules in total, 0 corresponds to the genes without cluster.
# # 
# # save(net,file='net.Rdata')
# 
# moduleColors <- labels2colors(net$colors)
# # Recalculate MEs with color labels
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# # Dissimilarity btw modules according to eigengenes
# MEDiss = 1 - cor(MEs0);
# # Cluster module eigengenes
# METree = hclust(as.dist(MEDiss), method = "average");
# plot(METree,
#      main = "Clustering of module eigengenes",
#      xlab = "",
#      sub = "")
# 
# 
# # threshold
# MEDissThres = 0.4
# abline(h = MEDissThres, col = "red")
# merge_modules = mergeCloseModules(datExpr, moduleColors, cutHeight = MEDissThres, verbose
#                                   = 3)
# 
# mergedColors = merge_modules$colors;
# mergedMEs = merge_modules$newMEs;
# plotDendroAndColors(net$dendrograms[[1]], cbind(moduleColors, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)










