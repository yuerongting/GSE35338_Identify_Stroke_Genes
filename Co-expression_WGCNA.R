# rm(list = ls())

library(mouse430a2.db)
library(Biobase)

getwd()
setwd("C:/Users/yrt05/Desktop/Systems biology project/GSE35338_RAW") # \\ abs path
# setwd("C:/Users/roy20001/Desktop/Systems Biology_Stroke-20210522T190533Z-001/Systems Biology_Stroke") # \\ abs path

options(stringsAsFactors = F)

### load data
load('GSE35338_eSet.Rdata')
load('exprSet_mydata_GSE35338_rma.Rdata')
# load('RMA_norm_data.Rdata')
exprSet = exprSet_mydata_GSE35338_rma

# eSet$pCox

a = eSet[[1]]


# exprSet = exprs(a)
# exprSet = exprSet_mydata_GSE35338_rma


# exprSet = exprSet_mydata_GSE35338_rma[,1:9]      # day1
# exprSet = exprSet_mydata_GSE35338_rma[,c(19:21,25:27)]   # day3
# exprSet = exprSet_mydata_GSE35338_rma[,c(22:24,28:30)]  # day 7
dat = exprSet
pdata = pData(eSet[[1]])[c(1:9,19:21,25:27,22:24,28:30),]

head(dat)

### Second: gene annotation
ann=read.delim(file='GPL1261-56135.txt',comment.char = '#',stringsAsFactors = F)
colnames(ann)

dat=as.data.frame(dat)

# dat=dat[,rownames(pdata)]
dat$ID=rownames(dat)


# annotation of probe
id_symbol_expr<-na.omit(merge(x=dat,y=ann[c('ID','Gene.Symbol')],by='ID',all.x=T))
symbol<-lapply(id_symbol_expr$Gene.Symbol,FUN = function(x){strsplit(x,'///')[[1]][1]})
id_symbol_expr$Gene.Symbol<-as.character(symbol)
# delete probes that not appear in GPL annotation file (few)

ids=id_symbol_expr[id_symbol_expr$Gene.Symbol != 'NA',]
ids=ids[ids$ID %in% rownames(dat),]
dat=dat[ids$ID,] 

dat=dat[,-22] # delete last column "ID"


#############################  
if(T){
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

# save(dat,ann,pdata, file = "all_21_days.Rdata")
load("all_21_days.Rdata")
dat_col = length(dat[1,])

for (i in 1:dat_col) {
  dat=dat[which(dat[,i]>0),]}
dim(dat) #21746  21


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
# five = dat
# five = dat[1:16000,]
five = dat[1:round(1/2*dim(dat)[1]),]
# five = dat[1:5000,]
# save(five,file = 'five.Rdata')


### Dendrogram
five=t(five)
five=as.data.frame(five)
# datExpr=five[-23,]
datExpr=five[-22,]
datExpr_1 = datExpr
# datExpr=five
sampleTree = hclust(dist(datExpr), method = "average");
plot(sampleTree)
datExpr_tree <- hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub = "", xlab = "",
     cex.axis = 0.9, cex.main = 1.2, cex.lab = 1, cex = 0.7)




#### data trait
if(F){
  # BiocManager::install("GEOquery")
  library(GEOquery)
  a = getGEO('GSE35338')
  look_at_a = pData(a[[1]])
  group_info = look_at_a[,46]
  # look_at_a$treatment = apply(look_at_a,1,group_info)
  # test = sapply(metadata, function(x) { x$b <- rep(8,10);return(x)})
  
  metadata=pData(a[[1]])[c(1:9,19:21,25:27,22:24,28:30),][,c(2,14,46)]
  datTraits = data.frame(gsm=rownames(datExpr),
                         treatment=trimws(sapply(as.character(metadata$characteristics_ch1.4),function(x) strsplit(x,":")[[1]][2]))
                         # treatment = rownames(datExpr)
                         # treatment = metadata[,2])
  )
  # save(datTraits, file = "datTraits.RData")
}

load("datTraits.RData")
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
library(WGCNA)

options(stringsAsFactors = FALSE)




dim(datExpr)

#### soft-thresholding
powers = c(seq(1,10, by=1), seq(12, 20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft[[2]]

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit (signed R^2)",type="n",main = paste("Scale independence"), ylim = c(0,0.9));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# soft thresholding power Beta
power = sft$powerEstimate

power = 8
###







###

if(T){
  net = blockwiseModules(
  datExpr_1,
  power = 8,
  # power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 100,
  reassignThreshold = 0, mergeCutHeight = 0.30,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F,
  verbose = 3
  )
  table(net$colors)
  }

Connectivity = softConnectivity(datExpr = datExpr, power = power) -1
par(mfrow = c(1,1))
scaleFreePlot(Connectivity, main = paste("soft threshold, power = ", power), truncated = F)


####
# ConnectivityCut =  16000  # number of most connected genes that will be considered
# 
# ConnectivityRank = rank(-Connectivity)
# restConnectivity = ConnectivityRank<= ConnectivityCut
# sum(restConnectivity)
# 
# 
# # Adjacency matrix for the selected connected genes
# ADJ = adjacency(datExpr[, restConnectivity], power = power)
ADJ = adjacency(datExpr, power = power)
gc()

# TOM based on adjacency matrix
TOM = TOMsimilarity(ADJ)
dissTOM = 1 - TOM
# dissTOM = TOMdist(ADJ)
gc()


### Hierarchical clustering with TOM matrix
hierTOM = hclust(as.dist(dissTOM), method = "average")
par(mfrow = c(1,1))
plot(hierTOM, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


### dynamic tree cutting algorithm
dynamicMods = cutreeDynamic(dendro = hierTOM, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 100);
table(dynamicMods)




#### plot show
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(hierTOM, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



# Calculate eigengenes
# MEList = moduleEigengenes(datExpr[, restConnectivity], colors = dynamicColors)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
dim(cor(MEs))
MEDiss = 1-cor(MEs);  # 21 samples * 20 modules
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");


par("mar")
par(mar=c(1,5,1,1))
#### Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.2   ###  fuse modules that has 80% correlation
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")


# Call an automatic merging function
# merge = mergeCloseModules(datExpr[, restConnectivity], dynamicColors, cutHeight = MEDissThres, verbose = 3)
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs;



####  4
if(T){
  mergedColors = labels2colors(net$colors)
  table(mergedColors)
  moduleColors=mergedColors

  plotDendroAndColors(net$dendrogram[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
}

if(F){
  #????????????????????????
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






#### 5

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
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  # dev.off()
  
  
  
  ### histogram of gene significance
  Stroke = as.data.frame(design[,1]);
  names(Stroke) = "Stroke"
  y=Stroke
  # GS1=as.numeric(cor(y,datExpr[, restConnectivity], use="p"))
  test = cor(y,datExpr, use="p")
  GS1=as.numeric(cor(y,datExpr, use="p"))
  GeneSignificance=abs(GS1)
  # Next module significance is defined as average gene significance.
  ModuleSignificance=tapply(GeneSignificance,
                            moduleColors, mean, na.rm=T)
  sizeGrWindow(8,7)
  par(mfrow = c(4,1))
  
  plotModuleSignificance(GeneSignificance,moduleColors)
  
}





####  6 



Stroke = as.data.frame(design[,1]);
names(Stroke) = "Stroke"
# module = "green"
module = "brown"
if(T){
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ## ????每??模??????????皮??森????系????
  ## MEs??每??模????每????????????
  ## datExpr??每????????每???????谋???量
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneModuleMembership[1:4,1:4]
  
  
  ## 只??连??????状????只?屑???
  ## ???????欠??? Stroke ??????????量0,1??????值??
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
  # dev.off()
  
}




## step 7 
# ??要?强??踊? TOM??????WGCNA?谋?准??图
# 然?????踊???同 模?? ???????? ??图
# ??同模???牟??尉???图
# ????模?????希???要?? intramodular connectivity
if(T){
  # nGenes = ncol(datExpr)
  # nSamples = nrow(datExpr)
  # geneTree = net$dendrograms[[1]]; 
  # dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6); 
  plotTOM = dissTOM^power;
  diag(plotTOM) = NA; 
  # TOMplot(plotTOM, hierTOM, moduleColors, main = "Network heatmap plot, all genes")
  nSelect = 400
  nSelect = 5000
  nSelect = 16000
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # There??s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
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
  # dev.off()
  
  
  
  #### relations btw modules and stroke
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  ## 只??连??????状????只?屑???
  ## ???????欠??? Luminal ??????????量0,1??????值??
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
  # dev.off()
  
  
  #### plot seperately
  
  # # Plot the dendrogram
  # sizeGrWindow(8,10);
  # par(cex = 1.0)
  # ## 模???慕?????
  # # png("step7-Eigengene-dendrogram-hclust.png",width = 800,height = 600)
  # plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
  #                       plotHeatmaps = FALSE)
  # # dev.off()
  # # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  # par(cex = 1.0)
  # ## ??状??模????
  # 
  # # png("step7-Eigengene-adjacency-heatmap.png",width = 800,height = 600)
  # plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
  #                       plotDendrograms = FALSE, xLabelsAngle = 90)
  # # dev.off()
  
}









## step 8 
# ??要?枪??木???某??模???诓??幕???
if(T){
  # Select module
  module = "brown";
  # module = "green";
  # Select module probes
  probes = colnames(datExpr) ## ??????????????probe???腔???
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  head(modProbes)
  
  # ????使??WGCNA???源?????图?秃艹???
  which.module="yellow";
  dat=datExpr[,moduleColors==which.module ] 
  plotMat(t(scale(dat)),nrgcols=30,rlabels=T,
          clabels=T,rcols=which.module,
          title=which.module )
  datExpr[1:4,1:4]
  dat=t(datExpr[,moduleColors==which.module ] )
  library(pheatmap)
  pheatmap(dat ,show_colnames =F,show_rownames = F) #????些??取??来??1000?????????诘?每一??取??????????来为一???碌谋???????
  n=t(scale(t(log(dat+1)))) # 'scale'???远?log-ratio??值???泄?一??
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  group_list=datTraits$treatment
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) 
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac )
  # ???院??????目????????械???状???氐?模??????
  # ??实未?鼐筒??遣?????????????
}




#### 9 
# export the selected module to Cytoscape
if(T){
  # Recalculate topological overlap
  TOM = TOMsimilarityFromExpr(datExpr, power = power); 
  # Select module
  module = "brown";
  # module = "green";
  # Select module probes
  probes = colnames(datExpr) ## ??????????????probe???腔???
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  ## 也????取指??模???幕?????
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  ## 模????应?幕?????系??
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
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, hierTOM, file = "NetworkConstruction.RData")











chooseTopHubInEachModule(
  datExpr, 
  mergedColors, 
  omitColors = "grey", 
  power = 8, 
  type = "signed", 
  )














modNames = substring(names(mergedMEs), 3)
module = "turquoise"
moduleGenes = moduleColors==module;



geneModule = as.data.frame(cor(datExpr[, restConnectivity], mergedMEs, use = "p"))

# selected_genes = geneModule[moduleGenes, match(module, modNames)]

selected_genes_name = rownames(geneModule[moduleGenes, ])


save(selected_genes_name,file = 'selected_genes_name_WGCNA.Rdata')

load('selected_genes_name_WGCNA.Rdata')



















###############################

cor <- WGCNA::cor
# minModuleSize = 50,  mergeCutHeight = 0.25 correspond to correlation of 0.75
net = blockwiseModules(datExpr, power = 8,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,saveTOMFileBase = "MouseTOM", 
                       verbose = 3)
cor<-stats::cor

class(net)
names(net)
table(net$colors)
# 0    1    2    3    4    5 
# 88 1792 1169 1034  642  275 

# clustering result: 5 modules in total, 0 corresponds to the genes without cluster.

save(net,file='net.Rdata')

moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Dissimilarity btw modules according to eigengenes
MEDiss = 1 - cor(MEs0);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")


# threshold
MEDissThres = 0.4
abline(h = MEDissThres, col = "red")
merge_modules = mergeCloseModules(datExpr, moduleColors, cutHeight = MEDissThres, verbose
                                  = 3)

mergedColors = merge_modules$colors;
mergedMEs = merge_modules$newMEs;
plotDendroAndColors(net$dendrograms[[1]], cbind(moduleColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)











par(mfrow=c(1,1))
verboseBarplot(GeneSignificance,mergedColors,main="Module Significance ",
               col=levels(factor(mergedColors)) ,xlab="Module" ) 







