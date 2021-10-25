getwd()
setwd('c:\\Users\\roy20001\\Desktop\\Systems biology project\\R script\\test')  # \\ abs path
# setwd("C:/Users/roy20001/Desktop/Systems Biology_Stroke-20210522T190533Z-001/Systems Biology_Stroke") # \\ abs path

#download data from GEO database, get file ".csv"
#function source: https://github.com/jmzeng1314/humanid/blob/master/R/downGSE.R
if(T){
  downGSE <- function(studyID = "GSE35338", destdir = ".") {
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
  
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
  
  
  save(eSet,file = 'GSE35338_eSet.Rdata')
  write.csv(exprSet, paste0(studyID, "_exprSet.csv"))  
  write.csv(pdata, paste0(studyID, "_metadata.csv"))  # classification of data
  return(eSet)
  }}

# BiocManager::install("GEOquery")
library(GEOquery)
eSet <- downGSE('GSE35338')

orig_data = exprs(eSet[[1]])

### load data
load('GSE35338_eSet.Rdata')


### download database
# BiocManager::install("mouse430a2.db")
library(mouse430a2.db)


### get group names
library(Biobase)
# library(GEOquery)
pdata = pData(eSet[[1]])




# pdata = pdata[1:30,]       #  all 21 days
# pdata = pdata[1:9,]       # day1
# pdata = pdata[c(19:21,25:27),]# day3
pdata = pdata[c(22:24,28:30),]  # day7
# table(ee)
dim(pdata)
# group_list=as.character(pdata[c(1:9,19:30),46])
# group_list=as.character(pdata[c(1:9),46])  # day 1
# group_list=as.character(pdata[c(1:6),46])  # day 3
group_list=as.character(pdata[c(1:6),46])  # day 7


### Normalization
# Import into R from local storage
# setwd("C:/Users/roy20001/Desktop/Systems Biology_Stroke-20210522T190533Z-001/data_cel")
setwd('C:/Users/roy20001/Desktop/Systems biology project/GSE35338_RAW')  # \\ abs path
# install.packages("BiocManager")
# BiocManager::install("affy")
library(affy)
mydata_GSE35338<-ReadAffy()

# data <- exprs(mydata_GSE35338[,c(1:9,19:21,25:27,22:24,28:30)][1])


### RMA normalization. Normalize to raw data (mydata_GSE35338)
if(T){
  
  # mydata_GSE35338_rma<-rma(mydata_GSE35338[,c(1:9,19:21,25:27,22:24,28:30)]) # all 21 days
  
  # mydata_GSE35338_rma<-rma(mydata_GSE35338[,1:9])      # day1
  # mydata_GSE35338_rma<-rma(mydata_GSE35338[,c(19:21,25:27)])   # day3
  mydata_GSE35338_rma<-rma(mydata_GSE35338[,c(22:24,28:30)])  # day 7
  
  exprSet_mydata_GSE35338_rma<-exprs(mydata_GSE35338_rma)
}

# save(exprSet_mydata_GSE35338_rma,file = 'exprSet_mydata_GSE35338_rma.Rdata')
# save(exprSet_mydata_GSE35338_rma,file = 'exprSet_mydata_GSE35338_rma_day1.Rdata')
# save(exprSet_mydata_GSE35338_rma,file = 'exprSet_mydata_GSE35338_rma_day3.Rdata')
save(exprSet_mydata_GSE35338_rma,file = 'exprSet_mydata_GSE35338_rma_day7.Rdata')

# load('exprSet_mydata_GSE35338_rma_day1.Rdata')
# load('exprSet_mydata_GSE35338_rma_day3.Rdata')
# load('exprSet_mydata_GSE35338_rma_day7.Rdata')

### Match probe ID with gene. (Preprocess to data)
exprSet = exprSet_mydata_GSE35338_rma

ids = toTable(mouse430a2SYMBOL) 
# ids = toTable(mouse430a2SYMBOL)   # regular gene for coding protein has size around 20k (obs).

length(unique(ids$symbol))      # number of genes
tail(sort(table(ids$symbol)))   # number of probes for each genes 
table(sort(table(ids$symbol)))  #  most genes are set with 1 probe (not original, is after processed by author)
plot(table(sort(table(ids$symbol))))

### filter probes
table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(ids)
dim(exprSet)

# same sequence of Series matrix and probe
ids=ids[match(rownames(exprSet),ids$probe_id ),]
head(ids)
exprSet[1:5,1:5] # match
dim(ids)

# classify according to 'Symbol': get single max probe for each gene 
# 'by' is classification
# each probe has mean, get the probe that has max mean value
if(F){
  tmp = by(exprSet,ids$symbol, 
           function(x) rownames(x)[which.max(rowMeans(x))] ) 
  dim(tmp)
  
  probes = as.character(tmp)
  dim(exprSet)
  
  # filter repeated probes
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  dim(exprSet)
  
  
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  exprSet[1:5,1:5]
  dim(exprSet)
}

identical(ids$probe_id,rownames(exprSet))
dat=exprSet


dim(ids)
dim(dat)


### Same function, replace probe names with genes names
if(T){
  ids$median=apply(dat,1,median)
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  ids=ids[!duplicated(ids$symbol),]
  dat=dat[ids$probe_id,]
  rownames(dat)=ids$symbol
  dat[1:4,1:4]
  dim(dat)
  }



### After matching, plot and check expression values
exprSet=dat
exprSet['Ighm',]
boxplot(exprSet[,1])
exprSet[1:5, 1:5]

### colorful boxplot
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(exprSet, col = cols,main="expression value",las=2)


### label the groups
library(reshape2)

dim(exprSet)
exprSet_L = melt(exprSet)


colnames(exprSet_L) = c('probe','sample','value')
exprSet_L$group = rep(group_list,each = nrow(exprSet))
head(exprSet_L)


### Clustering with sample distances
# colnames(exprSet) = paste(group_list,1:21,sep = '')    # all days

# colnames(exprSet) = paste(group_list,1:9,sep = '')    # day1
# colnames(exprSet) = paste(group_list,c(19:21,25:27),sep = '')  # day3
colnames(exprSet) = paste(group_list,c(22:24,28:30),sep = '')  # day7





### valid group name ###
library(stringr)
group_list_replace <- gsub("\\s+", "_", str_trim(group_list))
group_list_replace <- gsub("-", "_", str_trim(group_list_replace))
design_replace <- gsub("\\s+", "_", str_trim(group_list))
design_replace <- gsub("-", "_", str_trim(design_replace))

### label for grouping probes
# BiocManager::install("limma")
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list_replace))
colnames(design)=levels(factor(group_list_replace))
rownames(design)=colnames(exprSet)





## not change 'exprSet', use new variable 'exprSet_test'
exprSet_test = exprSet
colnames(exprSet_test) = gsub("\\s+", "_", str_trim(colnames(exprSet_test)))
colnames(exprSet_test) <- gsub("-", "_", str_trim(colnames(exprSet_test)))

rownames(design)=colnames(exprSet_test)


#################### Making rules for comparison (choose two of the groups) ####################
table(group_list_replace)


library('limma')
##### 'makecontrast' search DEG according to different labels
contrast.matrix<-makeContrasts("MCAO_induced_stroke - sham_surgery",levels = design)
contrast.matrix 
dim(contrast.matrix)
dim(exprSet_test)
head(exprSet_test)

# row1_mean = mean(exprSet_test[1,])


##### final procedure in DEG ########
##step1: (label) fits multiple linear models by weighted or generalized least squares
fit <- lmFit(exprSet_test,design)
dim(fit)
head(fit[[1]])
dim(design)
dim(exprSet_test)
colnames(fit)
rownames(fit)


##step2:  compute estimated coefficients and standard errors
fit2 <- contrasts.fit(fit, contrast.matrix)

fit2$coefficients  # estimated coefficients for each contrast for each probe.
dim(fit2$coefficients)
fit3 <- eBayes(fit2)
dim(fit3)
##step3
tempOutput = topTable(fit3, coef=1, n=Inf)
dim(tempOutput)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)


boxplot(nrDEG[,1])

### select DEG result ####
logFC_threshold = 0.30
find_DEG_up = nrDEG[which(nrDEG[1]> logFC_threshold & nrDEG[4]<0.05),]
find_DEG_down = nrDEG[which(nrDEG[1]< (-logFC_threshold) & nrDEG[4]<0.05),]

dim(find_DEG_up)
dim(find_DEG_down)
# find_DEG[,1]

# # ### save DEG
# write.table(nrDEG,"day1_DEG.txt")
# # # write.table(nrDEG,"day3_DEG.txt")
# # # write.table(nrDEG,"day7_DEG.txt")


### up-regulated
# save(find_DEG_up,file = 'day1_DEG_up.Rdata')
# save(find_DEG_up,file = 'day3_DEG_up.Rdata')
# save(find_DEG_up,file = 'day7_DEG_up.Rdata')
### down-regulated
# save(find_DEG_down,file = 'day1_DEG_down.Rdata')
# save(find_DEG_down,file = 'day3_DEG_down.Rdata')
# save(find_DEG_down,file = 'day7_DEG_down.Rdata')












### load data
load('day1_DEG_up.Rdata')
DEG_day1_up = find_DEG_up

load('day3_DEG_up.Rdata')
DEG_day3_up = find_DEG_up

load('day7_DEG_up.Rdata')
DEG_day7_up = find_DEG_up

load('day1_DEG_down.Rdata')
DEG_day1_down = find_DEG_down
load('day3_DEG_down.Rdata')
DEG_day3_down = find_DEG_down
load('day7_DEG_down.Rdata')
DEG_day7_down = find_DEG_down


# load('all_21_DEG_up_FC2.Rdata')
# load('all_21_DEG_down_FC2.Rdata')


### gene list (overlapped by Venn)
getname_day1_up = rownames(DEG_day1_up)
getname_day3_up = rownames(DEG_day3_up)
getname_day7_up = rownames(DEG_day7_up)

getname_day1_down = rownames(DEG_day1_down)
getname_day3_down = rownames(DEG_day3_down)
getname_day7_down = rownames(DEG_day7_down)
# 
# install.packages('tidyverse')
# install.packages('hrbrthemes')
# install.packages('tm')
# install.packages('proustr')
# install.packages('VennDiagram')
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
# 
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(getname_day1_up, getname_day3_up, getname_day7_up),
  category.names = c("Day 1" , "Day 3 " , "Day 7"),
  filename = 'Stroke_and_control_DEG_up.png',
  output=TRUE,
  
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)



venn.diagram(
  x = list(getname_day1_down, getname_day3_down, getname_day7_down),
  category.names = c("Day 1" , "Day 3 " , "Day 7"),
  filename = 'Stroke_and_control_DEG_down.png',
  output=TRUE,
  
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# install.packages('stringi')
library(stringi)


############################# find DEG 
DEG_up = rownames(find_DEG_up)
DEG_down = rownames(find_DEG_down)
DEG_total = c(DEG_up,DEG_down)

write(DEG_total,"DEG_total.txt",sep = '')  # save TXT data





###  DEG down together
# 919 genes, 277 genes and 25 genes
# gene_list_down <- c(getname_day1_down,getname_day3_down,getname_day7_down)  #  +  1221
getname_day_1_3_down <- c(getname_day1_down,getname_day3_down)   # 1196
day_1_3_down = getname_day_1_3_down[stri_duplicated( getname_day_1_3_down)]   # 157
getname_day_1_7_down <- c(getname_day1_down,getname_day7_down)   # 944
day_1_7_down = getname_day_1_7_down[stri_duplicated( getname_day_1_7_down)]
getname_day_3_7_down <- c(getname_day3_down,getname_day7_down)   # 176
day_3_7_down = getname_day_3_7_down[stri_duplicated( getname_day_3_7_down)]

getname_day_137_down <- c(day_1_3_down,day_1_7_down,day_3_7_down)
DEG_down <- getname_day_137_down[stri_duplicated( getname_day_137_down)]  # 6 overlap (three)
# write(DEG_down,"overlap_DEG_FC15_down.txt",sep = '')

###down

down_regulate <- c(day_1_3_down, getname_day7_down)
down_regulate_total <- down_regulate[stri_duplicated(down_regulate)]
write(down_regulate_total,"overlap_DEG_FC15_down.txt",sep = '')  # 3 genes



###  DEG up together
getname_day_1_3_up <- c(getname_day1_up,getname_day3_up)  # 1596
day_1_3_up = getname_day_1_3_up[stri_duplicated( getname_day_1_3_up)] 


getname_day_1_7_up <- c(getname_day1_up,getname_day7_up)  # 1544
day_1_7_up = getname_day_1_7_up[stri_duplicated( getname_day_1_7_up)]
getname_day_3_7_up <- c(getname_day3_up,getname_day7_up)  # 705
day_3_7_up = getname_day_3_7_up[stri_duplicated(getname_day_3_7_up)]

getname_day_137_up <- c(day_1_3_up,day_1_7_up,day_3_7_up)
DEG_up <- getname_day_137_up[stri_duplicated( getname_day_137_up)]  # 278 overlap (three)
# write(DEG_up,"overlap_DEG_FC15_up.txt",sep = '')

### up 
up_regulate <- c(day_1_3_up, getname_day7_up)
up_regulate_total <- up_regulate[stri_duplicated(up_regulate)]
write(up_regulate_total,"overlap_DEG_FC15_up.txt",sep = '')  # 139 genes





############################ Last step: overlap btw DEG and co-expression
# load('selected_genes_name_WGCNA.Rdata')
load('first_ten_record_greenyellow.Rdata.Rdata')
load('first_ten_record_lightyellow.Rdata')

library(stringr)
library(stringi)
two_module <- c(rownames(first_ten_record_greenyellow),rownames(first_ten_record_lightyellow))
# brown_and_green_repeat = brown_and_green[!stri_duplicated(brown_and_green)]
two_module_repeat = two_module[!stri_duplicated(two_module)]   # 324


DEG_total <- c(up_regulate_total,down_regulate_total)
# DEG_total <- c(DEG_down,DEG_up)

write(DEG_total,"DEG_total.txt",sep = '')  # total 142  (139+3) 



###  hubs from STRING
DEG_hub_name = c("Cxcl10","Cdkn1a","Psmb8","Atf3","Ccnd1","Usp18","Ifit3","Ifit2","Irf1","Psmb9",
                 "Ecm1","Anxa3","Anxa4","Il6","Pla2g4a","Gpx1","Lgals1","Fn1","Clcf1","Socs3","Ccl2",
                 "Timp1","S100a4","Vim","Nes","Anxa2","C3")



# total_gene = together_up[!stri_duplicated(together_up)]
library(VennDiagram)
# write(together_total,"together_total_21days.txt",sep = '')  # save TXT data

myCol2 <- c("antiquewhite","chartreuse")

venn.diagram(
  x = list(two_module_repeat, DEG_hub_name),
  category.names = c("Coexpress", "DEG"),
  filename = 'Two_method_together_all_21_days.png',
  output=TRUE,
  
  imagetype="png",
  height = 800 ,
  width = 1000 ,
  resolution = 400,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = myCol2,
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # 
  # Set names
  # cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  # rotation = 1,
)

repeated_genes <- intersect(two_module_repeat, DEG_hub_name)

length(repeated_genes)

write(repeated_genes,"repeated_genes.txt",sep = '')  # Cdkn1a









############### heatmap
library(pheatmap)
choose_gene = head(DEG_total,30)
choose_matrix = exprSet[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
pheatmap(choose_matrix, )
exprSet['Ptrh1',]


############### volcano plot
DEG=nrDEG
DEG$logFC
logFC_cutoff = 0.58
# logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'Up-regulated','Down-regulated'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nNumber of up-regulated gene is ',nrow(DEG[DEG$change =='Up-regulated',]) ,
                    '\nNumber of down-regulated gene is ',nrow(DEG[DEG$change =='Down-regulated',])
)

this_tile
# install.packages("ggplot2")
library("ggplot2")
head(DEG)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)





########################## online pathway analyses

### STRING
# BiocManager::install("STRINGdb")
library(STRINGdb)

data(diff_exp_example1)

gene_List = c(nrDEG$P.Value,nrDEG$logFC)
names(gene_List) = rownames(nrDEG)


string_db <- STRINGdb$new(version="11", species=10090,
                          score_threshold=700, input_directory= '') # 10090 Mus musculus
# deg_mapped <- string_db$map(together, "gene", removeUnmappedRows = TRUE )

example1_mapped = string_db$map(DEG_total, "gene", removeUnmappedRows = TRUE )
string_db$plot_ppi_enrichment( example1_mapped$STRING_id[1:1000] )

### mapping 
if(F){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GO.db")} # install GO database
# source("https://bioconductor.org/biocLite.R")
# BiocManager::install(version = "3.12")
# biocLite("STRINGdb")


# library(STRINGdb)
# string_db <- STRINGdb$new(version="11", species=10090,
#                           score_threshold=700, input_directory= '') # 10090 Mus musculus
# deg_mapped <- string_db$map( find_DEG, "gene", removeUnmappedRows = TRUE )



####################### enrichment KEGG
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
data(geneList, package="DOSE") # Disease Ontology Semantic and Enrichment analysis
# gene_getname <- names(geneList)[abs(geneList) > 2]
gene_getname.df <- bitr(DEG_total, fromType = "SYMBOL",
                        toType = c("ENSEMBL","ENTREZID"),
                        OrgDb = "mouse430a2.db")
head(gene_getname.df)

kk <- enrichKEGG(gene = gene_getname.df$ENTREZID,
                 organism = 'mmu',
                 pvalueCutoff = 0.05)
# kk$result
head(kk)[,-1]


### GSEA
# boxplot(nrDEG$logFC)
logFC_result = nrDEG$logFC

dim(kk$Description)
vignette('clusterProfiler')

