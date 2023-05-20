getwd()
setwd('H:\\Users\\yrt05\\Desktop\\Systems biology project\\R script\\test')  # \\ abs path

#download data from GEO database, get file ".csv"
#function source: https://github.com/jmzeng1314/humanid/blob/master/R/downGSE.R

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

if(T){
  downGSE <- function(studyID = "GSE137482", destdir = ".") {
    library(GEOquery)
    eSet <- getGEO(studyID, destdir = destdir, getGPL = F, GSEMatrix=T)
    
    # exprSet = exprs(eSet[[1]])
    pdata = pData(eSet[[1]])
    
    
    save(eSet,file = 'Seq_data.Rdata')
    # write.csv(exprSet, paste0(studyID, "_exprSet.csv"))  
    write.csv(pdata, paste0(studyID, "_metadata.csv"))  # classification of data
    return(eSet)
  }
  
  }

# getGEOSuppFiles('GSE137482')

# myColData <- phenoData(eSet) 

clindata <- eSet[["GSE137482_series_matrix.txt.gz"]]@phenoData@data

head(clindata)

# BiocManager::install("GEOquery")
library(GEOquery)
# eSet <- downGSE('GSE137482')

# orig_data = exprs(eSet[[1]])

### load data
load('Seq_data.Rdata')


### download database
# BiocManager::install("mouse430a2.db")
library(mouse430a2.db)


### get group names
library(Biobase)
# library(GEOquery)
pdata = pData(eSet[[1]])


### Normalization

setwd("H:\\Users\\yrt05\\Desktop\\Systems biology project\\GSE137482_RAW")
# install.packages("BiocManager")
# BiocManager::install("affy")
library(affy)
# mydata_GSE35338<-ReadAffy()
# data <- exprs(mydata_GSE35338[,c(1:9,19:21,25:27,22:24,28:30)][1])


setwd("H:\\Users\\yrt05\\Desktop\\Systems biology project\\GSE137482_RAW\\")
files <- list.files(pattern=".count.txt.gz")
## read data using loop
raw_counts <- read.delim(files[1], stringsAsFactors=FALSE, header = FALSE)
raw_counts <- as.matrix(raw_counts)

DF <- raw_counts
rownames(DF) <- raw_counts[,1]

library(stringr)


for (f in files) {
  dat <- read.delim(f, header=FALSE, stringsAsFactors=FALSE)
  
  col_nam <- gsub(".count.txt.gz", "", str_trim(f))
  colnames(dat)[2] <- c(col_nam)
  
  DF <- cbind(DF, dat[2])
}
DF <- DF[,c(-1,-2)]
DF <- DF[rowSums(DF[])>0,]  # delete all zero rows

### get group names
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biobase")

library(Biobase)

dim(pdata)
group_list <- pdata[,c(43)]
group_list=as.character(group_list)  

library(stringr)
group_list_replace <- gsub("\\s+", "_", str_trim(group_list))
group_list_replace <- gsub("-", "_", str_trim(group_list_replace))

### Normalization
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list_replace))
colnames(design)=levels(factor(group_list_replace))
rownames(design)=rownames(pdata[2])
# save(design,file='design.Rdata')

contrast.matrix<-makeContrasts("Control - MCAO",levels = design)

v <- voom( DF ,design,normalize="quantile")  ### normalization

data_sample_probe <- v[[1]]
data_sample_probe <- as.data.frame(data_sample_probe)

head(v[[1]])

fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tempOutput = topTable(fit3, coef=1, n=Inf)
DEG_voom = na.omit(tempOutput)
head(DEG_voom)

rownames(DEG_voom) <- gsub("\\..*","", str_trim(rownames(DEG_voom)))


### MAS5 normalization
# mydata_GSE35338_mas5<-mas5(mydata_GSE35338)
# # GC-RMA normalization
# biocLite("gcrma")
# mydata_GSE4536_gcrma<-gcrma(mydata_GSE4536)

### Farms normalization
# biocLite("farms")
# library(farms)
# mydata_GSE4536_farms<-qFarms(mydata_GSE4536)
# # dchip normalization
# mydata_GSE4536_dchip<-expresso(mydata_GSE4536,normalize.method="invariantset",bg.correct=FASE,pmcorrect.method="pmonly",summary.method="liwong")


# load('exprSet_mydata_GSE35338_rma.Rdata')


### Match probe ID with gene. (Preprocess to data)

library(dplyr)
DEG_voom %>% distinct()
DEG_voom <- DEG_voom[!duplicated(rownames(DEG_voom)), ] 

data_sample_probe <- data_sample_probe[!duplicated(rownames(data_sample_probe)), ] 

library(mouse430a2.db)
gene_probe_name <- mapIds(org.Mm.eg.db, keys = rownames(DEG_voom), keytype = "ENSEMBL", column = "SYMBOL", multiVals="first")
gene_probe_name <- gene_probe_name[!duplicated(gene_probe_name)] 

inds <- which(!is.na(gene_probe_name))
found_genes <- gene_probe_name[inds]
df2 <- DEG_voom[names(found_genes), ]
rownames(df2) <- found_genes

data_sample_probe <- data_sample_probe[names(found_genes), ]   # For use in "Co-expression"
rownames(data_sample_probe) <- found_genes
# save(data_sample_probe,file='data_sample_probe.Rdata')  # For use in "Co-expression"


### filter probes
ids = toTable(mouse430a2SYMBOL) 
table(rownames(df2) %in% ids$symbol)
exprSet <- df2
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$symbol,]
dim(exprSet)

# same sequence of Series matrix and probe
ids=ids[match(rownames(exprSet),ids$symbol),]
identical(ids$probe_id,rownames(exprSet)) ## Output "FALSE"
dat=exprSet

dim(ids)
dim(dat)


exprSet=dat
### colorful boxplot
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(exprSet, col = cols,main="expression value",las=2)
# boxplot(exprSet[,5])


### select DEG result ####
# nrDEG <- data_sample_probe
nrDEG <- dat
logFC_threshold = 1
find_DEG_up = nrDEG[which(nrDEG[1]> logFC_threshold & nrDEG[5]<0.01),]
find_DEG_down = nrDEG[which(nrDEG[1]< (-logFC_threshold) & nrDEG[5]<0.01),]


library(stringi)


############################# find DEG 
DEG_up = rownames(find_DEG_up)
DEG_down = rownames(find_DEG_down)
DEG_total = c(DEG_up,DEG_down)
length(DEG_total)
DEG_total_11_30 <- DEG_total


setwd("C:/Users/yrt05/Desktop/Systems biology project/GSE35338_RAW")
commen_gene <- read.delim("DEG_total.txt", header=FALSE)
DEG_two_data_common <- intersect(DEG_total_11_30, commen_gene$V1) ### Common DEGs from two data sets

# write(DEG_total,"DEG_new.txt",sep = '')  # save TXT data




############################ Last step: overlap btw DEG and co-expression
setwd("C:/Users/yrt05/Desktop/Systems biology project/GSE137482_RAW")
load('selected_genes_name_WGCNA_1130.Rdata')
# # together_up <- c(selected_genes_name,DEG_up)
# together_up <- c(selected_genes_name,DEG_up)
# together_down <- c(selected_genes_name,DEG_down)
# together_total <- c(selected_genes_name,DEG_up,DEG_down)

DEG_WGCNA_common <- intersect(selected_genes_name, DEG_total_11_30)
write(DEG_WGCNA_common,"DEG_WGCNA_common_11_30.txt",sep = '')

# 
library(VennDiagram)
# write(together_total,"Two_method_together_11_30.txt",sep = '')  # save TXT data

venn.diagram(
  x = list(selected_genes_name, DEG_total_11_30),
  category.names = c("Co-expres", "DEG"),
  filename = 'Two_method_together_11_30.png',
  # output=TRUE
)

intersect(DEG_two_data_common, DEG_WGCNA_common)





########################## online pathway analyses

### STRING
# BiocManager::install("STRINGdb")
library(STRINGdb)

# data(diff_exp_example1)
# head(DEG_voom)
rownames(exprSet)


nrDEG <- exprSet%>%filter(rownames(exprSet)%in%DEG_total)
# nrDEG <- exprSet%>%filter(rownames(exprSet)%in%DEG_WGCNA_common)

gene_List = c(nrDEG$adj.P.Val,nrDEG$logFC)
names(gene_List) = rownames(nrDEG)











### Pathway analysis plot
####################### enrichment KEGG


# BiocManager::install("clusterProfiler")
library(clusterProfiler)
data(geneList, package="DOSE") # Disease Ontology Semantic and Enrichment analysis
# gene_getname <- names(geneList)[abs(geneList) > 2]
gene_getname.df <- bitr(DEG_total_11_30, fromType = "SYMBOL",
                        toType = c("ENSEMBL","ENTREZID"),
                        OrgDb = "mouse430a2.db")
head(gene_getname.df)
# write(gene_getname.df$ENSEMBL, "DEG_ENSEMBL.txt", sep = '')


# install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method","auto")
# options(clusterProfiler.download.method = "wget")
R.utils::setOption("clusterProfiler.download.method","wget")
kk <- enrichKEGG(gene = gene_getname.df$ENTREZID,
                 organism = "hsa",
                 qvalueCutoff = 0.05)


font_size = 11
library(scales) 
head(kk@result)
kk@result[["Description"]][1:3]
kegg <- kk
p<- dotplot(kegg, showCategory=5, x = 'Count', font.size = font_size, orderBy="x")

p1 <- p + theme(
  plot.title = element_text(color="black", size=font_size, face="bold", hjust=0.5))



gene_id <- bitr(rownames(nrDEG), fromType = "SYMBOL",
                         toType = c("ENTREZID"),
                         OrgDb = "mouse430a2.db")
fold_vector <- nrDEG$logFC
names(fold_vector) <- gene_id$ENTREZID


### Pathway
p2 <- cnetplot(kk, showCategory=5, foldChange=fold_vector, node_label = "category", layout = "kk",
               colorEdge = TRUE, color_category ='steelblue', node_label_size = NULL,
               cex_category = 1, label_format = 5, cex_label_category = 1)
# min.value <- floor( min(p2$data$color, na.rm = TRUE) )
# max.value <- ceiling( max(p2$data$color, na.rm = TRUE) )
# p2 <- p2 + scale_color_gradientn(name = "Fold Change",
#                            colours = c("blue", "red"),
#                            values = rescale(c(min.value, 0, max.value)),
#                            limits= c(min.value, max.value), 
#                            breaks=c(min.value , 0, max.value) )
p2 <- p2 + theme(legend.position="right",
                 legend.text = element_text(size = font_size),
                 # legend.background = element_rect(fill="lightblue", size = 1),
                 legend.title = element_text(size = font_size, face = "bold"),
                 legend.justification = c("right", "top"),
                 # legend.margin = margin(6, 6, 6, 6)
                 )

plot_grid(p2, p1, labels = c("a","b"), label_size = 15,  nrow = 2, align = c("h"))




### Significant DEGs
significant_gene_id <- kk@result[["geneID"]][1]

significant_gene_id <- gsub("/"," ", str_trim(significant_gene_id))

significant_gene_id <- str_split(significant_gene_id, " ")[[1]]

significant_gene <- bitr(significant_gene_id, fromType = "ENTREZID",
                        toType = c("SYMBOL"),
                        OrgDb = "mouse430a2.db")
significant_gene_KEGG <- significant_gene$SYMBOL


intersect(significant_gene_KEGG, selected_genes_name)

intersect(intersect(significant_gene_KEGG, selected_genes_name), DEG_WGCNA_common)


up_DEG <- intersect(significant_gene_KEGG, rownames(find_DEG_up))
down_DEG <- intersect(significant_gene_KEGG, rownames(find_DEG_down))
up_DEG <- bitr(up_DEG, fromType = "SYMBOL",
                         toType = c("ENTREZID"),
                         OrgDb = "mouse430a2.db")
down_DEG <- bitr(down_DEG, fromType = "SYMBOL",
               toType = c("ENTREZID"),
               OrgDb = "mouse430a2.db")
rownames(nrDEG)
P_up <- dat[rownames(dat) %in% up_DEG$SYMBOL,4]
P_down <- dat[rownames(dat) %in% down_DEG$SYMBOL,4]



### GO enrichment
# 
# Genes<-c(significant_gene_id)
# # g <- goana(Genes)
# # BiocManager::install("org.mmu.eg.db")
# # library(org.mmu.eg.db)
# g <- goana(list(Up=up_DEG$ENTREZID, Down=down_DEG$ENTREZID, P.Up = P_up, P.Down = P_down), 
#            species="Mm", covariate = TRUE, plot = TRUE)
# GO_result <- topGO(g)
# GO_result$Term
# 
# GO_result$

### GO enrichment 2
library(clusterProfiler)
# de <- c(up_DEG$ENTREZID, down_DEG$ENTREZID)
de <- bitr(rownames(nrDEG), fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "mouse430a2.db")
de <- de$ENTREZID
go <- enrichGO(de, OrgDb = "org.Mm.eg.db", ont="all")
# ego <- enrichGO(de, OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
# go <- enrichDAVID(de, idType="ENTREZ_GENE_ID",  listType="Gene", annotation="KEGG_PATHWAY")
library(ggplot2)


library(GOplot)
data(EC)
test <- EC$david
setwd("C:/Users/yrt05/Desktop/Systems biology project/GSE137482_RAW")
readfile <- read.csv("gProfiler_mmusculus.csv")        ### web data through https://biit.cs.ut.ee/gprofiler/gost
readfile <- subset(readfile, select = c("source","term_name", "intersections", "adjusted_p_value", "term_id") )
colnames(readfile)[1:4] <-colnames(test)[c(1,3,4,5)]
readfile <- readfile[, c(1,5,2,3,4)]
colnames(readfile)[2] <-colnames(test)[2]
test3 <- as.factor(readfile$Genes)
length(test3)
test3 <- str_split(test3, ",")
# unlist(test3)
library("org.Mm.eg.db")
nrDEG_data <- cbind(rownames(nrDEG), nrDEG)
colnames(nrDEG_data) <- c("ID", colnames(nrDEG))
readfile$Category <- gsub("GO:","", str_trim(readfile$Category))
rownames(nrDEG_data) <- c(1:nrow(nrDEG_data))

# length(test3)
for (i in 1:length(test3)){
  map <- mapIds(org.Mm.eg.db, keys = unlist(str_split(unlist(test3[i]),",")) , keytype = "ENSEMBL", column = "SYMBOL", multiVals="first")
  # map <- mapIds(org.Mm.eg.db, keys = map , keytype = "SYMBOL", column = "PROBEID", multiVals="list")
  result <- data.frame(paste(toupper(map), collapse=', '))
  readfile$Genes[i] <- result
  # print(as.factor(paste(toupper(map), collapse=', ')))
}

readfile <- filter(readfile, Category %in% c("CC", "BP", "MF"))

circ <- circle_dat(readfile, nrDEG_data)


# GOBubble(circ, labels = 10)
# Reduce redundant terms with a gene overlap >= 0.75...
# reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# # ...and plot it
# GOBubble(reduced_circ, labels = 2.8)

# GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 10)  

# GOBar(subset(circ, category == 'BP'))







### GO plot



go2 <- dotplot(go, showCategory=5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

### Separate analyses
go_result <- summary(go)
go_result_All <- go_result$ID[order(go_result$Count, decreasing=T)]
# head(go_result_All)
go_result_All[1:5]
IDs <- c(go_result_All[1:5])

go1 <- GOCircle(circ, nsub = IDs,  rad1 = 2, rad2 = 3,  zsc.col = c('yellow', 'black', 'cyan'))

# go_result_BP <- filter(go_result, ONTOLOGY %in% c("BP"))
# go_result_BP <- go_result_BP$ID[order(go_result_BP$Count, decreasing=T)]
# go_result_MF <- filter(go_result, ONTOLOGY %in% c("MF"))
# go_result_MF <- go_result_MF$ID[order(go_result_MF$Count, decreasing=T)]
# go_result_CC <- filter(go_result, ONTOLOGY %in% c("CC"))
# go_result_CC <- go_result_CC$ID[order(go_result_CC$Count, decreasing=T)]


go2 <- go2 + theme(legend.position="right",
                 legend.text = element_text(size = font_size),
                 # legend.background = element_rect(fill="lightblue", size = 1),
                 legend.title = element_text(size = font_size, face = "bold"),
                 legend.justification = c("right", "top"),
                 # legend.margin = margin(6, 6, 6, 6)
)
# plot(go2)
library(cowplot)
# ,   rel_heights = c(0.3, 1.2)
# par(mfrow=c(1,2))
# dotplot(go, showCategory=5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
# GOCircle(circ, nsub = IDs, label.size = 6, rad1 = 3, rad2 = 4,  zsc.col = c('yellow', 'black', 'cyan'))

plot_grid(go2, go1, labels = c("a","b"), label_size = 22, ncol = 2, align = c("h"), greedy = TRUE, 
          rel_widths = c(1, 1.5),
          axis = "b", scale = c(1,1), cex= 2)






##### PPI analysis
first_ten_record_yellow # from WGCNA, single module, high connectivity

dim(nrDEG)

genes = rownames(nrDEG)
# interactions = get_string_ppi_matrix(genes, version = "11")
get.ppiNCBI <- function(g.n) {
  require(XML)
  ppi <- data.frame()
  for(i in 1:length(g.n)){
    o <- htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/?term=", g.n[i], sep=''))
    # check if interaction table exists
    exist <- length(getNodeSet(o, "//table//th[@id='inter-prod']"))>0
    if(exist){
      p <- getNodeSet(o, "//table")
      ## need to know which table is the good one
      for(j in 1:length(p)){
        int <- readHTMLTable(p[[j]])
        if(colnames(int)[2]=="Interactant"){break}
      }
      ppi <- rbind(ppi, data.frame(egID=g.n[i], intSymbol=int$`Other Gene`))
    }
    # play nice! and avoid being kicked out from NCBI servers
    Sys.sleep(1)
  }
  if(dim(ppi)[1]>0){
    ppi <- unique(ppi)
    print(paste(dim(ppi)[1], "interactions found"))
    return(ppi)
  } else{
    print("No interaction found")
  }
}

ppi <- get.ppiNCBI(genes)
# install.packages("remotes")
remotes::install_github("moosa-r/rbioapi")

library(rbioapi)
proteins_mapped <- rba_string_map_ids(ids = genes,
                                      species = 10090)
graph_1 <- rba_string_network_image(ids = proteins_mapped$stringId,
                                    image_format = "image",
                                    species = 10090,
                                    save_image = TRUE,
                                    required_score = 700,
                                    hide_disconnected_nodes = TRUE,
                                    network_flavor = "confidence")
# plot(graph_1)



### Read hub from DEG
# C:/Users/yrt05/Desktop/Systems biology project/GSE137482_RAW
DEG_sting <- read.csv("DEG_string.csv")
DEG_sting <- DEG_sting[, c("shared.name", "Degree", "Stress", "BetweennessCentrality")]
rownames(DEG_sting) <- DEG_sting$shared.name

percent <- 0.1
DEG_Stress <- rownames(DEG_sting[order(DEG_sting$Stress, decreasing = TRUE), ])[1:round(percent*(dim(DEG_sting)[1]))]
str_sort(DEG_Stress)

DEG_Betweenness <- rownames(DEG_sting[order(DEG_sting$BetweennessCentrality, decreasing = TRUE), ])[1:round(percent*(dim(DEG_sting)[1]))]
str_sort(DEG_Betweenness)

DEG_Degree <- rownames(DEG_sting[order(DEG_sting$Degree, decreasing = TRUE), ])[1:round(percent*(dim(DEG_sting)[1]))]

library("ggvenn")

B <-list('Stress'=DEG_Stress, 'Betweenness'=DEG_Betweenness, 'Degree' = DEG_Degree)
ggvenn(B,show_percentage=FALSE, stroke_size = 0.5, set_name_size = 8,
       text_color = "black",
       text_size = 8
       )


DEG_hub <- intersect(intersect(DEG_Stress, DEG_Betweenness), DEG_Degree)

hub_two_method <- intersect(DEG_hub, rownames(first_ten_record_yellow))


### Venn plot

# install.packages("ggvenn")
library("ggvenn")
A <-list('Co-expression'=rownames(first_ten_record_yellow),'DEG'=DEG_hub)
ggvenn(A,show_percentage=FALSE, stroke_size = 0.5, set_name_size = 8,
       text_color = "black",
       text_size = 8,
       fill_color=c("red","orange"))








### Hub GO
de_hub_dat <- bitr(hub_two_method, fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = "mouse430a2.db")
de_hub <- de_hub_dat$ENTREZID
hub_enrich <- enrichGO(de_hub, OrgDb = "org.Mm.eg.db", ont="all")

dotplot(hub_enrich, showCategory=5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

GO_hub_intersect <- intersect(rownames(hub_enrich@result), go_result$ID)
# counts(intersect(rownames(hub_enrich@result), go_result$ID))
frequency <- data.frame(t(stringdist::qgrams(freq = GO_hub_intersect, q = 2)))

freq_x <- sort(table(unlist(strsplit(GO_hub_intersect, " "))),      # Create frequency table
               decreasing = TRUE)
freq_x 



### Hub pathway
kk_hub <- enrichKEGG(gene = de_hub,
                 organism = "mmu",
                 qvalueCutoff = 0.05)

head(kk_hub@result)

### One particular pathway for a gene
# kk_hub@result[["Description"]][71]
# id_back_syml <- kk_hub@result[["geneID"]][71]
# de_hub_dat[de_hub_dat$ENTREZID == id_back_syml,]
# ### Two
# 
# kk_hub@result[["Description"]]

k=0
datalist = vector("list", length = length(kk_hub@result[["Description"]]))
for (j in kk_hub@result[["Description"]]){
  # pathways_select <- "Phagosome"
  ind <- which(kk_hub@result[["Description"]] == j)
  id_back_syml <- kk_hub@result[["geneID"]][ind]
  # id_back_syml <- gsub("/", " ", id_back_syml)
  
  id_back_syml <- str_split(id_back_syml, "/")
  
  # print(pathways_select)
  gene_name_list <- list()
  for (i in id_back_syml[[1]]){
    gene_select <- de_hub_dat[de_hub_dat$ENTREZID == i,]
    # print(gene_select$SYMBOL)
    gene_name_list <- append(gene_name_list, gene_select$SYMBOL)
  }
  k = k + 1
  datalist[[k]] <- append(j , gene_name_list[[1]])
  
}


pathways_select <- "Phagosome"
ind <- which(kk_hub@result[["Description"]] == pathways_select)
id_back_syml <- kk_hub@result[["geneID"]][ind]
# id_back_syml <- gsub("/", " ", id_back_syml)

id_back_syml <- str_split(id_back_syml, "/")

for (i in id_back_syml[[1]]){
  gene_select <- de_hub_dat[de_hub_dat$ENTREZID == i,]
  print(pathways_select)
  print(gene_select$SYMBOL)
}





gene_id_name <- kk_hub@result[["geneID"]]
gene_id_name <- gsub("/", " ", gene_id_name)
library(stringi)
hub_dat[,c(colnames(hub_dat),"Pathway")]=""
hub_dat[,c(1,2)] <- de_hub_dat
index = 0
for (i in de_hub_dat$ENTREZID){
  index <- index + 1
  give_gene_id <- i
  find_pathway <- stri_detect_fixed(gene_id_name, give_gene_id)
  hub_dat[index,3] <- paste(  list( kk_hub@result$Description[find_pathway])  ,collapse = ' ')
}
data.frame( list(kk_hub@result$Description[find_pathway]))



pathways_select <- "Phagosome"
stri_detect_fixed(hub_dat[,3], pathways_select)











