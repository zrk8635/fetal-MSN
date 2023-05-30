rm(list=ls())

# load packages

library(readxl)
library(Hmisc) # rcorr() return r and p , cor() only return r, cor.test() cannot compute matrix correlation
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# load data ---------------------------------------------------------------

setwd("fetus_individual_GM_network/roi/CorrelationAnalysis")
Labels <- read.table("CorticalRegionLabels_reorder_AAL.txt",sep = "\t")
colnames(Labels) <- c("rank","key","name","abbr","hemisphere","lobe","subregion")
Group <- read.table("GA.txt")
Subject <- read.table("subject.txt")
input <- read_xlsx("result.xlsx")
input <- as.matrix(input[,8:ncol(input)])

# set parameter -----------------------------------------------------------

n <- nrow(Subject); # number of subject
m <- 78; # number of cortical roi
l <- 8; # number of measure
features <- c("NumberofVertices","Thickness","SulcalDepth","MeanCurvature","GaussCurvature","Curvedness","Volume","SurfaceArea")
edge <- m*(m-1)/2
data_raw <- array(0,dim=c(l,m,n))
data <- array(NA,dim=c(l,m,n))

# extract measure matrix for subject ---------------------------------

for (i in 1:n)
{
  for (k in 1:l)
  {
    for (j in 1:m)
    {
      data_raw[k,j,i] <- input[i,((k-1)*m+j)]
    }
    # z-score normalization
    data[k,,i] <- scale(data_raw[k,,i],center=T,scale=T)
  }
  # reorder
  data[,,i] <- data[,c(Labels[,1]),i]
  write.table(data[,,i],paste0("fetus_individual_GM_network/roi/CorrelationAnalysis/data/subject/",Subject[i,],"/",Subject[i,],".txt"),
              sep = "\t", quote=FALSE, row.names = features, col.names = Labels[,4])
}

# correlation of subjects-------------------------------------------------------------

r <- array(1,dim=c(m,m,n));
p <- array(NA,dim=c(m,m,n));
for (i in 1:n)
{
  result <- rcorr(data[-8,,i],type='pearson') #delete surface area
  r[,,i]=result[["r"]];
  p[,,i]=result[["P"]];
  for (j in 1:m)
  {
    r[j,j,i]=0
  }
  rm("result")
  write.table(r[,,i],paste0("fetus_individual_GM_network/roi/CorrelationAnalysis/data/subject/",Subject[i,],"/",Subject[i,],"_correlation.txt"),
              sep = "\t", quote=FALSE,row.names = FALSE,col.names = FALSE)
}

# heatmap of correlation matrix of subject --------------------------------

for (i in 1:n)
{
  mat=as.matrix(r[,,i])
  rownames(mat)=Labels[,4]
  colnames(mat)=Labels[,4]
  pdf(paste0("fetus_individual_GM_network/roi/CorrelationAnalysis/data/subject/",Subject[i,],"/",Subject[i,],"_correaltion.pdf"),width=14,height=14)
  figure=ComplexHeatmap::Heatmap(mat,name='correlation',
                                 col = rev(colorRampPalette(brewer.pal(9, "Spectral"))(99)), 
                                 show_heatmap_legend = FALSE,
                                 # right_annotation = annot_row,
                                 # top_annotation = annot_col,
                                 cluster_rows=FALSE,
                                 show_row_dend = FALSE,
                                 show_row_names = FALSE,
                                 row_names_side = "left", 
                                 row_names_gp = gpar(fontsize = 6),
                                 cluster_columns=FALSE,
                                 show_column_dend = FALSE,
                                 show_column_names = FALSE,
                                 column_names_side = "bottom",
                                 column_names_gp = gpar(fontsize = 6)
  )
  print(figure)
  dev.off()
}

# threshold for correlation matrix of subject -------------------------------------------------

threshold <- array(1,dim=c(1,n))
sparsity <- array(1,dim=c(1,n))
for (i in 1:n)
{
  r_threshold=r[,,i]
  threshold[i]=min(apply(r_threshold,1,max))
  r_threshold=r_threshold[lower.tri(r_threshold)]
  r_threshold=sort(r_threshold,decreasing = TRUE)
  sparsity[i]=((which(r_threshold<threshold[i]))[1]-1)/edge
  sparsity[i]=ceiling(sparsity[i]*100)
}
write.table(sparsity,paste0("fetus_individual_GM_network/roi/CorrelationAnalysis/data/subject/","sparsity.txt"),
            sep = "\t", quote=FALSE, row.names = FALSE,col.names = FALSE)

# binary correlation matrix of subject ---------------------------------------------------

spars <- max(sparsity)
n_edge <- ceiling(edge*(spars/100))
r_bin <- array(0,dim=c(m,m,n));
for (i in 1:n)
{
  mat <- r[,,i]
  seq <- mat[lower.tri(mat)]
  seq <- sort(seq,decreasing = TRUE)
  thres <- seq[n_edge]
  mat[which(mat>=thres)] <- 1
  mat[which(mat<thres)] <- 0
  r_bin[,,i] <- mat
  rm("mat","seq")
  write.table(r_bin[,,i],paste0("fetus_individual_GM_network/roi/CorrelationAnalysis/data/subject/",Subject[i,],"/",Subject[i,],"_binary.txt"),
              sep = "\t", quote=FALSE, row.names = FALSE,col.names = FALSE)
}

# heatmap of binary correlation matrix of subject -------------------------

for (i in 1:n)
{
  mat=as.matrix(r_bin[,,i])
  rownames(mat)=Labels[,4]
  colnames(mat)=Labels[,4]
  pdf(paste0("fetus_individual_GM_network/roi/CorrelationAnalysis/data/subject/",Subject[i,],"/",Subject[i,],"_binary.pdf"),width=14,height=14)
  figure=ComplexHeatmap::Heatmap(mat,name='correlation',
                                 col = c("black","white"),
                                 show_heatmap_legend = FALSE,
                                 # right_annotation = annot_row,
                                 # top_annotation = annot_col,
                                 cluster_rows=FALSE,
                                 show_row_dend = FALSE,
                                 show_row_names = FALSE,
                                 row_names_side = "left",
                                 row_names_gp = gpar(fontsize = 6),
                                 cluster_columns=FALSE,
                                 show_column_dend = FALSE,
                                 show_column_names = FALSE,
                                 column_names_side = "bottom",
                                 column_names_gp = gpar(fontsize = 6)
  )
  print(figure)
  dev.off()
}
