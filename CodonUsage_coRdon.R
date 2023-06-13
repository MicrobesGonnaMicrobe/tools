#Codon usage bias analysis: Predicting gene expression

#coRdon: codon usage analysis in R
#R package for analysis of codone usage in unannotated or KEGG/COG annotated DNA sequences. 
#Calculates various measures of CU bias and CU-based predictors of gene expression, 
#and performs gene set enrichment analysis for annotated sequences. 
#Implements several methods for visualization of CU and enrichment analysis results.

#Instructions, vignette
#https://www.bioconductor.org/packages/release/bioc/vignettes/coRdon/inst/doc/coRdon.html

#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

#install coRdon
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("coRdon")

#installed other needed packages
BiocManager::install("KEGGREST")
BiocManager::install("ComplexHeatmap")

library(coRdon)
library(ggplot2)
library(Biobase)
library(KEGGREST)
library(grid)
library(ComplexHeatmap)

args = commandArgs(trailingOnly=TRUE)
# first flag is pathname where the genome files are stored, second flag is the project name
args[1] <- "C:/Users/phr001/Methods/Zetaproteobacteria/Analyses/Codon_bias/annotations_per_bin"
args[2] <- "Zetaproteobacteria"

# create a folder called after the project name to put the output files in
dir.create(paste0(args[1],args[2]))
# make a list of all the genome files in the project folder
files <- list.files(args[1],pattern = "\\.fa$")

for (fasta in files)
{
  pathtofile <- paste(args[1], fasta, sep = "/")
  dna <- readSet(file = pathtofile)
  codons <- codonTable(dna)
  milc <- MILC(codons, ribosomal = TRUE, filtering = "hard")
  
  if(!"ribosomal" %in% colnames(milc))
  {
    cat("no ribosomal genes found in",fasta);
    next
  }
  filename <- paste(pathtofile,"_enrichmentmatrix",sep = "")
  genomename <- tools::file_path_sans_ext(fasta)
  melp <- MELP(codons, ribosomal = TRUE, filtering = "hard", len.threshold = 100)
  genes <- getKO(codons)[getlen(codons) > 100]
  id <- getID(codons)
  l <- getlen(codons)
  # Creating a result table for each genome with unfiltered milc and "melp" values, as well as the original annotation
  milcunfiltered <- MILC(codons,self =FALSE, ribosomal = TRUE, filtering = "none")
  melpunfiltered <- MELP(codons, ribosomal = TRUE, filtering = "none")
  genesunfiltered <- getKO(codons)
  resulttable <- data.frame(milcunfiltered, melpunfiltered, genesunfiltered,id,l,dna)
  colnames(resulttable) <- c("MILC","MELP","KO","Annotation","Length","Sequence")
  #write.csv(resulttable,file= paste0(args[1],args[2],"/",genomename,"_resulttable.csv"),row.names=FALSE)
  write.csv2(resulttable,file= paste0(args[1],args[2],"/",genomename,"_resulttable.csv"),row.names=FALSE)
  
  if(!"ribosomal" %in% colnames(milc))
  {
    cat("no ribosomal genes found in",fasta);
    next
  }
  
  
  ct <- crossTab(genes, as.numeric(melp))
  #ct <- reduceCrossTab(ct, "module")
  enr<- enrichment(ct)
  assign(genomename,enr)
}
rm(enr)
dfs <- Filter(function(x) is(x, "AnnotatedDataFrame"), mget(ls()))
genomenames <-tools::file_path_sans_ext(files)
mat <- enrichMatrix(dfs, variable = "enrich")

paths <- names(KEGGREST::keggList("ko"))
paths <- regmatches(paths, regexpr("[[:alpha:]]{1}\\d{5}", paths))
pnames <- unname(KEGGREST::keggList("ko"))
ids <- match(rownames(mat), paths)
descriptions <- paste(rownames(mat),"|",pnames[ids])
rownames(mat) <- descriptions

#If the columns of the matrix should be ordered differently:
#matordered <- mat[,c(1, 2, 7, 8, 9, 3, 4, 10, 11, 12, 5,6,13, 14)]

#mat <- mat[apply(mat, 1, function(x) all(x!=0)), ]  <- this deletes all the KOs that are missing in some genomes

ht = ComplexHeatmap::Heatmap(
  mat, 
  name = "relative \nenrichment",
  col = circlize::colorRamp2( c(-100,0, 100), 
                              c("blue", "white", "red")),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  show_column_dend = FALSE, 
  show_row_dend = FALSE,
  width = ncol(mat)*unit(3, "mm"), 
  height = nrow(mat)*unit(3, "mm"),
  row_names_max_width = max_text_width(rownames(mat)),
  cluster_columns = FALSE,
  cluster_rows = TRUE)

ht = draw(ht)


w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)


pdf(file = paste0(args[1],args[2],"/",args[2],"_enrichment_heatmap_KO_1.pdf"),   
    width = w, # The width of the plot in inches
    height = h+1) # The height of the plot in inches  
ht
dev.off()

write.csv(mat,file= paste0(args[1],args[2],"/",args[2],"_enrichmentmatrix_KO.csv"))

write.table(mat, file = "Zetaproteobacteria_enrichmentmatrix_KO_simple.txt")