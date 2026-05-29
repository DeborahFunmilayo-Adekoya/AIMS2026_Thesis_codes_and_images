########################################
#importing the libraries
library(Seurat)
library(Matrix)
library(SeuratObject)
library(SeuratDisk)

#######################################
#setting the working directory as current work folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#confirming working directory
getwd()

##########################################
#loading in data
load("data_intestine.Rdata")
ls() #list all objects

#########################################
#proper extraction method to keep all information intact
class(seurat_obj_Intestine) #confirming the object type as Seurat object

#viewing it
seurat_obj_Intestine

#######################################
counts <- GetAssayData(seurat_obj_Intestine,assay="RNA" ,layer = "counts")
meta <- seurat_obj_Intestine@meta.data

#write the count to a matrix file
#method from matrix library
writeMM(counts,file="intestine_count.mtx")

#write the metadata to csv
write.csv(meta, "intestine_metadata.csv", row.names = TRUE)

#write the gene name to a table
write.table(rownames(seurat_obj_Intestine), "intestine_genes.tsv",col.names = FALSE,row.names = FALSE,quote=FALSE)

###############################################
###############################################
#The lung data
#loading in the lung data
load("data_sub_data_lung.Rdata")
ls() #list all objects

class(Seurat_Lung_sub)

Seurat_Lung_sub

#reading out the lung metadata
load("meta_data_sub_data_lung.Rdata")
ls()

#checking the datatype of the metadata lung file
class(meta_data_sub_lung)

#####################################saving to form for python
#writing to matrix. data in cell x gene format
writeMM(as(Seurat_Lung_sub,"sparseMatrix"), file="lung_matrix.mtx")

#export cell names from rows
write.table(rownames(Seurat_Lung_sub), file="lung_cell_names.tsv",
            sep = "\t", quote = FALSE, col.names= FALSE, row.names = FALSE)

#export gene names from columns
write.table(colnames(Seurat_Lung_sub), file="lung_gene_names.tsv",
            sep = "\t", quote = FALSE, col.names= FALSE, row.names = FALSE)


#export metadata
write.csv(meta_data_sub_lung, file = "lung_metadata.csv", row.names = TRUE)

dim(Seurat_Lung_sub)
