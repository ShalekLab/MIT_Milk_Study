library(DESeq2)

general_celltype_gene_exp <- read.csv("../../Data/processed_data/cell_subset_expression/epithelial_celltype_pseudobulk_counts.csv",row.names=1)
metadata <- read.csv("../../Data/processed_data/cell_subset_expression/epithelial_celltype_pseudobulk_metadata.csv", row.names=1)

celltype_gene_exp <- na.omit(celltype_gene_exp)
celltype_gene_exp <- t(celltype_gene_exp)

metadata <- metadata[colnames(celltype_gene_exp),]

# remove small samples

celltype_gene_exp <- celltype_gene_exp[,metadata$Cell_Number > 10]
metadata = metadata[colnames(celltype_gene_exp),]


for( i in unique(metadata$Epithelial.Cell.Subclusters)){
    metadata$is.thistype <- metadata$Epithelial.Cell.Subclusters==i
    dds_marker <- DESeqDataSetFromMatrix(countData = celltype_gene_exp, colData = metadata, 
                                     design = ~factor(donor)+ is.thistype)

    dds_thistype_wald <- DESeq(dds_marker,test="Wald")
    res <- results(dds_thistype_wald)
    write.csv(res, paste("../../Results/tables/Epithelial_Celltype_Pseudobulk_Marker_Genes/cluster_",i,"_marker_genes_pseudobulk.csv"),quote=FALSE)
}


# repeat for time varying genes
for( i in unique(metadata$Epithelial.Cell.Subclusters)){
    cepi <- celltype_gene_exp[,metadata$Epithelial.Cell.Subclusters == i]
    
    cepi <- na.omit(cepi)


    # remove very late samples and low cell number samples
    cepi_meta = metadata[colnames(cepi),]
    cepi <- cepi[,cepi_meta$time_post_partum_days < 400]
    cepi_meta = cepi_meta[colnames(cepi),]
    cepi <- cepi[,cepi_meta$Cell_Number > 10]
    cepi_meta = cepi_meta[colnames(cepi),]

    dim(cepi)

    dds <- DESeqDataSetFromMatrix(countData = cepi, colData = cepi_meta, design = ~ 0 +donor+ time_post_partum_days)
    dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ 0 +donor)
    write.csv(res,paste("../../Results/tables/Epithelial_Celltype_Pseudobulk_Time_Varying_Genes/cluster_",i,"_time_varying_genes_pseudobulk.csv"),quote=FALSE,row.names =TRUE)

    
    
 }
