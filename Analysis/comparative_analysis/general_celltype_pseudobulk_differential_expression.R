library(DESeq2)

general_celltype_gene_exp <- read.csv("../../Data/processed_data/cell_subset_expression/general_celltype_pseudobulk_counts.csv",row.names=1)

metadata <- read.csv("../../Data/processed_data/cell_subset_expression/general_celltype_pseudobulk_metadata.csv", row.names=1)

general_celltype_gene_exp <- na.omit(general_celltype_gene_exp)
general_celltype_gene_exp <- t(general_celltype_gene_exp)

metadata <- metadata[colnames(general_celltype_gene_exp),]

# remove small samples

general_celltype_gene_exp <- general_celltype_gene_exp[,metadata$Cell_Number > 10]
metadata = metadata[colnames(general_celltype_gene_exp),]



# marker genes for each celltype comparing expression within celltype vs rest
for( i in unique(metadata$simplified.celltypes)){
    print(i)
    metadata$is.thistype <- metadata$simplified.celltypes==i
    dds_marker <- DESeqDataSetFromMatrix(countData= general_celltype_gene_exp, colData = metadata, 
                                     design = ~factor(donor)+ is.thistype)

    dds_thistype_wald <- DESeq(dds_marker,test="Wald")
    res <- results(dds_thistype_wald)
    print("here")
    print(dim(res))
    write.csv(res, paste("../../Results/tables/General_Celltype_Pseudobulk_Marker_Genes/cluster_",i,"_marker_genes_pseudobulk.csv"),quote=FALSE)
}