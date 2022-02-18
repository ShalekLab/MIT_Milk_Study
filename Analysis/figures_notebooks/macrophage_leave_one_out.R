
library(Seurat)
library(tidyverse)
library(ggpubr)
macs<-readRDS("../../Data/processed_data/seuratobj_immune_macs_analysis_subset.rds")

loo_df = macs@meta.data[,c("RNA_snn_res.0.2","time_post_partum_weeks")]

for (d in unique(macs@meta.data$donor)){
one_sample_out <- subset(macs, subset = donor %in% c(d ), invert = TRUE)
one_sample_out <- RunPCA(one_sample_out, features = VariableFeatures(object = one_sample_out))
one_sample_out <- FindNeighbors( one_sample_out, dims = 1:20)
one_sample_out   <- FindClusters(one_sample_out, resolution = 0.2)

    loo_df[rownames(one_sample_out@meta.data), d]=one_sample_out@meta.data$seurat_clusters
    
    
    
    }


write.csv(loo_df, "../../Data/processed_data/macrophage_leave_one_donor_out_results.csv",quote=FALSE)


# dig deeper into BM05 left out
one_sample_out <- subset(macs, subset = donor %in% c("BM05" ), invert = TRUE)
one_sample_out <- RunPCA(one_sample_out, features = VariableFeatures(object = one_sample_out))
one_sample_out <- FindNeighbors( one_sample_out, dims = 1:20)
one_sample_out   <- FindClusters(one_sample_out, resolution = 0.2)


my_levels <- c(0,2,3,1,4)
one_sample_out@active.ident <- factor(x = yone_sample_out@active.ident, levels = my_levels)


# find markers with BM05 left out
Idents(object=one_sample_out)<-"seurat_clusters"
levels(one_sample_out) <- my_levels
seurat.markers.no.bm05 <- FindAllMarkers(one_sample_out, test.use="wilcox", only.pos=TRUE, min.pct=0.2,logfc.threshold=0.5)


top10.sub <- seurat.markers.no.bm05 %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DEmap <- DoHeatmap(one_sample_out, features = top10.sub$gene,raster=FALSE)  +theme(text = element_text(size = 8)) + scale_fill_gradientn(colors = colorRampPalette(c("#2C7BB6", "#FFFFBF", "#D7191C"))(256))
print(DEmap)
ggsave("../../Results/plots/figure_S6/heatmap_leave_macrophage_out.pdf",useDingbats=F)


# recalculate module scores with BM05 left out


meta.data.modulescores.M1 <- one_sample_out@meta.data[,c("seurat_clusters","M1_UP1")]
meta.data.modulescores.M1$module <- "M1"
colnames(meta.data.modulescores.M1) <- c("clusters","score","module")
meta.data.modulescores.M2 <- one_sample_out@meta.data[,c("seurat_clusters","M2_UP1")]
meta.data.modulescores.M2$module <- "M2"
colnames(meta.data.modulescores.M2) <- c("clusters","score","module")

mac.module.scores <- rbind(meta.data.modulescores.M1,meta.data.modulescores.M2)

p <- ggviolin(mac.module.scores, x = "module", y = "score", 
         fill = "module", palette = "jco",
         add = "boxplot", add.params = list(fill = "white")) +
         stat_compare_means(label.x = 0.5, method = "kruskal.test")

facet(p, facet.by = "clusters", nrow = 1, ncol = 5)
ggsave("../../Results/plots/figure_S6/M1M2_modules_leaveoutbm05.pdf", height = 4, width = 8, units = "in", useDingbats = F)
