## Trajectory analysis

library(Monocle3)

sce <- readRDS("sce.rds")

data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')

cell_metadata <- sce@meta.data

gene_annotation <- data.frame(gene_short_name = rownames(data))

rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "patient_loc")
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, preprocess_method = "Aligned")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds@int_colData$reducedDims$UMAP <- cds.embed

cds <- cluster_cells(cds)

cds <- learn_graph(cds,close_loop = T)

saveRDS(cds,file = "sce_monocle3.rds")
