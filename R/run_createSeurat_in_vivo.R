## Run create Seurat files for in vivo

## Fig 1; AA and healthy
aa_dirs        <- list.dirs("data/scRNAseq/aa/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
healthy_dirs   <- list.dirs("data/scRNAseq/healthy/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
folders        <- c(aa_dirs, healthy_dirs)
scrnaseq_files <- lapply(folders, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
aa_healthy_seurat  <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
aa_healthy_seurat  <- PercentageFeatureSet(aa_healthy_seurat, pattern = "^MT-", col.name = "percent.mt")
aa_healthy_seurat  <- PercentageFeatureSet(aa_healthy_seurat, pattern = "^RP", col.name = "percent.ribo")
aa_healthy_seurat  <- PercentageFeatureSet(aa_healthy_seurat, features = cycle.genes, col.name = "percent.cycle")
aa_healthy_seurat@meta.data$barcode   <- colnames(aa_healthy_seurat)
aa_healthy_seurat  <- aa_healthy_seurat %>% getQC()

## Put TCR
tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_barcode.txt")
aa_healthy_seurat <- mergeTCRtoSeurat(seurat_object = aa_healthy_seurat, tcr_df = tot_barcode)

## Get latents, use them for UMAPs, clustering
aa_healthy_seurat %>% getScviInput()
latents <- fread("results/scvi/latents_aa_healthy.txt")

## Now run python/run_scvi_example.py to obtain the latent dimensions
aa_healthy_seurat <- aa_healthy_seurat %>% putLatentsSeurat()
aa_healthy_seurat <- aa_healthy_seurat %>% getLatentUMAP() %>% getLatentClustering()

## Scale data
clonality_genes <- getClonalityGenes(aa_healthy_seurat)
unwanted_genes  <- getUnwantedGenes(aa_healthy_seurat)
aa_healthy_seurat <- aa_healthy_seurat %>% preprocessSeurat(cells.to.use = colnames(aa_healthy_seurat))

Idents(aa_healthy_seurat) <- aa_healthy_seurat$RNA_snn_res.0.2
aa_healthy_seurat$cluster <- Idents(aa_healthy_seurat)

## Get SingleR
aa_healthy_seurat <- aa_healthy_seurat %>% getSingler()

## Get DE-genes (pseudobulked)
aa_healthy_deg_markers     <- lapply(unique(aa_healthy_seurat$cluster), getDEGbyCluster)
aa_healthy_cluster_markers <- FindPseudoBulkMarkers(aa_healthy_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & abs(ave_logFC) > 1)

## Get pathways
aa_healthy_deg_pathways <- aa_healthy_deg_markers %>% getPathways()
aa_healthy_cluster_pathways <- aa_healthy_cluster_markers %>% getPathways()

## Fig2 onwards
scrnaseq_files <- lapply(public_dirs, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
public_seurat  <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
public_seurat  <- PercentageFeatureSet(public_seurat, pattern = "^MT-", col.name = "percent.mt")
public_seurat  <- PercentageFeatureSet(public_seurat, pattern = "^RP", col.name = "percent.ribo")
public_seurat  <- PercentageFeatureSet(public_seurat, features = cycle.genes, col.name = "percent.cycle")
public_seurat@meta.data$barcode   <- colnames(public_seurat)
public_seurat  <- public_seurat %>% getQC()

## Put TCR
public_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/public_barcode.txt")
public_seurat <- mergeTCRtoSeurat(seurat_object = public_seurat, tcr_df = public_barcode)

## Merge public_seurat
tot_seurat <- merge(aa_healthy_seurat, public_seurat)

## Get latents, scale
tot_seurat %>% getScviInput()
latents <- fread("results/scvi/latents_total.txt")
tot_seurat <- tot_seurat %>% putLatentsSeurat() %>% getLatentUMAP() %>% getLatentClustering()
tot_seurat <- tot_seurat %>% preprocessSeurat(cells.to.use = colnames(tot_seurat))

## Identify celltypes with over-clustering and SingleR
Idents(tot_seurat) <- tot_seurat$RNA_snn_res.3
tot_seurat$cluster <- Idents(tot_seurat)

cd4_seurat     <- cd4_seurat     %>% preprocessSeurat(cells.to.use = colnames(cd4_seurat))
cd8_seurat     <- cd8_seurat     %>% preprocessSeurat(cells.to.use = colnames(cd8_seurat))
nk_seurat      <- nk_seurat      %>% preprocessSeurat(cells.to.use = colnames(nk_seurat))
myeloid_seurat <- myeloid_seurat %>% preprocessSeurat(cells.to.use = colnames(myeloid_seurat))
b_seurat       <- b_seurat       %>% preprocessSeurat(cells.to.use = colnames(b_seurat))
hspc_seurat    <- hspc_seurat    %>% preprocessSeurat(cells.to.use = colnames(hspc_seurat))
