library(dplyr)
library(tidyr)
library(readr)
library(Seurat)
library(ggplot2)
library(cowplot)

# Run an integrated seurat analysis
set.seed(42)


# Read in the samples

# KWP (non resistant rat cell line)
KWP <- Read10X(data.dir = "/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/KWP/outs/filtered_feature_bc_matrix")
KWP <- CreateSeuratObject(counts = KWP , project = "KWP", min.cells = 3, min.features = 200)

# KWR (non resistant rat cell line)
KWR <- Read10X(data.dir = "/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/KWR/outs/filtered_feature_bc_matrix")
KWR <- CreateSeuratObject(counts = KWR , project = "KWR", min.cells = 3, min.features = 200)

# Merge the libraries
rat_all <-  merge(KWP, y = c( KWR), add.cell.ids = c("KWP", "KWR"), project = "Rat_all")

# Add a pct mt tag
rat_all[["percent.mt"]] <- PercentageFeatureSet(rat_all, pattern = "*MT*")

# Violin plot of features, counts and mitochondrial reads
VlnPlot(rat_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)

# Plot some more qc metrics
plot1 <- FeatureScatter(rat_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rat_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Subset down the data. Looks like our data is perhaps better quality than the tutorial or at 
# least has a lot more features per cell
rat_all_s <- subset(rat_all, subset = nFeature_RNA > 200 & percent.mt < 25)


# Plot some more qc metrics
plot1 <- FeatureScatter(rat_all_s, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rat_all_s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Normalise
rat_all_s <- NormalizeData(rat_all_s, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable genes
rat_all_s  <- FindVariableFeatures(rat_all_s, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(rat_all_s), 10)
plot1 <- VariableFeaturePlot(rat_all_s)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# Looks like we have some cell cycle genes as part of our variable genes
DimHeatmap(rat_all_s, dims = c(1, 2))

# Load up the seurat cell cycle genes

# Scale the data before PCA
all.genes <- rownames(rat_all_s)
rat_all_s <- ScaleData(rat_all_s, features = all.genes)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- data.frame(`Gene name` = cc.genes.updated.2019$s.genes, stringsAsFactors = F,check.names = F)
g2m.genes <- data.frame(`Gene name` = cc.genes.updated.2019$g2m.genes, stringsAsFactors = F,check.names = F)

# Try and match the Seurat cell cycle genes with rat cell cycle genes 

rat_to_human_genes <- read_tsv("/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/reference/Rat_human_homology.tsv")%>%
  select(`Rat gene name`, `Gene name`)

s.genes <- left_join(s.genes,rat_to_human_genes)%>%
  filter(!is.na(`Rat gene name`))%>%
  filter(!duplicated(`Gene name`))

g2m.genes <- left_join(g2m.genes,rat_to_human_genes)%>%
  filter(!is.na(`Rat gene name`))%>%
  filter(!duplicated(`Gene name`))

rat_all_s <- CellCycleScoring(rat_all_s, s.features = s.genes$`Rat gene name`, 
                              g2m.features = g2m.genes$`Rat gene name`, set.ident = TRUE)

# Phase has now been added to the metadata
head(rat_all_s[[]])

# Visualize the distribution of cell cycle markers across their cell cycle classification
# Expression of these genes is generally lower than the tutorial, but roughly matches the pattern
RidgePlot(rat_all_s, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)

# Runn a PCA on cell cycle genes 
rat_all_s <- RunPCA(rat_all_s, features = c(s.genes$`Rat gene name`, g2m.genes$`Rat gene name`))

# Looks like cell cycle state was the primary identifier
DimPlot(rat_all_s)

# Regress out the cell cycle genes
rat_all_s <- ScaleData(rat_all_s, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(rat_all_s))
saveRDS(rat_all_s, "/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/analysis_intermediate/cell_cycle_regressed_out.RDS")
# Save the RDS as this takes about 30 minutes
rat_all_s <- readRDS("/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/analysis_intermediate/cell_cycle_regressed_out.RDS")
# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
rat_all_s <- RunPCA(rat_all_s, features = c(s.genes$`Rat gene name`, g2m.genes$`Rat gene name`))
DimPlot(rat_all_s)

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
rat_all_s <- RunPCA(rat_all_s, features = VariableFeatures(object = rat_all_s), nfeatures.print = 10)

# Visualise results
print(rat_all_s[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(rat_all_s, dims = 1:2, reduction = "pca")
DimPlot(rat_all_s, reduction = "pca")
DimHeatmap(rat_all_s, dims = 1, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
rat_all_s <- JackStraw(rat_all_s, num.replicate = 100)
rat_all_s <- ScoreJackStraw(rat_all_s, dims = 1:20)

# From this plot I think we have a lot more signifant dims than the tutorial dataset
# I guess that means there is a lot going on consistently 
JackStrawPlot(rat_all_s, dims = 1:20)
ElbowPlot(rat_all_s)

# Clustering of the cells
rat_all_s <- FindNeighbors(rat_all_s, dims = 1:20)
rat_all_s <- FindClusters(rat_all_s, resolution = 1)

# Look at cluster IDs of the first 5 cells
head(Idents(rat_all_s), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
rat_all_s <- RunUMAP(rat_all_s, dims = 1:20)

p1 <- DimPlot(rat_all_s, label = TRUE, pt.size = 0.5,group.by = "orig.ident") + NoLegend()
p2 <- DimPlot(rat_all_s, label = TRUE, pt.size = 0.5) + NoLegend()
# Combine the plots 
CombinePlots(list(p1,p2))

# Look at mitochondiral genes contamination
FeaturePlot(rat_all_s, features = c("Mt-nd1", "Mt-co1", "Mt-cyb"),label = TRUE)


# Find markers for every cluster compared to all remaining cells, report only the positive ones
rat_all_s.markers <- FindAllMarkers(rat_all_s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Save a markers RDS as this takes a while to calculate as well
saveRDS(rat_all_s.markers , "/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/analysis_intermediate/cluster_markers.RDS")

top_markers <- rat_all_s.markers %>% group_by(cluster) %>% top_n(n = 300, wt = avg_logFC)
top_markers_small <- rat_all_s.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top10 <- top_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(rat_all_s, features = top10$gene) + NoLegend()

# Rename clusters
new.cluster.ids <- c("KWR 1", "KWP 1", "KWP 2", "KWR 2", "KWP 3", "KWR 3", 
                     "KWP 4", "Intermediate both", "Intermediate KWP", "KWP subcluster 1", "KWR 4",
                     "KWP subcluster 2" ,"Intermediate KWR", "KWR 5", "KWP subcluster 3")
names(new.cluster.ids) <- levels(rat_all_s)
rat_all_s <- RenameIdents(rat_all_s, new.cluster.ids)
# Plot renamed clusters
DimPlot(rat_all_s, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+
  ggsave("/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/Plots/Overall_seurat_clustering.pdf")

# Find differentially expressed features between main KWP and main KWR
main.de.markers <- FindMarkers(rat_all_s, ident.1 = "KWP 1", ident.2 = "KWR 1")

top_main_genes <- main.de.markers%>%
  tibble::rownames_to_column("Gene") %>%
  arrange(-abs(avg_logFC))%>%
  write_csv("/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/Tables/Main_cluster_gene_differences.csv")

# Plot the top 5 genes
FeaturePlot(rat_all_s, features = head(top_main_genes$Gene,1),label = TRUE)

# Plot Sstr2 since it is Rod's fave gene
p1 <- DimPlot(rat_all_s, label = TRUE, pt.size = 0.5,group.by = "orig.ident") + NoLegend()
p2 <- FeaturePlot(rat_all_s, features = "Sstr2",label = F)
# Combine the plots 
CombinePlots(list(p1,p2))+
  ggsave(paste0("/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/Plots/","Sstr2", ".pdf"))

# Loop through and plot top 20 genes

for(i in 1:20){
  
  gene <- top_main_genes$Gene[i]
  
  p1 <- DimPlot(rat_all_s, label = TRUE, pt.size = 0.5,group.by = "orig.ident") + NoLegend()
  p2 <- FeaturePlot(rat_all_s, features = gene,label = F)
  # Combine the plots 
  CombinePlots(list(p1,p2))+
    ggsave(paste0("/data/cephfs/punim0010/projects/Pattison_projects/SC/Single_nuclei_Kelly_rat/Plots/",gene, ".pdf"))
  
}

