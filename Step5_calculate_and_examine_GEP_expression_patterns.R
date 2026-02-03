## R script

# Load the SHOOT1 dataset and GEP information
library(Seurat)
shoot1 = readRDS('./data/shoot1.rds')
gep = read.table('AtEpidermisGGM.modules.txt', header = TRUE)
row.names(gep) = gep$Gene_ID

# Calculate GEP expression by averaging gene expression within each GEP
g = GetAssayData(shoot1)
e = aggregate(g, by = list(gep[row.names(g), 2]), mean)
row.names(e) = e[,1]
e = e[, -1]

# Add GEP expression as a new assay and visualize
shoot1[['g']] = CreateAssayObject(e)
FeaturePlot(shoot1, order = TRUE, raster = TRUE, 'g_M020')
FeaturePlot(shoot1, order = TRUE, raster = TRUE, 'g_M022')

