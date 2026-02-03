## R script

library(Seurat)
library(rhdf5)

# Integrate epidermal gene expression matrices (raw counts).
# DC1, DC2, DC3, DC4, MP1, NP1, NP2, NP3, NP4, NP5, and PC1 correspond to
# ATML, SHOOT1, SHOOT4, LEAF3, COTYLEDON, DARK SHOOT1, DARK SHOOT2,
# D2L 1h SHOOT, D2L 6h SHOOT, LIGHT SHOOT, and MATURE LEAF datasets, respectively.
dc1 = readRDS('./data/epidermis.DC1.rds') 
dc2 = readRDS('./data/epidermis.DC2.rds')
dc3 = readRDS('./data/epidermis.DC3.rds')
dc4 = readRDS('./data/epidermis.DC4.rds')
mp1 = readRDS('./data/epidermis.MP1.rds')
np1 = readRDS('./data/epidermis.NP1.rds')
np2 = readRDS('./data/epidermis.NP2.rds')
np3 = readRDS('./data/epidermis.NP3.rds')
np4 = readRDS('./data/epidermis.NP4.rds')
np5 = readRDS('./data/epidermis.NP5.rds')
pc1 = readRDS('./data/epidermis.PC1.rds')

# Merge datasets and normalize using SCTransform,
# regressing out batch effects by including original sample IDs as a covariate
merged = merge(dc1, c(dc2, dc3, dc4, mp1, np1, np2, np3, np4, np5, pc1))
merged = SCTransform(merged, vars.to.regress = "orig.ident")

# Save the normalized expression matrix in HDF5 format for MATLAB
expr_matrix = GetAssayData(merged)
h5createFile('at.epidermis.merged.h5')
h5createDataset('at.epidermis.merged.h5', 'x',
                dim(expr_matrix),
                chunk = c(1000, dim(expr_matrix)[2]))
h5write(as.matrix(expr_matrix), 'at.epidermis.merged.h5', 'x')
h5write(row.names(expr_matrix), 'at.epidermis.merged.h5', 'genenames')
h5write(colnames(expr_matrix), 'at.epidermis.merged.h5', 'cellnames')
savehistory('at.epidermis.merged.h5.rhistory')

