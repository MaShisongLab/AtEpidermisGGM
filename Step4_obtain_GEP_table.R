# R Script

# Read the MCL output file and convert it to a GEP table.
# Only gene modules (GEPs) with 15 or more genes are retained.
# Reference:
# https://github.com/MaShisongLab/SingleCellGGM_Network_Analysis_Tutorial

MCLresult <- readLines("out.ggm.epidermis.merged.edges.txt.I17")
gep_size_limit = 15
GGM_Modules <- as.data.frame(matrix(0, nrow = 0, ncol = 2))

for (i in 1:length(MCLresult)) {
  Module <- cbind(as.character(unlist(strsplit(MCLresult[i], '\t'))), i)
  if (nrow(Module) >= gep_size_limit) {
    GGM_Modules <- rbind(GGM_Modules, Module)
    j = i
  }
}

colnames(GGM_Modules) <- c("Gene_ID", "GEP_ID")
fmt = c("%01d","%02d","%03d","%04d","%05d")[min((floor(log(j, 10)) + 1), 5)]
GGM_Modules$Module_GEP_ID <- paste("M",
                                  sprintf(fmt, as.numeric(GGM_Modules$GEP_ID)),
                                  sep = "")

# Save the GEP table
write.table(GGM_Modules,
            "AtEpidermisGGM.modules.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)

