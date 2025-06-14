if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Hs.eg.db", "gplots", 
                       "RColorBrewer", "NMF", "BiasedUrn"), force=TRUE)
BiocManager::install("biomaRt")

library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db) #mouse
library(org.Hs.eg.db) #human 
library(biomaRt)
library(gplots)
library(RColorBrewer)
library(NMF)



mypath <- "/Users/siber/Desktop/rnaseq_DEGs/countdata"
setwd(mypath)

# List all text files
txt_files_ls <- list.files(path = mypath, pattern = "\\.txt$", full.names = TRUE)

# Read the first file
combined_set <- read.delim(txt_files_ls[1], stringsAsFactors = FALSE)

# Loop through the rest and merge
for (i in 2:length(txt_files_ls)) {
  next_file <- read.delim(txt_files_ls[i], stringsAsFactors = FALSE)
  combined_set <- merge(combined_set, next_file, by = "Gene")
}

#Creating Sample Info sheet
sampleInfo <- data.frame(
  FileName = paste0(rep(c("siNC", "siCTCF"), each = 3), 1:3),
  SampleName = paste0(rep(c("siNC", "siCTCF"), each = 3), 1:3),
  CellType = rep("HESC", 6),
  Status = rep(c("control", "treated"), each = 3),
  stringsAsFactors = FALSE
)
####
countdata <- combined_set[,-(1)]
rownames(countdata) <- combined_set[,1]
table(colnames(countdata)==sampleInfo$SampleName)


y <- DGEList(countdata)
y

# setting groups
group <- paste(sampleInfo$CellType, sampleInfo$Status, sep = ".")
group <- factor(group)

y$samples$group <- group
y$samples



# Set group and design matrix
group <- factor(paste(sampleInfo$CellType, sampleInfo$Status, sep = "."))
y$samples$group <- group

# Filter low-expressed genes
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize
y <- calcNormFactors(y)

# Design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Voom transformation
v <- voom(y, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design)

# Create contrast: treated vs control
cont.matrix <- makeContrasts(
  HESC.controlVsTreated = HESC.control - HESC.treated,
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Glimma MDS Plot
labels <- paste(sampleInfo$SampleName, sampleInfo$Status)
glMDSPlot(y, labels = labels, groups = group, folder = "mds_GSE299180")

# TopTable
topTable_DEG <- topTable(fit2, number = Inf, sort.by = "P")

# Plot volcano
with(top_all, {
  plot(logFC, -log10(P.Value), pch = 20, main = "Volcano Plot with Gene Names", xlim = c(-5, 5))
  abline(h = -log10(0.05), col = "red", lty = 2)
  abline(v = c(-1, 1), col = "blue", lty = 2)
  
  # Highlight DEGs with abs(logFC) > 1 & adj.P.Val < 0.05
  sig <- which(abs(logFC) > 1 & adj.P.Val < 0.05)
  text(logFC[sig], -log10(P.Value)[sig], labels = rownames(top_all)[sig], cex = 0.6, pos = 3, col = "darkgreen")
})
# Heatmap of top 100 DEGs
top100 <- topTable_DEG[1:100, ]
logCPM <- cpm(y, log = TRUE)
heat_data <- logCPM[rownames(top100), ]

# Optional color scaling
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)

heatmap.2(heat_data,
          scale = "row",
          trace = "none",
          col = rev(morecols(50)),
          ColSideColors = c("blue", "red")[group],
          margins = c(8, 6),
          main = "Top 100 DEGs Heatmap")

# Print top 50 DEGs (sorted by adjusted p-value)
top50 <- topTable(fit2, number = 50, sort.by = "P", adjust.method = "BH")
print(top50)

# For volcano plot with gene labels
top_all <- topTable(fit2, number = Inf, sort.by = "P", adjust.method = "BH")









