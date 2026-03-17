# RNA-seq Analysis Workflow
## GBA1 Deficiency in Parkinson's Disease — GSE315738

---

## Dataset

| Field | Value |
|-------|-------|
| Accession | GSE315738 |
| Samples | 54 (neurons + astrocytes from iPSC clones) |
| Organism | Homo sapiens |
| Platform | Illumina NovaSeq X Plus |
| Status | Public — no linked publication |

---

## Pipeline Overview

| Stage | Tool | Input | Output |
|-------|------|-------|--------|
| 1. Get raw data | SRA Toolkit | GEO accession | .fastq.gz files |
| 2. Quality control | FastQC / MultiQC | .fastq.gz | HTML QC reports |
| 3. Trimming | Trimmomatic | raw .fastq.gz | trimmed .fastq.gz |
| 4. Alignment | STAR | trimmed reads + GRCh38 | sorted .bam files |
| 5. Read counting | featureCounts | .bam + GTF | count matrix |
| 6. Diff. expression | DESeq2 (R) | count matrix | DEG tables |
| 7. Visualisation | ggplot2 / pheatmap | DESeq2 results | PCA, volcano, heatmap |
| 8. Pathway enrichment | clusterProfiler | DEG gene lists | GO/KEGG results |

---

## Three Comparisons

| # | Comparison | Design | Purpose |
|---|-----------|--------|---------|
| 1 | Neurons vs Astrocytes | `~ cell_type` | Cell-type specific GBA1 effects |
| 2 | GBA1 het vs hom vs ctrl (neurons only) | `~ genotype` | Mutation dose response |
| 3 | Patient vs healthy donor | `~ cell_type + genotype` | Disease vs health, controlling cell type |

---

## Stage 1 — Get Raw Data

Go to GSE315738 on GEO → scroll to bottom → click **SRA Run Selector** → select all 54 samples → download `SRR_Acc_List.txt`

```bash
conda install -c bioconda sra-tools

# Download
prefetch --option-file SRR_Acc_List.txt -O ./raw_data/

# Convert to fastq (2 files per sample)
fastq-dump --split-files --gzip ./raw_data/SRR* --outdir ./raw_data/
```

> `--split-files` produces `_1.fastq.gz` and `_2.fastq.gz` per sample.
> 108 files total (54 × 2). Expect ~200–400 GB disk space.

---

## Stage 2 — Quality Control

```bash
conda install -c bioconda fastqc multiqc

# Run FastQC on all files
fastqc ./raw_data/*.fastq.gz -o ./qc_reports/ -t 8

# Combine all reports into one HTML
multiqc ./qc_reports/ -o ./qc_summary/
```

Open `qc_summary/multiqc_report.html` and check:

| Check | Good | Action if bad |
|-------|------|---------------|
| Per base quality | Phred >28 (green) | Trim with Trimmomatic |
| Adapter content | Flat line at zero | Trim with Trimmomatic |
| Sequence duplication | High is OK for RNA-seq | No action needed |
| GC content | Bell curve ~50% | Investigate contamination |

**If adapters are flagged, trim:**

```bash
conda install -c bioconda trimmomatic

trimmomatic PE \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_trimmed.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_trimmed.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
```

> Run this as a loop over all 54 samples.

---

## Stage 3 — Alignment

```bash
conda install -c bioconda star

# Download genome reference (once)
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip *.gz

# Build genome index (once — ~30 min, needs 32 GB RAM)
STAR --runMode genomeGenerate \
     --genomeDir ./genome_index/ \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.109.gtf \
     --runThreadN 8

# Align each sample (run for all 54)
STAR --runThreadN 8 \
     --genomeDir ./genome_index/ \
     --readFilesIn sample_R1_trimmed.fastq.gz sample_R2_trimmed.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ./aligned/sample_name_
```

> Check `Log.final.out` per sample — alignment rate should be >80%.
> Use an HPC cluster or SLURM job array for all 54 samples in parallel.

---

## Stage 4 — Read Counting

```bash
conda install -c bioconda subread

featureCounts \
    -T 8 \
    -p \
    -a Homo_sapiens.GRCh38.109.gtf \
    -o ./counts/all_samples_counts.txt \
    ./aligned/*.bam
```

> `-p` is required for paired-end data.
> Output has 6 annotation columns then one count column per sample — strip the first 6 before loading into R.

---

## Stage 5 — Differential Expression (R)

**Install packages:**

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ggplot2", "pheatmap",
                       "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
```

**Load data and build metadata:**

```r
library(DESeq2)

counts_raw <- read.table('./counts/all_samples_counts.txt',
                          header=TRUE, skip=1, row.names=1)
counts <- counts_raw[, 6:ncol(counts_raw)]  # drop annotation columns

metadata <- data.frame(
  sample    = colnames(counts),
  cell_type = c(...),   # 'neuron' or 'astrocyte' for each sample
  genotype  = c(...),   # 'control_healthy','GBA_het','GBA_hom','isogenic_ctrl'
  clone     = c(...),
  row.names = colnames(counts)
)

# This MUST return TRUE before continuing
all(colnames(counts) == rownames(metadata))
```

### Comparison 1 — Neurons vs Astrocytes

```r
ctrl   <- metadata[metadata$genotype == 'control_healthy', ]
dds_ct <- DESeqDataSetFromMatrix(
    countData = counts[, rownames(ctrl)],
    colData   = ctrl,
    design    = ~ cell_type
)
dds_ct <- DESeq(dds_ct)
res_ct <- results(dds_ct, contrast=c('cell_type','neuron','astrocyte'), alpha=0.05)
```

### Comparison 2 — GBA1 dose response (neurons only)

```r
neuro    <- metadata[metadata$cell_type == 'neuron' &
                     metadata$genotype %in% c('GBA_het','GBA_hom','isogenic_ctrl'), ]
dds_dose <- DESeqDataSetFromMatrix(
    countData = counts[, rownames(neuro)],
    colData   = neuro,
    design    = ~ genotype
)
dds_dose$genotype <- relevel(dds_dose$genotype, ref='isogenic_ctrl')
dds_dose <- DESeq(dds_dose)

res_het <- results(dds_dose, contrast=c('genotype','GBA_het','isogenic_ctrl'), alpha=0.05)
res_hom <- results(dds_dose, contrast=c('genotype','GBA_hom','isogenic_ctrl'), alpha=0.05)
```

### Comparison 3 — Patient vs healthy donor

```r
pat     <- metadata[metadata$genotype %in% c('GBA_het','control_healthy'), ]
dds_pat <- DESeqDataSetFromMatrix(
    countData = counts[, rownames(pat)],
    colData   = pat,
    design    = ~ cell_type + genotype
)
dds_pat$genotype <- relevel(dds_pat$genotype, ref='control_healthy')
dds_pat <- DESeq(dds_pat)
res_pat <- results(dds_pat, contrast=c('genotype','GBA_het','control_healthy'), alpha=0.05)
```

**Extract and save significant genes:**

```r
get_sig <- function(res, lfc=1.0) {
    df      <- as.data.frame(res)
    df$gene <- rownames(df)
    sig     <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) >= lfc, ]
    sig[order(sig$padj), ]
}

sig_ct  <- get_sig(res_ct)
sig_het <- get_sig(res_het)
sig_hom <- get_sig(res_hom)
sig_pat <- get_sig(res_pat)

write.csv(sig_ct,  'DEGs_neuron_vs_astrocyte.csv')
write.csv(sig_het, 'DEGs_GBA_het_vs_ctrl.csv')
write.csv(sig_hom, 'DEGs_GBA_hom_vs_ctrl.csv')
write.csv(sig_pat, 'DEGs_patient_vs_healthy.csv')
```

> `padj < 0.05` = 5% false discovery rate. `log2FoldChange >= 1` = at least 2x change.
> Always use `padj`, never the raw `pvalue`.

---

## Stage 6 — Visualisation

**PCA:**

```r
library(ggplot2)

vsd  <- vst(dds_pat, blind=FALSE)
pca  <- plotPCA(vsd, intgroup=c('cell_type','genotype'), returnData=TRUE)
pctv <- round(100 * attr(pca, 'percentVar'))

ggplot(pca, aes(PC1, PC2, color=genotype, shape=cell_type)) +
    geom_point(size=4) +
    xlab(paste0('PC1: ', pctv[1], '% variance')) +
    ylab(paste0('PC2: ', pctv[2], '% variance')) +
    theme_bw()
ggsave('PCA_all_samples.png', dpi=150)
```

**Volcano plot:**

```r
library(EnhancedVolcano)

EnhancedVolcano(as.data.frame(res_het),
    lab      = rownames(res_het),
    x        = 'log2FoldChange',
    y        = 'padj',
    title    = 'GBA1 het vs isogenic control (neurons)',
    pCutoff  = 0.05,
    FCcutoff = 1.0)
ggsave('volcano_GBA_het_neurons.png', dpi=150)
```

**Heatmap — top 50 DEGs:**

```r
library(pheatmap)

top50 <- head(sig_het[order(sig_het$padj), ], 50)
mat   <- assay(vsd)[top50$gene, ]
mat   <- t(scale(t(mat)))  # z-score each gene

pheatmap(mat,
    annotation_col = data.frame(cell_type=pat$cell_type,
                                genotype=pat$genotype,
                                row.names=rownames(pat)),
    show_colnames  = FALSE,
    fontsize_row   = 7,
    filename       = 'heatmap_top50.png')
```

---

## Stage 7 — Pathway Enrichment

```r
library(clusterProfiler)
library(org.Hs.eg.db)

to_entrez <- function(symbols) {
    bitr(symbols, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
}

run_GO <- function(sig_df, label) {
    ids <- to_entrez(sig_df$gene)
    ego <- enrichGO(gene=ids$ENTREZID, OrgDb=org.Hs.eg.db, ont='BP',
                    pAdjustMethod='BH', pvalueCutoff=0.05, readable=TRUE)
    write.csv(as.data.frame(ego), paste0('GO_BP_', label, '.csv'))
    dotplot(ego, showCategory=20, title=paste('GO BP -', label))
    ggsave(paste0('GO_dotplot_', label, '.png'), dpi=150, height=10)
    return(ego)
}

go_het <- run_GO(sig_het, 'GBA_het_neurons')
go_hom <- run_GO(sig_hom, 'GBA_hom_neurons')
go_pat <- run_GO(sig_pat, 'patient_vs_healthy')

# Compare het vs hom pathways side by side
cmp <- compareCluster(
    list(het=to_entrez(sig_het$gene)$ENTREZID,
         hom=to_entrez(sig_hom$gene)$ENTREZID),
    fun='enrichGO', OrgDb=org.Hs.eg.db, ont='BP'
)
dotplot(cmp, title='Pathway comparison: GBA1 het vs hom')
ggsave('GO_comparison_het_hom.png', dpi=150, height=12)
```

> If `enrichGO` returns nothing, try lowering `pvalueCutoff` to `0.1` or switch `ont` to `'MF'` or `'CC'`.

---

## Expected Output Files

| File | Format | Description |
|------|--------|-------------|
| `DEGs_neuron_vs_astrocyte.csv` | CSV | Neuron vs astrocyte DEGs |
| `DEGs_GBA_het_vs_ctrl.csv` | CSV | GBA1 het vs isogenic control |
| `DEGs_GBA_hom_vs_ctrl.csv` | CSV | GBA1 hom vs isogenic control |
| `DEGs_patient_vs_healthy.csv` | CSV | Patient vs healthy donor |
| `GO_BP_GBA_het_neurons.csv` | CSV | GO Biological Process — GBA het |
| `GO_comparison_het_hom.png` | PNG | Side-by-side pathway dot plot |
| `PCA_all_samples.png` | PNG | PCA by genotype and cell type |
| `volcano_GBA_het_neurons.png` | PNG | Volcano plot for het vs ctrl |
| `heatmap_top50.png` | PNG | Top 50 DEGs heatmap |

---

## Key Thresholds

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `padj` cutoff | 0.05 | 5% false discovery rate |
| `log2FoldChange` | >= 1.0 | At least 2x change |
| VST transform | `blind=FALSE` | Uses design for accurate normalisation |
| Reference genotype | `isogenic_ctrl` / `control_healthy` | All fold changes relative to this |

---

> **Note:** Stages 1–4 (download through counting) should be run on a Linux HPC cluster.
> STAR alignment for 54 samples will take several hours on a laptop.
> Stages 5–7 run fine in R on a standard laptop with 16 GB RAM.

---

---

# Stage 9 — Hub Gene Analysis (WGCNA)

> This stage extends the DESeq2 results into a gene co-expression network,
> identifies the most biologically influential genes (hub genes), and validates
> them against protein-protein interaction databases. This is what separates
> a DEG list from a publishable mechanistic finding.

---

## What is WGCNA?

Weighted Gene Co-expression Network Analysis (WGCNA) groups genes into modules
based on how correlated their expression is across samples. Genes that go up and
down together are likely regulated together or part of the same pathway.

A **hub gene** is a gene with the highest connectivity within its module — it
correlates strongly with many other genes and is therefore likely a master
regulator or key driver of that module's biological function.

---

## Install WGCNA

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("WGCNA")

install.packages(c("igraph", "ggplot2", "dplyr", "reshape2"))
```

> WGCNA needs at least 16 GB RAM. If R crashes during network construction,
> reduce the number of input genes (see filtering step below).

---

## Step 1 — Prepare input matrix

WGCNA takes a normalised expression matrix as input, NOT raw counts.
Use the VST-transformed matrix from DESeq2. Rows = samples, columns = genes.

```r
library(WGCNA)
library(DESeq2)

# Allow multi-threading
allowWGCNAThreads()

# Use VST matrix from DESeq2 (already computed in Stage 6)
# vsd was created with: vsd <- vst(dds_pat, blind=FALSE)
expr_mat <- t(assay(vsd))   # transpose: rows=samples, cols=genes

# Filter to top 5000 most variable genes to reduce compute time
gene_var  <- apply(expr_mat, 2, var)
top_genes <- names(sort(gene_var, decreasing=TRUE))[1:5000]
expr_mat  <- expr_mat[, top_genes]

dim(expr_mat)   # should be: samples x 5000
```

> You can increase to 10,000 genes for a more comprehensive network if RAM allows.
> Below 3,000 genes the network becomes too sparse to be meaningful.

---

## Step 2 — Check for outlier samples

Bad samples will distort the entire network. Remove them before building.

```r
# Cluster samples to spot outliers
sample_tree <- hclust(dist(expr_mat), method="average")
plot(sample_tree, main="Sample clustering to detect outliers",
     sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

# Draw a cut line — any sample branching off far above the rest is an outlier
abline(h=80, col="red")   # adjust height based on your plot

# Remove outlier samples if any
# good_samples <- cutreeStatic(sample_tree, cutHeight=80, minSize=10)
# expr_mat <- expr_mat[good_samples == 1, ]
```

---

## Step 3 — Pick the soft threshold power

WGCNA requires a power parameter that determines how strongly correlations are
weighted. The goal is to find the lowest power where the network approximates
a scale-free topology (R² > 0.85).

```r
# Test a range of powers
powers <- c(1:10, seq(12, 20, by=2))

sft <- pickSoftThreshold(
    expr_mat,
    powerVector  = powers,
    verbose      = 5,
    networkType  = "signed"
)

# Plot to choose power
par(mfrow=c(1,2))

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit (R²)",
     type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.85, col="red")   # aim for R² > 0.85

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
```

> Pick the **lowest power** where the R² line crosses 0.85.
> For signed networks with RNA-seq data this is typically between 12 and 18.
> If no power reaches 0.85, use 12 and proceed — it just means noisier data.

---

## Step 4 — Build the co-expression network

This is the main WGCNA step. It builds the network and assigns genes to modules.
Each module gets a colour name (e.g. "turquoise", "blue", "brown").

```r
# Set your chosen power here
SOFT_POWER <- 14   # replace with value from Step 3

net <- blockwiseModules(
    expr_mat,
    power             = SOFT_POWER,
    networkType       = "signed",
    TOMType           = "signed",
    minModuleSize     = 30,
    mergeCutHeight    = 0.25,
    numericLabels     = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs          = TRUE,
    saveTOMFileBase   = "TOM",
    verbose           = 3
)

# How many modules were found and how big are they?
table(net$colors)

# Plot the dendrogram with module colours
plotDendroAndColors(
    net$dendrograms[[1]],
    net$colors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang         = 0.03,
    addGuide     = TRUE,
    guideHang    = 0.05
)
```

> `minModuleSize=30` means a module needs at least 30 genes.
> `mergeCutHeight=0.25` merges very similar modules — increase to 0.3 if
> you get too many small modules.

---

## Step 5 — Relate modules to traits

Find which modules are most associated with your conditions of interest
(genotype, cell type). These are the biologically relevant modules.

```r
# Build a numeric trait matrix from your metadata
# 1 = has the trait, 0 = does not
trait_mat <- data.frame(
    row.names  = rownames(expr_mat),
    GBA_het    = as.numeric(pat$genotype == "GBA_het"),
    GBA_hom    = as.numeric(pat$genotype == "GBA_hom"),
    is_neuron  = as.numeric(pat$cell_type == "neuron"),
    is_patient = as.numeric(pat$genotype != "control_healthy")
)

# Correlate module eigengenes with traits
module_eigengenes <- net$MEs

module_trait_cor <- cor(module_eigengenes, trait_mat, use="p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(expr_mat))

# Plot the correlation heatmap
textMatrix <- paste0(
    signif(module_trait_cor, 2), "\n(",
    signif(module_trait_pval, 1), ")"
)
dim(textMatrix) <- dim(module_trait_cor)

par(mar=c(6, 8.5, 3, 3))
labeledHeatmap(
    Matrix    = module_trait_cor,
    xLabels   = colnames(trait_mat),
    yLabels   = colnames(module_eigengenes),
    ySymbols  = colnames(module_eigengenes),
    colorLabels = FALSE,
    colors    = blueWhiteRed(50),
    textMatrix= textMatrix,
    setStdMargins = FALSE,
    cex.text  = 0.5,
    zlim      = c(-1, 1),
    main      = "Module-trait relationships"
)
```

> Modules with high correlation (close to 1 or -1) and low p-value
> are the ones to focus on. Write down their colour names.

---

## Step 6 — Identify hub genes

Hub genes are defined by two scores:
- **Module Membership (MM):** how correlated a gene is with its module's eigengene
- **Gene Significance (GS):** how correlated a gene is with the trait of interest

Hub genes have high MM (>0.8) AND high GS (>0.2).

```r
# Pick your module of interest — replace "turquoise" with your top module colour
MODULE_OF_INTEREST <- "turquoise"
TRAIT_OF_INTEREST  <- "GBA_het"   # column name from trait_mat

# Gene significance and module membership
gene_sig <- as.data.frame(cor(expr_mat, trait_mat[, TRAIT_OF_INTEREST], use="p"))
colnames(gene_sig) <- "GS"
gene_sig$GS_pval <- corPvalueStudent(as.numeric(gene_sig$GS), nrow(expr_mat))

mm <- as.data.frame(cor(expr_mat, module_eigengenes, use="p"))
mm_col <- paste0("MM", MODULE_OF_INTEREST)

# Filter to genes in the module of interest
module_genes <- names(net$colors[net$colors == MODULE_OF_INTEREST])
mm_module    <- mm[module_genes, mm_col, drop=FALSE]
gs_module    <- gene_sig[module_genes, , drop=FALSE]

hub_df <- data.frame(
    gene             = module_genes,
    ModuleMembership = mm_module[, 1],
    GeneSignificance = gs_module$GS,
    GS_pval          = gs_module$GS_pval
)

# Hub genes: high MM and high GS
hub_genes <- hub_df[abs(hub_df$ModuleMembership) > 0.8 &
                    abs(hub_df$GeneSignificance)  > 0.2, ]
hub_genes <- hub_genes[order(-abs(hub_genes$ModuleMembership)), ]

cat("Hub genes found:", nrow(hub_genes), "\n")
print(head(hub_genes, 20))

write.csv(hub_genes, paste0("hub_genes_", MODULE_OF_INTEREST, "_", TRAIT_OF_INTEREST, ".csv"), row.names=FALSE)

# Scatter plot of MM vs GS
plot(abs(mm_module[, 1]), abs(gs_module$GS),
     xlab=paste("Module Membership in", MODULE_OF_INTEREST, "module"),
     ylab=paste("Gene Significance for", TRAIT_OF_INTEREST),
     main="Module membership vs gene significance",
     pch=20, col="steelblue")
abline(v=0.8, h=0.2, col="red", lty=2)
```

---

## Step 7 — Validate hub genes with PPI network (STRING)

Protein-protein interaction networks confirm that your hub genes actually
interact biologically, not just statistically.

```r
install.packages("STRINGdb")
library(STRINGdb)

# Connect to STRING database (human = 9606)
string_db <- STRINGdb$new(
    version    = "11.5",
    species    = 9606,
    score_threshold = 400,   # medium confidence interactions
    input_directory = ""
)

# Map your hub gene symbols to STRING IDs
hub_mapped <- string_db$map(
    hub_genes,
    "gene",
    removeUnmappedRows = TRUE
)

# Plot the PPI network
string_db$plot_network(hub_mapped$STRING_id)

# Get interaction scores between hub genes
interactions <- string_db$get_interactions(hub_mapped$STRING_id)
write.csv(interactions, "hub_gene_PPI_interactions.csv", row.names=FALSE)
```

> Export the network and open it in **Cytoscape** (free software) for a
> publication-quality figure. Use the StringApp plugin inside Cytoscape
> to import directly from STRING.

---

## Step 8 — Cross-validate hub genes with DESeq2 results

Hub genes should also be differentially expressed. This step confirms your
hub genes are both statistically significant AND biologically central.

```r
# Check how many hub genes are also in your DEG lists
hub_gene_names <- hub_genes$gene

overlap_het <- intersect(hub_gene_names, sig_het$gene)
overlap_hom <- intersect(hub_gene_names, sig_hom$gene)
overlap_pat <- intersect(hub_gene_names, sig_pat$gene)

cat("Hub genes also DEGs in GBA het:    ", length(overlap_het), "\n")
cat("Hub genes also DEGs in GBA hom:    ", length(overlap_hom), "\n")
cat("Hub genes also DEGs in patient:    ", length(overlap_pat), "\n")

# These overlapping genes are your top candidates
final_candidates <- hub_genes[hub_genes$gene %in% union(overlap_het, overlap_pat), ]
final_candidates <- merge(final_candidates, sig_het[, c("gene","log2FoldChange","padj")],
                          by="gene", all.x=TRUE)

final_candidates <- final_candidates[order(-abs(final_candidates$ModuleMembership)), ]
write.csv(final_candidates, "final_hub_gene_candidates.csv", row.names=FALSE)
print(final_candidates)
```

---

## Updated Output Files

| File | Description |
|------|-------------|
| `hub_genes_turquoise_GBA_het.csv` | Hub genes in top module for GBA het trait |
| `hub_gene_PPI_interactions.csv` | STRING protein-protein interactions |
| `final_hub_gene_candidates.csv` | Hub genes that are also DEGs — top candidates |
| `TOM-block.1.RData` | Saved topological overlap matrix (reuse to avoid recomputing) |

---

## Updated Pipeline Overview

```
Raw FASTQ
    ↓ STAR alignment
BAM files
    ↓ featureCounts
Count matrix
    ↓ DESeq2
DEG tables (padj < 0.05, |log2FC| >= 1)
    ↓              ↓
Volcano        WGCNA
Heatmap        Co-expression modules
PCA            Module-trait correlation
               Hub genes (high MM + high GS)
                    ↓
               STRING PPI validation
                    ↓
               Cross-validate with DEGs
                    ↓
               Final hub gene candidates
```

---

## Key WGCNA Thresholds

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Soft power | 12–18 | Network scale-free topology fit |
| R² target | > 0.85 | Minimum acceptable fit |
| Min module size | 30 genes | Modules smaller than this are merged |
| Module membership | > 0.8 | Gene is strongly part of its module |
| Gene significance | > 0.2 | Gene correlates with the trait |
| STRING score | > 400 | Medium confidence PPI interaction |

---

> **Note:** Save the TOM (Topological Overlap Matrix) file after Step 4.
> It takes 30–60 minutes to compute and you do not want to rerun it.
> All subsequent steps can reload it with `load("TOM-block.1.RData")`.
