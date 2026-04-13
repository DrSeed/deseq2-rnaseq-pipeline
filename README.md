# DESeq2 RNA-seq Pipeline

> Your RNA-seq data is sitting there, untouched. You've been reading DESeq2 vignettes for a week. Stop. Run this pipeline and get your differentially expressed genes today.

## Why This Exists

Every RNA-seq experiment ends with the same question: which genes are changing? DESeq2 is the gold standard for answering that, but the jump from "I have a count matrix" to "here are my significant DEGs" trips people up more than it should.

This pipeline takes your raw count matrix and sample metadata, runs the full DESeq2 workflow, and hands you back a clean CSV of significant genes with publication-ready plots. No decision paralysis. No fiddling with parameters for three days.

## What's Actually Happening Under the Hood

1. **Low-count filtering** removes genes that are basically noise. If a gene isn't expressed in at least half your samples, it's not telling you anything useful.
2. **DESeq2 normalisation** accounts for library size differences between samples. This is critical because your "upregulated gene" might just be a sample with deeper sequencing.
3. **Log2 fold change shrinkage** (apeglm) pulls unreliable effect sizes toward zero. Without this, genes with two reads in one condition and four in another look like they have a massive fold change. They don't.
4. **Visualisation**: PCA to check if your samples actually cluster by condition (if they don't, you have bigger problems), volcano plot for the pretty picture reviewers expect, and a heatmap of top DEGs.

## The Practical Test

```bash
Rscript run_deseq2.R --counts data/counts.csv --metadata data/samples.csv --output results/
```

Your count matrix should have genes as rows, samples as columns. Your metadata needs a `condition` column with exactly two levels. That's it.

## When NOT to Use This

- **Fewer than 3 replicates per group?** DESeq2 will run but the statistics are questionable. Consider limma-voom instead.
- **Single-cell data?** This is for bulk RNA-seq. Use Scanpy or Seurat for scRNA-seq.
- **Time-series experiment?** You'll need to modify the design formula.

## The Bottom Line

The best differential expression analysis is the one you actually run. If your data has been sitting in a folder for more than a week, clone this repo and run it. You can always refine later. You can't refine nothing.
