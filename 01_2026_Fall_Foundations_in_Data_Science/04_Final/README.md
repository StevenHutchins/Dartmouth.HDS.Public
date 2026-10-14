# Final Project: COL1A1 Expression in TCGA Breast Invasive Carcinoma

Final project for HSE 711: Foundations in Data Science, Fall 2026, Dartmouth Health Data Science.

## Summary
What the project is about and the main question, in two or three sentences.
The data for this analysis were provided by the course instructor as two files: a matrix of RNA-sequencing
raw gene counts and a sample metadata file. The data originate from the Cancer Genome Atlas Breast
Invasive Carcinoma Project (TCGA-BRCA) [1, 2]. The count matrix contains RNA-sequencing raw gene
counts for 60,660 genes across 1,231 samples. Each sample is identified by its TCGA barcode and classified
by sample type as primary solid tumor, solid tissue normal, or metastatic. Sequencing depth varied widely
across samples, from about 19 to 114 million total reads (median about 58 million). Counts were not
normalized for library size, so comparisons between samples are approximate.
Visualizations include a histogram (log scale), scatter plot (log scale), violin plot (log scale), complex heatmap (log1p), and mean-difference (MA) plot.

## Data
- **Source:** Course-provided
- **Type:** public/synthetic

## Approach
- All analyses were performed in R (version 4.6.1) using RStudio and R Markdown. Data wrangling used
the tidyverse (2.0.0), mainly dplyr (1.2.1) for filtering, selecting, and summarizing and tidyr (1.3.2) for
reshaping data between wide and long formats. The histogram, scatter plot, violin plot, and mean-difference
plot were made with ggplot2 (4.0.3), with axis labels formatted by scales (1.4.0). The heatmap was
made with Bioconductor package ComplexHeatmap (2.28.0). Summary statistics (n, mean, median, standard
deviation, minimum and maximum) were calculated with a custom function using base R and exported as
a LaTeX table with knitr (1.52). Figures were saved as PDFs with ggsave and the base R pdf() device.
The final report was written in LaTeX using Overleaf.

## Key findings
- Across all analyses, the clearest pattern is that COL1A1 is expressed at higher levels in breast tumors than in normal breast tissue.
- COL1A1 is among the most highly expressed genes with one of the largest variances between tumor- and normal-tissue samples.
- FN1 shows a similar pattern and rises closely with COL1A1 across samples.
- Known limitations include lack of count normalization for sequencing depth (nearly six-fold variance across samples) and sample group size differences.

## Deliverables
- `report.pdf` or `analysis.Rmd`: the full write-up
- `output/`: charts and tables
- `slides.pdf`: presentation (if any)

## How to run it
1. Open `<project-name>.Rproj` in RStudio
2. Run `analysis.Rmd`

## Reflection
- I enjoyed learning about genetic data and how to use R to more deeply understand it. If I were to complete this project again, I would select a more unique secondary gene in hopes to find less alignment with the primary gene for more interesting visualizations. 