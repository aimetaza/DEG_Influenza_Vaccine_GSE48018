
# Transcriptomic Signature of Early Antiviral Response Following Influenza Vaccination  
## Differential Expression & Functional Enrichment Analysis (GSE48018)

## 📌 Overview
This project analyzes transcriptomic changes in peripheral blood samples following influenza vaccination using publicly available dataset GSE48018. The objective is to identify Differentially Expressed Genes (DEGs) between Day 1 (24 hours post-vaccination) and Day 0 (baseline), and to interpret their biological significance through enrichment analysis.

## 🧬 Dataset
- Accession: GSE48018
- Platform: Illumina HumanHT-12 v4.0 expression beadchip
- Samples: 431
- Comparison: Day 1 vs Day 0

## 🔬 Analysis Workflow
1. Data retrieval using GEOquery
2. Data preprocessing and quality inspection
3. Experimental design and contrast definition (Day1 – Day0)
4. Differential expression analysis using limma
5. DEG classification (UP/DOWN)
6. Visualization:
   - Volcano Plot
   - Heatmap (Top 50 DEGs)
   - UMAP (All genes, variable genes, DEGs)
7. Functional enrichment analysis:
   - Gene Ontology (GO)
   - KEGG Pathway

## 📊 Key Results
- Significant DEGs (FDR < 0.05):
  - 525 upregulated genes
  - 141 downregulated genes
- Upregulated genes are enriched in:
  - Response to virus
  - Type I interferon signaling
  - Antiviral defense mechanisms
- KEGG enrichment highlights:
  - Influenza A pathway
  - Interferon-related viral infection pathways
  - Antigen processing and presentation

These findings indicate strong activation of innate antiviral immune response 24 hours after vaccination.

## 🧠 Biological Interpretation
The transcriptomic profile at Day 1 is dominated by interferon-stimulated genes (ISGs), suggesting activation of early innate immune defense mechanisms prior to the development of adaptive immunity.

## 📂 Repository Structure
- `script_analysis.R` → Full R analysis script
- `Week3_Report.md / PDF` → Written report
- `plots/` → All generated figures

## ⚙️ Tools & Packages
- R
- limma
- GEOquery
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- pheatmap
- umap

## 👩‍💻 Author
Sarah Fahima Mumtaza  
Bioinformatics Undergraduate Student

