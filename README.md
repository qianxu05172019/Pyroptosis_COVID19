
# Single-Cell Transcriptome Analysis for Pyroptosis in COVID-19

## Project Overview

This repository hosts the R scripts and workflows used in our study on pyroptosis in COVID-19, as detailed in our published paper [View Article](https://www.sciencedirect.com/science/article/pii/S2667119022000052). These scripts cover data preprocessing, analysis, and visualization techniques that underpin our findings on the role of pyroptosis in COVID-19 severity.

## Tools and Techniques

### Data Handling and Analysis
- **Dimensionality Reduction**: Utilized **Seurat** for principal component analysis (PCA) and uniform manifold approximation and projection (UMAP).
- **Data Integration**: Employed **Harmony** to integrate datasets from multiple studies to minimize batch effects and enhance comparability.

### Visualization
- **Graphical Representations**: Leveraged **ggplot2** for generating expressive data visualizations.
- **Plot Arrangements**: Used **patchwork** to effectively arrange multiple plots in a cohesive layout.
![Pyroptosis in COVID-19](images/1-s2.0-S2667119022000052-gr1_lrg.jpg)


### Statistical and Enrichment Analysis
- **Statistical Testing**: Applied functions from **Seurat** and **dplyr** for subsetting, normalizing, and statistically analyzing single-cell datasets.
- **Enrichment Analysis**: Conducted gene set enrichment analysis using **fgsea** to identify crucial biological pathways involved in the condition.

## Repository Structure

- `Code/`: Contains all R scripts for analysis.
  - `calculate_pyroptosis.r`: Scripts for calculating pyroptosis scores.
  - `figure*.r`: Scripts for generating figures illustrating the analysis.
  - `prep_*.r`: Scripts for data preprocessing from various data sources.
  - `reply_to_reviewers/`: Responses and revisions based on peer review.
- `data/`: Directory for raw and processed data files (as referenced in the scripts).
- `results/`: Output directory for results including statistical summaries and figures.

## Getting Started

To run these analyses:
1. Ensure R and all required packages are installed.
2. Clone this repository.
3. Execute scripts within the `Code/` directory in sequential order as listed.

## Contributions

Contributions are welcome. Please fork the repository and submit pull requests with any enhancements. For major changes, open an issue first to discuss what you would like to change.

## Citation

If you utilize this repository, please cite the associated publication:
```
Xu, Q., Yang, Y., Zhang, X., & Cai, J. J. (2022). Association of pyroptosis and severeness of COVID-19 as revealed by integrated single-cell transcriptome data analysis. ImmunoInformatics, 6, 100013.
```

## License

This project is open-sourced under the MIT License.
