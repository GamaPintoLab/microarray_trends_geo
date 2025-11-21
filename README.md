# microarray_trends_geo
R code to retrieve platform- and assay-type–specific dataset counts from the NCBI Gene Expression Omnibus (GEO) and generate plots 

⸻

Microarray and HTS Trends in GEO (2000–2025)

This repository contains the R code used to retrieve platform- and assay-type–specific dataset counts from the NCBI Gene Expression Omnibus (GEO) and to generate the plots shown in Figure 1 of the corresponding book chapter. The analysis quantifies the historical use of major microarray technologies and the transition to high-throughput sequencing (HTS) between 2000 and 2025.

Overview

The scripts perform programmatic queries to GEO using the rentrez package, retrieving:
	•	Dataset counts by platform type (via the [PTYP] field).
	•	Dataset counts by assay type (via the [GTYP] field).
	•	Annual counts from 2000–2024 for each category.
The results are reshaped to long format and visualized using tidyverse tools.

The goal is to provide a reproducible workflow that reflects how platform usage and assay applications evolved over time.

Contents
	•	geo_microarray_query.R — Main script used to query GEO and aggregate counts.

Dependencies

The analysis was performed in R (version 4.3.1) using the following packages:

BiocManager
rentrez
tidyverse
ggbreak
scales

All packages are available from CRAN or Bioconductor.

Usage
	1.	Install dependencies:

install.packages("tidyverse")
install.packages("ggbreak")
install.packages("scales")
BiocManager::install("rentrez")


	2.	Run the script:

source("geo_microarray_query.R")


	3.	The script outputs:
	•	yearly counts for each platform or assay type
	•	long-format tables suitable for plotting
	•	ggplot2-based visualizations

Reproducibility

Queries rely on the NCBI E-utilities API through rentrez. Because GEO is updated continuously, counts may differ slightly if the script is run at a later date.

Citation

If you use this code or adapt the workflow, please cite:
  •	Gama-Carvalho and Sousa (2025). Microarrays and the dawn of Omics data analysis. in press
	•	Edgar et al. (2002). Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Research.
	•	Winter (2017). rentrez: an R package for the NCBI eUtils API. The R Journal.

DOI:

⸻

If you’d like, I can also produce a shorter or more formal version depending on the tone of the book chapter.
