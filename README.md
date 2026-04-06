# *Cryophytum* Biogeography in South Africa: coding workflow
<img width="191" height="20" alt="image" src="https://github.com/user-attachments/assets/9ecfa720-6d91-42a8-a4b5-a2425c9f0dbb" />

This repository contains the R code to reproduce the downstream analyses for the results in the paper:    

**Population genetics reveals insights into *Cryophytum* biogeography in South Africa**     

by: **van Steenderen, C.J.M.**, Sandenbergh, E., and Paterson, I.D. 2026

Centre for Biological Control, Department of Zoology and Entomology, Rhodes University, Makhanda/Grahamstown

The scripts:      

* **001_filtering_fulldata.R** filters the SNP output file (.vcf) from the Stacks workflow, using the `SNPfiltR` package. The **populations.snps.vcf** file is too large to host here, but is available from the CBC's Google Drive link: https://drive.google.com/file/d/1Bwo_qCxgxVFOgmslem3sMfiCX3D1Ez1k/view?usp=sharing
* **002_hybrid_analysis.R** runs population genetics statistics on the filtered SNP file, produces PCAs, maps, and fastSTRUCTURE plots
* **003_hybrids_climate_analysis.R** runs the analyses that investigate potential links between climate and edaphic variables and genetic distance between populations, and runs Mantel tests and correlation analyses

🌐The WordClim and SoilGRIDS variables need to be downloaded by the user using the R `geodata` package, as the data files are too large to provide in this repository.

All raw SNP data are available on the GenBank SRA under project ID **PRJNA1304995**.
This repository is linked to a permanent Zenodo doi: [10.5281/zenodo.19438438](https://doi.org/10.5281/zenodo.19438438)
