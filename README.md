# *Cryophytum* Biogeography in South Africa: coding workflow

This repository contains the R code to reproduce the downstream analyses for the results in the paper:    

**Population genetics reveals insights into *Cryophytum* biogeography in South Africa**       
by: van Steenderen, C.J.M., Sandenbergh, E., and Paterson, I.D. 2026

The scripts:      

* **001_filtering_fulldata.R** filters the SNP output file (.vcf) from the Stacks workflow, using the SNPfiltR package. The **populations.snps.vcf** file is too large to host here, but is available on the Google Drive link: https://drive.google.com/file/d/1Bwo_qCxgxVFOgmslem3sMfiCX3D1Ez1k/view?usp=sharing
* **002_hybrid_analysis.R** runs population genetics statistics on the filtered SNP file, produces PCAs, maps, and fastSTRUCTURE plots
* **003_hybrids_climate_analysis.R** runs the analyses that investigate potential links between climate and edaphic variables and genetic distance between populations, and runs Mantel tests and correlation analyses

🌐The WordClim and SoilGRIDS variables need to be downloaded by the user using the R geodata package, as the data files are too large to provide in this repository.

All raw SNP data are available on the GenBank SRA under project ID **PRJNA1304995**.
