# Project
1. This study is to understand how Jun and Ddit3 (transcription factors) mediate damaging effects upon a controlled optic nerve crash (CONC), a model of glaucoma
1. The purpose of this study is to find novel gene targets to treat optic nerved degeneration as is commonly seen in glaucoma

## Experimental design
1. Mouse: 
	- Adult mice ~3 month of age
	- Genotype: Wild type control (Control), Jun knockout (KO), Ddit3 KO, and Jun/Ddit3 double knockout 
1. Procedure: 
	- Left eye of the mouse undergone a controlled optic nerve crash - CONC
	- Right eye of the mouse remain intact or did not touch - DNT
	- 72hrs later, retina tissue both eyes was harvested and subjected to sequencing

## Data pre-processing: 
1. fastq files of RNA sequencing was processed using standard JAX in-house pipeline and raw transcript count was generated at the end of the pipeline. 

## Data Analysis workflow
1. Raw transcript counts were converted into counts per million (CPM)
2. Used standard *edgeR*'s quasi-likelihood pipeline to determine differentially expressed (DE) (q<0.05) genes comparing CONC and DNT for each genotype (Control, Ddit3, Jun, Ddit3/Jun)
3. Upload DE genes for each genotype onto Ingenuity Pathway Analysis (IPA) for comparison analysis in Canonical Pathway, Upstream Regulator, Disease and Function and Regulatory Effect
	- to understand the molecular consequences of CONC
	- to determine how Jun and Ddit3 knockout can modify the consequences of CONC
	
## Guidance for reviewing analysis
1. Quality control and generating DE genes using *edgeR* check:
	- "glaucoma_ONC.Rmd" file
	- "glaucoma_ONC.pdf" file
	- "glaucoma_ONC.pdf" file
1. The resulting of gene list and brief KEGG pathway overview  in each comparison from Rmd file, check: 
	- "analysis" folder 
1. RNA seq QC figures, check:
	- "figure" folder
1. IPA analysis and draft figures for manuscript, check: 
	- "report_9_9_18.pdf"
	
