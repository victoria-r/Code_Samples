# Code_Samples

This repository contains samples of relevant code I have worked on for assignments and class projects. 

These files contain information that may not be used to copy as they may still be part of relevant academic assignments. Please, if you are having trouble or are stuck in your code, reach out to your professor for guidence.

Thank you.

---------------------------------------------------------------------------------------------------------------------------------------------------------

# Cancer Genomics

## cancerGenomics.Rmd

RNnotebook file containing all the code and outputs from following the maftools vingette tutorial.

## cancerGenomics.sh

Bash script to trnder the RNotebook to an html file. 

## Getting Started

### Dependencies

* bash
* maftools 
* tcga_laml.maf.gz - TCGA LAML MAF file
* tcga_laml.annot.tsv - Metadata containing clinical information of survival history and histology
* LAML_sig_genes.txt.gz  
* APL_primary.maf.gz
* APL_relapse.maf.gz
* brca.maf.gz
* R
* RStudio - IDE
* Tutorial - https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html

## Results

View canerGenomics.md or cancerGenomics.html 

## Author

Victoria R. Liebsch-Aljawahiri

## Date Created

22 April 2022

## References 

Mayakonda, Anand, De-Chen Lin, Yassen Assenov, Christoph Plass, and H. Phillip Koeffler. 2018. “Maftools: Efficient and Comprehensive Analysis of Somatic Variants in Cancer.” *Genome Research* 28 (11): 1747–56..

Network, Cancer Genome Atlas Research, Timothy J. Ley, Christopher Miller, Li Ding, Benjamin J. Raphael, Andrew J. Mungall, A. Gordon Robertson, et al. 2013. “Genomic and Epigenomic Landscapes of Adult de Novo Acute Myeloid Leukemia.” *The New England Journal of Medicine* 368 (22): 2059–74..

---------------------------------------------------------------------------------------------------------------------------------------------------------

# translate_APOE.py

This script translates the APOE ortholog sequences from nucleotide sequences to amino acid sequences using BioPython. Creates a fasta file called apoe_aa.fasta that contains the translated orthologs. 

# align_APOE.sh

This script uses Clustal Omega to align the translated APOE orthologs from apoe_aa.fatsa. Creates a fasta file called aligned_apoe.fasta that contains the aligned orthologs.

## Getting Started

### Dependencies

* python3
* biopython
* clustalo

## Author

Victoria R. Liebsch

## Date Created - All programs

5 March 2022
