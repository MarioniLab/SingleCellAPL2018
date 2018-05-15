# scRNAseq

# To analyze pre-processed data

Create a directory called "data" within the current directory. Create a directory called "annotation_resources" within the current directory.

Download the idf, sdrf, and processed data from Array Express, accession number E-MTAB-6051.

Unzip the processed data and move both count tables to the data directory with the sdrf and idf files.

Obtain the necessary annotation files and put in the annotation_resources directory:
1. UCSC tables: ensGene and ensemblToGeneName tables for mm10 genome (I downloaded on 18/3/16); name them as follows: "ucsc_mm10_ensGene" and "ucsc_mm10_ensemblToGeneName"
2. List of mouse transcriptional regulatory genes from TF list in http://compbio.massey.ac.nz/apps/tcof/home/ (I downloaded 29 November 2016); name the file "BrowseTFTcoF-DB_Mouse.txt"

Run scripts in QC and analysis directories as described in their README files.

