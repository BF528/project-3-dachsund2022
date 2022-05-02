# Project Description

The transcriptomes of almost 12,000 individual pancreatic cells from four human donors and two mouse strains were determined using a droplet-based, single-cell RNA-seq approach developed by Baron et al. The inDrop method uses high-throughput droplet microfluidics to barcode the RNA from thousands of individual cells, implementing a sensitive linear amplification method for single-cell transcriptome library construction, and uses a systematic approach for capturing thousands of cells without pre-sorting. The cells were separated into 15 clusters that corresponded to previously identified cell types. Subpopulations within the ductal population, mechanisms of stellate cell activation, and heterogeneity in beta cell stress were discovered after a detailed investigation of each population independently.

The goal of this project was to use contemporary methodology and software to reproduce the conclusions of the publication. The objectives include analyzing the UMI counts to identify clusters and marker genes for distinct cell type populations, able to perform cell-by-gene quantification of UMI counts, performing quality control of the UMI counts matrix, which provide a biological meaning to the clustered cell types, and finally identifying novel genetic markers with them.

Workflow: Data Curator  >>>>>>  Programmer  >>>>>>>  Analyst  >>>>>>>>  Biologist

# Contributors

⧫ Vamshi Mallepalli (Data Curator)
⧫ Shreen Katyan (Programmer) 
⧫ Katherine Tu (Analyst) 
⧫ Sana Majid (Biologist)


# Repository Contents

Datacurator

A whitelist of barcodes is built using the raw reads in FASTQ file format to properly filter poor-quality samples by UMI count. The Salmon Alevin application, which is available on the command line, is then used to build a raw UMI count matrix.


Programmer

The Seurat conventional pre-processing workflow, is used to import and process the UMI count matrix generated above. Cells are grouped into subpopulations and low-quality readings are filtered.

Analyst

Using the marker genes supplied in Baron et alSupplementary .'s Data, the cell subpopulations are categorized into separate cell types. Each of these cell types' marker genes is kept, and a list of novel marker genes is exported for further investigation.

Biologist

The novel marker genes identified in the step above are examined for gene set enrichment, identifying key pathways that are over- or under-expressed, by cell type.


