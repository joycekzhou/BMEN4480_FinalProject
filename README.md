# BMEN4480_FinalProject

## microglia_clusters

Description:
Analysis of microglia scRNA-seq dataset to cluster and characterize microglial activation states

Input data:
Sousa, C., Golebiewska, A., Poovathingal, S.K. et al. Single-cell transcriptomics reveals distinct inflammation-induced microglia signatures. EMBO Reports (2018).

Output data:
Raw counts matrix and metadata csv files for subsequent RCTD analysis

## merge_microglia_all

Description: Merge microglia and primary motor cortex scRNA-seq datasets

Input data:
1) Sousa, C., Golebiewska, A., Poovathingal, S.K. et al. Single-cell transcriptomics reveals distinct inflammation-induced microglia signatures. EMBO Reports (2018).
2) Yao, Z., Liu, H., Xie, F. et al. An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types. bioRxiv (2020). Seurat object downloaded from from https://azimuth.hubmapconsortium.org/references/#Mouse%20-%20Motor%20Cortex


## microglialmarkerexpression
Located: In Spatial Analysis Folder

Description: Visualize top differentially expressed genes on spatial images using Seurat in RStudio

Input Data:
1) Hasel, P., Rose I.V.L., Sadick, J.S., Kim. R.D. et al. Neuroinflammatory astrocyte subtypes in the mouse brain. Nat Neurosci (2021)

## Running Deconvolution Algorithm with Custom Clusters
Primary script is called BMEN4480_Project_main.rmd, located in the Deconvolution_scripts_data folder

In order to run the deconvolution algorithm follow these steps:
1. Download the custom count matrices located here: https://drive.google.com/drive/folders/1ApGbanQ6f1Puu1QtiRCrkuIvmgzCCwh5?usp=sharing and move them into the Deconvolution_scripts_data folder
2. Make sure to install the necessary R dependencies (spacexr, Matrix, doParallel, and ggplot2)
3. Update the file locations to reflect the location this data is stored.
4. Run the markdown cells in the order they occur in the file.

General RCTD help can be found in the updated package, spacexr: https://github.com/dmcable/spacexr
