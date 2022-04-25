# Project: Deconvolving Rare Cell States from Spatial Transcriptomic Data
# Names: Anish Puligilla, Cody Slater, Joyce Zhou

rm(list = ls())

#Visualizing Microglial Gene Markers for Each Dataset
#LPS 1
setwd(".../spatial") #Set working directory where lowres image is contained. Should be within a folder called 'spatial'
LPS1_img <- Read10X_Image(image.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/LPS_1_Visium/spatial', image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
LPS1_spatdata <- Load10X_Spatial(data.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/LPS_1_Visium', filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE, image = NULL)
LPS1_spatdata <- SCTransform(LPS1_spatdata, assay = "Spatial", verbose = FALSE)

LPS1_plot1 <- SpatialFeaturePlot(LPS1_spatdata, features = c("P2ry12", "Cst3", "Siglech", "Gpr34","Ctsd"))
LPS1_plot2 <- SpatialFeaturePlot(LPS1_spatdata, features = c("Fth1","Ms4a6c","lgsf6","Msr1","Hspa5"))
LPS1_plot3 <- SpatialFeaturePlot(LPS1_spatdata, features = c("Malat1","Gm26924","Dst","Zeb2","Cacna1d"))

#LPS 2
setwd(".../spatial") #Set working directory where lowres image is contained. Should be within a folder called 'spatial'
LPS2_img <- Read10X_Image(image.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/LPS_2_Visium/spatial', image.name = "tissue_hires_image.png", filter.matrix = TRUE)
LPS2_spatdata <- Load10X_Spatial(data.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/LPS_2_Visium', filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE, image = LPS2_img)
LPS2_spatdata <- SCTransform(LPS2_spatdata, assay = "Spatial", verbose = FALSE)

LPS2_plot1 <- SpatialFeaturePlot(LPS2_spatdata, features = c("P2ry12", "Cst3", "Siglech", "Gpr34","Ctsd"))
LPS2_plot2 <- SpatialFeaturePlot(LPS2_spatdata, features = c("Fth1","Ms4a6c","lgsf6","Msr1","Hspa5"))
LPS2_plot3 <- SpatialFeaturePlot(LPS2_spatdata, features = c("Malat1","Gm26924","Dst","Zeb2","Cacna1d"))

#LPS 3
setwd(".../spatial") #Set working directory where lowres image is contained. Should be within a folder called 'spatial'
LPS3_img <- Read10X_Image(image.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/LPS_3_Visium/spatial', image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
LPS3_spatdata <- Load10X_Spatial(data.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/LPS_3_Visium', filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE, image = NULL)
LPS3_spatdata <- SCTransform(LPS3_spatdata, assay = "Spatial", verbose = FALSE)

LPS3_plot1 <- SpatialFeaturePlot(LPS3_spatdata, features = c("P2ry12", "Cst3", "Siglech", "Gpr34","Ctsd"))
LPS3_plot2 <- SpatialFeaturePlot(LPS3_spatdata, features = c("Fth1","Ms4a6c","lgsf6","Msr1","Hspa5"))
LPS3_plot3 <- SpatialFeaturePlot(LPS3_spatdata, features = c("Malat1","Gm26924","Dst","Zeb2","Cacna1d"))

#Sal 1
setwd(".../spatial") #Set working directory where lowres image is contained. Should be within a folder called 'spatial'
Sal1_img <- Read10X_Image(image.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/Saline_1_Visium/spatial', image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
Sal1_spatdata <- Load10X_Spatial(data.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/Saline_1_Visium', filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE, image = NULL)
Sal1_spatdata <- SCTransform(Sal1_spatdata, assay = "Spatial", verbose = FALSE)

Sal1_plot1 <- SpatialFeaturePlot(Sal1_spatdata, features = c("P2ry12", "Cst3", "Siglech", "Gpr34","Ctsd"))
Sal1_plot2 <- SpatialFeaturePlot(Sal1_spatdata, features = c("Fth1","Ms4a6c","lgsf6","Msr1","Hspa5"))
Sal1_plot3 <- SpatialFeaturePlot(Sal1_spatdata, features = c("Malat1","Gm26924","Dst","Zeb2","Cacna1d"))

#Sal 2
setwd(".../spatial") #Set working directory where lowres image is contained. Should be within a folder called 'spatial'
Sal2_img <- Read10X_Image(image.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/Saline_2_Visium/spatial', image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
Sal2_spatdata <- Load10X_Spatial(data.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/Saline_2_Visium', filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE, image = NULL)
Sal2_spatdata <- SCTransform(Sal2_spatdata, assay = "Spatial", verbose = FALSE)

Sal2_plot1 <- SpatialFeaturePlot(Sal2_spatdata, features = c("P2ry12", "Cst3", "Siglech", "Gpr34","Ctsd"))
Sal2_plot2 <- SpatialFeaturePlot(Sal2_spatdata, features = c("Fth1","Ms4a6c","lgsf6","Msr1","Hspa5"))
Sal2_plot3 <- SpatialFeaturePlot(Sal2_spatdata, features = c("Malat1","Gm26924","Dst","Zeb2","Cacna1d"))

#Sal 3
setwd(".../spatial") #Set working directory where lowres image is contained. Should be within a folder called 'spatial'
Sal3_img <- Read10X_Image(image.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/Saline_3_Visium/spatial', image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
Sal3_spatdata <- Load10X_Spatial(data.dir = 'C:/Users/Anish_Columbia.ANISH-PC/Desktop/Spatial_Data/Saline_3_Visium', filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE, image = NULL)
Sal3_spatdata <- SCTransform(Sal3_spatdata, assay = "Spatial", verbose = FALSE)

Sal3_plot1 <- SpatialFeaturePlot(Sal3_spatdata, features = c("P2ry12", "Cst3", "Siglech", "Gpr34","Ctsd"))
Sal3_plot2 <- SpatialFeaturePlot(Sal3_spatdata, features = c("Fth1","Ms4a6c","lgsf6","Msr1","Hspa5"))
Sal3_plot3 <- SpatialFeaturePlot(Sal3_spatdata, features = c("Malat1","Gm26924","Dst","Zeb2","Cacna1d"))
