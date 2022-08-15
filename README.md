# CellSegmentation

2D cell segmentation algorithms developed for human bronchial epithelial cells.

## Description

Human bronchial epithelial cells grown in air-liquid interface culture recapitulates pseudostratified, apicobasal structure found in vivo. The apical side includes differentiated cells (i.e. goblet and ciliated cells). The basal side includes mostly basal stem cells.

Apical or basal cells can be visualized from F-actin stained samples. From z-stack imaging through the apicobasal polarized epithelium, fields of view visualizing cell boundaries can be generated from maximum intensity projections.

The apical and basal cell segmentation use separate functions which both use marker-controlled watershed segmentation. The major differences include:

* `apicalseg` segments the cell boundaries followed by filtering incorrect segmentation (based on tortuosity, pixel intensities, and segment length) and filters by cell size
* `basalseg` segments the cell boundaries and filters by cell size

## Getting Started

### Dependencies

* Cell segmentation code was developed on a Windows 10 using MATLAB 2021a
* Images were pre-processed in Fiji/ImageJ 

### test_images

Samples images of cell boundaries for segmentation.

* Maximum intensity projection from F-actin z-stack images of pseudostratified epithelium
* Naming convention: Donor_Timepoint(hr)_Control/Pressure Treatment & Well #_Field of View_Region.tif
* Example: U13_72_P1_5_ACB.tif  
    * Donor: U13
    * Timepoint: 72 hour after treatment
    * Treatment: Pressure, Well #1
    * Field of View ID: #5
    * Region: Apical Cell Boundaries (ACB)

### lib

Library of segmentation functions

* `apicalseg` apical cell segmentation
* `basalseg` basal cell segmentation

### Executing program

The script `Shell_segment_cells.m` includes blocks of code using the segmentation and visualization functions.

## Authors

Thien-Khoi N. Phung
[@tkphung](https://twitter.com/tkphung)


## Acknowledgments

Inspiration, code snippets, etc.
* [Marker-controlled watershed segmentation](https://www.mathworks.com/help/images/marker-controlled-watershed-segmentation.html)