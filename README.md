# spatial-infer-cnv

To run inferCNV, it requires three different files, the filtered_feature_bc_matrix.h5 located inside of the spaceranger path, annotations as a csv file that contains normal annotations, and a gene order file.

The filtered_feature_bc_matrix.h5 should be located inside of the spaceranger output. Collect this file path.

The gene order file I mainly use is in `/diskmnt/Projects/Users/efang/SpatialInferCNV/siCNV_GeneOrderFile.tsv` (this needs to be changed inside of the script if you want to edit this)

The annotations require a reference of normal cells annotated as "normal." If your file does not contain this, you can either manually annotate a section of normal spots (recommended 16 of them) or you can use my find_normal script described later.

The output directory can be edited inside of the script itself. In the SpatialCNVcmdline4.R, the output directory points to:
`f"/diskmnt/Projects/Users/efang/SpatialInferCNV/InferCNVRunsOutputValidation/InferCNVrun_outputs-{LibraryName}"`

## Running find_normal
To use the script, you need to have the necessary environment located at: `/diskmnt/Projects/Users/efang/miniconda3/envs/FindNormal2`

Then using conda:
```
conda activate FindNormal2
```
For the files that you want to find normal annotations for:

Create a spreadsheet named specific_df.csv that has the columns of: ["LibraryName", "AnnotationFile",  "Spaceranger"]

Place the df csv inside of `/diskmnt/Projects/Users/efang/find_normal_squidpy`

And then
```
python /diskmnt/Projects/Users/efang/find_normal_squidpy/main_specific.py
```

# How it works:
The script using the find_normal.py file to find the normal annotations. 

The idea is to create a buffer (increase the radius of the polygon) around annotated and suspected "tumor" areas and annotate cells outside of that area as "normal." In addition to already marked regions, I also use the epcam marker and also find a buffer radius around those cells where there is a nonzero epcam value. 

The script runs to find areas of normal 6 different times before giving up.

1. 3.5 spot buffer around the epcam_buffer and 2.5 spot buffer around the marked tumor regions.
2. 2-3.5 spot buffer that is uniformally chosen for the epcam buffer and 2.5 spot buffer around the marked tumor regions.
3. 3 spot buffer around the epcam_buffer and 2.5 spot buffer around the marked tumor regions.
4. 1.5-3 spot buffer that is uniformally chosen for the epcam buffer and 2.5 spot buffer around the marked tumor regions.
5. 1-2.5 spot buffer that is uniformally chosen for the epcam buffer and 2.5 spot buffer around the marked tumor regions.
6. 0.5-2 spot buffer that is uniformally chosen for the epcam buffer and 2.5 spot buffer around the marked tumor regions.

## Running the script
The conda environment is located at: `/diskmnt/Projects/Users/efang/miniconda3/envs/SpatialInferCNV`

Activate the environment with:
```
conda activate SpatialInferCNV
```

Now with the annotation file, filtered_feature_bc_matrix.h5, and the library name: run the R file located in `/diskmnt/Projects/Users/efang/SpatialInferCNV/SpatialCNVcmdline4.R`

An example of how it would work is:
```
Rscript /diskmnt/Projects/Users/efang/SpatialInferCNV/SpatialCNVcmdline4.R /diskmnt/Projects/Users/efang/normalized_annotations/HT413C1-Th1K2A4U1Bp1-with_normal_annotations.csv /diskmnt/Datasets/Spatial_Transcriptomics/outputs_FFPE/Human/HT413C1/Th1K2/HT413C1-Th1K2A4U1Bp1/outs/filtered_feature_bc_matrix.h5 HT413C1-Th1K2A4U1Bp1
```
