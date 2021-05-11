# LIMMA-Python-implementation
This script is a python implementation of the Linear Models for Microarray Data (limma) package in R that helps perform differential gene expression analysis. Although limma was developed on microarray data, it's use is not limited to microarray data.  

INPUT - 2 files 
File 1 -- Expression data in a matrix, where each column represents an experiment or sample ID and the row represents a gene or probe.
File 2 -- Desgin matrix, where one column represents a status of the data (example:normal, cancer etc.) and other experiment or sample ID.

Have a look at the example files and make sure you have the right format.

One the files are setup, download the R packages to python. The script imports them to python, to do so specify the path where the  packages are downloaded.
For this analysis we need,rpy2 -- helps import packages from R to Python, limmma -- to run the analysis and writexl -- to export the output to an excel file.
