# LIMMA-Python-implementation
This script is a python implementation of the Linear Models for Microarray Data (limma) package in R that helps perform differential gene expression analysis. Although limma was developed on microarray data, it's use is not limited to microarray data.  

INPUT - 2 files.

File 1 -- Expression data in a matrix, where each column represents an experiment or sample ID and the row represents a gene or probe expression.

File 2 -- Desgin matrix, where one column represents a status of the data (example:normal, cancer etc.) and other experiment or sample ID. 
Have a look at the design file and make sure you have the right format. The 'Target'  column has 'zero' and 'one', which specify status of the data (example:normal, cancer etc).

Ones the files are setup, download the required packages. The script imports them to python, to do so specify the path where the  packages are downloaded.

Packages to download:

Python --> NumPy and Pandas.

R --> rpy2 -- helps import packages from R to Python, limma -- to run the analysis and writexl -- to export the output to an excel file.

If you run into any trouble do let me know.
