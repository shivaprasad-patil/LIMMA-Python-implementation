# -*- coding: utf-8 -*-
"""
@author: Shivaprasad
"""

import sys
import click
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.rinterface import RRuntimeError
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from statsmodels.stats.multitest import multipletests

data =  pd.read_excel('expression_file.xlsx') #replace your own data file
data = data.set_index('ID') #replace 'ID' with your own annotation if necessary
design = pd.read_excel('limma_design_file.xlsx') #replace with your own design file

#Import R libraries
base = importr('base')
stats = importr('stats')
limma = importr('limma', lib_loc='path to the downloaded R package limma')
writexl = importr('writexl', lib_loc='path to the downloaded R package writexl')
    
# Convert data and design pandas dataframes to R dataframes
with localconverter(ro.default_converter + pandas2ri.converter):
    r_data = ro.conversion.py2ri(data)
    r_design = ro.conversion.py2ri(design)
    # Use the genes index column from data as a R String Vector
    genes = ro.StrVector(
        [
            str(index)
            #added tovalues to convert to numpy array
            for index in data.index.tolist()
            #for index in data.index.tolist()
        ]
    )

# Create a model matrix using design's Target column using the R formula "~0 + f" to get all the unique factors as columns
f = base.factor(r_design.rx2('Target'), levels=base.unique(r_design.rx2('Target')))
form = Formula('~0 + f')
form.environment['f'] = f
r_design = stats.model_matrix(form)
r_design.colnames = base.levels(f)

# Fit the data to the design using lmFit from limma
fit = limma.lmFit(r_data, r_design)
# Make a contrasts matrix with the 1st and the last unique values
contrast_matrix = limma.makeContrasts(f"{r_design.colnames[0]}-{r_design.colnames[-1]}", levels=r_design)

# Fit the contrasts matrix to the lmFit data & calculate the bayesian fit
fit2 = limma.contrasts_fit(fit, contrast_matrix)
fit2 = limma.eBayes(fit2)

# topTreat the bayesian fit using the contrasts and add the genelist
r_output = limma.topTreat(fit2, coef=1, genelist=genes, number=np.Inf)
writexl.write_xlsx(r_output, "limma_output.xlsx")
