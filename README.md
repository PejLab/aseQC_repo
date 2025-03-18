# aseQC

**A**llele **S**pecific **E**xpression **Q**uality **C**ontrol.

This package has been authored by [Kaushik Ganapathy](kganapathy@scripps.edu) and [Eric Song](ejsong@ucsd.edu).

### Input Files

To install this package:

* Clone this [repository](https://github.com/Krganapa/aseQC_repo/) 
* Load in the devtools library in R using ```library(devtools)```
* Install the package using ```devtools::install("path_to_cloned_repo")```

### Input Files

Both refcounts and altcounts are R Dataframes with an index column and a column `name` containing the names of genes in addition to the columns for each sample. Ref counts will have the reference allele counts, and Alt counts will have the alternative allele counts from ASE Data. 

The format is outlined below:

|   | name   | sample 1 | sample 2 | ... | sample n |
| - | ------ | -------- | -------- | --- | -------- |
| 0 | Gene 1 | Count 1  | Count 2  | ... | Count n  |
| 1 | Gene 2 | Count 3  | Count 4  | ... | Count m  |
| ... | ... | ... | ... | ... | ... |

Both the index (unnamed) column and name column are **essential** for performing ASEQC fits.

### Output File

The output file is an R dataframe with several columns. The `std` column indicates the ASEQC score, and the `cutoff` indicates the cutoff score above which samples are considered to fail the test. 
